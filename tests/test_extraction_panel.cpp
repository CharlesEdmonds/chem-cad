// Headless render test for the Extraction Lab panel. Drives the real
// drawExtractionLab() through a null ImGui backend, so the cross-section,
// the suite import hand-off and the solute distribution panel are proven to
// work without a display.
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include <string>

#include "imgui.h"
#include "imgui_internal.h"

#include "chem/bridge.hpp"
#include "fluid/simulation.hpp"
#include "sol/funnel.hpp"
#include "sol/solubility.hpp"
#include "sol/solvent.hpp"
#include "ui/app_state.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"

using namespace chemcad;

namespace {

constexpr float kDisplayW = 1600.0f;
constexpr float kDisplayH = 1200.0f;

struct HeadlessImGui {
  explicit HeadlessImGui(float width = kDisplayW, float height = kDisplayH,
                         float fontScale = 1.0f) {
    IMGUI_CHECKVERSION();
    ctx = ImGui::CreateContext();
    ImGui::SetCurrentContext(ctx);
    ImGuiIO& io = ImGui::GetIO();
    io.DisplaySize = ImVec2(width, height);
    io.DeltaTime = 1.0f / 60.0f;
    io.IniFilename = nullptr;
    io.LogFilename = nullptr;
    io.FontGlobalScale = fontScale;
    unsigned char* pixels = nullptr;
    int w = 0, h = 0;
    io.Fonts->GetTexDataAsRGBA32(&pixels, &w, &h);
    io.Fonts->SetTexID(static_cast<ImTextureID>(1));
  }
  ~HeadlessImGui() { ImGui::DestroyContext(ctx); }
  HeadlessImGui(const HeadlessImGui&) = delete;
  HeadlessImGui& operator=(const HeadlessImGui&) = delete;
  ImGuiContext* ctx = nullptr;
};

int panelFrame(ui::AppState& st) {
  const ImGuiIO& io = ImGui::GetIO();
  ImGui::NewFrame();
  ImGui::SetNextWindowPos(ImVec2(0, 0));
  ImGui::SetNextWindowSize(io.DisplaySize);
  ImGui::Begin("Extraction Lab", nullptr,
               ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                   ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings);
  ui::drawExtractionLab(st);
  ImGui::End();
  ImGui::Render();
  return ImGui::GetDrawData()->TotalVtxCount;
}

int settledVertexCount(ui::AppState& st) {
  panelFrame(st);
  return panelFrame(st);
}

// What one frame of the panel actually produced: the text ImGui rendered (only
// text that reached the screen is logged, so this observes what the user reads)
// plus the geometry it emitted. `nonFinite` catches a NaN leaking out of a
// readout into vertex coordinates, which a vertex count cannot see; `overflow`
// is how far the page could be scrolled, which must always be zero.
struct FrameReport {
  std::string text;
  int vertices = 0;
  int nonFinite = 0;
  float overflow = 0.0f;

  bool says(const char* fragment) const {
    return text.find(fragment) != std::string::npos;
  }
};

FrameReport loggedFrame(ui::AppState& st) {
  FrameReport report;
  const ImGuiIO& io = ImGui::GetIO();
  ImGui::NewFrame();
  ImGui::SetNextWindowPos(ImVec2(0, 0));
  ImGui::SetNextWindowSize(io.DisplaySize);
  ImGui::Begin("Extraction Lab", nullptr,
               ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                   ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings);
  ImGui::LogToBuffer();
  ui::drawExtractionLab(st);
  // LogFinish() clears the buffer, so the copy has to happen first.
  report.text = ImGui::GetCurrentContext()->LogBuffer.c_str();
  ImGui::LogFinish();
  report.overflow = ImGui::GetScrollMaxY();
  ImGui::End();
  ImGui::Render();

  const ImDrawData* data = ImGui::GetDrawData();
  report.vertices = data->TotalVtxCount;
  for (int list = 0; list < data->CmdListsCount; ++list) {
    const ImDrawList* drawList = data->CmdLists[list];
    for (int vertex = 0; vertex < drawList->VtxBuffer.Size; ++vertex) {
      const ImVec2& position = drawList->VtxBuffer[vertex].pos;
      if (!std::isfinite(position.x) || !std::isfinite(position.y))
        ++report.nonFinite;
    }
  }
  return report;
}

// Two frames: the first builds the items the second measures and hovers.
FrameReport panelReport(ui::AppState& st) {
  panelFrame(st);
  return loggedFrame(st);
}

// ImGui can only hover an item that was submitted with the pointer already over
// it, so the position is established on one frame and read on the next.
FrameReport hoverAt(ui::AppState& st, ImVec2 point) {
  ImGui::GetIO().AddMousePosEvent(point.x, point.y);
  panelFrame(st);
  return loggedFrame(st);
}

// Full left click at `point`, reported on the release frame. Some widgets fire
// on press (IsItemClicked) and some on release (ImGui's combos); by the release
// frame both have run and anything they opened is on screen.
FrameReport clickAt(ui::AppState& st, ImVec2 point) {
  ImGuiIO& io = ImGui::GetIO();
  io.AddMousePosEvent(point.x, point.y);
  panelFrame(st);  // the item is submitted with the pointer already over it
  io.AddMouseButtonEvent(ImGuiMouseButton_Left, true);
  panelFrame(st);
  io.AddMouseButtonEvent(ImGuiMouseButton_Left, false);
  return loggedFrame(st);
}

// The charge editor's phase table, found in ImGui's persistent table pool by
// its column count. Synthetic clicks are then aimed at real cells rather than
// at guessed pixels.
const ImGuiTable* findTable(int columns) {
  ImGuiContext& g = *ImGui::GetCurrentContext();
  for (int i = 0; i < g.Tables.GetMapSize(); ++i)
    if (ImGuiTable* table = g.Tables.TryGetMapData(i))
      if (table->ColumnsCount == columns) return table;
  return nullptr;
}

// Click targets derived from the tables ImGui just drew, so a synthetic press
// lands on the real control. RowPosY1 is the last row submitted, which is the
// row these tests act on.
ImVec2 rowCellCentre(const ImGuiTable& table, int column) {
  const float glyph = ImGui::GetFontSize();
  return ImVec2(table.Columns[column].WorkMinX + glyph,
                table.RowPosY1 + table.RowCellPaddingY + glyph * 0.5f);
}

// The trash button sits one glyph plus a quarter gap past the colour swatch,
// which is drawn at the start of the table's last column.
ImVec2 removeButtonCentre(const ImGuiTable& table) {
  const float glyph = ImGui::GetFontSize();
  return ImVec2(table.Columns[6].WorkMinX + glyph * 1.5f +
                    ui::style::metrics().gap * 0.25f,
                table.RowPosY1 + table.RowCellPaddingY + glyph * 0.5f);
}

// The solver-quality combo, under the "Solver quality" caption in the last
// column of the transport grid.
ImVec2 qualityComboCentre(const ImGuiTable& table) {
  return ImVec2(table.Columns[3].WorkMinX + ImGui::GetFontSize() * 3.0f,
                table.RowPosY1 + table.RowCellPaddingY +
                    ImGui::GetTextLineHeightWithSpacing() +
                    ImGui::GetFrameHeight() * 0.5f);
}

constexpr double kPi = 3.14159265358979323846;

sol::Simulation physicsSimulation() {
  sol::Simulation sim;
  sim.vessel = sol::Vessel::SeparatoryFunnel;
  sim.vesselVolumeMl = 250.0;

  sol::Phase aqueous;
  aqueous.label = "Aqueous";
  aqueous.volumeMl = 125.0;
  aqueous.density = 1.10;
  aqueous.viscosity = 1.2;
  aqueous.interfacialTension = 30.0;
  aqueous.emulsionStability = 0.25;

  sol::Phase organic;
  organic.label = "Organic";
  organic.volumeMl = 125.0;
  organic.density = 0.82;
  organic.viscosity = 1.0;
  organic.interfacialTension = 30.0;
  organic.emulsionStability = 0.25;

  sim.phases = {aqueous, organic};
  sol::reset(sim);
  return sim;
}

double representedVolumeMl(const sol::Simulation& sim) {
  double total = 0.0;
  for (double settled : sim.settledMl) total += settled;
  for (const sol::Droplet& droplet : sim.droplets) total += double(droplet.parcelMl);
  return total;
}

double cloudRadiusM(const sol::Droplet& droplet) {
  return std::cbrt(3.0 * double(droplet.parcelMl) * 1e-6 / (4.0 * kPi));
}

double outlineMaxHalfWidth(sol::Vessel vessel, double height) {
  double maxWidth = 0.0;
  for (const core::Vec2& point : sol::vesselOutline(vessel, height)) {
    maxWidth = std::max(maxWidth, std::abs(double(point.x)));
  }
  return maxWidth;
}

}  // namespace

TEST_CASE("the vessel stage draws and responds to shaking") {
  HeadlessImGui gui;
  ui::AppState st;
  // The stage is permanently visible beside the console, and until the
  // particle solver has been built it draws the analytic vessel -- which is
  // exactly the state a headless test runs in.
  const int settled = settledVertexCount(st);

  // Two default phases are seeded on the first frame; five seconds into a
  // firm 6 s shake the column is full of droplets, every one of which is an
  // extra filled circle on the stage.
  REQUIRE(st.solubility.funnel.phases.size() == 2);
  sol::shake(st.solubility.funnel, sol::ShakeParams{6.0, 3.5, 0.06});
  for (int i = 0; i < 100; ++i) sol::step(st.solubility.funnel, 0.05);  // 5 s
  REQUIRE(st.solubility.funnel.shake.active);  // still shaking at 5 s of 6 s
  REQUIRE(sol::emulsifiedFraction(st.solubility.funnel) > 0.2);
  const int shaken = settledVertexCount(st);

  CHECK(shaken > settled);
}

TEST_CASE("a suite import replaces the charged phases, densest first") {
  HeadlessImGui gui;
  ui::AppState st;
  const sol::Solvent* water = sol::findSolvent("water");
  const sol::Solvent* dcm = sol::findSolvent("dcm");
  REQUIRE(water != nullptr);
  REQUIRE(dcm != nullptr);

  st.solubility.extractionImport.pending = true;
  st.solubility.extractionImport.solventIdA = "water";
  st.solubility.extractionImport.solventIdB = "dcm";
  st.solubility.extractionImport.volumeMlA = 60.0;
  st.solubility.extractionImport.volumeMlB = 40.0;
  panelFrame(st);

  CHECK(!st.solubility.extractionImport.pending);  // consumed exactly once
  const std::vector<sol::Phase>& phases = st.solubility.funnel.phases;
  REQUIRE(phases.size() == 2);
  // sol::reset sorts densest-first, so DCM leads even though water arrived
  // first -- the bottom of the funnel is index 0.
  CHECK(phases[0].label == dcm->name);
  CHECK(phases[1].label == water->name);
  CHECK(phases[0].density == doctest::Approx(dcm->density));
  CHECK(phases[1].density == doctest::Approx(water->density));
  CHECK(phases[0].volumeMl == doctest::Approx(40.0));
  CHECK(phases[1].volumeMl == doctest::Approx(60.0));
  CHECK(sol::totalVolumeMl(st.solubility.funnel) == doctest::Approx(100.0));

  // A second frame must not re-apply the consumed import.
  CHECK_NOTHROW(panelFrame(st));
  CHECK(st.solubility.funnel.phases.size() == 2);
}

TEST_CASE("a bad import id keeps the default phases") {
  HeadlessImGui gui;
  ui::AppState st;
  st.solubility.extractionImport.pending = true;
  st.solubility.extractionImport.solventIdA = "no-such-solvent";
  st.solubility.extractionImport.solventIdB = "water";
  CHECK_NOTHROW(panelFrame(st));

  CHECK(!st.solubility.extractionImport.pending);
  const std::vector<sol::Phase>& phases = st.solubility.funnel.phases;
  REQUIRE(phases.size() == 2);
  CHECK(phases[0].label != "no-such-solvent");
}

TEST_CASE("the solute distribution panel draws when a solute is available") {
  HeadlessImGui gui;

  // Partitioning is reported on the Charge tab beside the phases it splits.
  ui::AppState bare;
  bare.solubility.extractionTab = 1;
  const int bareVertices = settledVertexCount(bare);

  ui::AppState withSolute;
  withSolute.solubility.extractionTab = 1;
  withSolute.doc.molecules.push_back(chem::fromSmiles("OC(=O)c1ccccc1"));
  withSolute.touch();
  withSolute.solubility.solute =
      sol::describeSolute(withSolute.doc.molecules.front());
  withSolute.solubility.soluteValid = true;
  const int soluteVertices = settledVertexCount(withSolute);

  // The stacked partition bar and its labels are real geometry the bare
  // workspace does not draw.
  CHECK(soluteVertices > bareVertices);
}

TEST_CASE("parcel bookkeeping conserves volume through shake and long settling") {
  sol::Simulation sim = physicsSimulation();
  const double charged = sol::totalVolumeMl(sim);
  sol::shake(sim, sol::ShakeParams{5.0, 3.0, 0.05});

  for (int i = 0; i < 100; ++i) sol::step(sim, 0.05);
  REQUIRE(!sim.droplets.empty());
  CHECK(sim.droplets.size() >= 600);
  CHECK(sim.droplets.size() <= 1200);
  CHECK(sim.shake.sauterRadiusM >= 100e-6);
  CHECK(sim.shake.sauterRadiusM <= 400e-6);
  CHECK(std::abs(representedVolumeMl(sim) - charged) <= 1e-9);

  for (int i = 0; i < 1200; ++i) sol::step(sim, 0.05);
  CHECK(std::abs(representedVolumeMl(sim) - charged) <= 1e-9);
}

TEST_CASE("parcel clouds remain inside their analytic vessel over thirty seconds") {
  sol::Simulation sim = physicsSimulation();
  sol::shake(sim, sol::ShakeParams{5.0, 3.0, 0.05});
  const double height = sol::columnHeightM(sim);
  const double maxHalfWidth = outlineMaxHalfWidth(sim.vessel, height);
  constexpr double kGeometryTolerance = 2e-6;

  for (int frame = 0; frame < 600; ++frame) {
    sol::step(sim, 0.05);
    for (const sol::Droplet& droplet : sim.droplets) {
      const double cloud = cloudRadiusM(droplet);
      const double y = double(droplet.position.y);
      const double wall =
          maxHalfWidth * sol::vesselWidthAt(sim.vessel, y / height);
      // doctest 2.4.11 declares CAPTURE with a single parameter, so each value
      // gets its own call rather than relying on a variadic macro.
      CAPTURE(frame);
      CAPTURE(droplet.phase);
      CAPTURE(droplet.position.x);
      CAPTURE(droplet.position.y);
      CAPTURE(droplet.radius);
      CAPTURE(droplet.parcelMl);
      CHECK(droplet.radius >= 20e-6f);
      CHECK(droplet.radius <= 3e-3f);
      CHECK(y >= 0.0);
      CHECK(y <= height);
      CHECK(y - cloud >= -kGeometryTolerance);
      CHECK(y + cloud <= height + kGeometryTolerance);
      CHECK(std::abs(double(droplet.position.x)) + cloud <=
            wall + kGeometryTolerance);
    }
  }
}

TEST_CASE("column height and analytic profile reproduce rated vessel capacity") {
  constexpr std::array<sol::Vessel, 3> kVessels = {
      sol::Vessel::SeparatoryFunnel,
      sol::Vessel::DecantingFlask,
      sol::Vessel::GraduatedCylinder,
  };
  constexpr int kIntegrationSteps = 4096;

  for (sol::Vessel vessel : kVessels) {
    sol::Simulation sim;
    sim.vessel = vessel;
    sim.vesselVolumeMl = 250.0;
    const double height = sol::columnHeightM(sim);
    const double maxHalfWidth = outlineMaxHalfWidth(vessel, height);
    const std::vector<core::Vec2> outline = sol::vesselOutline(vessel, height);
    REQUIRE(outline.size() >= 515);
    REQUIRE(outline.front().x == outline.back().x);
    REQUIRE(outline.front().y == outline.back().y);
    double signedAreaTwice = 0.0;
    for (size_t i = 1; i < outline.size(); ++i) {
      signedAreaTwice += double(outline[i - 1].x) * double(outline[i].y) -
                         double(outline[i].x) * double(outline[i - 1].y);
    }
    CHECK(signedAreaTwice > 0.0);

    double integral = 0.0;
    for (int i = 0; i <= kIntegrationSteps; ++i) {
      const double t = double(i) / kIntegrationSteps;
      const double width = sol::vesselWidthAt(vessel, t);
      const double weight =
          (i == 0 || i == kIntegrationSteps) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);
      integral += weight * width * width;
    }
    integral /= 3.0 * kIntegrationSteps;
    const double integratedMl =
        kPi * maxHalfWidth * maxHalfWidth * height * integral * 1e6;
CAPTURE(static_cast<int>(vessel));
    CAPTURE(height);
    CAPTURE(maxHalfWidth);
    CAPTURE(integratedMl);
    CHECK(integratedMl == doctest::Approx(sim.vesselVolumeMl).epsilon(0.01));
  }
}

TEST_CASE("Squibb profile widens to its high shoulder and never forms a kite") {
  constexpr int kSamples = 720;
  double previous = sol::vesselWidthAt(sol::Vessel::SeparatoryFunnel, 0.0);
  for (int i = 1; i <= kSamples; ++i) {
    const double t = 0.72 * double(i) / kSamples;
    const double width = sol::vesselWidthAt(sol::Vessel::SeparatoryFunnel, t);
CAPTURE(t);
    CAPTURE(width);
    CAPTURE(previous);
    CHECK(width >= previous - 1e-12);
    previous = width;
  }
  CHECK(previous == doctest::Approx(1.0).epsilon(1e-12));

  for (int i = 1; i <= kSamples; ++i) {
    const double t = 0.72 + 0.14 * double(i) / kSamples;
    const double width = sol::vesselWidthAt(sol::Vessel::SeparatoryFunnel, t);
CAPTURE(t);
    CAPTURE(width);
    CAPTURE(previous);
    CHECK(width <= previous + 1e-12);
    previous = width;
  }
  CHECK(previous == doctest::Approx(0.20).epsilon(1e-12));
  CHECK(sol::vesselWidthAt(sol::Vessel::SeparatoryFunnel, 1.0) ==
        doctest::Approx(0.20));
}

// --------------------------------------------------------------- audit cases
// Everything below drives a real control or reads a real readout and asserts on
// what the panel put on screen, never on how it drew it.

TEST_CASE("the schematic hint quotes the slab the section actually cuts") {
  // drawLiquidSection keeps a particle when |y| < 0.75 * particleSpacingM, and
  // a Snapshot reports the SAME length for spacing and render radius: an SPH
  // free surface sits one dx out, not dx/2 (fluid/simulation.hpp). "1.5
  // particle radii" therefore advertised twice the slab that is drawn.
  const fluid::Snapshot snapshot;
  CHECK(snapshot.particleRadiusM == doctest::Approx(snapshot.particleSpacingM));

  HeadlessImGui gui;
  ui::AppState st;
  st.solubility.extractionRenderMode = ui::ExtractionRenderMode::Schematic2D;
  const FrameReport schematic = panelReport(st);
  CHECK(schematic.says("particle cut: |y| < 0.75 particle spacings"));
  CHECK_FALSE(schematic.says("1.5 particle radii"));

  // The solid stage has no section, so it must not claim one -- and the hint
  // has to change with the mode rather than being one fixed string.
  st.solubility.extractionRenderMode = ui::ExtractionRenderMode::Fluid3D;
  const FrameReport solid = panelReport(st);
  CHECK_FALSE(solid.says("particle cut"));
  CHECK(solid.says("right-drag orbit"));
}

TEST_CASE("the charge table leaves no cursor position ImGui rejects") {
  HeadlessImGui gui;
  ui::AppState st;
  st.solubility.extractionTab = 1;
  const FrameReport charge = panelReport(st);

  // Each cell activator used to restore the cursor one item-spacing past the
  // row's content extent and submit nothing there, which ImGui reports once per
  // cell every frame -- and asserts on in a debug build.
  INFO(charge.text);
  CHECK_FALSE(charge.says("MESSAGE FROM DEAR IMGUI"));
  CHECK_FALSE(charge.says("SetCursorPos"));
  CHECK(charge.vertices > 0);
  CHECK(charge.nonFinite == 0);
}

TEST_CASE("a click on a charge cell still opens the phase editor") {
  HeadlessImGui gui;
  ui::AppState st;
  st.solubility.extractionTab = 1;
  const FrameReport idle = panelReport(st);
  REQUIRE_FALSE(idle.says("Edit charged phase"));

  const ImGuiTable* table = findTable(7);
  REQUIRE(table != nullptr);
  const FrameReport opened = clickAt(st, rowCellCentre(*table, 0));
  CHECK(opened.says("Edit charged phase"));
  CHECK(opened.says("Density (g/mL)"));
  CHECK(opened.vertices > idle.vertices);
}

TEST_CASE("the last charged phase cannot be removed and the panel says why") {
  HeadlessImGui gui;
  ui::AppState st;
  st.solubility.extractionTab = 1;
  panelReport(st);
  REQUIRE(st.solubility.funnel.phases.size() == 2);

  const ImGuiTable* pair = findTable(7);
  REQUIRE(pair != nullptr);
  const ImVec2 removeSecond = removeButtonCentre(*pair);
  CHECK(hoverAt(st, removeSecond).says("Remove phase"));
  clickAt(st, removeSecond);
  REQUIRE(st.solubility.funnel.phases.size() == 1);
  const std::string survivor = st.solubility.funnel.phases.front().label;

  // drawExtractionLab re-seeds a default Aqueous + Dichloromethane charge
  // whenever the table empties, so removing the last phase would swap the
  // liquid for a different pair rather than remove it.
  panelReport(st);
  const ImGuiTable* single = findTable(7);
  REQUIRE(single != nullptr);
  const ImVec2 removeLast = removeButtonCentre(*single);
  CHECK(hoverAt(st, removeLast)
            .says("The vessel needs at least one charged phase."));
  clickAt(st, removeLast);
  REQUIRE(st.solubility.funnel.phases.size() == 1);
  CHECK(st.solubility.funnel.phases.front().label == survivor);
}

TEST_CASE("the wash plan agrees with the partition formula it is drawn from") {
  // Worked by hand, logP 2.00 with 100 mL aqueous against 50 mL organic and
  // 100 mg of solute:
  //   D          = 10^2.00 = 100
  //   D * V_org  = 100 * 50 = 5000 mL
  //   organic    = 100 mg * 5000 / (5000 + 100) = 98.039216 mg
  //   aqueous    = 100 - 98.039216            =  1.960784 mg
  //   q          = V_aq / (D V_org + V_aq) = 100 / 5100 = 0.01960784
  //   1 wash     = 1 - q   = 0.9803922  -> short of 99%
  //   2 washes   = 1 - q^2 = 0.9996155  -> clears it
  const sol::Partition worked = sol::partition(100.0, 2.0, 100.0, 50.0);
  CHECK(worked.mgOrganic == doctest::Approx(98.039216).epsilon(1e-7));
  CHECK(worked.mgAqueous == doctest::Approx(1.960784).epsilon(1e-6));
  CHECK(worked.fractionOrganic == doctest::Approx(0.9803922).epsilon(1e-7));
  const double q = worked.mgAqueous / 100.0;
  CHECK(q == doctest::Approx(100.0 / 5100.0));
  CHECK(1.0 - q < 0.99);
  CHECK(1.0 - q * q >= 0.99);

  HeadlessImGui gui;
  ui::AppState st;
  st.solubility.extractionTab = 2;
  st.solubility.soluteValid = true;
  st.solubility.solute.logP = 2.0;
  panelReport(st);
  REQUIRE(st.solubility.funnel.phases.size() == 2);
  // sol::reset sorts densest-first, so index 0 is the organic layer and index 1
  // the lighter aqueous one the panel picks by itself.
  REQUIRE(st.solubility.funnel.phases[0].label == "Dichloromethane");
  st.solubility.funnel.phases[0].volumeMl = 50.0;
  st.solubility.funnel.phases[1].volumeMl = 100.0;
  CHECK(panelReport(st).says(
      "2 washes with fresh Dichloromethane recovers >= 99% of the solute."));

  // logP 8 is past the six decades sol::partition clamps the distribution ratio
  // to. The wash card used to re-derive D as an unclamped 10^logP, so with a
  // large volume ratio it recommended one wash where the split bar's own
  // distribution needs two.
  st.solubility.solute.logP = 8.0;
  st.solubility.funnel.phases[0].volumeMl = 0.4;
  st.solubility.funnel.phases[1].volumeMl = 5000.0;
  const double clamped = sol::partition(1.0, 8.0, 5000.0, 0.4).mgAqueous;
  CHECK(clamped == doctest::Approx(5000.0 / (1.0e6 * 0.4 + 5000.0)));
  CHECK(1.0 - clamped < 0.99);
  CHECK(1.0 - clamped * clamped >= 0.99);
  CHECK(panelReport(st).says("2 washes with fresh Dichloromethane"));
}

TEST_CASE("the organic segment is named for every layer it pools") {
  HeadlessImGui gui;
  ui::AppState st;
  st.solubility.extractionTab = 2;
  st.solubility.soluteValid = true;
  st.solubility.solute.logP = 1.0;
  panelReport(st);

  REQUIRE(st.solubility.funnel.phases.size() == 2);
  CHECK(panelReport(st).says("with fresh Dichloromethane"));

  // A third layer is pooled into the same organic volume, so naming that pool
  // after the first of them would credit one layer with the whole extraction.
  sol::Phase third = st.solubility.funnel.phases[0];
  third.label = "Toluene";
  third.density = 0.87;
  st.solubility.funnel.phases.push_back(third);
  const FrameReport pooled = panelReport(st);
  CHECK(pooled.says("with fresh Organic (2 phases)"));
  CHECK_FALSE(pooled.says("with fresh Dichloromethane"));
}

TEST_CASE("an unusable partition explains itself instead of printing nan") {
  auto expectNoGarbage = [](ui::AppState& state, const char* what) {
    for (int tab = 1; tab <= 2; ++tab) {
      state.solubility.extractionTab = tab;
      const FrameReport report = panelReport(state);
      INFO(what << " on tab " << tab << " ->\n" << report.text);
      CHECK_FALSE(report.says("nan"));
      CHECK_FALSE(report.says("inf"));
      CHECK(report.nonFinite == 0);
      // The card still draws its own explanation; a blank card is the failure
      // this is guarding against just as much as a printed NaN is.
      CHECK(report.vertices > 0);
      CHECK(report.overflow == 0.0f);
    }
  };

  HeadlessImGui gui;
  {
    ui::AppState noSolute;
    expectNoGarbage(noSolute, "no solute");
  }
  {
    // A logP the property model never resolved used to reach the wash chart as
    // q = nan and print "(q = nan)" across the recommendation.
    ui::AppState brokenLogP;
    brokenLogP.solubility.soluteValid = true;
    brokenLogP.solubility.solute.logP =
        std::numeric_limits<double>::quiet_NaN();
    expectNoGarbage(brokenLogP, "unresolved logP");
  }
  {
    ui::AppState infiniteLogP;
    infiniteLogP.solubility.soluteValid = true;
    infiniteLogP.solubility.solute.logP =
        std::numeric_limits<double>::infinity();
    expectNoGarbage(infiniteLogP, "infinite logP");
  }
  {
    ui::AppState emptyOrganic;
    emptyOrganic.solubility.soluteValid = true;
    emptyOrganic.solubility.solute.logP = 1.5;
    panelFrame(emptyOrganic);
    REQUIRE(emptyOrganic.solubility.funnel.phases.size() == 2);
    emptyOrganic.solubility.funnel.phases[0].volumeMl = 0.0;
    expectNoGarbage(emptyOrganic, "an emptied organic layer");
  }
  {
    ui::AppState singlePhase;
    singlePhase.solubility.soluteValid = true;
    singlePhase.solubility.solute.logP = 1.5;
    panelFrame(singlePhase);
    singlePhase.solubility.funnel.phases.resize(1);
    expectNoGarbage(singlePhase, "a single charged phase");
  }
}

TEST_CASE("a shake that cannot run names the reason, and does not cry fault") {
  HeadlessImGui gui;
  // Headless: the particle solver is never built, which is exactly the state
  // the panel used to present as a failure.
  ui::AppState waiting;
  const FrameReport deferred = panelReport(waiting);
  CHECK(deferred.says("Shake is unavailable"));
  CHECK(deferred.says(
      "The particle solver starts once the Extraction workspace is on top."));
  // That line replaces one about the physics RATE, which is a property of a
  // solver that exists.
  CHECK_FALSE(
      deferred.says("Physics rate appears after the first completed step."));

  ui::AppState building;
  building.solubility.fluidBuildPending = true;
  const FrameReport pending = panelReport(building);
  CHECK(pending.says("Shake is unavailable: The particle solver is still being built."));

  // The Live physics card draws the same reason. Its headline and body are
  // custom-drawn rather than ImGui text, so the proof that it reads the reason
  // is that its geometry changes with it.
  ui::AppState idleCard;
  idleCard.solubility.extractionTab = 2;
  const int idleVertices = panelReport(idleCard).vertices;
  ui::AppState pendingCard;
  pendingCard.solubility.extractionTab = 2;
  pendingCard.solubility.fluidBuildPending = true;
  const int pendingVertices = panelReport(pendingCard).vertices;
  CHECK(idleVertices > 0);
  CHECK(pendingVertices != idleVertices);
}

TEST_CASE("the solver-quality picker shows the throughput it measured") {
  HeadlessImGui gui;
  ui::AppState st;
  st.solubility.fluidPresetRealTimeFactor[1] = 0.37;
  st.solubility.fluidPresetRealTimeFactorValid[1] = true;
  st.solubility.fluidPresetMeasuredParticles[1] = 4242;
  panelReport(st);

  const ImGuiTable* transport = findTable(4);
  REQUIRE(transport != nullptr);
  const FrameReport open = clickAt(st, qualityComboCentre(*transport));
  INFO(open.text);
  REQUIRE(open.says("compression limit"));  // the picker really opened
  CHECK(open.says("measured 0.37x at ~4242"));

  // Only the preset that was actually measured carries a measurement; the other
  // two stay honest estimates.
  const size_t first = open.text.find("measured ");
  REQUIRE(first != std::string::npos);
  CHECK(open.text.find("measured ", first + 1) == std::string::npos);
}

TEST_CASE("the extraction panel fits narrow and short surfaces") {
  struct Surface {
    const char* name;
    float width;
    float height;
  };
  // Narrow enough that the stage stacks above the console instead of sitting
  // beside it, and short enough that every card competes for the same rows.
  const Surface surfaces[] = {
      {"900x700", 900.0f, 700.0f},
      {"820x1180", 820.0f, 1180.0f},
      {"1120x620", 1120.0f, 620.0f},
      {"1600x1200", 1600.0f, 1200.0f},
      {"2560x1080", 2560.0f, 1080.0f},
  };

  for (const Surface& surface : surfaces) {
    for (const float scale : {1.0f, 1.5f, 1.75f}) {
      for (int tab = 0; tab <= 2; ++tab) {
        HeadlessImGui gui(surface.width, surface.height, scale);
        ui::AppState st;
        st.solubility.extractionTab = tab;
        st.solubility.soluteValid = true;
        st.solubility.solute.logP = 1.5;
        const FrameReport report = panelReport(st);
        INFO(surface.name << " scale " << scale << " tab " << tab << ":\n"
                          << report.text);
        CHECK(report.overflow == 0.0f);
        CHECK(report.vertices > 0);
        CHECK(report.nonFinite == 0);
        // ImGui's own layout complaints are how clipped and out-of-bounds
        // content announces itself; none of them may survive a resize.
        CHECK_FALSE(report.says("MESSAGE FROM DEAR IMGUI"));
      }
    }
  }
}
// Headless render test for the Solubility Suite panel. Drives the real
// drawSolubilitySuite() through a null ImGui backend, so the ratio plot and
// the solvent screen are proven to emit geometry without a display. The
// separatory funnel lives in the Extraction Lab now -- see
// test_extraction_panel.cpp for its coverage.
//
// The composition charts are also INPUTS, and the cases at the bottom drive
// them with synthetic clicks: charts3d::ternary() returns -2 with the clicked
// barycentric composition, charts::linePlot() returns the hovered sample, and
// both are meant to travel all the way to sb.prediction. The whole chain
// (pointer -> sb.ratios -> prediction -> on-chart marker) is asserted link by
// link, because a chart that moves its marker while the headline number stays
// put is worse than a chart that ignores the click.
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

#include "imgui.h"
#include "imgui_internal.h"

#include "chem/bridge.hpp"
#include "sol/solubility.hpp"
#include "sol/solvent.hpp"
#include "ui/app_state.hpp"
#include "ui/charts3d.hpp"
#include "ui/solubility_state.hpp"
#include "ui/ui.hpp"

using namespace chemcad;

namespace {

constexpr float kDisplayW = 1600.0f;
constexpr float kDisplayH = 1200.0f;

struct HeadlessImGui {
  HeadlessImGui() {
    IMGUI_CHECKVERSION();
    ctx = ImGui::CreateContext();
    ImGui::SetCurrentContext(ctx);
    ImGuiIO& io = ImGui::GetIO();
    io.DisplaySize = ImVec2(kDisplayW, kDisplayH);
    io.DeltaTime = 1.0f / 60.0f;
    io.IniFilename = nullptr;
    io.LogFilename = nullptr;
    unsigned char* pixels = nullptr;
    int w = 0, h = 0;
    io.Fonts->GetTexDataAsRGBA32(&pixels, &w, &h);
    io.Fonts->SetTexID(static_cast<ImTextureID>(1));
  }
  ~HeadlessImGui() { ImGui::DestroyContext(ctx); }
  ImGuiContext* ctx = nullptr;
};

int panelFrame(ui::AppState& st) {
  ImGui::NewFrame();
  ImGui::SetNextWindowPos(ImVec2(0, 0));
  ImGui::SetNextWindowSize(ImVec2(kDisplayW, kDisplayH));
  ImGui::Begin("Solubility Suite", nullptr,
               ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                   ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings);
  ui::drawSolubilitySuite(st);
  ImGui::End();
  ImGui::Render();
  return ImGui::GetDrawData()->TotalVtxCount;
}

// Two frames: the panel caches its solute, sweep and screening on the first
// one, so the second frame is the one that actually has curves and tables to
// draw.
int settledVertexCount(ui::AppState& st) {
  panelFrame(st);
  return panelFrame(st);
}

void loadSolute(ui::AppState& st, const char* smiles) {
  st.doc.molecules.push_back(chem::fromSmiles(smiles));
  st.touch();
}

void loadBinaryBlend(ui::AppState& st) {
  loadSolute(st, "OC(=O)c1ccccc1");
  st.solubility.solventIds[0] = "water";
  st.solubility.solventIds[1] = "ethanol";
  st.solubility.ratios[0] = 1.0f;
  st.solubility.ratios[1] = 1.0f;
}

void loadTernaryBlend(ui::AppState& st) {
  loadBinaryBlend(st);
  st.solubility.solventCount = 3;
  st.solubility.solventIds[2] = "toluene";
  st.solubility.ratios[2] = 1.0f;
  // The coarsest sample grid: it changes how densely the simplex is sampled,
  // never where the triangle is drawn, and a click probe renders dozens of
  // frames.
  st.solubility.sweepSteps = 2;
}

// ----------------------------------------------------------- input injection
// ImGui can only hover an item that already existed on an earlier frame, so one
// synthetic click costs three frames: settle the pointer, deliver the press,
// deliver the release so the next probe starts from a clean button state.
void moveMouse(ImVec2 at) {
  ImGuiIO& io = ImGui::GetIO();
  io.MousePos = at;
  io.AddMousePosEvent(at.x, at.y);
}

void clickAt(ui::AppState& st, ImVec2 at) {
  ImGuiIO& io = ImGui::GetIO();
  moveMouse(at);
  panelFrame(st);
  // Deliberately no second position event: ImGui's event trickling stops
  // draining the queue when a button change follows a *changed* pointer
  // position, which would push this press onto a later frame than the release.
  io.AddMouseButtonEvent(ImGuiMouseButton_Left, true);
  panelFrame(st);
  io.AddMouseButtonEvent(ImGuiMouseButton_Left, false);
  panelFrame(st);
}

// ImGui floors the pointer position when it drains the event queue, so every
// probe is placed on whole pixels and the expected composition is derived from
// the same whole pixels.
ImVec2 floorPos(ImVec2 p) { return ImVec2(std::floor(p.x), std::floor(p.y)); }

// A card is a child window, so ImGui itself knows where its content region is.
// That is the one piece of layout the probes take on trust; the triangle inside
// it is then measured from the chart's own answers.
ImGuiWindow* findChildWindow(const char* idFragment) {
  ImGuiContext& g = *ImGui::GetCurrentContext();
  for (ImGuiWindow* window : g.Windows) {
    if (window->Name && std::strstr(window->Name, idFragment) != nullptr) return window;
  }
  return nullptr;
}

// A composition no click can produce: three parts that do not sum to one, so
// "the chart wrote something" and "the chart wrote nothing" are distinguishable
// without knowing what the answer should be.
constexpr std::array<float, 3> kUntouched{2.0f, 3.0f, 4.0f};

struct Probe {
  bool accepted = false;  // the chart turned this pointer position into a blend
  double a = 0.0;
  double b = 0.0;
  double c = 0.0;
};

Probe probeClick(ui::AppState& st, ImVec2 at) {
  st.solubility.ratios = kUntouched;
  clickAt(st, at);
  const std::array<float, 3>& ratios = st.solubility.ratios;
  return {ratios != kUntouched, ratios[0], ratios[1], ratios[2]};
}

// The composition triangle as charts3d.cpp constructs it: apex centred on the
// chart's top inset, base one equilateral height below it.
struct Triangle {
  float centreX = 0.0f;
  float topY = 0.0f;
  float height = 0.0f;
  float side = 0.0f;
};

// Barycentric weights the chart owes a pointer at `at`, inverting that same
// construction: apex weight falls with depth, the base weights split by how far
// the pointer sits from the centre line.
std::array<double, 3> weightsAt(const Triangle& tri, ImVec2 at) {
  const double depth = (at.y - tri.topY) / tri.height;
  const double lateral = (at.x - tri.centreX) / tri.side;
  return {1.0 - depth, depth * 0.5 - lateral, depth * 0.5 + lateral};
}

ImVec2 pointFor(const Triangle& tri, double a, double b, double c) {
  return ImVec2(tri.centreX + static_cast<float>((c - b) * 0.5) * tri.side,
                tri.topY + static_cast<float>(1.0 - a) * tri.height);
}

// Measures the live triangle from the chart's own replies rather than replaying
// the panel's layout arithmetic. Three probes are enough: the apex weight is
// linear in y and the base weights split linearly in x, so a vertical pair
// fixes the apex and the height and a horizontal pair fixes the axis and the
// side. Nothing here assumes where the triangle sits -- the equilateral
// relation between height and side is then checked, not relied on.
Triangle measureTriangle(ui::AppState& st) {
  ImGuiWindow* card = findChildWindow("blend_response_card");
  REQUIRE(card != nullptr);
  const ImRect content = card->ContentRegionRect;
  const float probeX = std::floor((content.Min.x + content.Max.x) * 0.5f);

  // Walk up from the bottom of the card's content: the first pointer position
  // the chart turns into a composition is just inside the triangle's base.
  float baseY = 0.0f;
  Probe base;
  for (int step = 1; step <= 140 && !base.accepted; ++step) {
    baseY = std::floor(content.Max.y) - static_cast<float>(step) * 3.0f;
    REQUIRE(baseY > content.Min.y);
    base = probeClick(st, ImVec2(probeX, baseY));
  }
  REQUIRE(base.accepted);

  Probe upper;
  float upperY = baseY;
  for (const float rise : {160.0f, 120.0f, 80.0f, 40.0f, 20.0f, 10.0f}) {
    upperY = baseY - rise;
    if (upperY <= content.Min.y) continue;
    upper = probeClick(st, ImVec2(probeX, upperY));
    if (upper.accepted) break;
  }
  REQUIRE(upper.accepted);

  Probe lateral;
  float offset = 0.0f;
  for (const float candidate : {96.0f, 64.0f, 32.0f, 16.0f, 8.0f}) {
    offset = candidate;
    lateral = probeClick(st, ImVec2(probeX + offset, baseY));
    if (lateral.accepted) break;
  }
  REQUIRE(lateral.accepted);

  Triangle tri;
  // Depth into the triangle is 1 - apexWeight, so the apex weight FALLS as y
  // grows: the rise and the weight difference are both negative here.
  tri.height = static_cast<float>((upperY - baseY) / (base.a - upper.a));
  tri.topY = baseY - static_cast<float>(1.0 - base.a) * tri.height;
  const double lateralBase = (base.c - base.b) * 0.5;
  const double lateralOff = (lateral.c - lateral.b) * 0.5;
  tri.side = static_cast<float>(offset / (lateralOff - lateralBase));
  tri.centreX = probeX - static_cast<float>(lateralBase) * tri.side;

  // Height and side were measured on independent axes; an equilateral triangle
  // ties them together, and the axis must be the centre of the chart, which
  // spans the card's full content width.
  CHECK(tri.side == doctest::Approx(tri.height * 2.0 / std::sqrt(3.0)).epsilon(5e-3));
  CHECK(std::fabs(tri.centreX - (content.Min.x + content.Max.x) * 0.5f) < 1.5f);
  return tri;
}

// The prediction the panel owes the current sliders, recomputed the way
// resolveBlend() + drawSolubilitySuite() do: unset slots and zero parts drop
// out, pH is self-buffered and there is no background electrolyte by default.
double predictionForRatios(const ui::SolubilityState& sb) {
  std::vector<sol::Component> components;
  for (int i = 0; i < sb.solventCount; ++i) {
    const sol::Solvent* solvent = sol::findSolvent(sb.solventIds[static_cast<size_t>(i)]);
    const float parts = sb.ratios[static_cast<size_t>(i)];
    if (solvent && parts > 0.0f) components.push_back({solvent, static_cast<double>(parts)});
  }
  return sol::predict(sb.solute, components, static_cast<double>(sb.temperatureC),
                      nullptr, 0.0, sol::kAutoPH)
      .gramsPerMillilitre;
}

double purePrediction(const ui::SolubilityState& sb, const char* solventId) {
  const sol::Solvent* solvent = sol::findSolvent(solventId);
  REQUIRE(solvent != nullptr);
  const std::vector<sol::Component> components{{solvent, 1.0}};
  return sol::predict(sb.solute, components, static_cast<double>(sb.temperatureC),
                      nullptr, 0.0, sol::kAutoPH)
      .gramsPerMillilitre;
}

// ------------------------------------------------ the triangle on its own
// A degenerate response -- every sampled composition predicting the same
// solubility, which is what an insoluble solute gives -- is reproduced by
// driving the chart directly rather than by hunting for a solute the model
// happens to flatten. Origin and size are the test's, so the geometry is
// arithmetic: at this aspect the triangle is height-limited, which puts the
// dead centre of the rect on its axis at half depth, i.e. one half apex and one
// quarter of each base vertex.
const ImVec2 kTernaryOrigin(40.0f, 40.0f);
const ImVec2 kTernarySize(600.0f, 300.0f);

struct TernaryPick {
  int result = 0;
  double a = 0.0;
  double b = 0.0;
  double c = 0.0;
  int vertices = 0;
};

TernaryPick ternaryFrame(const std::vector<double>& points,
                         const std::vector<double>& values, bool press) {
  ImGuiIO& io = ImGui::GetIO();
  const ImVec2 centre(kTernaryOrigin.x + kTernarySize.x * 0.5f,
                      kTernaryOrigin.y + kTernarySize.y * 0.5f);
  io.MousePos = centre;
  io.AddMousePosEvent(centre.x, centre.y);
  if (press) io.AddMouseButtonEvent(ImGuiMouseButton_Left, true);
  ImGui::NewFrame();
  ImGui::SetNextWindowPos(ImVec2(0, 0));
  ImGui::SetNextWindowSize(ImVec2(kDisplayW, kDisplayH));
  ImGui::Begin("Ternary", nullptr,
               ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                   ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings);
  ImGui::SetCursorScreenPos(kTernaryOrigin);
  const ui::charts3d::TernaryStyle style;
  TernaryPick pick;
  pick.result = ui::charts3d::ternary("##bare_ternary", points.data(), values.data(),
                                      static_cast<int>(values.size()), kTernarySize,
                                      style, &pick.a, &pick.b, &pick.c);
  ImGui::End();
  ImGui::Render();
  pick.vertices = ImGui::GetDrawData()->TotalVtxCount;
  if (press) io.AddMouseButtonEvent(ImGuiMouseButton_Left, false);
  return pick;
}

std::vector<double> simplexGrid(int divisions) {
  std::vector<double> points;
  for (int i = 0; i <= divisions; ++i) {
    for (int j = 0; i + j <= divisions; ++j) {
      points.push_back(static_cast<double>(divisions - i - j) / divisions);
      points.push_back(static_cast<double>(i) / divisions);
      points.push_back(static_cast<double>(j) / divisions);
    }
  }
  return points;
}

}  // namespace

TEST_CASE("the panel opens on a binary blend so the ratio plot is reachable") {
  ui::AppState st;
  CHECK(st.solubility.solventCount == 2);
  CHECK(st.solubility.extractionImport.pending == false);
}

TEST_CASE("the ratio plot draws geometry once a solute and two solvents exist") {
  HeadlessImGui gui;

  // Baseline: a solute but no solvents, so there is nothing to sweep.
  ui::AppState bare;
  loadSolute(bare, "OC(=O)c1ccccc1");
  const int bareVertices = settledVertexCount(bare);

  ui::AppState plotted;
  loadBinaryBlend(plotted);
  const int plottedVertices = settledVertexCount(plotted);

  // The curve, its shaded area, the axis frame, the tick labels and the peak
  // marker are all real draw commands, so the plotted panel cannot cost the
  // same as the empty one.
  CHECK(plottedVertices > bareVertices);
}

TEST_CASE("a ternary blend draws more than a binary one") {
  HeadlessImGui gui;

  ui::AppState binary;
  loadBinaryBlend(binary);
  const int binaryVertices = settledVertexCount(binary);

  ui::AppState ternary;
  loadBinaryBlend(ternary);
  ternary.solubility.solventCount = 3;
  ternary.solubility.solventIds[2] = "toluene";
  ternary.solubility.ratios[2] = 1.0f;
  const int ternaryVertices = settledVertexCount(ternary);

  // The shaded simplex is one filled triangle per cell, so it dwarfs a polyline.
  CHECK(ternaryVertices > binaryVertices);
}

TEST_CASE("the solvent screen fills and stays sorted for a valid solute") {
  HeadlessImGui gui;
  ui::AppState st;
  loadBinaryBlend(st);
  settledVertexCount(st);

  const std::vector<sol::ScreenRow>& rows = st.solubility.screening;
  REQUIRE(rows.size() == sol::solvents().size());
  for (size_t i = 1; i < rows.size(); ++i) {
    CHECK(rows[i - 1].prediction.gramsPerMillilitre >=
          rows[i].prediction.gramsPerMillilitre);
  }

  // The screen invalidates when the temperature moves.
  st.solubility.temperatureC = 60.0f;
  settledVertexCount(st);
  CHECK(st.solubility.screeningSignature !=
        std::string(""));
}

TEST_CASE("the panel survives a missing solvent id and an unparseable solute") {
  HeadlessImGui gui;
  ui::AppState st;
  st.solubility.solventIds[0] = "this-solvent-does-not-exist";
  st.solubility.solventIds[1] = "water";
  CHECK_NOTHROW(settledVertexCount(st));

  st.solubility.useSketch = false;
  st.solubility.manualSmiles = "not-a-valid-smiles";
  CHECK_NOTHROW(settledVertexCount(st));
}

TEST_CASE("clicking the composition triangle sets the blend it was clicked at") {
  HeadlessImGui gui;
  ui::AppState st;
  loadTernaryBlend(st);
  settledVertexCount(st);

  ui::SolubilityState& sb = st.solubility;
  REQUIRE(sb.soluteValid);
  REQUIRE(sb.solventCount == 3);
  REQUIRE_FALSE(sb.ternarySurface.nodes.empty());
  // The triangle, not the 3D landscape, is the view a three-solvent blend opens
  // on -- and it is the only one of the two that takes a composition.
  REQUIRE(sb.ternaryView == 0);

  const Triangle tri = measureTriangle(st);

  // Apex is Solvent A (water), bottom-left Solvent B (ethanol), bottom-right
  // Solvent C (toluene): one probe near a vertex, one at the centroid, one
  // asymmetric interior point.
  const double kCases[][3] = {{0.90, 0.05, 0.05},
                              {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
                              {0.30, 0.20, 0.50},
                              {0.06, 0.05, 0.89}};
  for (const double(&want)[3] : kCases) {
    const ImVec2 at = floorPos(pointFor(tri, want[0], want[1], want[2]));
    const std::array<double, 3> expected = weightsAt(tri, at);
    const Probe got = probeClick(st, at);
    REQUIRE(got.accepted);
    CHECK(got.a == doctest::Approx(expected[0]).epsilon(1e-3));
    CHECK(got.b == doctest::Approx(expected[1]).epsilon(1e-3));
    CHECK(got.c == doctest::Approx(expected[2]).epsilon(1e-3));

    // Unit consistency. The chart writes normalised fractions into sliders
    // labelled "Volume parts" on a 0..10 range. That is coherent rather than
    // confusing: every consumer -- resolveBlend(), the marker, Apply peak, Send
    // to Extraction -- renormalises, so 0.33/0.33/0.33 parts is the same blend
    // as 1/1/1 parts and still reads legibly at the slider's "%.2f parts". What
    // must not happen is a click leaving a slider off its own scale.
    for (const float part : sb.ratios) {
      CHECK(part >= 0.0f);
      CHECK(part <= 10.0f);
    }
    // The marker renormalises the three parts, so a click writing fractions
    // that sum to one has to put the marker back under the pointer. Sub-pixel,
    // because the round trip is screen -> composition -> screen.
    const double total = got.a + got.b + got.c;
    CHECK(total == doctest::Approx(1.0).epsilon(1e-5));
    const ImVec2 marker = pointFor(tri, got.a / total, got.b / total, got.c / total);
    CHECK(std::fabs(marker.x - at.x) < 1.0f);
    CHECK(std::fabs(marker.y - at.y) < 1.0f);
  }

  // Inside the widget rect but outside the triangle: the bottom-left corner of
  // the chart's bounding box sits left of the b-a edge and below the base.
  ImGuiWindow* card = findChildWindow("blend_response_card");
  REQUIRE(card != nullptr);
  const ImVec2 corner(std::floor(card->ContentRegionRect.Min.x) + 2.0f,
                      std::floor(card->ContentRegionRect.Max.y) - 2.0f);
  const std::array<double, 3> cornerWeights = weightsAt(tri, corner);
  REQUIRE(std::min({cornerWeights[0], cornerWeights[1], cornerWeights[2]}) < 0.0);
  CHECK_FALSE(probeClick(st, corner).accepted);

  // And the landscape view is genuinely a different chart: the same pointer
  // position that just set a composition sets nothing once the 3D surface is up.
  sb.ternaryView = 1;
  const ImVec2 centroid = floorPos(pointFor(tri, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0));
  CHECK_FALSE(probeClick(st, centroid).accepted);
  sb.ternaryView = 0;
}

TEST_CASE("a composition clicked on the triangle drives the reported prediction") {
  HeadlessImGui gui;
  ui::AppState st;
  loadTernaryBlend(st);
  settledVertexCount(st);
  ui::SolubilityState& sb = st.solubility;
  const Triangle tri = measureTriangle(st);

  const ImVec2 nearWater = floorPos(pointFor(tri, 0.94, 0.03, 0.03));
  REQUIRE(probeClick(st, nearWater).accepted);
  const double atWater = sb.prediction.gramsPerMillilitre;
  // The link that matters: the headline number is the prediction OF the clicked
  // composition, not of whatever the sliders held before the click.
  CHECK(atWater == doctest::Approx(predictionForRatios(sb)).epsilon(1e-9));

  const ImVec2 nearEthanol = floorPos(pointFor(tri, 0.03, 0.94, 0.03));
  REQUIRE(probeClick(st, nearEthanol).accepted);
  const double atEthanol = sb.prediction.gramsPerMillilitre;
  CHECK(atEthanol == doctest::Approx(predictionForRatios(sb)).epsilon(1e-9));

  // Two materially different corners of the simplex must not report the same
  // solubility, and which way the number moves must agree with what the same
  // model says about the two pure solvents.
  const double scale = std::max(std::fabs(atWater), std::fabs(atEthanol));
  CHECK(std::fabs(atEthanol - atWater) > 0.01 * scale);
  CHECK((atEthanol > atWater) ==
        (purePrediction(sb, "ethanol") > purePrediction(sb, "water")));

  // The latency is one frame, and one frame is correct: drawSolubilitySuite()
  // predicts from sb.ratios before it draws the chart that reads the click, so
  // the headline catches up on the frame after the press, never on the press
  // itself. Pinned down here because two frames of lag, or none at all, would
  // both mean the chain had been rewired.
  ImGuiIO& io = ImGui::GetIO();
  sb.ratios = kUntouched;
  moveMouse(nearWater);
  panelFrame(st);
  const double stale = sb.prediction.gramsPerMillilitre;
  io.AddMouseButtonEvent(ImGuiMouseButton_Left, true);
  panelFrame(st);
  CHECK(sb.ratios != kUntouched);
  CHECK(sb.prediction.gramsPerMillilitre == doctest::Approx(stale).epsilon(1e-12));
  io.AddMouseButtonEvent(ImGuiMouseButton_Left, false);
  panelFrame(st);
  CHECK(sb.prediction.gramsPerMillilitre ==
        doctest::Approx(predictionForRatios(sb)).epsilon(1e-9));
  CHECK_FALSE(sb.prediction.gramsPerMillilitre == doctest::Approx(stale));
}

TEST_CASE("clicking the binary response curve sets the working point and the prediction follows") {
  HeadlessImGui gui;
  ui::AppState st;
  loadBinaryBlend(st);
  settledVertexCount(st);
  ui::SolubilityState& sb = st.solubility;
  REQUIRE(sb.sweep.size() > 2);

  ImGuiWindow* card = findChildWindow("blend_response_card");
  REQUIRE(card != nullptr);
  const ImRect content = card->ContentRegionRect;
  // The curve fills the card's remaining content, and its item rect reaches the
  // bottom of that region, so a click just above the bottom edge is inside the
  // plot and clear of the sample-count and units controls at the top.
  const float y = std::floor(content.Max.y) - 4.0f;
  const float leftX = std::floor(content.Min.x + content.GetWidth() * 0.32f);
  const float rightX = std::floor(content.Min.x + content.GetWidth() * 0.84f);

  const Probe left = probeClick(st, ImVec2(leftX, y));
  REQUIRE(left.accepted);
  // The plot snaps to a swept sample, so the working point is always a
  // composition the panel actually has a prediction for.
  CHECK(left.a + left.b == doctest::Approx(1.0).epsilon(1e-5));
  bool snapped = false;
  for (const sol::SweepPoint& point : sb.sweep) {
    if (std::fabs(point.fractions[0] - left.a) < 1e-5) snapped = true;
  }
  CHECK(snapped);
  const double atLeft = sb.prediction.gramsPerMillilitre;
  CHECK(atLeft == doctest::Approx(predictionForRatios(sb)).epsilon(1e-9));

  const Probe right = probeClick(st, ImVec2(rightX, y));
  REQUIRE(right.accepted);
  // The x axis is Solvent A volume percent, so clicking further right must ask
  // for more of Solvent A.
  CHECK(right.a > left.a);
  CHECK(right.a + right.b == doctest::Approx(1.0).epsilon(1e-5));
  const double atRight = sb.prediction.gramsPerMillilitre;
  CHECK(atRight == doctest::Approx(predictionForRatios(sb)).epsilon(1e-9));

  const double scale = std::max(std::fabs(atLeft), std::fabs(atRight));
  CHECK(std::fabs(atRight - atLeft) > 0.01 * scale);

  // Sweep the pointer across the plot: the working point must track it
  // monotonically, always land on a swept sample, and reach both pure ends.
  // A mapping that collapsed to a couple of buckets would satisfy the two
  // probes above and fail here.
  double previous = -1.0;
  int distinct = 0;
  double lowest = 1.0;
  double highest = 0.0;
  // ImRect::Contains is exclusive on the far edge, so the sweep stays two
  // pixels inside the chart on both sides.
  const float sweepLeft = std::floor(content.Min.x) + 2.0f;
  const float sweepRight = std::floor(content.Max.x) - 2.0f;
  for (int step = 0; step <= 24; ++step) {
    const float x = std::floor(sweepLeft + (sweepRight - sweepLeft) *
                                               static_cast<float>(step) / 24.0f);
    const Probe swept = probeClick(st, ImVec2(x, y));
    CAPTURE(step);
    CAPTURE(swept.a);
    REQUIRE(swept.accepted);
    bool onSample = false;
    for (const sol::SweepPoint& point : sb.sweep) {
      if (std::fabs(point.fractions[0] - swept.a) < 1e-5) onSample = true;
    }
    CHECK(onSample);
    CHECK(swept.a >= previous);
    if (swept.a > previous) ++distinct;
    previous = swept.a;
    lowest = std::min(lowest, swept.a);
    highest = std::max(highest, swept.a);
  }
  CHECK(lowest <= 0.05);
  CHECK(highest >= 0.95);
  CHECK(distinct >= 8);

  // A click outside the curve's rect leaves the working point alone.
  CHECK_FALSE(probeClick(st, ImVec2(std::floor(content.Min.x) + 2.0f,
                                    std::floor(content.Min.y) + 2.0f))
                  .accepted);
}

TEST_CASE("the composition triangle still takes a click when the response is flat") {
  HeadlessImGui gui;
  const std::vector<double> points = simplexGrid(6);
  const size_t count = points.size() / 3;
  const std::vector<double> flat(count, 4.2);
  std::vector<double> varying(count);
  for (size_t i = 0; i < count; ++i) varying[i] = points[i * 3 + 1];

  // First frame establishes the item, second delivers the press.
  ternaryFrame(points, varying, false);
  const TernaryPick varyingPick = ternaryFrame(points, varying, true);
  REQUIRE(varyingPick.result == -2);
  CHECK(varyingPick.a == doctest::Approx(0.5).epsilon(1e-3));
  CHECK(varyingPick.b == doctest::Approx(0.25).epsilon(1e-3));
  CHECK(varyingPick.c == doctest::Approx(0.25).epsilon(1e-3));

  // A response that reads the same everywhere -- an insoluble solute -- has
  // nothing to colour, but it must not cost the user the composition click the
  // blend card advertises.
  ternaryFrame(points, flat, false);
  const TernaryPick flatPick = ternaryFrame(points, flat, true);
  CHECK(flatPick.result == -2);
  CHECK(flatPick.a == doctest::Approx(varyingPick.a).epsilon(1e-6));
  CHECK(flatPick.b == doctest::Approx(varyingPick.b).epsilon(1e-6));
  CHECK(flatPick.c == doctest::Approx(varyingPick.c).epsilon(1e-6));

  // It also still draws itself -- outline, labels and a uniform fill -- while
  // the varying response adds isolines on top, which proves the values are read
  // rather than ignored.
  CHECK(flatPick.vertices > 0);
  CHECK(varyingPick.vertices > flatPick.vertices);
}

// Headless render test for the Extraction Lab panel. Drives the real
// drawExtractionLab() through a null ImGui backend, so the cross-section,
// the suite import hand-off and the solute distribution panel are proven to
// work without a display.
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <algorithm>
#include <array>
#include <cmath>

#include <string>

#include "imgui.h"

#include "chem/bridge.hpp"
#include "sol/funnel.hpp"
#include "sol/solubility.hpp"
#include "sol/solvent.hpp"
#include "ui/app_state.hpp"
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

TEST_CASE("the cross-section draws and responds to shaking") {
  HeadlessImGui gui;
  ui::AppState st;
  const int settled = settledVertexCount(st);

  // Two default phases are seeded on the first frame; five seconds into a
  // firm 6 s shake the column is full of droplets, every one of which is an
  // extra filled circle.
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

  ui::AppState bare;
  const int bareVertices = settledVertexCount(bare);

  ui::AppState withSolute;
  withSolute.doc.molecules.push_back(chem::fromSmiles("OC(=O)c1ccccc1"));
  withSolute.touch();
  withSolute.solubility.solute =
      sol::describeSolute(withSolute.doc.molecules.front());
  withSolute.solubility.soluteValid = true;
  const int soluteVertices = settledVertexCount(withSolute);

  // The stacked partition bar and its labels are real geometry on top of the
  // bare vessel.
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
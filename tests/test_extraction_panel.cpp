// Headless render test for the Extraction Lab panel. Drives the real
// drawExtractionLab() through a null ImGui backend, so the cross-section,
// the suite import hand-off and the solute distribution panel are proven to
// work without a display.
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

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

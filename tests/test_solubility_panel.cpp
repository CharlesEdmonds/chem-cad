// Headless render test for the Solubility Suite panel. Drives the real
// drawSolubilitySuite() through a null ImGui backend, so the ratio plot and the
// funnel cross-section are proven to emit geometry without a display.
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <string>

#include "imgui.h"

#include "chem/bridge.hpp"
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
  ImGui::Begin("Solubility Suite", nullptr,
               ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                   ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings);
  ui::drawSolubilitySuite(st);
  ImGui::End();
  ImGui::Render();
  return ImGui::GetDrawData()->TotalVtxCount;
}

// Two frames: the panel caches its solute and sweep on the first one, so the
// second frame is the one that actually has a curve to draw.
int settledVertexCount(ui::AppState& st) {
  panelFrame(st);
  return panelFrame(st);
}

void loadSolute(ui::AppState& st, const char* smiles) {
  st.doc.molecules.push_back(chem::fromSmiles(smiles));
  st.touch();
}

}  // namespace

TEST_CASE("the panel opens on a binary blend so the ratio plot is reachable") {
  ui::AppState st;
  CHECK(st.solubility.solventCount == 2);
}

TEST_CASE("the ratio plot draws geometry once a solute and two solvents exist") {
  HeadlessImGui gui;

  // Baseline: a solute but no solvents, so there is nothing to sweep.
  ui::AppState bare;
  loadSolute(bare, "OC(=O)c1ccccc1");
  const int bareVertices = settledVertexCount(bare);

  ui::AppState plotted;
  loadSolute(plotted, "OC(=O)c1ccccc1");
  REQUIRE(sol::findSolvent("water") != nullptr);
  REQUIRE(sol::findSolvent("ethanol") != nullptr);
  plotted.solubility.solventIds[0] = "water";
  plotted.solubility.solventIds[1] = "ethanol";
  plotted.solubility.ratios[0] = 1.0f;
  plotted.solubility.ratios[1] = 1.0f;
  const int plottedVertices = settledVertexCount(plotted);

  // The curve, its shaded area, the axis frame, the tick labels and the peak
  // marker are all real draw commands, so the plotted panel cannot cost the
  // same as the empty one.
  CHECK(plottedVertices > bareVertices);
}

TEST_CASE("a ternary blend draws more than a binary one") {
  HeadlessImGui gui;

  ui::AppState binary;
  loadSolute(binary, "OC(=O)c1ccccc1");
  binary.solubility.solventCount = 2;
  binary.solubility.solventIds[0] = "water";
  binary.solubility.solventIds[1] = "ethanol";
  binary.solubility.ratios[0] = 1.0f;
  binary.solubility.ratios[1] = 1.0f;
  const int binaryVertices = settledVertexCount(binary);

  ui::AppState ternary;
  loadSolute(ternary, "OC(=O)c1ccccc1");
  ternary.solubility.solventCount = 3;
  ternary.solubility.solventIds[0] = "water";
  ternary.solubility.solventIds[1] = "ethanol";
  ternary.solubility.solventIds[2] = "toluene";
  ternary.solubility.ratios[0] = 1.0f;
  ternary.solubility.ratios[1] = 1.0f;
  ternary.solubility.ratios[2] = 1.0f;
  const int ternaryVertices = settledVertexCount(ternary);

  // The shaded simplex is one filled triangle per cell, so it dwarfs a polyline.
  CHECK(ternaryVertices > binaryVertices);
}

TEST_CASE("the funnel cross-section draws and responds to shaking") {
  HeadlessImGui gui;
  ui::AppState st;
  loadSolute(st, "OC(=O)c1ccccc1");
  const int settled = settledVertexCount(st);

  // The phase editor seeds two phases on the first frame; shaking fills the
  // column with droplets, every one of which is an extra filled circle.
  REQUIRE(st.solubility.funnel.phases.size() == 2);
  sol::shake(st.solubility.funnel, 1.0);
  REQUIRE(sol::emulsifiedFraction(st.solubility.funnel) > 0.2);
  const int shaken = settledVertexCount(st);

  CHECK(shaken > settled);
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

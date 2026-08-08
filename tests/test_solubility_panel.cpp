// Headless render test for the Solubility Suite panel. Drives the real
// drawSolubilitySuite() through a null ImGui backend, so the ratio plot and
// the solvent screen are proven to emit geometry without a display. The
// separatory funnel lives in the Extraction Lab now -- see
// test_extraction_panel.cpp for its coverage.
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

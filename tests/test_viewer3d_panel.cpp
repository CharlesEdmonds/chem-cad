// Headless render test for the 3D molecule viewer: drives the real
// drawViewer3D() through a null ImGui backend, proving the turntable
// renderer emits geometry and follows the sketch.
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <string>

#include "imgui.h"

#include "chem/bridge.hpp"
#include "ui/app_state.hpp"
#include "ui/ui.hpp"

using namespace chemcad;

namespace {

constexpr float kDisplayW = 1400.0f;
constexpr float kDisplayH = 1000.0f;

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
  ImGui::Begin("3D View", nullptr,
               ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                   ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings);
  ui::drawViewer3D(st);
  ImGui::End();
  ImGui::Render();
  return ImGui::GetDrawData()->TotalVtxCount;
}

}  // namespace

TEST_CASE("the viewer renders an empty state without a sketch") {
  HeadlessImGui gui;
  ui::AppState st;
  CHECK_NOTHROW(panelFrame(st));
  CHECK(!st.viewer3d.hasModel);
}

TEST_CASE("the viewer embeds the sketch and draws real geometry") {
  HeadlessImGui gui;

  ui::AppState bare;
  const int bareVertices = panelFrame(bare);

  ui::AppState st;
  st.doc.molecules.push_back(chem::fromSmiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C"));
  st.touch();
  const int vertices = panelFrame(st);

  REQUIRE(st.viewer3d.hasModel);
  CHECK(st.viewer3d.formula == "C8H10N4O2");
  CHECK(vertices > bareVertices);
}

TEST_CASE("the model re-embeds when the sketch changes") {
  HeadlessImGui gui;
  ui::AppState st;
  st.doc.molecules.push_back(chem::fromSmiles("CCO"));
  const uint64_t firstRevision = st.docRevision;
  st.touch();
  panelFrame(st);
  REQUIRE(st.viewer3d.hasModel);
  const size_t ethanolAtoms = st.viewer3d.model.atoms.size();

  // A new molecule at a new revision must replace the cached conformer.
  st.doc.molecules.clear();
  st.doc.molecules.push_back(chem::fromSmiles("c1ccccc1"));
  st.touch();
  REQUIRE(st.docRevision != firstRevision);
  panelFrame(st);
  CHECK(st.viewer3d.model.atoms.size() != ethanolAtoms);
  CHECK(st.viewer3d.formula == "C6H6");
}

TEST_CASE("every render style draws without throwing") {
  HeadlessImGui gui;
  ui::AppState st;
  st.doc.molecules.push_back(chem::fromSmiles("c1ccccc1"));
  st.touch();
  for (int style = 0; style < 3; ++style) {
    st.viewer3d.style = style;
    CHECK_NOTHROW(panelFrame(st));
  }
}

// The "one page per tab" contract, enforced mechanically.
//
// Every workspace is supposed to measure the surface it was given and divide
// it, so that no panel ever grows a scrollbar it did not ask for. That property
// is invisible in a screenshot of one window size and easy to break by adding a
// single fixed-height card, so it is checked here across the display shapes the
// application actually meets: a small laptop panel, 1080p, an ultrawide, a
// nearly-square 4:3 dock and a 4K surface, each at several UI scales.
//
// A panel that legitimately scrolls INSIDE a budgeted child (the Toolbox
// results list, the Planner's route list) still passes: the assertion is about
// the panel's own window, which is what the user sees as "the page".
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <string>
#include <vector>

#include "imgui.h"

#include "ui/app_state.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"

using namespace chemcad;

namespace {

struct HeadlessImGui {
  explicit HeadlessImGui(float width, float height, float fontScale) {
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

using PanelFn = void (*)(ui::AppState&);

struct Panel {
  const char* name;
  PanelFn draw;
};

struct Surface {
  const char* name;
  float width;
  float height;
};

// Draws the panel filling the display and reports how far it could be scrolled.
// A positive value means the content overflowed the page.
float overflowOf(ui::AppState& st, PanelFn draw, const char* title) {
  float worst = 0.0f;
  // Two frames: docking, tab state and every animated value settle on the
  // second, and a layout that only overflows once settled would otherwise pass.
  for (int frame = 0; frame < 2; ++frame) {
    ImGui::NewFrame();
    const ImGuiIO& io = ImGui::GetIO();
    ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f));
    ImGui::SetNextWindowSize(io.DisplaySize);
    if (ImGui::Begin(title, nullptr,
                     ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                         ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings)) {
      draw(st);
      worst = ImGui::GetScrollMaxY();
    }
    ImGui::End();
    ImGui::Render();
  }
  return worst;
}

}  // namespace

TEST_CASE("every workspace fits its page across display shapes and scales") {
  const Panel panels[] = {
      {"Extraction", &ui::drawExtractionLab},
      {"Solubility", &ui::drawSolubilitySuite},
      {"Planner", &ui::drawReactionPlanner},
      {"Toolbox", &ui::drawToolbox},
      {"Periodic", &ui::drawPeriodicTable},
      {"Properties", &ui::drawPropertiesPanel},
      {"Tools", &ui::drawToolPalette},
      {"Viewer", &ui::drawViewer3D},
  };

  // Real shapes, not round numbers: a 16:9 laptop, 1080p, a 21:9 ultrawide, a
  // 5:4 dock and a 4K surface. The ultrawide and the 5:4 are the ones that
  // catch a layout branching on width alone.
  const Surface surfaces[] = {
      {"1366x768", 1366.0f, 768.0f},   {"1920x1080", 1920.0f, 1080.0f},
      {"2560x1080", 2560.0f, 1080.0f}, {"1280x1024", 1280.0f, 1024.0f},
      {"3840x2160", 3840.0f, 2160.0f},
  };

  // Approximates the OS display scale and the app's own CHEMCAD_UI_SCALE: both
  // land on the font size, which is the unit every layout is written against.
  // 2.5 is not academic -- a 4K TV at 200% zoom with the default 1.25 preference
  // resolves to exactly that.
  const float scales[] = {1.0f, 1.25f, 1.75f, 2.5f};

  for (const Surface& surface : surfaces) {
    for (const float scale : scales) {
      for (const Panel& panel : panels) {
        HeadlessImGui gui(surface.width, surface.height, scale);
        ui::AppState st;
        const float overflow = overflowOf(st, panel.draw, panel.name);
        INFO(panel.name << " at " << surface.name << " scale " << scale
                        << " overflowed its page by " << overflow << " px");
        CHECK(overflow == 0.0f);
      }
    }
  }
}

TEST_CASE("the extraction console fits on every one of its tabs") {
  // The stage stays put while the console swaps, so each tab is a separate
  // layout and each has to fit on its own.
  for (int tab = 0; tab <= 2; ++tab) {
    for (const float scale : {1.0f, 1.5f}) {
      HeadlessImGui gui(1366.0f, 768.0f, scale);
      ui::AppState st;
      st.solubility.extractionTab = tab;
      const float overflow = overflowOf(st, &ui::drawExtractionLab, "Extraction");
      INFO("extraction tab " << tab << " at scale " << scale << " overflowed by "
                             << overflow << " px");
      CHECK(overflow == 0.0f);
    }
  }
}

#include "imgui.h"

#include <filesystem>

#include "core/paths.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"

namespace chemcad::ui {

namespace {

style::Metrics gMetrics;

ImFont* gBody = nullptr;
ImFont* gSemibold = nullptr;
ImFont* gMono = nullptr;

ImFont* addFont(const std::filesystem::path& path) {
  if (!std::filesystem::exists(path)) return nullptr;
  return ImGui::GetIO().Fonts->AddFontFromFileTTF(path.string().c_str(), 13.0f);
}

}  // namespace

namespace style {

ImU32 u32(ImVec4 c, float alphaMul) {
  c.w *= alphaMul;
  return ImGui::ColorConvertFloat4ToU32(c);
}

ImU32 mix(ImVec4 a, ImVec4 b, float t, float alphaMul) {
  const ImVec4 c{a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t,
                 a.z + (b.z - a.z) * t, (a.w + (b.w - a.w) * t) * alphaMul};
  return ImGui::ColorConvertFloat4ToU32(c);
}

const Metrics& metrics() { return gMetrics; }

namespace fonts {

void load() {
  const std::filesystem::path dir = core::assetsDir() / "fonts";
  gBody = addFont(dir / "Inter-Regular.ttf");
  gSemibold = addFont(dir / "Inter-SemiBold.ttf");
  gMono = addFont(dir / "JetBrainsMono-Regular.ttf");
  if (!gBody) gBody = ImGui::GetIO().Fonts->AddFontDefault();
  if (!gSemibold) gSemibold = gBody;
  if (!gMono) gMono = gBody;
}

ImFont* body() { return gBody; }
ImFont* semibold() { return gSemibold ? gSemibold : gBody; }
ImFont* mono() { return gMono ? gMono : (gBody ? gBody : gSemibold); }

}  // namespace fonts

bool pushFont(ImFont* font) {
  if (!font) return false;
  ImGui::PushFont(font);
  return true;
}

void popFont(bool pushed) {
  if (pushed) ImGui::PopFont();
}

}  // namespace style

void applyTheme(float uiScale) {
  ImGui::StyleColorsDark();
  ImGuiStyle& s = ImGui::GetStyle();
  using namespace style::col;

  // -------------------------------------------------------------- geometry
  s.WindowRounding = 6.0f;
  s.ChildRounding = 6.0f;
  s.FrameRounding = 5.0f;
  s.GrabRounding = 5.0f;
  s.PopupRounding = 8.0f;
  s.ScrollbarRounding = 8.0f;
  s.TabRounding = 5.0f;
  s.WindowPadding = ImVec2(10, 10);
  s.FramePadding = ImVec2(8, 5);
  s.ItemSpacing = ImVec2(8, 6);
  s.ItemInnerSpacing = ImVec2(6, 5);
  s.IndentSpacing = 18.0f;
  s.ScrollbarSize = 11.0f;
  s.GrabMinSize = 12.0f;
  s.WindowBorderSize = 1.0f;
  s.ChildBorderSize = 1.0f;
  s.FrameBorderSize = 1.0f;  // hairline outlines on every framed widget
  s.PopupBorderSize = 1.0f;
  s.TabBorderSize = 0.0f;
  s.TabBarOverlineSize = 2.0f;
  s.WindowTitleAlign = ImVec2(0.5f, 0.5f);
  // Docked panels don't expose the stock "Hide tab bar" menu: ChemCAD owns
  // the layout and keeps panel titles consistently visible.
  s.WindowMenuButtonPosition = ImGuiDir_None;
  s.SeparatorTextBorderSize = 1.0f;
  s.SeparatorTextPadding = ImVec2(16, 2);
  s.DockingSeparatorSize = 2.0f;

  // --------------------------------------------------------------- palette
  ImVec4* c = s.Colors;
  c[ImGuiCol_Text] = Text;
  c[ImGuiCol_TextDisabled] = TextFaint;
  c[ImGuiCol_TextLink] = Teal;
  c[ImGuiCol_TextSelectedBg] = {Accent.x, Accent.y, Accent.z, 0.30f};

  c[ImGuiCol_WindowBg] = BgPanel;
  c[ImGuiCol_ChildBg] = {0, 0, 0, 0};  // transparent; cards opt into surfaces
  c[ImGuiCol_PopupBg] = {BgRaised.x, BgRaised.y, BgRaised.z, 0.98f};
  c[ImGuiCol_Border] = Border;
  c[ImGuiCol_BorderShadow] = {0, 0, 0, 0};

  c[ImGuiCol_FrameBg] = BgSurface;
  c[ImGuiCol_FrameBgHovered] = BgRaised;
  c[ImGuiCol_FrameBgActive] = {BorderStrong.x, BorderStrong.y, BorderStrong.z, 0.55f};

  c[ImGuiCol_TitleBg] = BgDeep;
  c[ImGuiCol_TitleBgActive] = BgSurface;
  c[ImGuiCol_TitleBgCollapsed] = BgDeep;
  c[ImGuiCol_MenuBarBg] = BgDeep;

  c[ImGuiCol_ScrollbarBg] = {0, 0, 0, 0};
  c[ImGuiCol_ScrollbarGrab] = {BorderStrong.x, BorderStrong.y, BorderStrong.z, 0.70f};
  c[ImGuiCol_ScrollbarGrabHovered] = BorderStrong;
  c[ImGuiCol_ScrollbarGrabActive] = AccentActive;

  c[ImGuiCol_CheckMark] = Accent;
  c[ImGuiCol_SliderGrab] = Accent;
  c[ImGuiCol_SliderGrabActive] = AccentHover;

  c[ImGuiCol_Button] = BgRaised;
  c[ImGuiCol_ButtonHovered] = {BorderStrong.x, BorderStrong.y, BorderStrong.z, 0.75f};
  c[ImGuiCol_ButtonActive] = {AccentActive.x, AccentActive.y, AccentActive.z, 0.85f};

  c[ImGuiCol_Header] = {BgRaised.x, BgRaised.y, BgRaised.z, 0.70f};
  c[ImGuiCol_HeaderHovered] = {BorderStrong.x, BorderStrong.y, BorderStrong.z, 0.55f};
  c[ImGuiCol_HeaderActive] = {Accent.x, Accent.y, Accent.z, 0.30f};

  c[ImGuiCol_Separator] = Border;
  c[ImGuiCol_SeparatorHovered] = AccentActive;
  c[ImGuiCol_SeparatorActive] = Accent;

  c[ImGuiCol_ResizeGrip] = {Accent.x, Accent.y, Accent.z, 0.15f};
  c[ImGuiCol_ResizeGripHovered] = {Accent.x, Accent.y, Accent.z, 0.45f};
  c[ImGuiCol_ResizeGripActive] = Accent;

  c[ImGuiCol_InputTextCursor] = Accent;  // ImGuiCol_InputTextCursor (1.92)
  c[ImGuiCol_TabHovered] = {BgRaised.x, BgRaised.y, BgRaised.z, 0.90f};
  c[ImGuiCol_Tab] = {BgDeep.x, BgDeep.y, BgDeep.z, 0.90f};
  c[ImGuiCol_TabSelected] = BgPanel;  // blends the tab into its panel
  c[ImGuiCol_TabSelectedOverline] = Accent;
  c[ImGuiCol_TabDimmed] = {BgDeep.x, BgDeep.y, BgDeep.z, 0.80f};
  c[ImGuiCol_TabDimmedSelected] = {BgPanel.x, BgPanel.y, BgPanel.z, 0.80f};
  c[ImGuiCol_TabDimmedSelectedOverline] = {Accent.x, Accent.y, Accent.z, 0.40f};

  c[ImGuiCol_DockingPreview] = {Accent.x, Accent.y, Accent.z, 0.45f};
  c[ImGuiCol_DockingEmptyBg] = BgDeep;

  c[ImGuiCol_TableHeaderBg] = BgSurface;
  c[ImGuiCol_TableBorderStrong] = BorderStrong;
  c[ImGuiCol_TableBorderLight] = Border;
  c[ImGuiCol_TableRowBg] = {0, 0, 0, 0};
  c[ImGuiCol_TableRowBgAlt] = {1, 1, 1, 0.02f};

  c[ImGuiCol_NavCursor] = Accent;
  c[ImGuiCol_ModalWindowDimBg] = {0, 0, 0, 0.55f};

  // --------------------------------------------------------------- metrics
  gMetrics.radiusSm = 4.0f * uiScale;
  gMetrics.radiusMd = 6.0f * uiScale;
  gMetrics.radiusLg = 10.0f * uiScale;
  gMetrics.hairline = 1.0f;
  gMetrics.gap = 8.0f * uiScale;
  gMetrics.iconSize = 22.0f * uiScale;

  s.ScaleAllSizes(uiScale);
  s.FontScaleDpi = uiScale;
}

}  // namespace chemcad::ui

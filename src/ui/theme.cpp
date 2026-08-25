#include "imgui.h"

#include <filesystem>

#include "core/paths.hpp"
#include "ui/display_scale.hpp"
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

ImVec4 withAlpha(ImVec4 color, float alphaMul) {
  color.w *= alphaMul;
  return color;
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
  if (gBody) return;

  const std::filesystem::path dir = core::assetsDir() / "fonts";
  gBody = addFont(dir / "Inter-Regular.ttf");
  gSemibold = addFont(dir / "Inter-SemiBold.ttf");
  gMono = addFont(dir / "JetBrainsMono-Regular.ttf");
  if (!gBody) gBody = ImGui::GetIO().Fonts->AddFontDefault();
  if (!gSemibold) gSemibold = gBody;
  if (!gMono) gMono = gBody;
}

// These accessors are documented never to return null: callers measure with
// `font->CalcTextSizeA(...)`, which would dereference it. When the TTFs are
// absent -- headless tests, or a stripped install -- the atlas still holds
// ImGui's built-in font, so fall back to whatever the frame is currently using
// rather than handing back a pointer every call site would have to check.
ImFont* body() { return gBody ? gBody : ImGui::GetFont(); }
ImFont* semibold() { return gSemibold ? gSemibold : body(); }
ImFont* mono() { return gMono ? gMono : body(); }

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
  // One factor for the whole app: every non-text length goes through dp(), and
  // the sketch camera's pixels-per-bond rides on it too.
  setDisplayScale(uiScale);
  uiScale = displayScale();

  gMetrics.radiusSm = 2.0f * uiScale;
  gMetrics.radiusMd = 4.0f * uiScale;
  gMetrics.radiusLg = 10.0f * uiScale;
  gMetrics.hairline = 1.0f * uiScale;
  gMetrics.gap = 8.0f * uiScale;
  gMetrics.iconSize = 22.0f * uiScale;

  ImGuiStyle& s = ImGui::GetStyle();
  s = ImGuiStyle{};
  ImGui::StyleColorsDark(&s);
  s.ScaleAllSizes(uiScale);
  s.FontScaleDpi = uiScale;

  style::fonts::load();
  const float fontSize = style::fonts::body()->LegacySize * uiScale;
  const float framePaddingY = fontSize * 0.25f;
  using namespace style::col;

  // -------------------------------------------------------------- geometry
  s.WindowRounding = gMetrics.radiusMd;
  s.ChildRounding = gMetrics.radiusMd;
  s.FrameRounding = gMetrics.radiusSm;
  s.GrabRounding = gMetrics.radiusSm;
  s.PopupRounding = gMetrics.radiusLg;
  s.ScrollbarRounding = gMetrics.radiusSm;
  s.TabRounding = gMetrics.radiusSm;
  s.ImageRounding = gMetrics.radiusSm;
  s.TreeLinesRounding = gMetrics.radiusSm;
  s.MenuItemRounding = gMetrics.radiusSm;
  s.SelectableRounding = gMetrics.radiusSm;
  s.DragDropTargetRounding = gMetrics.radiusSm;

  s.WindowPadding = ImVec2(gMetrics.gap * 1.25f, gMetrics.gap);
  s.FramePadding = ImVec2(gMetrics.gap, framePaddingY);
  s.ItemSpacing = ImVec2(gMetrics.gap, gMetrics.gap * 0.625f);
  s.ItemInnerSpacing = ImVec2(gMetrics.gap * 0.75f, gMetrics.gap * 0.5f);
  s.CellPadding = ImVec2(gMetrics.gap * 0.5f, fontSize * 0.125f);
  s.IndentSpacing = fontSize + gMetrics.gap * 0.625f;
  s.ScrollbarSize = fontSize * 0.55f;
  s.ScrollbarPadding = gMetrics.hairline;
  s.GrabMinSize = fontSize;

  s.WindowBorderSize = gMetrics.hairline;
  s.ChildBorderSize = gMetrics.hairline;
  s.FrameBorderSize = gMetrics.hairline;
  s.TabBorderSize = gMetrics.hairline;
  s.PopupBorderSize = gMetrics.hairline;
  s.TabBarBorderSize = gMetrics.hairline;
  s.TabBarOverlineSize = gMetrics.hairline * 2.0f;
  s.ImageBorderSize = gMetrics.hairline;
  s.TreeLinesSize = gMetrics.hairline;
  s.InputTextCursorSize = gMetrics.hairline;
  s.SeparatorSize = gMetrics.hairline;
  s.SeparatorTextBorderSize = gMetrics.hairline;
  s.DragDropTargetBorderSize = gMetrics.hairline * 2.0f;
  s.DragDropTargetPadding = gMetrics.gap * 0.5f;

  s.WindowTitleAlign = ImVec2(0.5f, 0.5f);
  // Docked panels don't expose the stock "Hide tab bar" menu: ChemCAD owns
  // the layout and keeps panel titles consistently visible.
  s.WindowMenuButtonPosition = ImGuiDir_None;
  s.SeparatorTextPadding = ImVec2(gMetrics.gap * 2.0f, framePaddingY * 0.5f);
  s.DockingSeparatorSize = gMetrics.hairline * 2.0f;

  // --------------------------------------------------------------- palette
  ImVec4* c = s.Colors;
  c[ImGuiCol_Text] = Text;
  c[ImGuiCol_TextDisabled] = TextFaint;
  c[ImGuiCol_TextLink] = AccentHover;
  c[ImGuiCol_TextSelectedBg] = withAlpha(Accent, 0.30f);

  c[ImGuiCol_WindowBg] = BgPanel;
  c[ImGuiCol_ChildBg] = withAlpha(BgPanel, 0.0f);
  c[ImGuiCol_PopupBg] = withAlpha(BgRaised, 0.98f);
  c[ImGuiCol_Border] = Border;
  c[ImGuiCol_BorderShadow] = withAlpha(BgDeep, 0.0f);

  c[ImGuiCol_FrameBg] = BgSurface;
  c[ImGuiCol_FrameBgHovered] = BgRaised;
  c[ImGuiCol_FrameBgActive] = withAlpha(BgRaised, 0.88f);

  c[ImGuiCol_TitleBg] = BgDeep;
  c[ImGuiCol_TitleBgActive] = BgSurface;
  c[ImGuiCol_TitleBgCollapsed] = BgDeep;
  c[ImGuiCol_MenuBarBg] = BgDeep;

  c[ImGuiCol_ScrollbarBg] = withAlpha(BgDeep, 0.0f);
  c[ImGuiCol_ScrollbarGrab] = withAlpha(Border, 0.52f);
  c[ImGuiCol_ScrollbarGrabHovered] = withAlpha(BorderStrong, 0.82f);
  c[ImGuiCol_ScrollbarGrabActive] = AccentActive;

  c[ImGuiCol_CheckMark] = Accent;
  c[ImGuiCol_CheckboxSelectedBg] = withAlpha(Accent, 0.22f);
  c[ImGuiCol_SliderGrab] = Accent;
  c[ImGuiCol_SliderGrabActive] = AccentHover;

  c[ImGuiCol_Button] = withAlpha(Accent, 0.14f);
  c[ImGuiCol_ButtonHovered] = withAlpha(AccentHover, 0.28f);
  c[ImGuiCol_ButtonActive] = withAlpha(AccentActive, 0.42f);

  c[ImGuiCol_Header] = withAlpha(Accent, 0.12f);
  c[ImGuiCol_HeaderHovered] = withAlpha(AccentHover, 0.22f);
  c[ImGuiCol_HeaderActive] = withAlpha(AccentActive, 0.32f);

  c[ImGuiCol_Separator] = GridLine;
  c[ImGuiCol_SeparatorHovered] = AccentHover;
  c[ImGuiCol_SeparatorActive] = Accent;

  c[ImGuiCol_ResizeGrip] = withAlpha(Accent, 0.12f);
  c[ImGuiCol_ResizeGripHovered] = withAlpha(AccentHover, 0.42f);
  c[ImGuiCol_ResizeGripActive] = Accent;

  c[ImGuiCol_InputTextCursor] = Accent;
  c[ImGuiCol_TabHovered] = withAlpha(AccentHover, 0.16f);
  c[ImGuiCol_Tab] = withAlpha(BgDeep, 0.62f);
  c[ImGuiCol_TabSelected] = withAlpha(Accent, 0.12f);
  c[ImGuiCol_TabSelectedOverline] = AccentHover;
  c[ImGuiCol_TabDimmed] = withAlpha(BgDeep, 0.38f);
  c[ImGuiCol_TabDimmedSelected] = withAlpha(Accent, 0.07f);
  c[ImGuiCol_TabDimmedSelectedOverline] = withAlpha(Accent, 0.46f);

  c[ImGuiCol_DockingPreview] = withAlpha(Accent, 0.45f);
  c[ImGuiCol_DockingEmptyBg] = BgDeep;

  c[ImGuiCol_PlotLines] = Data;
  c[ImGuiCol_PlotLinesHovered] = DataBright;
  c[ImGuiCol_PlotHistogram] = DataDim;
  c[ImGuiCol_PlotHistogramHovered] = DataBright;

  c[ImGuiCol_TableHeaderBg] = BgSurface;
  c[ImGuiCol_TableBorderStrong] = DataDim;
  c[ImGuiCol_TableBorderLight] = GridLine;
  c[ImGuiCol_TableRowBg] = withAlpha(BgPanel, 0.0f);
  c[ImGuiCol_TableRowBgAlt] = withAlpha(BgRaised, 0.14f);

  c[ImGuiCol_TreeLines] = GridLine;
  c[ImGuiCol_DragDropTarget] = AccentHover;
  c[ImGuiCol_DragDropTargetBg] = withAlpha(Accent, 0.12f);
  c[ImGuiCol_UnsavedMarker] = Danger;
  c[ImGuiCol_NavCursor] = Data;
  c[ImGuiCol_NavWindowingHighlight] = DataBright;
  c[ImGuiCol_NavWindowingDimBg] = withAlpha(BgDeep, 0.62f);
  c[ImGuiCol_ModalWindowDimBg] = withAlpha(BgDeep, 0.72f);
}

}  // namespace chemcad::ui

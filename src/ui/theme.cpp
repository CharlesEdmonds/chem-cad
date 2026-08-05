#include "imgui.h"

#include "ui/ui.hpp"

namespace chemcad::ui {

void applyTheme(float uiScale) {
  ImGui::StyleColorsDark();
  ImGuiStyle& s = ImGui::GetStyle();

  s.WindowRounding = 4.0f;
  s.ChildRounding = 4.0f;
  s.FrameRounding = 3.0f;
  s.GrabRounding = 3.0f;
  s.PopupRounding = 4.0f;
  s.ScrollbarRounding = 6.0f;
  s.TabRounding = 4.0f;
  s.WindowPadding = ImVec2(8, 8);
  s.FramePadding = ImVec2(7, 4);
  s.ItemSpacing = ImVec2(7, 6);
  s.WindowBorderSize = 1.0f;
  s.FrameBorderSize = 0.0f;
  s.WindowTitleAlign = ImVec2(0.5f, 0.5f);

  ImVec4* c = s.Colors;
  c[ImGuiCol_WindowBg] = ImVec4(0.11f, 0.12f, 0.14f, 1.00f);
  c[ImGuiCol_ChildBg] = ImVec4(0.13f, 0.14f, 0.17f, 1.00f);
  c[ImGuiCol_PopupBg] = ImVec4(0.10f, 0.11f, 0.13f, 0.98f);
  c[ImGuiCol_Border] = ImVec4(0.24f, 0.26f, 0.30f, 1.00f);
  c[ImGuiCol_FrameBg] = ImVec4(0.18f, 0.20f, 0.24f, 1.00f);
  c[ImGuiCol_FrameBgHovered] = ImVec4(0.24f, 0.27f, 0.33f, 1.00f);
  c[ImGuiCol_FrameBgActive] = ImVec4(0.28f, 0.32f, 0.40f, 1.00f);
  c[ImGuiCol_TitleBg] = ImVec4(0.09f, 0.10f, 0.12f, 1.00f);
  c[ImGuiCol_TitleBgActive] = ImVec4(0.14f, 0.17f, 0.22f, 1.00f);
  c[ImGuiCol_MenuBarBg] = ImVec4(0.13f, 0.14f, 0.17f, 1.00f);
  c[ImGuiCol_Header] = ImVec4(0.20f, 0.31f, 0.45f, 0.70f);
  c[ImGuiCol_HeaderHovered] = ImVec4(0.26f, 0.40f, 0.58f, 0.85f);
  c[ImGuiCol_HeaderActive] = ImVec4(0.30f, 0.47f, 0.68f, 1.00f);
  c[ImGuiCol_Button] = ImVec4(0.20f, 0.23f, 0.28f, 1.00f);
  c[ImGuiCol_ButtonHovered] = ImVec4(0.28f, 0.36f, 0.48f, 1.00f);
  c[ImGuiCol_ButtonActive] = ImVec4(0.34f, 0.47f, 0.66f, 1.00f);
  c[ImGuiCol_Tab] = ImVec4(0.15f, 0.17f, 0.21f, 1.00f);
  c[ImGuiCol_TabHovered] = ImVec4(0.28f, 0.40f, 0.56f, 1.00f);
  c[ImGuiCol_TabSelected] = ImVec4(0.22f, 0.33f, 0.48f, 1.00f);
  c[ImGuiCol_Separator] = ImVec4(0.24f, 0.26f, 0.30f, 1.00f);
  c[ImGuiCol_Text] = ImVec4(0.90f, 0.92f, 0.95f, 1.00f);
  c[ImGuiCol_TextDisabled] = ImVec4(0.48f, 0.51f, 0.56f, 1.00f);
  c[ImGuiCol_CheckMark] = ImVec4(0.45f, 0.72f, 1.00f, 1.00f);
  c[ImGuiCol_SliderGrab] = ImVec4(0.36f, 0.55f, 0.80f, 1.00f);
  c[ImGuiCol_DockingPreview] = ImVec4(0.26f, 0.45f, 0.70f, 0.70f);

  s.ScaleAllSizes(uiScale);
  ImGui::GetStyle().FontScaleDpi = uiScale;
}

}  // namespace chemcad::ui

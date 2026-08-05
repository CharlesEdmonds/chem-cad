#include <algorithm>

#include "imgui.h"

#include "ui/canvas_internal.hpp"
#include "ui/ui.hpp"

namespace chemcad::ui {

void drawCanvas(AppState& st) {
  static canvas::Runtime runtime;

  ImVec2 size = ImGui::GetContentRegionAvail();
  size.x = std::max(size.x, 1.0f);
  size.y = std::max(size.y, 1.0f);
  const ImVec2 origin = ImGui::GetCursorScreenPos();
  ImGui::InvisibleButton("##canvas", size,
                         ImGuiButtonFlags_MouseButtonLeft |
                             ImGuiButtonFlags_MouseButtonRight |
                             ImGuiButtonFlags_MouseButtonMiddle);
  const ImVec2 rectMin = ImGui::GetItemRectMin();
  const ImVec2 rectMax = ImGui::GetItemRectMax();
  const bool hovered = ImGui::IsItemHovered();
  const bool active = ImGui::IsItemActive();

  st.canvasOrigin = {rectMin.x, rectMin.y};
  st.canvasSize = {size.x, size.y};
  const canvas::CanvasRect rect{{rectMin.x, rectMin.y},
                                {size.x, size.y},
                                rectMin,
                                rectMax};

  canvas::hitTest(st, rect, hovered);
  canvas::handleInput(st, runtime, rect, hovered, active);
  canvas::hitTest(st, rect, hovered);

  ImDrawList* draw = ImGui::GetWindowDrawList();
  draw->PushClipRect(rectMin, rectMax, true);
  canvas::render(st, runtime, rect);
  draw->PopClipRect();
}

}  // namespace chemcad::ui

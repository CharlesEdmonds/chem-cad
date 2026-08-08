#include <algorithm>
#include <cstdio>

#include "imgui.h"

#include "ui/canvas_internal.hpp"
#include "ui/icons.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"
#include "ui/widgets.hpp"

namespace chemcad::ui {
namespace {

// One ghost button of the zoom overlay.
bool overlayButton(const char* id, ImVec2 pos, ImVec2 size, const char* tooltip) {
  ImGui::SetCursorScreenPos(pos);
  ImGui::InvisibleButton(id, size);
  const bool clicked = ImGui::IsItemClicked();
  if (ImGui::IsItemHovered() && tooltip) ImGui::SetTooltip("%s", tooltip);
  return clicked;
}

// Bottom-right canvas chrome: zoom out / zoom % / zoom in / fit. Submitted
// after the canvas's own button, so they win hover where they overlap.
void drawZoomOverlay(AppState& st, ImVec2 rectMin, ImVec2 rectMax) {
  const style::Metrics& m = style::metrics();
  const float h = ImGui::GetFontSize() * 1.7f;
  const float gap = m.gap * 0.5f;
  const float margin = m.gap * 1.2f;
  const core::Vec2 centre{(rectMin.x + rectMax.x) * 0.5f, (rectMin.y + rectMax.y) * 0.5f};
  const core::Vec2 origin{rectMin.x, rectMin.y};

  char pct[8];
  std::snprintf(pct, sizeof(pct), "%.0f%%", std::round(st.cam.zoom * 100.0f));
  const float pctW = ImGui::CalcTextSize(pct).x + m.gap * 1.2f;

  const float totalW = h * 3.0f + pctW + gap * 3.0f;
  const float y = rectMax.y - margin - h;
  float x = rectMax.x - margin - totalW;

  ImDrawList* dl = ImGui::GetWindowDrawList();
  const auto chrome = [&](ImVec2 p, ImVec2 s, bool hovered) {
    dl->AddRectFilled(p, ImVec2(p.x + s.x, p.y + s.y),
                      style::u32(style::col::BgRaised, hovered ? 0.95f : 0.72f), m.radiusSm);
    dl->AddRect(p, ImVec2(p.x + s.x, p.y + s.y), style::u32(style::col::BorderStrong, 0.8f),
                m.radiusSm, 0, m.hairline);
  };

  if (overlayButton("##zoom_out", ImVec2(x, y), ImVec2(h, h), "Zoom out"))
    st.cam.zoomAt(1.0f / 1.2f, centre, origin);
  chrome(ImVec2(x, y), ImVec2(h, h), ImGui::IsItemHovered());
  icons::draw(dl, icons::Icon::Minus, ImVec2(x + h * 0.5f, y + h * 0.5f), h * 0.34f,
              style::u32(style::col::Text));
  x += h + gap;

  if (overlayButton("##zoom_reset", ImVec2(x, y), ImVec2(pctW, h), "Reset zoom (Ctrl+0)"))
    st.cam.zoom = 1.0f;
  chrome(ImVec2(x, y), ImVec2(pctW, h), ImGui::IsItemHovered());
  {
    const ImVec2 w = ImGui::CalcTextSize(pct);
    dl->AddText(ImVec2(x + (pctW - w.x) * 0.5f, y + (h - w.y) * 0.5f), style::u32(style::col::Text),
                pct);
  }
  x += pctW + gap;

  if (overlayButton("##zoom_in", ImVec2(x, y), ImVec2(h, h), "Zoom in (mouse wheel)"))
    st.cam.zoomAt(1.2f, centre, origin);
  chrome(ImVec2(x, y), ImVec2(h, h), ImGui::IsItemHovered());
  icons::draw(dl, icons::Icon::Plus, ImVec2(x + h * 0.5f, y + h * 0.5f), h * 0.34f,
              style::u32(style::col::Text));
  x += h + gap;

  if (overlayButton("##zoom_fit", ImVec2(x, y), ImVec2(h, h), "Fit to window (Ctrl+F)"))
    st.cam.fit(st.doc, st.canvasSize);
  chrome(ImVec2(x, y), ImVec2(h, h), ImGui::IsItemHovered());
  icons::draw(dl, icons::Icon::ZoomFit, ImVec2(x + h * 0.5f, y + h * 0.5f), h * 0.34f,
              style::u32(style::col::Text));
}

}  // namespace

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

  drawZoomOverlay(st, rectMin, rectMax);
}

}  // namespace chemcad::ui

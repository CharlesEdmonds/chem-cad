#include "ui/widgets.hpp"

#include <algorithm>

#include "imgui_internal.h"

namespace chemcad::ui::widgets {

float hoverT(ImGuiID id, bool hovered) {
  ImGuiStorage* storage = ImGui::GetStateStorage();
  float t = storage->GetFloat(id, 0.0f);
  const float target = hovered ? 1.0f : 0.0f;
  const float step = style::metrics().animSpeed * ImGui::GetIO().DeltaTime;
  t += std::clamp(target - t, -step, step);
  storage->SetFloat(id, t);
  return t;
}

void sectionHeader(const char* label, ImVec4 accent) {
  const style::Metrics& m = style::metrics();
  ImGui::Spacing();
  const float fs = ImGui::GetFontSize();
  const ImVec2 pos = ImGui::GetCursorScreenPos();
  const float avail = ImGui::GetContentRegionAvail().x;

  // Reserve the header line, then draw into it.
  ImGui::Dummy(ImVec2(avail, fs + m.gap * 0.9f));
  ImDrawList* dl = ImGui::GetWindowDrawList();

  const float tickW = std::max(2.0f, fs * 0.16f);
  const float tickH = fs * 0.85f;
  const float tickY = pos.y + (fs - tickH) * 0.5f + 1.0f;
  dl->AddRectFilled(ImVec2(pos.x, tickY), ImVec2(pos.x + tickW, tickY + tickH),
                    style::u32(accent), tickW * 0.5f);

  const ImVec2 textSize = ImGui::CalcTextSize(label);
  const bool pushed = style::pushFont(style::fonts::semibold());
  dl->AddText(ImVec2(pos.x + tickW + m.gap * 0.75f, pos.y), style::u32(style::col::TextDim),
              label);
  style::popFont(pushed);

  const float lineX = pos.x + tickW + m.gap * 0.75f + textSize.x + m.gap;
  const float lineY = pos.y + fs * 0.5f + 1.0f;
  if (lineX < pos.x + avail - m.gap) {
    dl->AddLine(ImVec2(lineX, lineY), ImVec2(pos.x + avail, lineY),
                style::u32(style::col::Border), m.hairline);
  }
  ImGui::Spacing();
}

bool primaryButton(const char* label, ImVec2 size) {
  const bool pushed = style::pushFont(style::fonts::semibold());
  ImGui::PushStyleColor(ImGuiCol_Button, style::col::Accent);
  ImGui::PushStyleColor(ImGuiCol_ButtonHovered, style::col::AccentHover);
  ImGui::PushStyleColor(ImGuiCol_ButtonActive, style::col::AccentActive);
  ImGui::PushStyleColor(ImGuiCol_Text, style::col::OnAccent);
  ImGui::PushStyleColor(ImGuiCol_Border, style::col::AccentActive);
  const bool clicked = ImGui::Button(label, size);
  ImGui::PopStyleColor(5);
  style::popFont(pushed);
  return clicked;
}

bool ghostButton(const char* label, ImVec2 size) {
  if (size.x <= 0.0f) {
    const ImGuiStyle& s = ImGui::GetStyle();
    size.x = ImGui::CalcTextSize(label).x + s.FramePadding.x * 2.0f;
  }
  if (size.y <= 0.0f) size.y = ImGui::GetFrameHeight();

  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::InvisibleButton(label, size);
  const bool clicked = ImGui::IsItemClicked();
  const float t = hoverT(ImGui::GetItemID(), ImGui::IsItemHovered());
  const ImVec2 max(min.x + size.x, min.y + size.y);
  const style::Metrics& m = style::metrics();

  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddRectFilled(min, max, style::mix(style::col::BgSurface, style::col::BgRaised, t,
                                         0.35f + 0.65f * t),
                    m.radiusMd);
  dl->AddRect(min, max, style::mix(style::col::Border, style::col::BorderStrong, t),
              m.radiusMd, 0, m.hairline);

  const ImVec2 textSize = ImGui::CalcTextSize(label);
  dl->AddText(ImVec2((min.x + max.x - textSize.x) * 0.5f,
                     (min.y + max.y - textSize.y) * 0.5f),
              style::mix(style::col::TextDim, style::col::Text, 0.35f + 0.65f * t), label);
  return clicked;
}

bool iconButton(const char* id, icons::Icon icon, ImVec2 size, bool selected,
                const char* tooltip) {
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::InvisibleButton(id, size);
  const bool clicked = ImGui::IsItemClicked();
  const bool hovered = ImGui::IsItemHovered();
  const float t = hoverT(ImGui::GetItemID(), hovered);
  const ImVec2 max(min.x + size.x, min.y + size.y);
  const style::Metrics& m = style::metrics();
  const ImVec2 centre((min.x + max.x) * 0.5f, (min.y + max.y) * 0.5f);
  const float glyph = std::min(m.iconSize, std::min(size.x, size.y) * 0.62f);

  ImDrawList* dl = ImGui::GetWindowDrawList();
  ImU32 glyphColor;
  if (selected) {
    dl->AddRectFilled(min, max, style::mix(style::col::Accent, style::col::AccentHover, t),
                      m.radiusMd);
    dl->AddRect(min, max, style::u32(style::col::AccentActive), m.radiusMd, 0,
                m.hairline);
    glyphColor = style::u32(style::col::OnAccent);
  } else {
    dl->AddRectFilled(min, max, style::mix(style::col::BgSurface, style::col::BgRaised, t),
                      m.radiusMd);
    dl->AddRect(min, max, style::mix(style::col::Border, style::col::BorderStrong, t),
                m.radiusMd, 0, m.hairline);
    glyphColor = style::mix(style::col::TextDim, style::col::Text, 0.25f + 0.75f * t);
  }
  icons::draw(dl, icon, centre, glyph, glyphColor);

  if (hovered && tooltip) ImGui::SetTooltip("%s", tooltip);
  return clicked;
}

void badge(const char* text, ImVec4 color) {
  const style::Metrics& m = style::metrics();
  const float fs = ImGui::GetFontSize();
  const float padX = m.gap * 0.8f;
  const ImVec2 textSize = ImGui::CalcTextSize(text);
  const ImVec2 size(textSize.x + padX * 2.0f, fs + m.gap * 0.55f);

  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::Dummy(size);
  const ImVec2 max(min.x + size.x, min.y + size.y);

  ImDrawList* dl = ImGui::GetWindowDrawList();
  const float radius = size.y * 0.5f;
  dl->AddRectFilled(min, max, style::u32(color, 0.16f), radius);
  dl->AddRect(min, max, style::u32(color, 0.55f), radius, 0, m.hairline);
  dl->AddText(ImVec2(min.x + padX, min.y + (size.y - fs) * 0.5f), style::u32(color), text);
}

void statCard(const char* label, const char* value, ImVec2 size) {
  const style::Metrics& m = style::metrics();
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::Dummy(size);
  const ImVec2 max(min.x + size.x, min.y + size.y);

  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddRectFilled(min, max, style::u32(style::col::BgSurface), m.radiusMd);
  dl->AddRect(min, max, style::u32(style::col::Border), m.radiusMd, 0, m.hairline);

  const float pad = m.gap * 0.9f;
  const float labelSize = ImGui::GetFontSize() * 0.82f;
  const bool mono = style::pushFont(style::fonts::mono());
  const ImVec2 valueSize = ImGui::CalcTextSize(value);
  dl->AddText(ImVec2(min.x + pad, min.y + pad), style::u32(style::col::Text), value);
  style::popFont(mono);
  dl->AddText(nullptr, labelSize, ImVec2(min.x + pad, min.y + size.y - pad - labelSize),
              style::u32(style::col::TextFaint), label);
  (void)valueSize;
}

bool beginCard(const char* id, ImVec2 size, ImVec4 bg, ImGuiWindowFlags extraFlags) {
  ImGui::PushStyleColor(ImGuiCol_ChildBg, bg);
  ImGui::PushStyleVar(ImGuiStyleVar_ChildRounding, style::metrics().radiusMd);
  const ImGuiChildFlags childFlags =
      ImGuiChildFlags_Borders | (size.y == 0.0f ? ImGuiChildFlags_AutoResizeY : 0);
  const bool open = ImGui::BeginChild(id, size, childFlags, extraFlags);
  if (!open) endCard();
  return open;
}

void endCard() {
  ImGui::EndChild();
  ImGui::PopStyleVar();
  ImGui::PopStyleColor();
}

}  // namespace chemcad::ui::widgets

#include "ui/widgets.hpp"

#include <algorithm>
#include <cctype>
#include <cfloat>
#include <cstddef>
#include <cstdio>
#include <vector>

#include "sol/solvent.hpp"

#include "imgui_internal.h"

namespace chemcad::ui::widgets {

namespace {

int resizeStringInput(ImGuiInputTextCallbackData* data) {
  if (data->EventFlag != ImGuiInputTextFlags_CallbackResize) return 0;
  auto* value = static_cast<std::string*>(data->UserData);
  value->resize(static_cast<std::size_t>(data->BufTextLen));
  data->Buf = value->data();
  return 0;
}

std::vector<ImVec2> toolbarStarts;

}  // namespace

bool containsCaseInsensitive(std::string_view haystack, std::string_view needle) {
  if (needle.empty()) return true;
  return std::search(haystack.begin(), haystack.end(), needle.begin(), needle.end(),
                     [](unsigned char left, unsigned char right) {
                       return std::tolower(left) == std::tolower(right);
                     }) != haystack.end();
}

bool stringInputWithHint(const char* id, const char* hint, std::string& value) {
  return stringInputWithHint(id, hint, value, 0);
}

bool stringInputWithHint(const char* id, const char* hint, std::string& value,
                         ImGuiInputTextFlags flags, bool mono) {
  flags |= ImGuiInputTextFlags_CallbackResize;
  const bool pushed = mono ? style::pushFont(style::fonts::mono()) : false;
  const bool changed = ImGui::InputTextWithHint(
      id, hint, value.data(), value.capacity() + 1, flags, resizeStringInput, &value);
  style::popFont(pushed);
  return changed;
}

bool solventCombo(const char* id, std::string& solventId, std::string& query) {
  const sol::Solvent* current = sol::findSolvent(solventId);
  const char* preview = current ? current->name.c_str() : "Select solvent";
  bool changed = false;

  // ImGui caps a combo popup at roughly eight rows. Scrolling to reach water
  // would defeat the point of curating, so the popup is sized to hold the whole
  // curated list plus the search field and its hint; a wide search result set
  // then scrolls inside the same box.
  const float rowHeight = ImGui::GetTextLineHeightWithSpacing();
  const int curatedRows = static_cast<int>(sol::commonSolvents().size());
  const float resultsHeight =
      rowHeight * static_cast<float>(std::clamp(curatedRows, 6, 16));
  const float popupChrome =
      ImGui::GetFrameHeightWithSpacing() + rowHeight * 2.0f + ImGui::GetStyle().WindowPadding.y * 2.0f;
  ImGui::SetNextWindowSizeConstraints(ImVec2(ImGui::GetFontSize() * 13.0f, 0.0f),
                                      ImVec2(FLT_MAX, resultsHeight + popupChrome));

  if (ImGui::BeginCombo(id, preview)) {
    if (ImGui::IsWindowAppearing()) {
      query.clear();
      ImGui::SetKeyboardFocusHere();
    }

    char hint[64];
    std::snprintf(hint, sizeof(hint), "Search all %zu solvents...", sol::solvents().size());
    ImGui::SetNextItemWidth(-FLT_MIN);
    stringInputWithHint("##solvent_search", hint, query);
    ImGui::Separator();

    ImGui::BeginChild("##solvent_results", ImVec2(0.0f, resultsHeight),
                      ImGuiChildFlags_Borders);
    bool matched = false;
    const auto drawResult = [&](const sol::Solvent& solvent) {
      if (!query.empty() && !containsCaseInsensitive(solvent.name, query) &&
          !containsCaseInsensitive(solvent.id, query)) {
        return;
      }
      matched = true;
      const bool selected = solvent.id == solventId;
      if (ImGui::Selectable(solvent.name.c_str(), selected)) {
        solventId = solvent.id;
        changed = true;
      }
      if (selected) ImGui::SetItemDefaultFocus();
    };

    if (query.empty()) {
      for (const sol::Solvent* solvent : sol::commonSolvents()) drawResult(*solvent);
    } else {
      for (const sol::Solvent& solvent : sol::solvents()) drawResult(solvent);
    }

    if (!matched) {
      ImGui::PushStyleColor(ImGuiCol_Text, style::col::TextDim);
      ImGui::Text("No solvent matches \"%s\"", query.c_str());
      ImGui::PopStyleColor();
    }
    ImGui::EndChild();

    if (query.empty()) {
      ImGui::PushStyleColor(ImGuiCol_Text, style::col::TextDim);
      ImGui::TextUnformatted("Type above to search every solvent");
      ImGui::PopStyleColor();
    }
    ImGui::EndCombo();
  }
  return changed;
}

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
  const float glyph =
      std::min(m.iconSize * 1.25f, std::min(size.x, size.y) * 0.80f);

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

bool disclosure(const char* id, const char* label, const char* summary, bool defaultOpen,
                icons::Icon icon, ImVec4 accent) {
  const style::Metrics& m = style::metrics();
  const ImVec2 min = ImGui::GetCursorScreenPos();
  const ImVec2 size(ImGui::GetContentRegionAvail().x, ImGui::GetFrameHeight());
  ImGuiStorage* storage = ImGui::GetStateStorage();
  const ImGuiID key = ImGui::GetID(id);
  bool open = storage->GetBool(key, defaultOpen);

  ImGui::PushID(id);
  ImGui::InvisibleButton("##disclosure_row", size);
  const ImGuiID itemId = ImGui::GetItemID();
  if (ImGui::IsItemClicked()) {
    open = !open;
    storage->SetBool(key, open);
  }
  const bool hovered = ImGui::IsItemHovered();
  ImGui::PopID();
  const float hover = hoverT(itemId, hovered);
  const ImGuiID animKey = ImHashStr("##disclosure_open", 0, key);
  float openT = storage->GetFloat(animKey, defaultOpen ? 1.0f : 0.0f);
  const float step = m.animSpeed * ImGui::GetIO().DeltaTime;
  openT += std::clamp((open ? 1.0f : 0.0f) - openT, -step, step);
  storage->SetFloat(animKey, openT);

  const ImVec2 max(min.x + size.x, min.y + size.y);
  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddRectFilled(min, max,
                    style::mix(style::col::BgSurface, style::col::BgRaised,
                               std::max(hover, openT * 0.35f)),
                    m.radiusSm);

  float x = min.x + m.gap * 0.55f;
  const ImVec2 centre(x + m.iconSize * 0.32f, (min.y + max.y) * 0.5f);
  icons::draw(dl, openT > 0.5f ? icons::Icon::ChevronDown : icons::Icon::ChevronRight,
              centre, m.iconSize * 0.64f,
              style::mix(style::col::TextFaint, style::col::TextDim, hover));
  x += m.iconSize * 0.72f + m.gap * 0.35f;

  if (icon != icons::Icon::None) {
    icons::draw(dl, icon, ImVec2(x + m.iconSize * 0.5f, centre.y), m.iconSize,
                style::u32(accent));
    x += m.iconSize + m.gap * 0.55f;
  }

  const bool semibold = style::pushFont(style::fonts::semibold());
  const ImVec2 labelSize = ImGui::CalcTextSize(label);
  dl->AddText(ImVec2(x, min.y + (size.y - labelSize.y) * 0.5f),
              style::u32(style::col::Text), label);
  style::popFont(semibold);

  if (summary && summary[0] != '\0') {
    const ImVec2 summarySize = ImGui::CalcTextSize(summary);
    const float summaryX = max.x - m.gap * 0.65f - summarySize.x;
    if (summaryX > x + labelSize.x + m.gap) {
      dl->AddText(ImVec2(summaryX, min.y + (size.y - summarySize.y) * 0.5f),
                  style::u32(style::col::TextDim), summary);
    }
  }
  return open;
}

void cardHeader(icons::Icon icon, const char* title, const char* subtitle, ImVec4 accent) {
  const style::Metrics& m = style::metrics();
  const float fs = ImGui::GetFontSize();
  const ImVec2 min = ImGui::GetCursorScreenPos();
  const float avail = ImGui::GetContentRegionAvail().x;
  const float height = std::max(fs, m.iconSize) + m.gap * 0.35f;
  ImGui::Dummy(ImVec2(avail, height));

  ImDrawList* dl = ImGui::GetWindowDrawList();
  const ImVec2 glyphCentre(min.x + m.iconSize * 0.5f, min.y + height * 0.5f);
  icons::draw(dl, icon, glyphCentre, m.iconSize, style::u32(accent));

  const float titleX = min.x + m.iconSize + m.gap * 0.7f;
  const bool semibold = style::pushFont(style::fonts::semibold());
  const ImVec2 titleSize = ImGui::CalcTextSize(title);
  dl->AddText(ImVec2(titleX, min.y + (height - titleSize.y) * 0.5f),
              style::u32(style::col::Text), title);
  style::popFont(semibold);

  float lineX = titleX + titleSize.x + m.gap;
  if (subtitle && subtitle[0] != '\0') {
    const ImVec2 subtitleSize = ImGui::CalcTextSize(subtitle);
    dl->AddText(ImVec2(lineX, min.y + (height - subtitleSize.y) * 0.5f),
                style::u32(style::col::TextFaint), subtitle);
    lineX += subtitleSize.x + m.gap;
  }
  if (lineX < min.x + avail) {
    const float lineY = min.y + height * 0.5f;
    dl->AddLine(ImVec2(lineX, lineY), ImVec2(min.x + avail, lineY),
                style::u32(style::col::Border), m.hairline);
  }
}

void beginToolbar(const char* id) {
  const style::Metrics& m = style::metrics();
  const ImVec2 min = ImGui::GetCursorScreenPos();
  const float width = ImGui::GetContentRegionAvail().x;
  const float height = ImGui::GetFrameHeight() + m.gap;
  toolbarStarts.push_back(min);

  // The toolbar has a standard height and consumes the available row, so its
  // backing can be painted before child controls without draw-list channels.
  ImDrawList* dl = ImGui::GetWindowDrawList();
  const ImVec2 max(min.x + width, min.y + height);
  dl->AddRectFilled(min, max, style::u32(style::col::BgRaised), m.radiusMd);
  dl->AddRect(min, max, style::u32(style::col::Border), m.radiusMd, 0, m.hairline);
  ImGui::Dummy(ImVec2(width, height));
  ImGui::SetCursorScreenPos(ImVec2(min.x + m.gap * 0.5f, min.y + m.gap * 0.5f));
  ImGui::PushID(id);
  ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(m.gap * 0.5f, m.gap * 0.5f));
  ImGui::BeginGroup();
}

void endToolbar() {
  IM_ASSERT(!toolbarStarts.empty());
  ImGui::EndGroup();
  ImGui::PopStyleVar();
  ImGui::PopID();

  const ImVec2 min = toolbarStarts.back();
  toolbarStarts.pop_back();
  const float bottom = min.y + ImGui::GetFrameHeight() + style::metrics().gap;
  const ImVec2 cursor = ImGui::GetCursorScreenPos();
  ImGui::SetCursorScreenPos(ImVec2(min.x, std::max(cursor.y, bottom)));
}

void toolbarSeparator() {
  const style::Metrics& m = style::metrics();
  const ImVec2 min = ImGui::GetCursorScreenPos();
  const ImVec2 size(m.gap, ImGui::GetFrameHeight());
  ImGui::Dummy(size);
  const float x = min.x + size.x * 0.5f;
  ImGui::GetWindowDrawList()->AddLine(
      ImVec2(x, min.y), ImVec2(x, min.y + size.y), style::u32(style::col::Border),
      m.hairline);
}

bool segmented(const char* id, const char* const* labels, int count, int& index,
               float width) {
  if (count <= 0) return false;

  const style::Metrics& m = style::metrics();
  if (width <= 0.0f) width = ImGui::GetContentRegionAvail().x;
  const ImVec2 min = ImGui::GetCursorScreenPos();
  const ImVec2 size(width, ImGui::GetFrameHeight());
  ImGui::InvisibleButton(id, size);
  const ImGuiID itemId = ImGui::GetItemID();
  const bool hovered = ImGui::IsItemHovered();
  const float hover = hoverT(itemId, hovered);
  const float cellWidth = size.x / static_cast<float>(count);
  int hoveredCell = -1;
  if (hovered) {
    hoveredCell = std::clamp(
        static_cast<int>((ImGui::GetIO().MousePos.x - min.x) / cellWidth), 0, count - 1);
  }

  bool changed = false;
  if (ImGui::IsItemClicked()) {
    const int clickedCell = std::clamp(
        static_cast<int>((ImGui::GetIO().MousePos.x - min.x) / cellWidth), 0, count - 1);
    if (clickedCell != index) {
      index = clickedCell;
      changed = true;
    }
  }

  ImGuiStorage* storage = ImGui::GetStateStorage();
  const ImGuiID animKey = ImHashStr("##segmented_position", 0, itemId);
  const int selected = std::clamp(index, 0, count - 1);
  float position = storage->GetFloat(animKey, static_cast<float>(selected));
  const float blend = std::clamp(m.animSpeed * ImGui::GetIO().DeltaTime, 0.0f, 1.0f);
  position += (static_cast<float>(selected) - position) * blend;
  storage->SetFloat(animKey, position);

  const ImVec2 max(min.x + size.x, min.y + size.y);
  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddRectFilled(min, max,
                    style::mix(style::col::BgSurface, style::col::BgRaised, hover, 0.85f),
                    m.radiusMd);
  dl->AddRect(min, max, style::u32(style::col::Border), m.radiusMd, 0, m.hairline);
  const ImVec2 indicatorMin(min.x + cellWidth * position, min.y);
  const ImVec2 indicatorMax(indicatorMin.x + cellWidth, max.y);
  dl->AddRectFilled(indicatorMin, indicatorMax, style::u32(style::col::Accent), m.radiusMd);

  for (int i = 0; i < count; ++i) {
    const char* text = labels[i] ? labels[i] : "";
    const ImVec2 textSize = ImGui::CalcTextSize(text);
    const ImVec2 textPos(min.x + cellWidth * (static_cast<float>(i) + 0.5f) -
                             textSize.x * 0.5f,
                         min.y + (size.y - textSize.y) * 0.5f);
    const ImU32 color =
        i == index
            ? style::u32(style::col::OnAccent)
            : style::mix(style::col::TextDim, style::col::Text,
                         i == hoveredCell ? 1.0f : 0.0f);
    dl->AddText(textPos, color, text);
  }
  return changed;
}

bool segmentedIcons(const char* id, const icons::Icon* glyphs,
                    const char* const* tooltips, int count, int& index, float width) {
  if (count <= 0) return false;

  const style::Metrics& m = style::metrics();
  if (width <= 0.0f) width = ImGui::GetContentRegionAvail().x;
  const ImVec2 min = ImGui::GetCursorScreenPos();
  const ImVec2 size(width, ImGui::GetFrameHeight());
  ImGui::InvisibleButton(id, size);
  const ImGuiID itemId = ImGui::GetItemID();
  const bool hovered = ImGui::IsItemHovered();
  const float hover = hoverT(itemId, hovered);
  const float cellWidth = size.x / static_cast<float>(count);
  int hoveredCell = -1;
  if (hovered) {
    hoveredCell = std::clamp(
        static_cast<int>((ImGui::GetIO().MousePos.x - min.x) / cellWidth), 0, count - 1);
  }

  bool changed = false;
  if (ImGui::IsItemClicked()) {
    const int clickedCell = std::clamp(
        static_cast<int>((ImGui::GetIO().MousePos.x - min.x) / cellWidth), 0, count - 1);
    if (clickedCell != index) {
      index = clickedCell;
      changed = true;
    }
  }

  ImGuiStorage* storage = ImGui::GetStateStorage();
  const ImGuiID animKey = ImHashStr("##segmented_icon_position", 0, itemId);
  const int selected = std::clamp(index, 0, count - 1);
  float position = storage->GetFloat(animKey, static_cast<float>(selected));
  const float blend = std::clamp(m.animSpeed * ImGui::GetIO().DeltaTime, 0.0f, 1.0f);
  position += (static_cast<float>(selected) - position) * blend;
  storage->SetFloat(animKey, position);

  const ImVec2 max(min.x + size.x, min.y + size.y);
  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddRectFilled(min, max,
                    style::mix(style::col::BgSurface, style::col::BgRaised, hover, 0.85f),
                    m.radiusMd);
  dl->AddRect(min, max, style::u32(style::col::Border), m.radiusMd, 0, m.hairline);
  const ImVec2 indicatorMin(min.x + cellWidth * position, min.y);
  const ImVec2 indicatorMax(indicatorMin.x + cellWidth, max.y);
  dl->AddRectFilled(indicatorMin, indicatorMax, style::u32(style::col::Accent), m.radiusMd);

  const float glyphSize = size.y * 0.62f;
  for (int i = 0; i < count; ++i) {
    const ImU32 color =
        i == index
            ? style::u32(style::col::OnAccent)
            : style::mix(style::col::TextDim, style::col::Text,
                         i == hoveredCell ? 1.0f : 0.0f);
    icons::draw(dl, glyphs[i],
                ImVec2(min.x + cellWidth * (static_cast<float>(i) + 0.5f),
                       min.y + size.y * 0.5f),
                glyphSize, color);
  }
  if (hoveredCell >= 0 && tooltips && tooltips[hoveredCell] &&
      tooltips[hoveredCell][0] != '\0') {
    ImGui::SetTooltip("%s", tooltips[hoveredCell]);
  }
  return changed;
}

bool toggle(const char* id, const char* label, bool& value, const char* tooltip) {
  const style::Metrics& m = style::metrics();
  const float height = ImGui::GetFrameHeight();
  const float trackWidth = height * 1.8f;
  const ImVec2 labelSize = ImGui::CalcTextSize(label);
  const ImVec2 size(trackWidth + m.gap + labelSize.x, height);
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::InvisibleButton(id, size);
  const bool clicked = ImGui::IsItemClicked();
  if (clicked) value = !value;
  const bool hovered = ImGui::IsItemHovered();
  const float hover = hoverT(ImGui::GetItemID(), hovered);

  ImGuiStorage* storage = ImGui::GetStateStorage();
  const ImGuiID animKey = ImHashStr("##toggle_position", 0, ImGui::GetItemID());
  float position = storage->GetFloat(animKey, value ? 1.0f : 0.0f);
  const float blend = std::clamp(m.animSpeed * ImGui::GetIO().DeltaTime, 0.0f, 1.0f);
  position += ((value ? 1.0f : 0.0f) - position) * blend;
  storage->SetFloat(animKey, position);

  const float trackHeight = height * 0.62f;
  const float trackY = min.y + (height - trackHeight) * 0.5f;
  const ImVec2 trackMin(min.x, trackY);
  const ImVec2 trackMax(min.x + trackWidth, trackY + trackHeight);
  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddRectFilled(trackMin, trackMax,
                    style::mix(style::col::BgRaised, style::col::Accent, position),
                    m.radiusLg);
  dl->AddRect(trackMin, trackMax,
              style::mix(style::col::Border, style::col::AccentActive, position),
              m.radiusLg, 0, m.hairline);

  const float knobDiameter = trackHeight - m.hairline * 2.0f;
  const float knobMinX = trackMin.x + m.hairline;
  const float knobMaxX = trackMax.x - m.hairline - knobDiameter;
  const float knobX = knobMinX + (knobMaxX - knobMinX) * position;
  dl->AddCircleFilled(ImVec2(knobX + knobDiameter * 0.5f, trackY + trackHeight * 0.5f),
                      knobDiameter * 0.5f,
                      style::mix(style::col::TextDim, style::col::OnAccent, position));
  dl->AddText(ImVec2(min.x + trackWidth + m.gap,
                     min.y + (height - labelSize.y) * 0.5f),
              style::mix(style::col::TextDim, style::col::Text, hover), label);

  if (hovered && tooltip) ImGui::SetTooltip("%s", tooltip);
  return clicked;
}

bool chip(const char* label, bool selected, ImVec4 accent) {
  const style::Metrics& m = style::metrics();
  const float height = ImGui::GetFontSize() + m.gap * 0.8f;
  const ImVec2 textSize = ImGui::CalcTextSize(label);
  const ImVec2 size(textSize.x + m.gap * 2.0f, height);
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::InvisibleButton(label, size);
  const bool clicked = ImGui::IsItemClicked();
  const float hover = hoverT(ImGui::GetItemID(), ImGui::IsItemHovered());
  const ImVec2 max(min.x + size.x, min.y + size.y);

  ImDrawList* dl = ImGui::GetWindowDrawList();
  if (selected) {
    dl->AddRectFilled(min, max, style::u32(accent, 0.22f), m.radiusLg);
    dl->AddRect(min, max, style::u32(accent), m.radiusLg, 0, m.hairline);
    dl->AddText(ImVec2(min.x + m.gap, min.y + (height - textSize.y) * 0.5f),
                style::u32(style::col::Text), label);
  } else {
    dl->AddRectFilled(min, max,
                      style::mix(style::col::BgSurface, style::col::BgRaised, hover),
                      m.radiusLg);
    dl->AddRect(min, max,
                style::mix(style::col::Border, style::col::BorderStrong, hover),
                m.radiusLg, 0, m.hairline);
    dl->AddText(ImVec2(min.x + m.gap, min.y + (height - textSize.y) * 0.5f),
                style::mix(style::col::TextDim, style::col::Text, hover), label);
  }
  return clicked;
}

float actionButtonWidth(icons::Icon icon, const char* label, bool primary) {
  const style::Metrics& m = style::metrics();
  const bool hasIcon = icon != icons::Icon::None;
  const bool semibold = primary ? style::pushFont(style::fonts::semibold()) : false;
  const float textWidth = ImGui::CalcTextSize(label).x;
  style::popFont(semibold);
  return (hasIcon ? m.iconSize + m.gap * 0.6f : 0.0f) + textWidth + m.gap * 2.0f;
}

bool actionButton(const char* id, icons::Icon icon, const char* label, ImVec2 size,
                  bool primary, const char* tooltip) {
  const style::Metrics& m = style::metrics();
  const float glyphSize = m.iconSize;
  const bool hasIcon = icon != icons::Icon::None;
  const float pairGap = hasIcon ? m.gap * 0.6f : 0.0f;
  const bool semibold = primary ? style::pushFont(style::fonts::semibold()) : false;
  const ImVec2 textSize = ImGui::CalcTextSize(label);
  if (size.x <= 0.0f) {
    size.x = (hasIcon ? glyphSize : 0.0f) + pairGap + textSize.x + m.gap * 2.0f;
  }
  if (size.y <= 0.0f) size.y = ImGui::GetFrameHeight();

  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::InvisibleButton(id, size);
  const bool clicked = ImGui::IsItemClicked();
  const bool hovered = ImGui::IsItemHovered();
  const float hover = hoverT(ImGui::GetItemID(), hovered);
  const ImVec2 max(min.x + size.x, min.y + size.y);
  ImDrawList* dl = ImGui::GetWindowDrawList();
  if (primary) {
    dl->AddRectFilled(min, max,
                      style::mix(style::col::Accent, style::col::AccentHover, hover),
                      m.radiusMd);
    dl->AddRect(min, max, style::u32(style::col::AccentActive), m.radiusMd, 0,
                m.hairline);
  } else {
    dl->AddRectFilled(min, max,
                      style::mix(style::col::BgSurface, style::col::BgRaised, hover,
                                 0.35f + 0.65f * hover),
                      m.radiusMd);
    dl->AddRect(min, max,
                style::mix(style::col::Border, style::col::BorderStrong, hover),
                m.radiusMd, 0, m.hairline);
  }

  const float pairWidth = (hasIcon ? glyphSize : 0.0f) + pairGap + textSize.x;
  float x = min.x + (size.x - pairWidth) * 0.5f;
  if (hasIcon) {
    icons::draw(dl, icon, ImVec2(x + glyphSize * 0.5f, min.y + size.y * 0.5f),
                glyphSize,
                primary
                    ? style::u32(style::col::OnAccent)
                    : style::mix(style::col::TextDim, style::col::Text, hover));
    x += glyphSize + pairGap;
  }
  dl->AddText(ImVec2(x, min.y + (size.y - textSize.y) * 0.5f),
              primary
                  ? style::u32(style::col::OnAccent)
                  : style::mix(style::col::TextDim, style::col::Text, hover),
              label);
  style::popFont(semibold);

  if (hovered && tooltip) ImGui::SetTooltip("%s", tooltip);
  return clicked;
}

bool glyphSlider(const char* id, icons::Icon icon, const char* label, float& value,
                 float minValue, float maxValue, const char* format,
                 const char* tooltip) {
  const style::Metrics& m = style::metrics();
  const float frameHeight = ImGui::GetFrameHeight();
  const ImVec2 glyphMin = ImGui::GetCursorScreenPos();
  ImGui::Dummy(ImVec2(frameHeight, frameHeight));
  const bool glyphHovered = ImGui::IsItemHovered();
  icons::draw(ImGui::GetWindowDrawList(), icon,
              ImVec2(glyphMin.x + frameHeight * 0.5f, glyphMin.y + frameHeight * 0.5f),
              std::min(m.iconSize, frameHeight * 0.72f), style::u32(style::col::TextDim));

  ImGui::SameLine(0.0f, m.gap * 0.5f);
  ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x);
  ImGui::PushID(id);
  ImGui::PushStyleColor(ImGuiCol_FrameBg, style::col::BgSurface);
  ImGui::PushStyleColor(ImGuiCol_SliderGrab, style::col::Accent);
  const bool changed = ImGui::SliderFloat("##slider", &value, minValue, maxValue, format);
  ImGui::PopStyleColor(2);
  ImGui::PopID();

  const bool sliderHovered = ImGui::IsItemHovered();
  const ImVec2 frameMin = ImGui::GetItemRectMin();
  const ImVec2 frameMax = ImGui::GetItemRectMax();
  const float labelSize = ImGui::GetFontSize() * 0.82f;
  ImGui::GetWindowDrawList()->AddText(
      nullptr, labelSize,
      ImVec2(frameMin.x + ImGui::GetStyle().FramePadding.x,
             frameMin.y + (frameMax.y - frameMin.y - labelSize) * 0.5f),
      style::u32(style::col::TextFaint), label);

  if ((glyphHovered || sliderHovered) && tooltip) ImGui::SetTooltip("%s", tooltip);
  return changed;
}

void metric(const char* caption, const char* value, const char* unit, const char* delta,
            ImVec4 accent) {
  const style::Metrics& m = style::metrics();
  const float fs = ImGui::GetFontSize();
  const float small = fs * 0.82f;
  const float lineGap = m.gap * 0.35f;
  const float height = small + lineGap + fs;
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::Dummy(ImVec2(ImGui::GetContentRegionAvail().x, height));

  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddText(nullptr, small, min, style::u32(style::col::TextFaint), caption);
  const float valueY = min.y + small + lineGap;
  const bool mono = style::pushFont(style::fonts::mono());
  const ImVec2 valueSize = ImGui::CalcTextSize(value);
  dl->AddText(ImVec2(min.x, valueY), style::u32(accent), value);
  style::popFont(mono);

  float x = min.x + valueSize.x;
  if (unit && unit[0] != '\0') {
    dl->AddText(nullptr, small, ImVec2(x, valueY + fs - small),
                style::u32(style::col::TextDim), unit);
    x += ImGui::CalcTextSize(unit).x * (small / fs);
  }
  if (delta && delta[0] != '\0') {
    x += m.gap * 0.6f;
    const ImVec4 deltaColor =
        delta[0] == '+' ? style::col::Success
                        : (delta[0] == '-' ? style::col::Danger : style::col::TextDim);
    dl->AddText(nullptr, small, ImVec2(x, valueY + fs - small),
                style::u32(deltaColor), delta);
  }
}

void keyValue(const char* key, const char* value, ImVec4 accent) {
  const style::Metrics& m = style::metrics();
  const float width = ImGui::GetContentRegionAvail().x;
  const float height = ImGui::GetTextLineHeightWithSpacing();
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::Dummy(ImVec2(width, height));
  const bool hovered = ImGui::IsItemHovered();

  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddText(min, style::u32(style::col::TextDim), key);

  const ImVec2 keySize = ImGui::CalcTextSize(key);
  const bool mono = style::pushFont(style::fonts::mono());
  const ImVec2 valueSize = ImGui::CalcTextSize(value);
  // A long value would otherwise be drawn straight through the key. Clip it to
  // the space the key leaves and let hover reveal the rest, because silently
  // overlapping two strings is worse than truncating one.
  const float room = std::max(width - keySize.x - m.gap, 0.0f);
  const float x = min.x + width - std::min(valueSize.x, room);
  const ImVec4 clip(min.x + keySize.x + m.gap, min.y, min.x + width, min.y + height);
  dl->AddText(nullptr, 0.0f, ImVec2(x, min.y), style::u32(accent), value, nullptr, 0.0f,
              &clip);
  style::popFont(mono);

  if (hovered && valueSize.x > room) ImGui::SetTooltip("%s: %s", key, value);
}

void progressRow(const char* caption, float fraction, const char* value, ImVec4 accent) {
  const style::Metrics& m = style::metrics();
  const float width = ImGui::GetContentRegionAvail().x;
  const float textHeight = ImGui::GetFontSize();
  const float trackHeight = m.hairline * 3.0f;
  const float trackGap = m.gap * 0.4f;
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::Dummy(ImVec2(width, textHeight + trackGap + trackHeight));

  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddText(min, style::u32(style::col::TextDim), caption);
  const bool mono = style::pushFont(style::fonts::mono());
  const ImVec2 valueSize = ImGui::CalcTextSize(value);
  dl->AddText(ImVec2(min.x + width - valueSize.x, min.y), style::u32(style::col::Text),
              value);
  style::popFont(mono);

  const float trackY = min.y + textHeight + trackGap;
  const ImVec2 trackMin(min.x, trackY);
  const ImVec2 trackMax(min.x + width, trackY + trackHeight);
  dl->AddRectFilled(trackMin, trackMax, style::u32(style::col::BgRaised), m.radiusSm);
  const float fill = std::clamp(fraction, 0.0f, 1.0f);
  if (fill > 0.0f) {
    dl->AddRectFilled(trackMin, ImVec2(min.x + width * fill, trackMax.y),
                      style::u32(accent), m.radiusSm);
  }
}

void helpMarker(const char* text) {
  const float size = ImGui::GetFontSize();
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::PushID(text);
  ImGui::InvisibleButton("##help_marker", ImVec2(size, size));
  const bool hovered = ImGui::IsItemHovered();
  const float hover = hoverT(ImGui::GetItemID(), hovered);
  icons::draw(ImGui::GetWindowDrawList(), icons::Icon::Info,
              ImVec2(min.x + size * 0.5f, min.y + size * 0.5f), size,
              style::mix(style::col::TextFaint, style::col::Text, hover));
  ImGui::PopID();

  if (hovered && text) {
    ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
    ImGui::SetTooltip("%s", text);
    ImGui::PopTextWrapPos();
  }
}

void emptyState(icons::Icon icon, const char* headline, const char* body) {
  const style::Metrics& m = style::metrics();
  const float fs = ImGui::GetFontSize();
  const float glyphSize = fs * 2.6f;
  const ImVec2 min = ImGui::GetCursorScreenPos();
  const float width = std::max(ImGui::GetContentRegionAvail().x, fs);
  const float bodyWidth = width * 0.7f;

  const bool semibold = style::pushFont(style::fonts::semibold());
  const ImVec2 headlineSize = ImGui::CalcTextSize(headline);
  style::popFont(semibold);
  const ImVec2 bodySize = ImGui::CalcTextSize(body, nullptr, false, bodyWidth);
  const float contentHeight =
      glyphSize + m.gap + headlineSize.y + m.gap * 0.7f + bodySize.y;

  // The block reserves its own content height before centring in it. Measuring
  // only the available region collapses to a single line inside an auto-height
  // card, which clipped everything below the glyph.
  const float height =
      std::max(ImGui::GetContentRegionAvail().y, contentHeight + m.gap * 2.0f);
  ImGui::Dummy(ImVec2(width, height));

  const float top = min.y + std::max((height - contentHeight) * 0.5f, 0.0f);
  const float centreX = min.x + width * 0.5f;

  ImDrawList* dl = ImGui::GetWindowDrawList();
  icons::draw(dl, icon, ImVec2(centreX, top + glyphSize * 0.5f), glyphSize,
              style::u32(style::col::TextFaint));
  const float headlineY = top + glyphSize + m.gap;
  dl->AddText(style::fonts::semibold(), fs,
              ImVec2(centreX - headlineSize.x * 0.5f, headlineY),
              style::u32(style::col::TextDim), headline);
  const ImVec2 bodyPos(centreX - bodyWidth * 0.5f,
                       headlineY + headlineSize.y + m.gap * 0.7f);
  dl->AddText(nullptr, fs, bodyPos, style::u32(style::col::TextFaint), body, nullptr,
              bodyWidth);
}

void notice(icons::Icon icon, const char* text, ImVec4 accent) {
  const style::Metrics& m = style::metrics();
  const float fs = ImGui::GetFontSize();
  const float width = ImGui::GetContentRegionAvail().x;
  const float padding = m.gap;
  const float glyphSize = m.iconSize;
  const float textWidth =
      std::max(width - padding * 2.0f - glyphSize - m.gap, fs);
  const ImVec2 textSize = ImGui::CalcTextSize(text, nullptr, false, textWidth);
  const float height = std::max(glyphSize, textSize.y) + padding * 2.0f;
  const ImVec2 min = ImGui::GetCursorScreenPos();
  const ImVec2 max(min.x + width, min.y + height);
  ImGui::Dummy(ImVec2(width, height));

  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddRectFilled(min, max, style::u32(accent, 0.12f), m.radiusMd);
  dl->AddRect(min, max, style::u32(accent, 0.40f), m.radiusMd, 0, m.hairline);
  const ImVec2 glyphCentre(min.x + padding + glyphSize * 0.5f,
                           min.y + height * 0.5f);
  icons::draw(dl, icon, glyphCentre, glyphSize, style::u32(accent));
  const ImVec2 textPos(min.x + padding + glyphSize + m.gap,
                       min.y + (height - textSize.y) * 0.5f);
  dl->AddText(nullptr, fs, textPos, style::u32(style::col::Text), text, nullptr,
              textWidth);
}

}  // namespace chemcad::ui::widgets

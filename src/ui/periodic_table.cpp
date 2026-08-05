#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "imgui.h"

#include "ui/element_data.hpp"
#include "ui/icons.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"
#include "ui/widgets.hpp"

namespace chemcad::ui {
namespace {

std::string lower(std::string text) {
  std::transform(text.begin(), text.end(), text.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return text;
}

bool matches(const ElementData& element, const std::string& query) {
  if (query.empty()) return true;
  const std::string symbol = lower(element.symbol);
  const std::string name = lower(element.name);
  return symbol.starts_with(query) || name.find(query) != std::string::npos;
}

// Category hues. These read as small stripes and tinted symbols on dark
// tiles, so they stay saturated.
ImVec4 categoryColor(const std::string& category) {
  if (category == "alkali-metal") return {0.90f, 0.42f, 0.36f, 1.0f};
  if (category == "alkaline-earth") return {0.90f, 0.63f, 0.34f, 1.0f};
  if (category == "transition-metal") return {0.51f, 0.63f, 0.84f, 1.0f};
  if (category == "post-transition") return {0.45f, 0.73f, 0.75f, 1.0f};
  if (category == "metalloid") return {0.49f, 0.76f, 0.53f, 1.0f};
  if (category == "nonmetal") return {0.42f, 0.71f, 0.88f, 1.0f};
  if (category == "halogen") return {0.73f, 0.55f, 0.86f, 1.0f};
  if (category == "noble-gas") return {0.80f, 0.50f, 0.74f, 1.0f};
  if (category == "lanthanide") return {0.82f, 0.60f, 0.47f, 1.0f};
  if (category == "actinide") return {0.72f, 0.53f, 0.57f, 1.0f};
  return {0.55f, 0.60f, 0.65f, 1.0f};
}

void selectElement(AppState& st, const ElementData& element) {
  st.currentElement = element.z;
  st.tool = Tool::Atom;
  st.statusMessage = "Element: " + element.name + " (" + element.symbol + ")";
}

void elementTooltip(const ElementData& element) {
  ImGui::BeginTooltip();
  const bool pushed = style::pushFont(style::fonts::semibold());
  ImGui::Text("%s", element.name.c_str());
  style::popFont(pushed);
  ImGui::Separator();
  ImGui::Text("Symbol: %s", element.symbol.c_str());
  ImGui::Text("Atomic number: %u", static_cast<unsigned>(element.z));
  ImGui::Text("Standard atomic weight: %.6g", element.mass);
  ImGui::TextColored(categoryColor(element.category), "Category: %s",
                     element.category.c_str());
  ImGui::EndTooltip();
}

// One element tile: behaviour via an invisible button, everything else drawn.
void drawTile(ImDrawList* dl, ImVec2 min, float cell, const ElementData& element,
              bool selected, bool match, float hoverT) {
  const style::Metrics& m = style::metrics();
  const float alpha = match ? 1.0f : 0.20f;
  const float radius = std::min(m.radiusMd, cell * 0.22f);
  const ImVec2 max(min.x + cell, min.y + cell);
  const ImVec4 cat = categoryColor(element.category);

  const ImU32 fill = style::mix(style::col::BgSurface, cat, 0.05f + 0.13f * hoverT, alpha);
  dl->AddRectFilled(min, max, fill, radius);

  if (selected) {
    // Soft outer glow, then the solid amber frame.
    dl->AddRect(ImVec2(min.x - 2.0f, min.y - 2.0f), ImVec2(max.x + 2.0f, max.y + 2.0f),
                style::u32(style::col::Accent, 0.28f * alpha), radius + 2.0f, 0, 4.0f);
    dl->AddRect(min, max, style::u32(style::col::Accent, alpha), radius, 0,
                m.hairline + 1.0f);
  } else {
    dl->AddRect(min, max, style::mix(style::col::Border, cat, 0.25f + 0.60f * hoverT, alpha),
                radius, 0, m.hairline);
  }

  // Category stripe along the bottom edge.
  const float stripeH = std::max(2.0f, cell * 0.10f);
  const float inset = cell * 0.14f;
  dl->AddRectFilled(ImVec2(min.x + inset, max.y - inset * 0.55f - stripeH),
                    ImVec2(max.x - inset, max.y - inset * 0.55f),
                    style::u32(cat, alpha * (0.75f + 0.25f * hoverT)), stripeH * 0.5f);

  ImFont* font = ImGui::GetFont();
  const float fs = ImGui::GetFontSize();
  const float symbolSize = std::min(fs, cell * 0.46f);
  const ImU32 symbolCol =
      selected ? style::u32(style::col::Accent, alpha)
               : style::mix(cat, style::col::Text, 0.42f + 0.25f * hoverT, alpha);
  const ImVec2 symbolExtent = font->CalcTextSizeA(symbolSize, 10000.0f, 0.0f,
                                                  element.symbol.c_str());
  dl->AddText(font, symbolSize,
              ImVec2(min.x + (cell - symbolExtent.x) * 0.5f,
                     min.y + (cell - symbolExtent.y) * 0.5f - cell * 0.04f),
              symbolCol, element.symbol.c_str());

  // Atomic number, top-left.
  if (cell >= fs * 1.25f) {
    char zbuf[8];
    std::snprintf(zbuf, sizeof(zbuf), "%u", static_cast<unsigned>(element.z));
    const float zSize = std::min(fs * 0.62f, cell * 0.26f);
    dl->AddText(font, zSize, ImVec2(min.x + cell * 0.10f, min.y + cell * 0.07f),
                style::u32(style::col::TextFaint, alpha), zbuf);
  }
}

// Ln/An spacer tile: dashed outline, dim label.
void drawPlaceholder(ImDrawList* dl, ImVec2 min, float cell, const char* label) {
  const style::Metrics& m = style::metrics();
  const ImVec2 max(min.x + cell, min.y + cell);
  const ImU32 col = style::u32(style::col::TextFaint, 0.55f);
  // Dashed square: four dashed edges.
  const float dash = cell * 0.14f;
  const float gap = cell * 0.10f;
  auto dashed = [&](ImVec2 a, ImVec2 b) {
    const float dx = b.x - a.x, dy = b.y - a.y;
    const float len = std::sqrt(dx * dx + dy * dy);
    if (len < 0.01f) return;
    const float ux = dx / len, uy = dy / len;
    for (float s = 0.0f; s < len; s += dash + gap) {
      const float e = std::min(s + dash, len);
      dl->AddLine(ImVec2(a.x + ux * s, a.y + uy * s), ImVec2(a.x + ux * e, a.y + uy * e),
                  col, m.hairline);
    }
  };
  ImVec2 cmin(min.x + 1.0f, min.y + 1.0f), cmax(max.x - 1.0f, max.y - 1.0f);
  dashed(cmin, ImVec2(cmax.x, cmin.y));
  dashed(ImVec2(cmax.x, cmin.y), cmax);
  dashed(cmax, ImVec2(cmin.x, cmax.y));
  dashed(ImVec2(cmin.x, cmax.y), cmin);

  ImFont* font = ImGui::GetFont();
  const float size = std::min(ImGui::GetFontSize() * 0.8f, cell * 0.4f);
  const ImVec2 extent = font->CalcTextSizeA(size, 10000.0f, 0.0f, label);
  dl->AddText(font, size,
              ImVec2(min.x + (cell - extent.x) * 0.5f, min.y + (cell - extent.y) * 0.5f),
              col, label);
}

}  // namespace

void drawPeriodicTable(AppState& st) {
  static char search[96]{};
  static bool loadErrorSurfaced = false;

  const auto& elements = elementTable();
  if (!elementTableLoadError().empty() && !loadErrorSurfaced) {
    st.statusMessage = elementTableLoadError();
    loadErrorSurfaced = true;
  }

  ImGui::SetNextItemWidth(-1.0f);
  const bool enter = ImGui::InputTextWithHint(
      "##element_search", "Search element...", search, sizeof(search),
      ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_AutoSelectAll);
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "Type a name or symbol, Enter picks a unique match");
  // Magnifier decoration inside the frame's right edge.
  {
    const ImVec2 rmin = ImGui::GetItemRectMin();
    const ImVec2 rmax = ImGui::GetItemRectMax();
    const float g = (rmax.y - rmin.y) * 0.44f;
    icons::draw(ImGui::GetWindowDrawList(), icons::Icon::Search,
                ImVec2(rmax.x - g * 0.9f, (rmin.y + rmax.y) * 0.5f), g,
                style::u32(style::col::TextFaint));
  }

  const std::string query = lower(search);
  std::vector<const ElementData*> matched;
  if (!query.empty()) {
    for (const auto& element : elements) {
      if (matches(element, query)) matched.push_back(&element);
    }
  }
  if (enter && matched.size() == 1) selectElement(st, *matched.front());

  // All 18 groups must be reachable without horizontal scrolling: C, N and O
  // are the most-used tiles in the app, and they sit in the right-hand
  // groups. Size the cells from the space actually available. The content
  // width is exactly 18 cells + 17 gaps — an 18th gap would always overflow.
  const float gap = std::max(2.0f, ImGui::GetStyle().ItemSpacing.x * 0.35f);
  const float avail = ImGui::GetContentRegionAvail().x;
  const float preferred = std::max(30.0f, ImGui::GetFontSize() * 2.05f);
  const float fitted = (avail - 17.0f * gap) / 18.0f;
  const float cell = std::clamp(fitted, 18.0f, preferred);
  const float step = cell + gap;
  const float width = 18.0f * cell + 17.0f * gap;
  const float height = 10.0f * cell + 9.0f * gap;

  ImGui::SetNextWindowContentSize(ImVec2(width, height));
  if (!ImGui::BeginChild("##ptable_grid", ImVec2(0, 0), ImGuiChildFlags_None,
                         ImGuiWindowFlags_HorizontalScrollbar)) {
    ImGui::EndChild();
    return;
  }

  ImDrawList* dl = ImGui::GetWindowDrawList();
  const ImVec2 origin = ImGui::GetCursorPos();
  auto placeholderAt = [&](int group, int period, const char* label, const char* tooltip) {
    const ImVec2 pos(origin.x + (group - 1) * step, origin.y + (period - 1) * step);
    ImGui::SetCursorPos(pos);
    ImGui::InvisibleButton("##placeholder", ImVec2(cell, cell));
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", tooltip);
    const ImVec2 screen = ImGui::GetItemRectMin();
    drawPlaceholder(dl, screen, cell, label);
  };
  placeholderAt(3, 6, "Ln", "Lanthanides (57-71), shown below");
  placeholderAt(3, 7, "An", "Actinides (89-103), shown below");

  for (const auto& element : elements) {
    int column = element.group;
    int row = element.period;
    if (element.z >= 57 && element.z <= 71) {
      column = 3 + static_cast<int>(element.z) - 57;
      row = 9;
    } else if (element.z >= 89 && element.z <= 103) {
      column = 3 + static_cast<int>(element.z) - 89;
      row = 10;
    } else if (column == 0) {
      continue;
    }

    const bool match = matches(element, query);
    const ImVec2 pos(origin.x + (column - 1) * step, origin.y + (row - 1) * step);
    ImGui::SetCursorPos(pos);
    ImGui::PushID(static_cast<int>(element.z));
    ImGui::InvisibleButton("##tile", ImVec2(cell, cell));
    const bool clicked = ImGui::IsItemClicked();
    const bool hovered = ImGui::IsItemHovered();
    const float t = widgets::hoverT(ImGui::GetItemID(), hovered);
    const ImVec2 screen = ImGui::GetItemRectMin();

    const bool selected = st.currentElement == element.z;
    drawTile(dl, screen, cell, element, selected, match, t);

    if (clicked) selectElement(st, element);
    if (hovered) elementTooltip(element);
    ImGui::PopID();
  }

  ImGui::EndChild();
}

}  // namespace chemcad::ui

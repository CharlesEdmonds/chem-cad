#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "imgui.h"

#include "ui/atom_view.hpp"
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

const char* categoryDisplayName(const std::string& category) {
  if (category == "alkali-metal") return "Alkali metal";
  if (category == "alkaline-earth") return "Alkaline earth metal";
  if (category == "transition-metal") return "Transition metal";
  if (category == "post-transition") return "Post-transition metal";
  if (category == "metalloid") return "Metalloid";
  if (category == "nonmetal") return "Nonmetal";
  if (category == "halogen") return "Halogen";
  if (category == "noble-gas") return "Noble gas";
  if (category == "lanthanide") return "Lanthanide";
  if (category == "actinide") return "Actinide";
  return "Unknown category";
}

void selectElement(AppState& st, const ElementData& element) {
  st.currentElement = element.z;
  st.tool = Tool::Atom;
  st.statusMessage = "Element: " + element.name + " (" + element.symbol + ")";
}

void elementTooltip(const ElementData& element) {
  const style::Metrics& metrics = style::metrics();
  const float fs = ImGui::GetFontSize();
  const ImVec4 category = categoryColor(element.category);
  const ImVec2 cardSize(fs * 21.0f, fs * 33.5f);

  // Keep the hover card close to the pointer without letting it run off the
  // current viewport. BeginTooltip remains backend-agnostic for headless use.
  const ImGuiViewport* viewport = ImGui::GetMainViewport();
  const ImVec2 mouse = ImGui::GetMousePos();
  ImVec2 cardPos(mouse.x + metrics.gap * 2.0f, mouse.y + metrics.gap * 1.5f);
  const float viewportRight = viewport->WorkPos.x + viewport->WorkSize.x;
  const float viewportBottom = viewport->WorkPos.y + viewport->WorkSize.y;
  if (cardPos.x + cardSize.x > viewportRight - metrics.gap)
    cardPos.x = mouse.x - cardSize.x - metrics.gap * 2.0f;
  if (cardPos.y + cardSize.y > viewportBottom - metrics.gap)
    cardPos.y = mouse.y - cardSize.y - metrics.gap * 1.5f;
  cardPos.x = std::max(viewport->WorkPos.x + metrics.gap, cardPos.x);
  cardPos.y = std::max(viewport->WorkPos.y + metrics.gap, cardPos.y);

  ImGui::SetNextWindowPos(cardPos, ImGuiCond_Always);
  ImGui::SetNextWindowSize(cardSize, ImGuiCond_Always);
  ImGui::SetNextWindowBgAlpha(1.0f);
  ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding,
                      ImVec2(metrics.gap * 1.25f, metrics.gap * 1.15f));
  ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, metrics.radiusLg);
  ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, metrics.hairline);
  ImGui::PushStyleColor(ImGuiCol_PopupBg, style::col::BgRaised);
  ImGui::PushStyleColor(ImGuiCol_Border, style::col::BorderStrong);
  ImGui::BeginTooltip();

  ImDrawList* dl = ImGui::GetWindowDrawList();
  const ImVec2 windowMin = ImGui::GetWindowPos();
  const ImVec2 windowMax(windowMin.x + ImGui::GetWindowSize().x,
                         windowMin.y + ImGui::GetWindowSize().y);
  const float stripH = std::max(metrics.hairline * 3.0f, metrics.gap * 0.65f);
  dl->AddRectFilled(windowMin, ImVec2(windowMax.x, windowMin.y + stripH),
                    style::u32(category), metrics.radiusLg,
                    ImDrawFlags_RoundCornersTop);
  ImGui::Dummy(ImVec2(0.0f, stripH * 0.45f));

  // Large symbol tile and identity block.
  const float tile = fs * 3.35f;
  const ImVec2 tileMin = ImGui::GetCursorScreenPos();
  ImGui::Dummy(ImVec2(tile, tile));
  const ImVec2 tileMax(tileMin.x + tile, tileMin.y + tile);
  dl->AddRectFilled(tileMin, tileMax,
                    style::mix(style::col::BgSurface, category, 0.20f),
                    metrics.radiusMd);
  dl->AddRect(tileMin, tileMax, style::u32(category, 0.85f),
              metrics.radiusMd, 0, metrics.hairline);
  ImFont* symbolFont =
      style::fonts::semibold() ? style::fonts::semibold() : ImGui::GetFont();
  const float symbolSize = fs * 1.85f;
  const ImVec2 symbolExtent =
      symbolFont->CalcTextSizeA(symbolSize, 10000.0f, 0.0f,
                                element.symbol.c_str());
  dl->AddText(symbolFont, symbolSize,
              ImVec2(tileMin.x + (tile - symbolExtent.x) * 0.5f,
                     tileMin.y + (tile - symbolExtent.y) * 0.5f),
              style::mix(category, style::col::Text, 0.42f),
              element.symbol.c_str());

  ImGui::SameLine(0.0f, metrics.gap);
  ImGui::BeginGroup();
  const bool nameFont = style::pushFont(style::fonts::semibold());
  ImGui::TextUnformatted(element.name.c_str());
  style::popFont(nameFont);
  const bool zFont = style::pushFont(style::fonts::mono());
  ImGui::TextColored(style::col::TextDim, "Z = %u",
                     static_cast<unsigned>(element.z));
  style::popFont(zFont);
  widgets::badge(categoryDisplayName(element.category), category);
  ImGui::EndGroup();

  ImGui::Spacing();
  ImGui::Separator();
  ImGui::Spacing();

  // Animated atom stage: orbital clouds by default, Bohr shells as backup.
  static bool bohrMode = false;
  const float stage = ImGui::GetContentRegionAvail().x;
  const ImVec2 stageMin = ImGui::GetCursorScreenPos();
  ImGui::Dummy(ImVec2(stage, stage));
  const ImVec2 stageMax(stageMin.x + stage, stageMin.y + stage);
  drawAtomModel(element, stageMin, stageMax, bohrMode);
  {
    // Mode toggle rides the stage's top-right corner.
    const char* modeLabel = bohrMode ? "BOHR" : "ORBITAL";
    ImGui::SetCursorScreenPos(ImVec2(stageMax.x - ImGui::CalcTextSize(modeLabel).x -
                                         metrics.gap * 2.2f,
                                     stageMin.y + metrics.gap * 0.7f));
    if (ImGui::SmallButton(modeLabel)) bohrMode = !bohrMode;
    if (ImGui::IsItemHovered())
      ImGui::SetTooltip("Switch between orbital-cloud and Bohr-shell representations");
    ImGui::SetCursorScreenPos(ImVec2(stageMin.x, stageMax.y));
  }
  {
    const int neutrons =
        std::max(0, static_cast<int>(std::lround(element.mass)) - element.z);
    char buf[96];
    std::snprintf(buf, sizeof(buf), "%d p+   %d n0   %d e-", element.z, neutrons,
                  element.z);
    const ImVec2 labelSize = ImGui::CalcTextSize(buf);
    dl->AddText(style::fonts::mono(), fs * 0.9f,
                ImVec2((stageMin.x + stageMax.x - labelSize.x) * 0.5f,
                       stageMax.y - labelSize.y - metrics.gap * 0.6f),
                style::u32(style::col::TextDim), buf);
  }

  ImGui::Spacing();

  // ElementData supplies number, mass, group and period; the configuration is
  // derived from Z. Keep the values aligned and monospaced.
  constexpr ImGuiTableFlags tableFlags =
      ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_NoSavedSettings |
      ImGuiTableFlags_RowBg;
  if (ImGui::BeginTable("##element_stats", 2, tableFlags)) {
    ImGui::TableSetupColumn("##label", ImGuiTableColumnFlags_WidthStretch, 0.62f);
    ImGui::TableSetupColumn("##value", ImGuiTableColumnFlags_WidthStretch, 0.38f);

    const auto row = [](const char* label, const char* value) {
      ImGui::TableNextRow();
      ImGui::TableNextColumn();
      ImGui::TextDisabled("%s", label);
      ImGui::TableNextColumn();
      const bool mono = style::pushFont(style::fonts::mono());
      ImGui::TextUnformatted(value);
      style::popFont(mono);
    };

    char atomicNumber[8];
    std::snprintf(atomicNumber, sizeof(atomicNumber), "%u",
                  static_cast<unsigned>(element.z));
    char atomicMass[32];
    std::snprintf(atomicMass, sizeof(atomicMass), "%.6g", element.mass);
    char group[16];
    std::snprintf(group, sizeof(group), "%d", element.group);
    char period[16];
    std::snprintf(period, sizeof(period), "%d", element.period);
    row("Atomic number", atomicNumber);
    row("Standard atomic weight", atomicMass);
    row("Group", element.group == 0 ? "f-block" : group);
    row("Period", period);
    row("Configuration", elementConfigString(element).c_str());
    ImGui::EndTable();
  }

  ImGui::Spacing();
  ImGui::TextColored(style::col::Accent, "CLICK");
  ImGui::SameLine(0.0f, metrics.gap * 0.5f);
  ImGui::TextDisabled("Select element and activate the Atom tool");

  ImGui::EndTooltip();
  ImGui::PopStyleColor(2);
  ImGui::PopStyleVar(3);
}

// One element tile: behaviour via an invisible button, everything else drawn.
void drawTile(ImDrawList* dl, ImVec2 min, ImVec2 size, const ElementData& element,
              bool selected, bool match, float hoverT) {
  const style::Metrics& m = style::metrics();
  const float alpha = match ? 1.0f : 0.20f;
  const float cell = std::min(size.x, size.y);
  const float radius = std::min(m.radiusMd, cell * 0.22f);
  const ImVec2 max(min.x + size.x, min.y + size.y);
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
  const float symbolSize = std::min(fs, size.x * 0.56f);
  const ImU32 symbolCol =
      selected ? style::u32(style::col::Accent, alpha)
               : style::mix(cat, style::col::Text, 0.42f + 0.25f * hoverT, alpha);
  const ImVec2 symbolExtent = font->CalcTextSizeA(symbolSize, 10000.0f, 0.0f,
                                                  element.symbol.c_str());
  dl->AddText(font, symbolSize,
              ImVec2(min.x + (size.x - symbolExtent.x) * 0.5f,
                     min.y + (size.y - symbolExtent.y) * 0.5f - size.y * 0.035f),
              symbolCol, element.symbol.c_str());

  // Atomic number, top-left.
  if (size.x >= fs * 1.25f) {
    char zbuf[8];
    std::snprintf(zbuf, sizeof(zbuf), "%u", static_cast<unsigned>(element.z));
    const float zSize = std::min(fs * 0.62f, size.x * 0.26f);
    dl->AddText(font, zSize, ImVec2(min.x + size.x * 0.10f, min.y + size.y * 0.07f),
                style::u32(style::col::TextFaint, alpha), zbuf);
  }
}

// Ln/An spacer tile: dashed outline, dim label.
void drawPlaceholder(ImDrawList* dl, ImVec2 min, ImVec2 size, const char* label) {
  const style::Metrics& m = style::metrics();
  const float cell = std::min(size.x, size.y);
  const ImVec2 max(min.x + size.x, min.y + size.y);
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
  const float fontSize = std::min(ImGui::GetFontSize() * 0.8f, cell * 0.4f);
  const ImVec2 extent = font->CalcTextSizeA(fontSize, 10000.0f, 0.0f, label);
  dl->AddText(font, fontSize,
              ImVec2(min.x + (size.x - extent.x) * 0.5f,
                     min.y + (size.y - extent.y) * 0.5f),
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

  // Reserve compact gutters for group and period labels, then fit all 18
  // groups to the remaining width. The draw-list labels never advance the
  // cursor, so every behaviour item stays inside the declared content rect.
  const ImVec2 avail = ImGui::GetContentRegionAvail();
  const float fs = ImGui::GetFontSize();
  const float gapX = std::max(1.0f, ImGui::GetStyle().ItemSpacing.x * 0.18f);
  const float gapY = std::max(2.0f, ImGui::GetStyle().ItemSpacing.y * 0.32f);
  const float labelGutter = fs * 1.35f;
  const float groupHeader = fs * 1.25f;
  const float preferredW = std::max(fs * 2.10f, 32.0f);
  const float fittedW = (avail.x - labelGutter - 17.0f * gapX) / 18.0f;
  const float cellW = std::clamp(fittedW, fs * 0.92f, preferredW);
  const float fittedH = (avail.y - groupHeader - 9.0f * gapY) / 10.0f;
  const float cellH = std::clamp(fittedH, cellW, cellW * 1.65f);
  const float stepX = cellW + gapX;
  const float stepY = cellH + gapY;
  const float width = labelGutter + 18.0f * cellW + 17.0f * gapX;
  const float height = groupHeader + 10.0f * cellH + 9.0f * gapY;

  ImGui::SetNextWindowContentSize(ImVec2(width, height));
  ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
  const bool visible = ImGui::BeginChild(
      "##ptable_grid", ImVec2(0, 0), ImGuiChildFlags_None,
      ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
  ImGui::PopStyleVar();
  if (!visible) {
    ImGui::EndChild();
    return;
  }

  ImDrawList* dl = ImGui::GetWindowDrawList();
  const ImVec2 contentOrigin = ImGui::GetCursorPos();
  const ImVec2 origin(contentOrigin.x + labelGutter, contentOrigin.y + groupHeader);
  const ImVec2 contentScreen = ImGui::GetCursorScreenPos();
  const ImVec2 originScreen(contentScreen.x + labelGutter,
                            contentScreen.y + groupHeader);

  for (int groupIndex = 1; groupIndex <= 18; ++groupIndex) {
    char label[4];
    std::snprintf(label, sizeof(label), "%d", groupIndex);
    const ImVec2 extent = ImGui::CalcTextSize(label);
    dl->AddText(ImVec2(originScreen.x + (groupIndex - 1) * stepX +
                               (cellW - extent.x) * 0.5f,
                       contentScreen.y + (groupHeader - extent.y) * 0.5f),
                style::u32(style::col::TextFaint), label);
  }
  for (int periodIndex = 1; periodIndex <= 7; ++periodIndex) {
    char label[4];
    std::snprintf(label, sizeof(label), "%d", periodIndex);
    const ImVec2 extent = ImGui::CalcTextSize(label);
    dl->AddText(ImVec2(contentScreen.x + (labelGutter - extent.x) * 0.5f,
                       originScreen.y + (periodIndex - 1) * stepY +
                           (cellH - extent.y) * 0.5f),
                style::u32(style::col::TextFaint), label);
  }
  const char* detachedLabels[] = {"Ln", "An"};
  for (int i = 0; i < 2; ++i) {
    const ImVec2 extent = ImGui::CalcTextSize(detachedLabels[i]);
    dl->AddText(ImVec2(contentScreen.x + (labelGutter - extent.x) * 0.5f,
                       originScreen.y + (8 + i) * stepY +
                           (cellH - extent.y) * 0.5f),
                style::u32(style::col::TextFaint), detachedLabels[i]);
  }
  auto placeholderAt = [&](int group, int period, const char* label, const char* tooltip) {
    const ImVec2 pos(origin.x + (group - 1) * stepX, origin.y + (period - 1) * stepY);
    ImGui::SetCursorPos(pos);
    ImGui::PushID(label);
    ImGui::InvisibleButton("##placeholder", ImVec2(cellW, cellH));
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", tooltip);
    const ImVec2 screen = ImGui::GetItemRectMin();
    drawPlaceholder(dl, screen, ImVec2(cellW, cellH), label);
    ImGui::PopID();
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
    const ImVec2 pos(origin.x + (column - 1) * stepX, origin.y + (row - 1) * stepY);
    ImGui::SetCursorPos(pos);
    ImGui::PushID(static_cast<int>(element.z));
    ImGui::InvisibleButton("##tile", ImVec2(cellW, cellH));
    const bool clicked = ImGui::IsItemClicked();
    const bool hovered = ImGui::IsItemHovered();
    const float t = widgets::hoverT(ImGui::GetItemID(), hovered);
    const ImVec2 screen = ImGui::GetItemRectMin();

    const bool selected = st.currentElement == element.z;
    drawTile(dl, screen, ImVec2(cellW, cellH), element, selected, match, t);

    if (clicked) selectElement(st, element);
    if (hovered) elementTooltip(element);
    ImGui::PopID();
  }

  ImGui::EndChild();
}

}  // namespace chemcad::ui

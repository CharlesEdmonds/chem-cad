#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <limits>
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

enum class ColourMode { Category, AtomicMass };

const char* categoryDisplayName(const std::string& category);

std::string lower(std::string text) {
  std::transform(text.begin(), text.end(), text.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return text;
}

bool matches(const ElementData& element, const std::string& query) {
  if (query.empty()) return true;
  const std::string symbol = lower(element.symbol);
  const std::string name = lower(element.name);
  const std::string category = lower(categoryDisplayName(element.category));
  const std::string number = std::to_string(element.z);
  return symbol.starts_with(query) || name.find(query) != std::string::npos ||
         number.starts_with(query) || category.find(query) != std::string::npos;
}

bool exactMatch(const ElementData& element, const std::string& query) {
  return lower(element.symbol) == query || lower(element.name) == query ||
         std::to_string(element.z) == query;
}

// Category identity is deliberately constrained to the shared palette so the
// table remains coherent with every other panel under alternate UI scaling.
ImVec4 categoryColor(const std::string& category) {
  if (category == "alkali-metal") return style::col::Danger;
  if (category == "alkaline-earth") return style::col::Accent;
  if (category == "transition-metal") return style::col::Violet;
  if (category == "post-transition") return style::col::Teal;
  if (category == "metalloid") return style::col::Success;
  if (category == "nonmetal") return style::col::Teal;
  if (category == "halogen") return style::col::Violet;
  if (category == "noble-gas") return style::col::AccentHover;
  if (category == "lanthanide") return style::col::AccentActive;
  if (category == "actinide") return style::col::Danger;
  return style::col::TextDim;
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

ImU32 propertyRamp(float value, float alpha = 1.0f) {
  const float t = std::clamp(value, 0.0f, 1.0f);
  if (t <= 0.5f) {
    return style::mix(style::col::BgSurface, style::col::Teal, t * 2.0f, alpha);
  }
  return style::mix(style::col::Teal, style::col::Accent, (t - 0.5f) * 2.0f, alpha);
}

float normaliseMass(const ElementData& element, double minMass, double maxMass) {
  if (maxMass <= minMass) return 0.0f;
  return static_cast<float>((element.mass - minMass) / (maxMass - minMass));
}

std::string fitText(ImFont* font, float fontSize, const std::string& text, float width) {
  if (width <= 0.0f) return {};
  const float unlimited = std::numeric_limits<float>::max();
  if (font->CalcTextSizeA(fontSize, unlimited, 0.0f, text.c_str()).x <= width)
    return text;

  constexpr const char* suffix = "...";
  const float suffixWidth = font->CalcTextSizeA(fontSize, unlimited, 0.0f, suffix).x;
  if (suffixWidth > width) return {};

  std::size_t end = text.size();
  while (end > 0) {
    --end;
    while (end > 0 && (static_cast<unsigned char>(text[end]) & 0xC0u) == 0x80u) --end;
    const float candidateWidth =
        font->CalcTextSizeA(fontSize, unlimited, 0.0f, text.data(), text.data() + end).x;
    if (candidateWidth + suffixWidth <= width) return text.substr(0, end) + suffix;
  }
  return {};
}

void elementTooltip(const ElementData& element) {
  const float fs = ImGui::GetFontSize();
  ImGui::SetNextWindowSizeConstraints(ImVec2(fs * 18.0f, 0.0f),
                                      ImVec2(fs * 27.0f, fs * 24.0f));
  ImGui::PushStyleColor(ImGuiCol_PopupBg, style::col::BgRaised);
  ImGui::PushStyleColor(ImGuiCol_Border, style::col::BorderStrong);
  ImGui::BeginTooltip();

  widgets::cardHeader(icons::Icon::Atom, element.name.c_str(), element.symbol.c_str(),
                      categoryColor(element.category));

  char number[16];
  std::snprintf(number, sizeof(number), "%u", static_cast<unsigned>(element.z));
  char mass[32];
  std::snprintf(mass, sizeof(mass), "%.6g u", element.mass);
  char group[16];
  std::snprintf(group, sizeof(group), "%d", element.group);
  char period[16];
  std::snprintf(period, sizeof(period), "%d", element.period);

  widgets::keyValue("Symbol", element.symbol.c_str(), categoryColor(element.category));
  widgets::keyValue("Name", element.name.c_str());
  widgets::keyValue("Atomic number", number);
  widgets::keyValue("Atomic mass", mass);
  widgets::keyValue("Group", element.group == 0 ? "f-block" : group);
  widgets::keyValue("Period", period);
  widgets::keyValue("Category", categoryDisplayName(element.category),
                    categoryColor(element.category));

  ImGui::Spacing();
  widgets::notice(icons::Icon::Atom, "Click to select this element and activate the Atom tool.",
                  style::col::Accent);

  ImGui::EndTooltip();
  ImGui::PopStyleColor(2);
}

// Each tile keeps three readable levels of hierarchy even when the table must
// scroll. The name is omitted only if its measured ellipsis cannot fit.
void drawTile(ImDrawList* dl, ImVec2 min, ImVec2 size, const ElementData& element,
              bool selected, bool queryActive, bool searchMatch, float hover,
              ColourMode mode, double minMass, double maxMass) {
  const style::Metrics& m = style::metrics();
  const float fs = ImGui::GetFontSize();
  const float alpha = !queryActive || searchMatch ? 1.0f : 0.22f;
  const float radius = std::min(m.radiusMd, fs * 0.45f);
  const ImVec2 max(min.x + size.x, min.y + size.y);
  const ImVec4 category = categoryColor(element.category);
  const float massT = normaliseMass(element, minMass, maxMass);

  if (mode == ColourMode::Category) {
    dl->AddRectFilled(min, max,
                      style::mix(style::col::BgSurface, category, 0.10f + hover * 0.12f,
                                 alpha),
                      radius);
  } else {
    dl->AddRectFilled(min, max, propertyRamp(massT, alpha * (0.72f + hover * 0.18f)),
                      radius);
  }

  const ImU32 modeBorder = mode == ColourMode::Category
                               ? style::mix(style::col::Border, category,
                                            0.30f + hover * 0.48f, alpha)
                               : propertyRamp(massT, alpha * (0.65f + hover * 0.25f));
  dl->AddRect(min, max, modeBorder, radius, 0, m.hairline);

  if (queryActive && searchMatch) {
    const float inset = m.hairline * 2.0f;
    dl->AddRect(ImVec2(min.x + inset, min.y + inset),
                ImVec2(max.x - inset, max.y - inset), style::u32(style::col::Teal, alpha),
                std::max(radius - inset, 0.0f), 0, m.hairline * 2.0f);
  }
  if (selected) {
    const float glow = m.hairline * 2.0f;
    dl->AddRect(ImVec2(min.x - glow, min.y - glow),
                ImVec2(max.x + glow, max.y + glow),
                style::u32(style::col::Accent, 0.28f * alpha), radius + glow, 0,
                m.hairline * 3.0f);
    dl->AddRect(min, max, style::u32(style::col::Accent, alpha), radius, 0,
                m.hairline * 2.0f);
  }

  const float padX = fs * 0.34f;
  const float numberSize = fs * 0.62f;
  ImFont* body = style::fonts::body() ? style::fonts::body() : ImGui::GetFont();
  ImFont* semibold =
      style::fonts::semibold() ? style::fonts::semibold() : ImGui::GetFont();

  char zbuf[8];
  std::snprintf(zbuf, sizeof(zbuf), "%u", static_cast<unsigned>(element.z));
  dl->AddText(body, numberSize, ImVec2(min.x + padX, min.y + fs * 0.20f),
              style::u32(style::col::TextFaint, alpha), zbuf);

  const float symbolSize = fs * 1.42f;
  const ImVec2 symbolExtent =
      semibold->CalcTextSizeA(symbolSize, size.x, 0.0f, element.symbol.c_str());
  const float symbolY = min.y + fs * 1.14f;
  dl->AddText(semibold, symbolSize,
              ImVec2(min.x + (size.x - symbolExtent.x) * 0.5f, symbolY),
              selected ? style::u32(style::col::Accent, alpha)
                       : style::u32(style::col::Text, alpha),
              element.symbol.c_str());

  const float stripeHeight = m.hairline * 2.0f;
  const float stripeY = max.y - fs * 0.28f;
  const float stripeInset = fs * 0.38f;
  const ImU32 stripe = mode == ColourMode::Category
                           ? style::u32(category, alpha * (0.72f + hover * 0.28f))
                           : propertyRamp(massT, alpha);
  dl->AddRectFilled(ImVec2(min.x + stripeInset, stripeY - stripeHeight),
                    ImVec2(max.x - stripeInset, stripeY), stripe,
                    stripeHeight * 0.5f);

  const float nameSize = fs * 0.66f;
  const float nameWidth = size.x - padX * 2.0f;
  const std::string fittedName = fitText(body, nameSize, element.name, nameWidth);
  const ImVec2 nameExtent =
      fittedName.empty()
          ? ImVec2(0.0f, 0.0f)
          : body->CalcTextSizeA(nameSize, nameWidth, 0.0f, fittedName.c_str());
  const float nameY = stripeY - stripeHeight - fs * 0.24f - nameExtent.y;
  const float availableNameHeight = nameY - (symbolY + symbolExtent.y);
  if (!fittedName.empty() && availableNameHeight >= fs * 0.08f) {
    dl->AddText(body, nameSize,
                ImVec2(min.x + (size.x - nameExtent.x) * 0.5f, nameY),
                style::u32(style::col::TextDim, alpha), fittedName.c_str());
  }
}

void drawPlaceholder(ImDrawList* dl, ImVec2 min, ImVec2 size, const char* label) {
  const style::Metrics& m = style::metrics();
  const float fs = ImGui::GetFontSize();
  const ImVec2 max(min.x + size.x, min.y + size.y);
  const ImU32 color = style::u32(style::col::TextFaint, 0.55f);
  const float dash = fs * 0.48f;
  const float dashGap = fs * 0.32f;
  auto dashed = [&](ImVec2 a, ImVec2 b) {
    const float dx = b.x - a.x;
    const float dy = b.y - a.y;
    const float length = std::sqrt(dx * dx + dy * dy);
    if (length <= m.hairline) return;
    const float ux = dx / length;
    const float uy = dy / length;
    for (float start = 0.0f; start < length; start += dash + dashGap) {
      const float end = std::min(start + dash, length);
      dl->AddLine(ImVec2(a.x + ux * start, a.y + uy * start),
                  ImVec2(a.x + ux * end, a.y + uy * end), color, m.hairline);
    }
  };

  const float inset = m.hairline;
  const ImVec2 innerMin(min.x + inset, min.y + inset);
  const ImVec2 innerMax(max.x - inset, max.y - inset);
  dashed(innerMin, ImVec2(innerMax.x, innerMin.y));
  dashed(ImVec2(innerMax.x, innerMin.y), innerMax);
  dashed(innerMax, ImVec2(innerMin.x, innerMax.y));
  dashed(ImVec2(innerMin.x, innerMax.y), innerMin);

  ImFont* font = style::fonts::semibold() ? style::fonts::semibold() : ImGui::GetFont();
  const float fontSize = fs * 0.78f;
  const ImVec2 extent = font->CalcTextSizeA(fontSize, size.x, 0.0f, label);
  dl->AddText(font, fontSize,
              ImVec2(min.x + (size.x - extent.x) * 0.5f,
                     min.y + (size.y - extent.y) * 0.5f),
              color, label);
}

void drawMassLegend(double minMass, double maxMass) {
  const style::Metrics& m = style::metrics();
  const float fs = ImGui::GetFontSize();
  const float width = std::max(ImGui::GetContentRegionAvail().x, fs);
  const float labelHeight = ImGui::GetTextLineHeight();
  const float stripHeight = std::max(m.hairline * 4.0f, fs * 0.30f);
  const float totalHeight = labelHeight + m.gap * 0.45f + stripHeight;
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::Dummy(ImVec2(width, totalHeight));
  ImDrawList* dl = ImGui::GetWindowDrawList();

  char low[32];
  char high[32];
  std::snprintf(low, sizeof(low), "%.4g u", minMass);
  std::snprintf(high, sizeof(high), "%.4g u", maxMass);
  dl->AddText(min, style::u32(style::col::TextDim), low);
  const ImVec2 highExtent = ImGui::CalcTextSize(high);
  dl->AddText(ImVec2(min.x + width - highExtent.x, min.y),
              style::u32(style::col::TextDim), high);

  const float stripY = min.y + labelHeight + m.gap * 0.45f;
  constexpr int segments = 48;
  for (int i = 0; i < segments; ++i) {
    const float x0 = min.x + width * static_cast<float>(i) / segments;
    const float x1 = min.x + width * static_cast<float>(i + 1) / segments;
    dl->AddRectFilled(ImVec2(x0, stripY), ImVec2(x1, stripY + stripHeight),
                      propertyRamp((static_cast<float>(i) + 0.5f) / segments));
  }
  dl->AddRect(ImVec2(min.x, stripY), ImVec2(min.x + width, stripY + stripHeight),
              style::u32(style::col::BorderStrong), m.radiusSm, 0, m.hairline);
}

}  // namespace

void drawPeriodicTable(AppState& st) {
  static std::string search;
  static std::string previousQuery;
  static uint8_t pendingJump = 0;
  static int colourModeIndex = 0;
  static bool loadErrorSurfaced = false;

  const auto& elements = elementTable();
  if (!elementTableLoadError().empty() && !loadErrorSurfaced) {
    st.statusMessage = elementTableLoadError();
    loadErrorSurfaced = true;
  }

  const style::Metrics& m = style::metrics();
  const float fs = ImGui::GetFontSize();
  widgets::beginToolbar("##periodic_search_toolbar");
  const float glyphFrame = ImGui::GetFrameHeight();
  const ImVec2 glyphMin = ImGui::GetCursorScreenPos();
  ImGui::Dummy(ImVec2(glyphFrame, glyphFrame));
  icons::draw(ImGui::GetWindowDrawList(), icons::Icon::Search,
              ImVec2(glyphMin.x + glyphFrame * 0.5f, glyphMin.y + glyphFrame * 0.5f),
              std::min(m.iconSize, glyphFrame * 0.66f), style::u32(style::col::TextDim));
  ImGui::SameLine();
  ImGui::SetNextItemWidth(-1.0f);
  const bool submitted = widgets::stringInputWithHint(
      "##element_search", "Symbol, name or atomic number", search,
      ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_AutoSelectAll);
  widgets::endToolbar();

  const std::string query = lower(search);
  std::vector<const ElementData*> matched;
  if (!query.empty()) {
    for (const auto& element : elements) {
      if (matches(element, query)) matched.push_back(&element);
    }
  }

  if (query != previousQuery) {
    pendingJump = 0;
    if (!matched.empty()) {
      const auto exact = std::find_if(matched.begin(), matched.end(), [&](const ElementData* e) {
        return exactMatch(*e, query);
      });
      pendingJump = static_cast<uint8_t>((exact != matched.end() ? *exact : matched.front())->z);
    }
    previousQuery = query;
  }
  if (submitted && matched.size() == 1) selectElement(st, *matched.front());

  ImGui::Spacing();
  static constexpr const char* modeLabels[] = {"Category", "Atomic mass"};
  widgets::segmented("##element_colour_mode", modeLabels,
                     2, colourModeIndex);

  double minMass = 0.0;
  double maxMass = 0.0;
  if (!elements.empty()) {
    minMass = elements.front().mass;
    maxMass = elements.front().mass;
    for (const auto& element : elements) {
      minMass = std::min(minMass, element.mass);
      maxMass = std::max(maxMass, element.mass);
    }
  }
  const ColourMode colourMode = colourModeIndex == 0 ? ColourMode::Category
                                                      : ColourMode::AtomicMass;
  if (colourMode == ColourMode::AtomicMass) {
    ImGui::Spacing();
    drawMassLegend(minMass, maxMass);
  }
  ImGui::Spacing();

  // A font-derived minimum protects legibility. Narrow panels scroll instead of
  // shrinking chemical symbols into an unreadable mosaic.
  const ImVec2 avail = ImGui::GetContentRegionAvail();
  const float gapX = std::max(m.hairline, ImGui::GetStyle().ItemSpacing.x * 0.24f);
  const float gapY = std::max(m.hairline * 2.0f, ImGui::GetStyle().ItemSpacing.y * 0.38f);
  const float labelGutter = fs * 1.45f;
  const float groupHeader = ImGui::GetTextLineHeightWithSpacing();
  const float minCellWidth = fs * 3.70f;
  const float maxCellWidth = fs * 5.20f;
  const float fittedWidth = (avail.x - labelGutter - gapX * 17.0f) / 18.0f;
  const float cellWidth = std::clamp(fittedWidth, minCellWidth, maxCellWidth);
  const float cellHeight = std::max(fs * 4.45f, cellWidth * 1.08f);
  const float stepX = cellWidth + gapX;
  const float stepY = cellHeight + gapY;
  const float width = labelGutter + cellWidth * 18.0f + gapX * 17.0f;
  const float height = groupHeader + cellHeight * 10.0f + gapY * 9.0f;

  ImGui::SetNextWindowContentSize(ImVec2(width, height));
  ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
  const bool visible = ImGui::BeginChild("##ptable_grid", ImVec2(0.0f, 0.0f),
                                         ImGuiChildFlags_None,
                                         ImGuiWindowFlags_HorizontalScrollbar);
  ImGui::PopStyleVar();
  if (!visible) {
    ImGui::EndChild();
    return;
  }

  ImDrawList* dl = ImGui::GetWindowDrawList();
  const ImVec2 contentOrigin = ImGui::GetCursorPos();
  const ImVec2 origin(contentOrigin.x + labelGutter, contentOrigin.y + groupHeader);
  const ImVec2 contentScreen = ImGui::GetCursorScreenPos();
  const ImVec2 originScreen(contentScreen.x + labelGutter, contentScreen.y + groupHeader);

  for (int groupIndex = 1; groupIndex <= 18; ++groupIndex) {
    char label[4];
    std::snprintf(label, sizeof(label), "%d", groupIndex);
    const ImVec2 extent = ImGui::CalcTextSize(label);
    dl->AddText(ImVec2(originScreen.x + (groupIndex - 1) * stepX +
                                   (cellWidth - extent.x) * 0.5f,
                           contentScreen.y + (groupHeader - extent.y) * 0.5f),
                style::u32(style::col::TextFaint), label);
  }
  for (int periodIndex = 1; periodIndex <= 7; ++periodIndex) {
    char label[4];
    std::snprintf(label, sizeof(label), "%d", periodIndex);
    const ImVec2 extent = ImGui::CalcTextSize(label);
    dl->AddText(ImVec2(contentScreen.x + (labelGutter - extent.x) * 0.5f,
                           originScreen.y + (periodIndex - 1) * stepY +
                               (cellHeight - extent.y) * 0.5f),
                style::u32(style::col::TextFaint), label);
  }
  static constexpr const char* detachedLabels[] = {"Ln", "An"};
  for (int i = 0; i < 2; ++i) {
    const ImVec2 extent = ImGui::CalcTextSize(detachedLabels[i]);
    dl->AddText(ImVec2(contentScreen.x + (labelGutter - extent.x) * 0.5f,
                           originScreen.y + (8 + i) * stepY +
                               (cellHeight - extent.y) * 0.5f),
                style::u32(style::col::TextFaint), detachedLabels[i]);
  }

  auto placeholderAt = [&](int group, int period, const char* label, const char* tooltip) {
    const ImVec2 pos(origin.x + (group - 1) * stepX, origin.y + (period - 1) * stepY);
    ImGui::SetCursorPos(pos);
    ImGui::PushID(label);
    ImGui::InvisibleButton("##placeholder", ImVec2(cellWidth, cellHeight));
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", tooltip);
    drawPlaceholder(dl, ImGui::GetItemRectMin(), ImVec2(cellWidth, cellHeight), label);
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

    const bool searchMatch = matches(element, query);
    const ImVec2 pos(origin.x + (column - 1) * stepX, origin.y + (row - 1) * stepY);
    ImGui::SetCursorPos(pos);
    ImGui::PushID(static_cast<int>(element.z));
    ImGui::InvisibleButton("##tile", ImVec2(cellWidth, cellHeight));
    const bool clicked = ImGui::IsItemClicked();
    const bool hovered = ImGui::IsItemHovered();
    const float hover = widgets::hoverT(ImGui::GetItemID(), hovered);

    drawTile(dl, ImGui::GetItemRectMin(), ImVec2(cellWidth, cellHeight), element,
             st.currentElement == element.z, !query.empty(), searchMatch, hover, colourMode,
             minMass, maxMass);

    if (pendingJump == element.z) {
      ImGui::SetScrollHereX(0.5f);
      ImGui::SetScrollHereY(0.5f);
      pendingJump = 0;
    }
    if (clicked) selectElement(st, element);
    if (hovered) elementTooltip(element);
    ImGui::PopID();
  }

  ImGui::EndChild();
}

}  // namespace chemcad::ui

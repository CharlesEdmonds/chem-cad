#include <algorithm>
#include <cctype>
#include <cstdio>
#include <string>
#include <vector>

#include "imgui.h"

#include "ui/charts.hpp"
#include "ui/element_data.hpp"
#include "ui/icons.hpp"
#include "ui/layout.hpp"
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


// Categories remain distinguishable without spending amber, which is reserved
// for the selected tile and other actions.
ImVec4 categoryColor(const std::string& category) {
  if (category == "alkali-metal") return style::col::DataBright;
  if (category == "alkaline-earth") return style::col::Data;
  if (category == "transition-metal") return style::col::DataDim;
  if (category == "post-transition") return style::col::Teal;
  if (category == "metalloid") return style::col::DataBright;
  if (category == "nonmetal") return style::col::Data;
  if (category == "halogen") return style::col::DataDim;
  if (category == "noble-gas") return style::col::Teal;
  if (category == "lanthanide") return style::col::DataBright;
  if (category == "actinide") return style::col::DataDim;
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
    return style::mix(style::col::BgSurface, style::col::DataDim, t * 2.0f, alpha);
  }
  return style::mix(style::col::DataDim, style::col::DataBright,
                    (t - 0.5f) * 2.0f, alpha);
}

float normaliseMass(const ElementData& element, const charts::Axis& axis) {
  return static_cast<float>(axis.normalise(element.mass));
}

void elementTooltip(const ElementData& element) {
  const layout::Frame frame = layout::measure();
  ImGui::SetNextWindowSizeConstraints(ImVec2(frame.em * 18.0f, 0.0f),
                                      ImVec2(frame.em * 27.0f, frame.em * 24.0f));
  ImGui::PushStyleColor(ImGuiCol_PopupBg, style::col::BgRaised);
  ImGui::PushStyleColor(ImGuiCol_Border, style::col::BorderStrong);
  ImGui::BeginTooltip();

  widgets::cardHeader(icons::Icon::Atom, element.name.c_str(), element.symbol.c_str(),
                      categoryColor(element.category));

  char number[16];
  char mass[32];
  char group[16];
  char period[16];
  std::snprintf(number, sizeof(number), "%u", static_cast<unsigned>(element.z));
  std::snprintf(mass, sizeof(mass), "%.6g u", element.mass);
  std::snprintf(group, sizeof(group), "%d", element.group);
  std::snprintf(period, sizeof(period), "%d", element.period);

  widgets::keyValue("Symbol", element.symbol.c_str(), categoryColor(element.category));
  widgets::keyValue("Name", element.name.c_str());
  widgets::keyValue("Atomic number", number);
  widgets::keyValue("Atomic mass", mass, style::col::DataBright);
  widgets::keyValue("Group", element.group == 0 ? "f-block" : group);
  widgets::keyValue("Period", period);
  widgets::keyValue("Category", categoryDisplayName(element.category),
                    categoryColor(element.category));
  ImGui::Spacing();
  widgets::notice(icons::Icon::Atom,
                  "Click to select this element and activate the Atom tool.",
                  style::col::DataDim);

  ImGui::EndTooltip();
  ImGui::PopStyleColor(2);
}

void drawTile(ImDrawList* dl, ImVec2 min, ImVec2 size, const ElementData& element,
              bool selected, bool queryActive, bool searchMatch, float hover,
              ColourMode mode, const charts::Axis& massAxis,
              const layout::Frame& frame) {
  const style::Metrics& m = style::metrics();
  const float alpha = !queryActive || searchMatch ? 1.0f : 0.22f;
  const ImVec2 max(min.x + size.x, min.y + size.y);
  const ImVec4 category = categoryColor(element.category);
  const float massT = normaliseMass(element, massAxis);
  const ImU32 fill = mode == ColourMode::Category
                         ? style::mix(style::col::BgSurface, category,
                                      0.10f + hover * 0.12f, alpha)
                         : propertyRamp(massT, alpha * (0.72f + hover * 0.18f));
  dl->AddRectFilled(min, max, fill, std::min(m.radiusMd, frame.em * 0.45f));

  const ImU32 border = mode == ColourMode::Category
                           ? style::mix(style::col::Border, category,
                                        0.30f + hover * 0.48f, alpha)
                           : propertyRamp(massT, alpha * (0.65f + hover * 0.25f));
  dl->AddRect(min, max, border, std::min(m.radiusMd, frame.em * 0.45f), 0,
              m.hairline);

  if (queryActive && searchMatch) {
    const float inset = m.hairline * 2.0f;
    dl->AddRect(ImVec2(min.x + inset, min.y + inset),
                ImVec2(max.x - inset, max.y - inset),
                style::u32(style::col::Teal, alpha), 0.0f, 0, m.hairline * 2.0f);
  }
  if (selected) {
    dl->AddRect(min, max, style::u32(style::col::Accent, alpha),
                std::min(m.radiusMd, frame.em * 0.45f), 0, m.hairline * 2.0f);
  }

  ImFont* body = style::fonts::body() ? style::fonts::body() : ImGui::GetFont();
  ImFont* semibold =
      style::fonts::semibold() ? style::fonts::semibold() : ImGui::GetFont();
  const float textWidth = std::max(size.x - frame.gap * 0.5f, 0.0f);
  char number[8];
  std::snprintf(number, sizeof(number), "%u", static_cast<unsigned>(element.z));

  // Width and height decide the information level; the font itself never falls
  // beneath the shared readability floor.
  const bool numberFits = layout::fits(number, textWidth);
  const bool enoughForNumber = size.y >= frame.row * 1.75f;
  const bool showNumber = numberFits && enoughForNumber;
  const bool nameFits = layout::fits(element.name.c_str(), textWidth);
  const bool enoughForName = size.y >= frame.row * 2.75f;
  const bool showName = showNumber && nameFits && enoughForName;

  const float symbolSize = layout::labelFont(
      std::min(frame.em, std::min(size.x * 0.45f, size.y * 0.42f)));
  const ImVec2 symbolExtent =
      semibold->CalcTextSizeA(symbolSize, size.x, 0.0f, element.symbol.c_str());
  float symbolY = min.y + (size.y - symbolExtent.y) * 0.5f;
  if (showNumber) symbolY += frame.row * (showName ? -0.05f : 0.12f);

  dl->PushClipRect(min, max, true);
  if (showNumber) {
    const float numberSize = layout::labelFont(frame.em * 0.72f);
    dl->AddText(body, numberSize,
                ImVec2(min.x + frame.gap * 0.24f, min.y + frame.gap * 0.12f),
                style::u32(style::col::TextFaint, alpha), number);
  }
  dl->AddText(semibold, symbolSize,
              ImVec2(min.x + (size.x - symbolExtent.x) * 0.5f, symbolY),
              selected ? style::u32(style::col::Accent, alpha)
                       : style::u32(style::col::Text, alpha),
              element.symbol.c_str());
  if (showName) {
    const float nameSize = layout::labelFont(frame.em * 0.72f);
    const ImVec2 nameExtent =
        body->CalcTextSizeA(nameSize, textWidth, 0.0f, element.name.c_str());
    dl->AddText(body, nameSize,
                ImVec2(min.x + (size.x - nameExtent.x) * 0.5f,
                       max.y - nameExtent.y - frame.gap * 0.18f),
                style::u32(style::col::TextDim, alpha), element.name.c_str());
  }
  dl->PopClipRect();
}

void drawPlaceholder(ImDrawList* dl, ImVec2 min, ImVec2 size, const char* label,
                     const layout::Frame& frame) {
  const style::Metrics& m = style::metrics();
  const ImVec2 max(min.x + size.x, min.y + size.y);
  dl->AddRect(min, max, style::u32(style::col::GridLine, 0.7f),
              std::min(m.radiusSm, frame.em * 0.3f), 0, m.hairline);
  const float fontSize = layout::labelFont(frame.em * 0.78f);
  ImFont* font = style::fonts::semibold() ? style::fonts::semibold() : ImGui::GetFont();
  const ImVec2 extent = font->CalcTextSizeA(fontSize, size.x, 0.0f, label);
  dl->AddText(font, fontSize,
              ImVec2(min.x + (size.x - extent.x) * 0.5f,
                     min.y + (size.y - extent.y) * 0.5f),
              style::u32(style::col::TextFaint), label);
}

void drawMassLegend(const charts::Axis& axis, ImVec2 size,
                    const layout::Frame& frame) {
  const style::Metrics& m = style::metrics();
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::Dummy(size);
  if (size.x <= 0.0f || size.y <= 0.0f) return;

  ImDrawList* dl = ImGui::GetWindowDrawList();
  const float labelSize = layout::labelFont(frame.em * 0.72f);
  const float stripHeight = std::max(m.hairline * 3.0f, size.y - frame.row);
  const float stripY = min.y + size.y - stripHeight;
  const int segments = std::max(16, static_cast<int>(frame.ems() * 2.0f));
  for (int i = 0; i < segments; ++i) {
    const float x0 = min.x + size.x * static_cast<float>(i) / segments;
    const float x1 = min.x + size.x * static_cast<float>(i + 1) / segments;
    dl->AddRectFilled(ImVec2(x0, stripY), ImVec2(x1, stripY + stripHeight),
                      propertyRamp((static_cast<float>(i) + 0.5f) / segments));
  }
  dl->AddRect(ImVec2(min.x, stripY), ImVec2(min.x + size.x, stripY + stripHeight),
              style::u32(style::col::GridLine), m.radiusSm, 0, m.hairline);

  for (int tick = 0; tick < axis.ticks; ++tick) {
    const double value = axis.min + axis.step * static_cast<double>(tick);
    if (value > axis.max + axis.step * 0.25) break;
    const float x = min.x + static_cast<float>(axis.normalise(value)) * size.x;
    char label[32];
    std::snprintf(label, sizeof(label), "%.*f u", axis.decimals, value);
    const ImVec2 extent = ImGui::CalcTextSize(label);
    const float labelX = std::clamp(x - extent.x * 0.5f, min.x,
                                    min.x + size.x - extent.x);
    dl->AddText(nullptr, labelSize, ImVec2(labelX, min.y),
                style::u32(style::col::DataDim), label);
    dl->AddLine(ImVec2(x, stripY), ImVec2(x, stripY + stripHeight),
                style::u32(style::col::GridLine), m.hairline);
  }
}

void advanceVerticalGap(float gap) {
  layout::nextRow(ImGui::GetCursorPosY() + gap -
                       ImGui::GetStyle().ItemSpacing.y);
}

}  // namespace

void drawPeriodicTable(AppState& st) {
  static std::string search;
  static int colourModeIndex = 0;
  static bool loadErrorSurfaced = false;

  const auto& elements = elementTable();
  if (!elementTableLoadError().empty() && !loadErrorSurfaced) {
    st.statusMessage = elementTableLoadError();
    loadErrorSurfaced = true;
  }

  const layout::Frame page = layout::measure();
  const style::Metrics& m = style::metrics();
  const float budget = std::min(page.size.y, layout::pageHeight());
  const bool massMode = colourModeIndex == static_cast<int>(ColourMode::AtomicMass);
  float rows[4]{};
  if (massMode) {
    const float weights[] = {0.0f, 0.0f, 0.0f, 1.0f};
    const float minimums[] = {
        page.control, page.control, page.row * 1.6f, 0.0f};
    layout::distribute(budget, weights, minimums, 4, page.gap, rows);
  } else {
    const float weights[] = {0.0f, 0.0f, 1.0f};
    const float minimums[] = {page.control, page.control, 0.0f};
    float compactRows[3]{};
    layout::distribute(budget, weights, minimums, 3, page.gap,
                       compactRows);
    rows[0] = compactRows[0];
    rows[1] = compactRows[1];
    rows[3] = compactRows[2];
  }

  double minMass = 0.0;
  double maxMass = 0.0;
  if (!elements.empty()) {
    minMass = maxMass = elements.front().mass;
    for (const auto& element : elements) {
      minMass = std::min(minMass, element.mass);
      maxMass = std::max(maxMass, element.mass);
    }
  }
  const int targetTicks = std::clamp(static_cast<int>(page.ems() / 12.0f), 2, 7);
  const charts::Axis massAxis = charts::niceAxis(minMass, maxMass, targetTicks);

  if (ImGui::BeginChild("##periodic_search", ImVec2(page.size.x, rows[0]),
                        ImGuiChildFlags_None,
                        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    widgets::beginToolbar("##periodic_search_toolbar");
    const float glyphFrame = ImGui::GetFrameHeight();
    const ImVec2 glyphMin = ImGui::GetCursorScreenPos();
    ImGui::Dummy(ImVec2(glyphFrame, glyphFrame));
    icons::draw(ImGui::GetWindowDrawList(), icons::Icon::Search,
                ImVec2(glyphMin.x + glyphFrame * 0.5f,
                       glyphMin.y + glyphFrame * 0.5f),
                std::min(m.iconSize, glyphFrame * 0.66f),
                style::u32(style::col::TextDim));
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
    if (submitted && matched.size() == 1) selectElement(st, *matched.front());
  }
  ImGui::EndChild();
  advanceVerticalGap(page.gap);

  if (ImGui::BeginChild("##periodic_mode", ImVec2(page.size.x, rows[1]),
                        ImGuiChildFlags_None,
                        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    static constexpr const char* labels[] = {"Category", "Atomic mass"};
    widgets::segmented("##element_colour_mode", labels, 2, colourModeIndex);
  }
  ImGui::EndChild();
  advanceVerticalGap(page.gap);

  if (rows[2] > 0.0f) {
    if (ImGui::BeginChild("##periodic_legend", ImVec2(page.size.x, rows[2]),
                          ImGuiChildFlags_None,
                          ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
      drawMassLegend(massAxis, ImVec2(ImGui::GetContentRegionAvail().x, rows[2]), page);
    }
    ImGui::EndChild();
    advanceVerticalGap(page.gap);
  }

  ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
  ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0.0f, 0.0f));
  const bool gridVisible = ImGui::BeginChild(
      "##ptable_grid", ImVec2(page.size.x, rows[3]), ImGuiChildFlags_None,
      ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
  ImGui::PopStyleVar(2);
  if (!gridVisible) {
    ImGui::EndChild();
    return;
  }

  const layout::Frame grid = layout::measure(ImVec2(page.size.x, rows[3]));
  if (grid.size.x <= 0.0f || grid.size.y <= 0.0f) {
    ImGui::EndChild();
    return;
  }
  const float gapX =
      std::min(std::max(m.hairline, grid.gap * 0.12f),
               grid.size.x / 36.0f);
  const float gapY =
      std::min(std::max(m.hairline, grid.gap * 0.12f),
               grid.size.y / 20.0f);
  const float rawCellWidth = (grid.size.x - gapX * 17.0f) / 18.0f;
  const float rawCellHeight = (grid.size.y - gapY * 9.0f) / 10.0f;
  const bool showColumnLabels = rawCellHeight >= grid.row * 1.6f;
  const bool showRowLabels = rawCellWidth >= grid.em * 1.6f;
  const float groupHeader = showColumnLabels ? grid.row : 0.0f;
  const float labelGutter = showRowLabels ? grid.em * 1.4f : 0.0f;
  const float cellWidth =
      std::max((grid.size.x - labelGutter - gapX * 17.0f) / 18.0f, 0.0f);
  const float cellHeight =
      std::max((grid.size.y - groupHeader - gapY * 9.0f) / 10.0f, 0.0f);
  const float stepX = cellWidth + gapX;
  const float stepY = cellHeight + gapY;
  const ImVec2 cursor = ImGui::GetCursorPos();
  const ImVec2 screen = ImGui::GetCursorScreenPos();
  const ImVec2 origin(cursor.x + labelGutter, cursor.y + groupHeader);
  const ImVec2 originScreen(screen.x + labelGutter, screen.y + groupHeader);
  ImDrawList* dl = ImGui::GetWindowDrawList();
  const float gridLabelSize = layout::labelFont(grid.em * 0.72f);

  if (showColumnLabels) {
    for (int group = 1; group <= 18; ++group) {
      char label[4];
      std::snprintf(label, sizeof(label), "%d", group);
      const ImVec2 extent = ImGui::CalcTextSize(label);
      dl->AddText(nullptr, gridLabelSize,
                  ImVec2(originScreen.x + (group - 1) * stepX +
                             (cellWidth - extent.x) * 0.5f,
                         screen.y + (groupHeader - extent.y) * 0.5f),
                  style::u32(style::col::TextFaint), label);
    }
  }
  if (showRowLabels) {
    for (int period = 1; period <= 7; ++period) {
      char label[4];
      std::snprintf(label, sizeof(label), "%d", period);
      const ImVec2 extent = ImGui::CalcTextSize(label);
      dl->AddText(nullptr, gridLabelSize,
                  ImVec2(screen.x + (labelGutter - extent.x) * 0.5f,
                         originScreen.y + (period - 1) * stepY +
                             (cellHeight - extent.y) * 0.5f),
                  style::u32(style::col::TextFaint), label);
    }
    static constexpr const char* detachedLabels[] = {"Ln", "An"};
    for (int index = 0; index < 2; ++index) {
      const ImVec2 extent = ImGui::CalcTextSize(detachedLabels[index]);
      dl->AddText(nullptr, gridLabelSize,
                  ImVec2(screen.x + (labelGutter - extent.x) * 0.5f,
                         originScreen.y + (8 + index) * stepY +
                             (cellHeight - extent.y) * 0.5f),
                  style::u32(style::col::TextFaint), detachedLabels[index]);
    }
  }

  auto placeholderAt = [&](int group, int period, const char* label,
                           const char* tooltip) {
    ImGui::SetCursorPos(
        ImVec2(origin.x + (group - 1) * stepX, origin.y + (period - 1) * stepY));
    ImGui::PushID(label);
    ImGui::InvisibleButton("##placeholder", ImVec2(cellWidth, cellHeight));
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", tooltip);
    drawPlaceholder(dl, ImGui::GetItemRectMin(), ImVec2(cellWidth, cellHeight), label,
                    grid);
    ImGui::PopID();
  };
  placeholderAt(3, 6, "Ln", "Lanthanides (57-71), shown below");
  placeholderAt(3, 7, "An", "Actinides (89-103), shown below");

  const std::string query = lower(search);
  const ColourMode colourMode = colourModeIndex == 0 ? ColourMode::Category
                                                      : ColourMode::AtomicMass;
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

    ImGui::SetCursorPos(
        ImVec2(origin.x + (column - 1) * stepX, origin.y + (row - 1) * stepY));
    ImGui::PushID(static_cast<int>(element.z));
    ImGui::InvisibleButton("##tile", ImVec2(cellWidth, cellHeight));
    const bool clicked = ImGui::IsItemClicked();
    const bool hovered = ImGui::IsItemHovered();
    const float hover = widgets::hoverT(ImGui::GetItemID(), hovered);
    drawTile(dl, ImGui::GetItemRectMin(), ImVec2(cellWidth, cellHeight), element,
             st.currentElement == element.z, !query.empty(), matches(element, query), hover,
             colourMode, massAxis, grid);
    if (clicked) selectElement(st, element);
    if (hovered) elementTooltip(element);
    ImGui::PopID();
  }

  ImGui::EndChild();
}

}  // namespace chemcad::ui

#include <algorithm>
#include <cctype>
#include <cstring>
#include <string>
#include <vector>

#include "imgui.h"

#include "ui/element_data.hpp"
#include "ui/ui.hpp"

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

ImVec4 categoryColor(const std::string& category) {
  if (category == "alkali-metal") return {0.86f, 0.38f, 0.31f, 1.0f};
  if (category == "alkaline-earth") return {0.86f, 0.60f, 0.29f, 1.0f};
  if (category == "transition-metal") return {0.46f, 0.58f, 0.78f, 1.0f};
  if (category == "post-transition") return {0.42f, 0.69f, 0.70f, 1.0f};
  if (category == "metalloid") return {0.45f, 0.72f, 0.49f, 1.0f};
  if (category == "nonmetal") return {0.35f, 0.66f, 0.83f, 1.0f};
  if (category == "halogen") return {0.69f, 0.51f, 0.82f, 1.0f};
  if (category == "noble-gas") return {0.75f, 0.45f, 0.70f, 1.0f};
  if (category == "lanthanide") return {0.77f, 0.55f, 0.43f, 1.0f};
  if (category == "actinide") return {0.67f, 0.48f, 0.52f, 1.0f};
  return {0.50f, 0.55f, 0.60f, 1.0f};
}

ImVec4 adjusted(ImVec4 color, float amount) {
  color.x = std::clamp(color.x + amount, 0.0f, 1.0f);
  color.y = std::clamp(color.y + amount, 0.0f, 1.0f);
  color.z = std::clamp(color.z + amount, 0.0f, 1.0f);
  return color;
}

void selectElement(AppState& st, const ElementData& element) {
  st.currentElement = element.z;
  st.tool = Tool::Atom;
  st.statusMessage = "Element: " + element.name + " (" + element.symbol + ")";
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

  const std::string query = lower(search);
  std::vector<const ElementData*> matched;
  if (!query.empty()) {
    for (const auto& element : elements) {
      if (matches(element, query)) matched.push_back(&element);
    }
  }
  if (enter && matched.size() == 1) selectElement(st, *matched.front());

  // All 18 groups must be reachable without horizontal scrolling: C, N and O
  // are the most-used buttons in the app, and they sit in the right-hand
  // groups. Size the cells from the space actually available, shrinking the
  // symbol font to match rather than pushing columns off-screen.
  const float gap = std::max(2.0f, ImGui::GetStyle().ItemSpacing.x * 0.35f);
  const float avail = ImGui::GetContentRegionAvail().x;
  const float preferred = std::max(30.0f, ImGui::GetFontSize() * 2.05f);
  const float fitted = (avail - 17.0f * gap) / 18.0f;
  const float cell = std::clamp(fitted, 18.0f, preferred);
  const float step = cell + gap;
  const float width = 18.0f * step;
  const float height = 10.0f * step;
  // Two-letter symbols must still fit inside a shrunken cell.
  const float symbolScale = std::min(1.0f, cell / (ImGui::GetFontSize() * 1.9f));

  ImGui::SetNextWindowContentSize(ImVec2(width, height));
  if (!ImGui::BeginChild("##ptable_grid", ImVec2(0, 0), ImGuiChildFlags_None,
                         ImGuiWindowFlags_HorizontalScrollbar)) {
    ImGui::EndChild();
    return;
  }

  const ImVec2 origin = ImGui::GetCursorPos();
  auto drawPlaceholder = [&](int group, int period, const char* label, const char* tooltip) {
    ImGui::SetCursorPos(ImVec2(origin.x + (group - 1) * step, origin.y + (period - 1) * step));
    ImGui::BeginDisabled();
    if (symbolScale < 1.0f) ImGui::SetWindowFontScale(symbolScale);
    ImGui::Button(label, ImVec2(cell, cell));
    if (symbolScale < 1.0f) ImGui::SetWindowFontScale(1.0f);
    ImGui::EndDisabled();
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) ImGui::SetTooltip("%s", tooltip);
  };
  drawPlaceholder(3, 6, "Ln", "Lanthanides (57-71), shown below");
  drawPlaceholder(3, 7, "An", "Actinides (89-103), shown below");

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
    ImGui::SetCursorPos(ImVec2(origin.x + (column - 1) * step, origin.y + (row - 1) * step));
    ImGui::PushID(static_cast<int>(element.z));
    if (!match) ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.25f);

    const ImVec4 base = categoryColor(element.category);
    const float luminance = 0.2126f * base.x + 0.7152f * base.y + 0.0722f * base.z;
    ImGui::PushStyleColor(ImGuiCol_Button, base);
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, adjusted(base, 0.10f));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, adjusted(base, -0.08f));
    ImGui::PushStyleColor(ImGuiCol_Text,
                          luminance > 0.58f ? ImVec4(0.08f, 0.09f, 0.10f, 1.0f)
                                           : ImVec4(0.97f, 0.98f, 1.0f, 1.0f));

    const bool selected = st.currentElement == element.z;
    if (selected) {
      ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 2.5f);
      ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(1.0f, 0.95f, 0.35f, 1.0f));
    }
    if (symbolScale < 1.0f) ImGui::SetWindowFontScale(symbolScale);
    if (ImGui::Button(element.symbol.c_str(), ImVec2(cell, cell))) selectElement(st, element);
    if (symbolScale < 1.0f) ImGui::SetWindowFontScale(1.0f);
    if (selected) {
      ImGui::PopStyleColor();
      ImGui::PopStyleVar();
    }
    ImGui::PopStyleColor(4);

    if (ImGui::IsItemHovered()) {
      ImGui::BeginTooltip();
      ImGui::TextUnformatted(element.name.c_str());
      ImGui::Separator();
      ImGui::Text("Symbol: %s", element.symbol.c_str());
      ImGui::Text("Atomic number: %u", static_cast<unsigned>(element.z));
      ImGui::Text("Standard atomic weight: %.6g", element.mass);
      ImGui::Text("Category: %s", element.category.c_str());
      ImGui::EndTooltip();
    }
    if (!match) ImGui::PopStyleVar();
    ImGui::PopID();
  }

  ImGui::EndChild();
}

}  // namespace chemcad::ui

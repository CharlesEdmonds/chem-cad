#include "imgui.h"

#include <algorithm>
#include <cctype>
#include <exception>
#include <set>
#include <string>
#include <vector>

#include "rxn/kb.hpp"
#include "ui/charts.hpp"
#include "ui/icons.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"
#include "ui/widgets.hpp"
namespace chemcad::ui {
namespace {

struct ToolboxState {
  std::string query;
  std::string selectedType;
  std::string selectedSubstrate;
  std::vector<std::string> types;
  std::vector<std::string> substrates;
  std::vector<int> typeCounts;
  int layout = 0;
  bool indexed = false;
};


std::string upperCopy(const std::string& value) {
  std::string upper;
  upper.reserve(value.size());
  for (unsigned char ch : value) upper.push_back(static_cast<char>(std::toupper(ch)));
  return upper;
}

std::string join(const std::vector<std::string>& values) {
  std::string result;
  for (const std::string& value : values) {
    if (!result.empty()) result += ", ";
    result += value;
  }
  return result;
}

bool matchesQuery(const rxn::ReactionTemplate& reaction, const std::string& query) {
  if (query.empty()) return true;
  if (widgets::containsCaseInsensitive(reaction.name, query) ||
      widgets::containsCaseInsensitive(rxn::reactionType(reaction), query) ||
      widgets::containsCaseInsensitive(reaction.substrate, query) ||
      widgets::containsCaseInsensitive(reaction.conditions, query) ||
      widgets::containsCaseInsensitive(reaction.outcome, query) ||
      widgets::containsCaseInsensitive(reaction.notes, query) ||
      widgets::containsCaseInsensitive(reaction.source, query)) {
    return true;
  }
  for (const std::string& reagent : reaction.reagents) {
    if (widgets::containsCaseInsensitive(reagent, query)) return true;
  }
  for (const std::string& tag : reaction.tags) {
    if (widgets::containsCaseInsensitive(tag, query)) return true;
  }
  return false;
}


void drawChipRow(const char* id, const std::vector<std::string>& values,
                 std::string& selectedValue, ImVec4 accent) {
  const float gap = style::metrics().gap * 0.5f;
  const float availableWidth = ImGui::GetContentRegionAvail().x;
  float usedWidth = 0.0f;

  ImGui::PushID(id);
  for (const std::string& value : values) {
    const float chipWidth = ImGui::CalcTextSize(value.c_str()).x +
                            ImGui::GetStyle().FramePadding.x * 2.0f;
    if (usedWidth > 0.0f && usedWidth + gap + chipWidth <= availableWidth) {
      ImGui::SameLine(0.0f, gap);
      usedWidth += gap;
    } else if (usedWidth > 0.0f) {
      usedWidth = 0.0f;
    }
    if (widgets::chip(value.c_str(), selectedValue == value, accent)) {
      selectedValue = selectedValue == value ? std::string() : value;
    }
    usedWidth += chipWidth;
  }
  ImGui::PopID();
}

void drawReactionCard(const rxn::ReactionTemplate& reaction, const std::string& type) {
  ImGui::PushID(reaction.id.c_str());
  const float width = ImGui::GetContentRegionAvail().x;
  if (widgets::beginCard("##reaction_card", ImVec2(width, 0.0f), style::col::BgSurface)) {
    widgets::cardHeader(icons::Icon::Reaction, reaction.name.c_str(), type.c_str(),
                        style::col::Accent);

    widgets::keyValue("Substrate class",
                      reaction.substrate.empty() ? "—" : reaction.substrate.c_str(),
                      style::col::Teal);
    const char* keyReagent = reaction.reagents.empty() ? "—" : reaction.reagents.front().c_str();
    widgets::keyValue("Key reagent", keyReagent, style::col::Text);

    const std::string& description = reaction.notes.empty() ? reaction.outcome : reaction.notes;
    ImGui::PushStyleColor(ImGuiCol_Text, style::col::TextDim);
    ImGui::TextUnformatted(description.empty() ? "No summary recorded." : description.c_str());
    ImGui::PopStyleColor();
    if (ImGui::IsItemHovered() && !description.empty()) {
      ImGui::SetTooltip("%s", description.c_str());
    }

    if (widgets::disclosure("##reaction_details", "Details",
                            "Conditions, scope and provenance", false, icons::Icon::Book,
                            style::col::Teal)) {
      ImGui::Indent();
      const std::string reagents = join(reaction.reagents);
      widgets::keyValue("Reagents", reagents.empty() ? "—" : reagents.c_str());
      widgets::keyValue("Conditions",
                        reaction.conditions.empty() ? "—" : reaction.conditions.c_str());
      widgets::keyValue("Outcome", reaction.outcome.empty() ? "—" : reaction.outcome.c_str(),
                        style::col::Accent);
      widgets::keyValue("Scope / limitations",
                        reaction.notes.empty() ? "—" : reaction.notes.c_str());

      if (!reaction.source.empty()) {
        ImGui::TextColored(style::col::TextDim, "Provenance");
        ImGui::SameLine();
        const std::string displayedSource =
            reaction.source == "standard practice"
                ? "standard practice (no procedure citation)"
                : reaction.source;
        widgets::badge(displayedSource.c_str(), style::col::Teal);
      }

      if (!reaction.tags.empty()) {
        const std::string tags = join(reaction.tags);
        widgets::keyValue("Tags", tags.c_str(), style::col::Violet);
      }
      if (!reaction.byproducts.empty()) {
        const std::string byproducts = join(reaction.byproducts);
        widgets::keyValue("Byproducts", byproducts.c_str());
      }
      widgets::keyValue("Reaction SMARTS",
                        reaction.smarts.empty() ? "—" : reaction.smarts.c_str());
      const std::string priority = std::to_string(reaction.priority);
      widgets::keyValue("Library priority", priority.c_str());
      ImGui::Unindent();
    }
    widgets::endCard();
  }
  ImGui::PopID();
}

}  // namespace

void drawToolbox(AppState& st) {
  (void)st;
  static ToolboxState state;

  const std::vector<rxn::ReactionTemplate>* knowledge = nullptr;
  try {
    knowledge = &rxn::knowledgeBase();
  } catch (const std::exception& error) {
    widgets::sectionHeader("REACTION TOOLBOX", style::col::Danger);
    std::string message = "Reaction library unavailable: ";
    message += error.what();
    widgets::notice(icons::Icon::Warning, message.c_str(), style::col::Danger);
    return;
  }

  if (!state.indexed) {
    std::set<std::string> typeSet;
    std::set<std::string> substrateSet;
    for (const rxn::ReactionTemplate& reaction : *knowledge) {
      typeSet.insert(rxn::reactionType(reaction));
      if (!reaction.substrate.empty()) substrateSet.insert(reaction.substrate);
    }
    state.types.assign(typeSet.begin(), typeSet.end());
    state.substrates.assign(substrateSet.begin(), substrateSet.end());
    state.typeCounts.assign(state.types.size(), 0);
    for (const rxn::ReactionTemplate& reaction : *knowledge) {
      const std::string type = rxn::reactionType(reaction);
      const auto found = std::lower_bound(state.types.begin(), state.types.end(), type);
      if (found != state.types.end() && *found == type) {
        ++state.typeCounts[static_cast<std::size_t>(found - state.types.begin())];
      }
    }
    state.indexed = true;
  }

  widgets::sectionHeader("REACTION TOOLBOX", style::col::Accent);

  std::vector<std::string> histogramAnnotations(state.types.size());
  std::vector<charts::BarRow> histogramRows;
  histogramRows.reserve(state.types.size());
  for (std::size_t i = 0; i < state.types.size(); ++i) {
    histogramAnnotations[i] = std::to_string(state.typeCounts[i]);
    charts::BarRow row;
    row.label = state.types[i].c_str();
    row.value = static_cast<double>(state.typeCounts[i]);
    row.annotation = histogramAnnotations[i].c_str();
    row.accent = style::col::Teal;
    row.selected = state.selectedType == state.types[i];
    histogramRows.push_back(row);
  }

  if (widgets::beginCard("##library_overview",
                         ImVec2(ImGui::GetContentRegionAvail().x, 0.0f),
                         style::col::BgSurface)) {
    widgets::cardHeader(icons::Icon::ChartBars, "Library composition",
                        "Reaction families", style::col::Teal);
    const float chartHeight =
        ImGui::GetTextLineHeightWithSpacing() * static_cast<float>(histogramRows.size());
    const int clickedType =
        charts::rankedBars("##reaction_type_histogram", histogramRows.data(),
                           static_cast<int>(histogramRows.size()),
                           ImVec2(ImGui::GetContentRegionAvail().x, chartHeight));
    if (clickedType >= 0 && clickedType < static_cast<int>(state.types.size())) {
      state.selectedType = state.types[static_cast<std::size_t>(clickedType)];
    }
    widgets::endCard();
  }

  const float frameHeight = ImGui::GetFrameHeight();
  const float glyphSlot = frameHeight * 0.75f;
  const float layoutWidth = frameHeight * 2.25f;
  const float resetWidth =
      frameHeight + ImGui::CalcTextSize("Reset").x + style::metrics().gap * 1.5f;
  const icons::Icon layouts[] = {icons::Icon::Grid, icons::Icon::List};
  const char* layoutTooltips[] = {"Card grid", "Dense list"};

  widgets::beginToolbar("##toolbox_toolbar");
  const float toolbarWidth = ImGui::GetContentRegionAvail().x;
  const ImVec2 glyphMin = ImGui::GetCursorScreenPos();
  icons::draw(ImGui::GetWindowDrawList(), icons::Icon::Search,
              ImVec2(glyphMin.x + glyphSlot * 0.5f, glyphMin.y + frameHeight * 0.5f),
              ImGui::GetFontSize(), style::u32(style::col::TextDim));
  ImGui::Dummy(ImVec2(glyphSlot, frameHeight));
  ImGui::SameLine(0.0f, style::metrics().gap * 0.5f);
  const float searchWidth =
      std::max(ImGui::GetFontSize() * 10.0f,
               toolbarWidth - glyphSlot - layoutWidth - resetWidth -
                   style::metrics().gap * 6.0f);
  ImGui::SetNextItemWidth(searchWidth);
  widgets::stringInputWithHint("##query", "Search reactions, substrates, reagents",
                               state.query);
  ImGui::SameLine();
  widgets::toolbarSeparator();
  ImGui::SameLine();
  widgets::segmentedIcons("##layout", layouts, layoutTooltips, 2, state.layout, layoutWidth);
  ImGui::SameLine();
  widgets::toolbarSeparator();
  ImGui::SameLine();
  if (widgets::actionButton("##reset_filters", icons::Icon::Close, "Reset",
                            ImVec2(resetWidth, frameHeight))) {
    state.query.clear();
    state.selectedType.clear();
    state.selectedSubstrate.clear();
  }
  widgets::endToolbar();

  auto matchesFilters = [&](const rxn::ReactionTemplate& reaction) {
    return (state.selectedType.empty() || rxn::reactionType(reaction) == state.selectedType) &&
           (state.selectedSubstrate.empty() ||
            reaction.substrate == state.selectedSubstrate) &&
           matchesQuery(reaction, state.query);
  };

  std::size_t previewHits = 0;
  for (const rxn::ReactionTemplate& reaction : *knowledge) {
    if (matchesFilters(reaction)) ++previewHits;
  }
  const int activeFilterCount = static_cast<int>(!state.query.empty()) +
                                static_cast<int>(!state.selectedType.empty()) +
                                static_cast<int>(!state.selectedSubstrate.empty());
  const std::string filterSummary =
      std::to_string(activeFilterCount) + " active - " + std::to_string(previewHits) +
      " of " + std::to_string(knowledge->size());

  if (widgets::disclosure("##filters", "Filters", filterSummary.c_str(), false,
                          icons::Icon::Filter, style::col::Accent)) {
    ImGui::Indent();
    widgets::sectionHeader("REACTION TYPE", style::col::Violet);
    drawChipRow("reaction_types", state.types, state.selectedType, style::col::Violet);
    widgets::sectionHeader("SUBSTRATE", style::col::Teal);
    // A bounded chip well keeps the full taxonomy available without displacing results.
    const float substrateHeight = ImGui::GetFrameHeight() * 3.4f;
    if (ImGui::BeginChild("##substrate_chips", ImVec2(0.0f, substrateHeight))) {
      drawChipRow("reaction_substrates", state.substrates, state.selectedSubstrate,
                  style::col::Teal);
    }
    ImGui::EndChild();
    ImGui::Unindent();
  }

  std::vector<const rxn::ReactionTemplate*> visible;
  visible.reserve(knowledge->size());
  for (const rxn::ReactionTemplate& reaction : *knowledge) {
    if (matchesFilters(reaction)) visible.push_back(&reaction);
  }
  std::sort(visible.begin(), visible.end(), [](const auto* lhs, const auto* rhs) {
    const std::string lhsType = rxn::reactionType(*lhs);
    const std::string rhsType = rxn::reactionType(*rhs);
    if (lhsType != rhsType) return lhsType < rhsType;
    return lhs->name < rhs->name;
  });

  std::set<std::string> visibleTypes;
  for (const rxn::ReactionTemplate* reaction : visible) {
    visibleTypes.insert(rxn::reactionType(*reaction));
  }

  ImGui::Spacing();
  const std::string resultSummary =
      std::to_string(visible.size()) + " reactions - " +
      std::to_string(visibleTypes.size()) + " types";
  widgets::badge(resultSummary.c_str(), style::col::Accent);

  if (visible.empty()) {
    std::vector<std::string> activeNames;
    if (!state.query.empty()) activeNames.push_back("search \"" + state.query + "\"");
    if (!state.selectedType.empty()) activeNames.push_back("type " + state.selectedType);
    if (!state.selectedSubstrate.empty()) {
      activeNames.push_back("substrate " + state.selectedSubstrate);
    }
    const std::string emptyBody =
        activeNames.empty()
            ? "The loaded reaction library contains no templates."
            : "Active filters: " + join(activeNames) + ". Clear them to browse the library.";
    widgets::emptyState(icons::Icon::Search, "No reactions match", emptyBody.c_str());
    if (!activeNames.empty()) {
      const float clearWidth = frameHeight + ImGui::CalcTextSize("Clear filters").x +
                               style::metrics().gap * 2.0f;
      const float cursorX = ImGui::GetCursorPosX();
      const float centredX =
          cursorX + std::max(0.0f, (ImGui::GetContentRegionAvail().x - clearWidth) * 0.5f);
      ImGui::SetCursorPosX(centredX);
      if (widgets::actionButton("##clear_empty_filters", icons::Icon::Close,
                                "Clear filters", ImVec2(clearWidth, frameHeight), true)) {
        state.query.clear();
        state.selectedType.clear();
        state.selectedSubstrate.clear();
      }
    }
    return;
  }

  if (state.layout == 1) {
    const ImGuiTableFlags flags =
        ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerH |
        ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_NoSavedSettings;
    if (ImGui::BeginTable("##reaction_list", 3, flags,
                          ImVec2(ImGui::GetContentRegionAvail().x, 0.0f))) {
      ImGui::TableSetupColumn("Reaction", ImGuiTableColumnFlags_WidthStretch, 0.34f);
      ImGui::TableSetupColumn("Type", ImGuiTableColumnFlags_WidthStretch, 0.18f);
      ImGui::TableSetupColumn("Summary", ImGuiTableColumnFlags_WidthStretch, 0.48f);
      ImGui::TableHeadersRow();
      for (const rxn::ReactionTemplate* reaction : visible) {
        ImGui::PushID(reaction->id.c_str());
        ImGui::TableNextRow(0, ImGui::GetFrameHeight());
        ImGui::TableSetColumnIndex(0);
        const bool headingFont = style::pushFont(style::fonts::semibold());
        ImGui::TextUnformatted(reaction->name.c_str());
        style::popFont(headingFont);
        ImGui::TableSetColumnIndex(1);
        const std::string type = rxn::reactionType(*reaction);
        widgets::badge(type.c_str(), style::col::Teal);
        ImGui::TableSetColumnIndex(2);
        const std::string& summary =
            reaction->notes.empty() ? reaction->outcome : reaction->notes;
        ImGui::TextUnformatted(summary.empty() ? "No summary recorded." : summary.c_str());
        if (ImGui::IsItemHovered() && !summary.empty()) {
          ImGui::SetTooltip("%s", summary.c_str());
        }
        ImGui::PopID();
      }
      ImGui::EndTable();
    }
    return;
  }

  std::size_t groupBegin = 0;
  while (groupBegin < visible.size()) {
    const std::string type = rxn::reactionType(*visible[groupBegin]);
    std::size_t groupEnd = groupBegin + 1;
    while (groupEnd < visible.size() && rxn::reactionType(*visible[groupEnd]) == type) {
      ++groupEnd;
    }

    const std::string heading = upperCopy(type);
    widgets::sectionHeader(heading.c_str(), style::col::Accent);
    const float availableWidth = ImGui::GetContentRegionAvail().x;
    const float minimumCardWidth = ImGui::GetFontSize() * 24.0f;
    const int columnCount =
        std::clamp(static_cast<int>(availableWidth / minimumCardWidth), 1, 3);
    ImGui::PushID(type.c_str());
    if (ImGui::BeginTable("##reaction_grid", columnCount,
                          ImGuiTableFlags_SizingStretchSame |
                              ImGuiTableFlags_NoSavedSettings)) {
      for (std::size_t index = groupBegin; index < groupEnd; ++index) {
        if ((index - groupBegin) % static_cast<std::size_t>(columnCount) == 0) {
          ImGui::TableNextRow();
        }
        ImGui::TableSetColumnIndex(
            static_cast<int>((index - groupBegin) % static_cast<std::size_t>(columnCount)));
        drawReactionCard(*visible[index], type);
      }
      ImGui::EndTable();
    }
    ImGui::PopID();
    groupBegin = groupEnd;
  }
}

}  // namespace chemcad::ui

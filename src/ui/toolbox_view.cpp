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
#include "ui/layout.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"
#include "ui/widgets.hpp"

namespace chemcad::ui {
namespace {

struct ToolboxState {
  std::string query;
  std::string selectedType;
  std::string selectedSubstrate;
  std::string selectedReactionId;
  std::vector<std::string> types;
  std::vector<std::string> substrates;
  std::vector<int> typeCounts;
  int layout = 1;
  bool indexed = false;
};

std::string sentenceCase(const std::string& value) {
  std::string result = value;
  if (!result.empty()) {
    result.front() = static_cast<char>(
        std::toupper(static_cast<unsigned char>(result.front())));
  }
  return result;
}

std::string join(const std::vector<std::string>& values) {
  std::string result;
  for (const std::string& value : values) {
    if (!result.empty()) result += ", ";
    result += value;
  }
  return result;
}

const char* valueOrDash(const std::string& value) {
  return value.empty() ? "—" : value.c_str();
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

void drawReactionDetails(const rxn::ReactionTemplate& reaction) {
  const std::string reagents = join(reaction.reagents);
  const std::string tags = join(reaction.tags);
  const std::string byproducts = join(reaction.byproducts);
  const std::string arity = std::to_string(reaction.arity);
  const std::string priority = std::to_string(reaction.priority);
  const std::string displayedSource =
      reaction.source == "standard practice"
          ? "standard practice (no procedure citation)"
          : reaction.source;

  widgets::keyValue("Template ID", valueOrDash(reaction.id), style::col::Data);
  widgets::keyValue("Substrate class", valueOrDash(reaction.substrate),
                    style::col::Data);
  widgets::keyValue("Reactant count", arity.c_str(), style::col::Data);
  widgets::keyValue("Reagents", reagents.empty() ? "—" : reagents.c_str(),
                    style::col::Data);
  widgets::keyValue("Conditions", valueOrDash(reaction.conditions), style::col::Data);
  widgets::keyValue("Outcome", valueOrDash(reaction.outcome), style::col::Data);
  widgets::keyValue("Notes", valueOrDash(reaction.notes), style::col::Data);
  widgets::keyValue("Source", displayedSource.empty() ? "—" : displayedSource.c_str(),
                    style::col::Data);
  widgets::keyValue("Tags", tags.empty() ? "—" : tags.c_str(), style::col::Data);
  widgets::keyValue("Byproducts", byproducts.empty() ? "—" : byproducts.c_str(),
                    style::col::Data);
  widgets::keyValue("Reaction SMARTS", valueOrDash(reaction.smarts), style::col::Data);
  widgets::keyValue("Library priority", priority.c_str(), style::col::Data);
}

void drawReactionCard(const rxn::ReactionTemplate& reaction, const std::string& type,
                      float width) {
  ImGui::PushID(reaction.id.c_str());
  if (widgets::beginCard("##reaction_card", ImVec2(width, 0.0f),
                         style::col::BgSurface)) {
    widgets::cardHeader(icons::Icon::Reaction, reaction.name.c_str(), type.c_str(),
                        style::col::Data);
    widgets::keyValue("Substrate class", valueOrDash(reaction.substrate),
                      style::col::Data);
    const char* keyReagent =
        reaction.reagents.empty() ? "—" : reaction.reagents.front().c_str();
    widgets::keyValue("Key reagent", keyReagent, style::col::Data);

    const std::string& description =
        reaction.notes.empty() ? reaction.outcome : reaction.notes;
    ImGui::PushStyleColor(ImGuiCol_Text, style::col::TextDim);
    ImGui::TextUnformatted(description.empty() ? "No summary recorded."
                                               : description.c_str());
    ImGui::PopStyleColor();
    if (ImGui::IsItemHovered() && !description.empty()) {
      ImGui::SetTooltip("%s", description.c_str());
    }

    if (widgets::disclosure("##reaction_details", "Details",
                            "Conditions, notes and provenance", false,
                            icons::Icon::Book, style::col::Accent)) {
      ImGui::Indent();
      drawReactionDetails(reaction);
      ImGui::Unindent();
    }
    widgets::endCard();
  }
  ImGui::PopID();
}

void drawSubstrateCombo(ToolboxState& state, float width) {
  ImGui::SetNextItemWidth(width);
  const char* preview = state.selectedSubstrate.empty()
                            ? "All substrate classes"
                            : state.selectedSubstrate.c_str();
  if (!ImGui::BeginCombo("##substrate_filter", preview)) return;

  const bool allSelected = state.selectedSubstrate.empty();
  if (ImGui::Selectable("All substrate classes", allSelected)) {
    state.selectedSubstrate.clear();
  }
  if (allSelected) ImGui::SetItemDefaultFocus();
  for (const std::string& substrate : state.substrates) {
    const bool selected = state.selectedSubstrate == substrate;
    if (ImGui::Selectable(substrate.c_str(), selected)) {
      state.selectedSubstrate = substrate;
    }
    if (selected) ImGui::SetItemDefaultFocus();
  }
  ImGui::EndCombo();
}

void resetFilters(ToolboxState& state) {
  state.query.clear();
  state.selectedType.clear();
  state.selectedSubstrate.clear();
}

void drawSearchToolbar(ToolboxState& state, const layout::Frame& frame) {
  widgets::beginToolbar("##toolbox_search_toolbar");
  const float width = std::max(ImGui::GetContentRegionAvail().x - frame.gap,
                               frame.em);
  ImGui::SetNextItemWidth(width);
  widgets::stringInputWithHint("##query", "Search reactions, substrates, reagents",
                               state.query);
  widgets::endToolbar();
}

void drawFilterToolbar(ToolboxState& state, const layout::Frame& frame,
                       bool compact) {
  const icons::Icon layouts[] = {icons::Icon::Grid, icons::Icon::List};
  const char* layoutTooltips[] = {"Card grid", "Dense list"};
  const float layoutWidth = frame.control * 2.25f;
  const float resetWidth = compact
                               ? frame.control
                               : widgets::actionButtonWidth(icons::Icon::Close, "Reset");

  widgets::beginToolbar("##toolbox_filter_toolbar");
  // What the toolbar's own furniture costs. beginToolbar insets both ends by
  // half a gap and sets ItemSpacing to half a gap, and each toolbarSeparator()
  // is a full gap wide -- all of them in style::metrics().gap, NOT the frame's
  // density-scaled gap. Leave any of that out and the last control is pushed
  // off the end, the more so the higher the display scale.
  const float unit = style::metrics().gap;
  const float furniture = unit * 0.5f          // right inset
                          + unit * 0.5f * 4.0f  // four SameLine spacings
                          + unit * 2.0f;        // two separators
  const float available = ImGui::GetContentRegionAvail().x;
  const float substrateWidth =
      std::max(frame.em * (compact ? 7.0f : 10.0f),
               available - layoutWidth - resetWidth - furniture);
  drawSubstrateCombo(state, substrateWidth);
  ImGui::SameLine();
  widgets::toolbarSeparator();
  ImGui::SameLine();
  widgets::segmentedIcons("##layout", layouts, layoutTooltips, 2, state.layout,
                          layoutWidth);
  ImGui::SameLine();
  widgets::toolbarSeparator();
  ImGui::SameLine();
  if (compact) {
    if (widgets::iconButton("##reset_filters", icons::Icon::Close,
                            ImVec2(resetWidth, frame.control), false,
                            "Reset filters")) {
      resetFilters(state);
    }
  } else if (widgets::actionButton("##reset_filters", icons::Icon::Close, "Reset",
                                   ImVec2(resetWidth, frame.control))) {
    resetFilters(state);
  }
  widgets::endToolbar();
}

void drawCombinedToolbar(ToolboxState& state, const layout::Frame& frame) {
  const icons::Icon layouts[] = {icons::Icon::Grid, icons::Icon::List};
  const char* layoutTooltips[] = {"Card grid", "Dense list"};
  const float layoutWidth = frame.control * 2.25f;
  const float resetWidth = widgets::actionButtonWidth(icons::Icon::Close, "Reset");
  const float substrateWidth = frame.em * 12.0f;

  widgets::beginToolbar("##toolbox_toolbar");
  // See drawFilterToolbar: five SameLine spacings here, plus the right inset
  // and the two separators.
  const float unit = style::metrics().gap;
  const float furniture = unit * 0.5f + unit * 0.5f * 5.0f + unit * 2.0f;
  const float available = ImGui::GetContentRegionAvail().x;
  const float searchWidth =
      std::max(frame.em * 8.0f,
               available - substrateWidth - layoutWidth - resetWidth - furniture);
  ImGui::SetNextItemWidth(searchWidth);
  widgets::stringInputWithHint("##query", "Search reactions, substrates, reagents",
                               state.query);
  ImGui::SameLine();
  drawSubstrateCombo(state, substrateWidth);
  ImGui::SameLine();
  widgets::toolbarSeparator();
  ImGui::SameLine();
  widgets::segmentedIcons("##layout", layouts, layoutTooltips, 2, state.layout,
                          layoutWidth);
  ImGui::SameLine();
  widgets::toolbarSeparator();
  ImGui::SameLine();
  if (widgets::actionButton("##reset_filters", icons::Icon::Close, "Reset",
                            ImVec2(resetWidth, frame.control))) {
    resetFilters(state);
  }
  widgets::endToolbar();
}

// Does the one-row toolbar actually fit? Its items have fixed widths, so on a
// narrow page -- or a normal page at a high display scale, which is the same
// thing measured in ems -- the run overflows and the last control is pushed off
// the end. When it does not fit, the caller splits it into two rows.
bool combinedToolbarFits(const layout::Frame& frame) {
  const float unit = style::metrics().gap;
  const float furniture = unit * 0.5f * 6.0f + unit * 2.0f;
  const float fixed = frame.em * 8.0f      // search, at its floor
                      + frame.em * 12.0f   // substrate combo
                      + frame.control * 2.25f  // layout toggle
                      + widgets::actionButtonWidth(icons::Icon::Close, "Reset");
  return fixed + furniture <= frame.size.x;
}

void drawRankedBarsCard(ToolboxState& state,
                        const std::vector<charts::BarRow>& rows, ImVec2 size) {
  if (!widgets::beginCard("##library_composition", size, style::col::BgSurface,
                          ImGuiWindowFlags_NoScrollbar |
                              ImGuiWindowFlags_NoScrollWithMouse)) {
    return;
  }
  widgets::cardHeader(icons::Icon::ChartBars, "Library composition",
                      "Reaction families", style::col::Data);
  const layout::Frame chartFrame = layout::measure();
  const float weight[] = {1.0f};
  float chartHeight[1] = {};
  layout::distribute(chartFrame.size.y, weight, nullptr, 1, 0.0f, chartHeight);
  const int clicked = charts::rankedBars(
      "##reaction_type_histogram", rows.data(), static_cast<int>(rows.size()),
      ImVec2(chartFrame.size.x, chartHeight[0]));
  if (clicked >= 0 && clicked < static_cast<int>(state.types.size())) {
    const std::string& type = state.types[static_cast<std::size_t>(clicked)];
    state.selectedType = state.selectedType == type ? std::string() : type;
  }
  widgets::endCard();
}

void drawArityCard(const charts::StackSegment* segments, int segmentCount,
                   const std::string& total, ImVec2 size) {
  if (!widgets::beginCard("##reactant_count", size, style::col::BgSurface,
                          ImGuiWindowFlags_NoScrollbar |
                              ImGuiWindowFlags_NoScrollWithMouse)) {
    return;
  }
  widgets::cardHeader(icons::Icon::Layers, "Reactant count",
                      "Inputs per template", style::col::Data);
  const layout::Frame chartFrame = layout::measure();
  const float weight[] = {1.0f};
  float chartHeight[1] = {};
  layout::distribute(chartFrame.size.y, weight, nullptr, 1, 0.0f, chartHeight);
  charts::donut("##reaction_arity", segments, segmentCount, total.c_str(),
                "templates", ImVec2(chartFrame.size.x, chartHeight[0]));
  widgets::endCard();
}

void drawOverview(ToolboxState& state, const std::vector<charts::BarRow>& rows,
                  const charts::StackSegment* segments, int segmentCount,
                  const std::string& total, ImVec2 size) {
  const layout::Frame frame = layout::measure(size);
  const int columns = layout::columnsThatFit(frame, 18.0f);
  if (columns >= 2) {
    const float width = layout::columnWidth(frame, 2);
    drawRankedBarsCard(state, rows, ImVec2(width, size.y));
    ImGui::SameLine(0.0f, frame.gap);
    drawArityCard(segments, segmentCount, total, ImVec2(width, size.y));
    return;
  }

  const float weights[] = {3.0f, 1.0f};
  const float minimums[] = {frame.row * 4.0f, frame.row * 3.0f};
  float heights[2] = {};
  layout::distribute(size.y, weights, minimums, 2, frame.gap, heights);
  const float firstRowY = ImGui::GetCursorPosY();
  drawRankedBarsCard(state, rows, ImVec2(size.x, heights[0]));
  layout::nextRow(firstRowY + heights[0] + frame.gap);
  drawArityCard(segments, segmentCount, total, ImVec2(size.x, heights[1]));
}

void drawReactionTable(const std::vector<const rxn::ReactionTemplate*>& reactions,
                       ToolboxState& state) {
  const widgets::Column columns[] = {
      {"Reaction", false, true, nullptr, 14.0f},
      {"Type", false, false, nullptr, 8.0f},
      {"Substrate class", false, false, nullptr, 10.0f},
      {"Key reagent", false, false, nullptr, 10.0f},
  };
  if (!widgets::beginDataTable("##reaction_list", columns, 4,
                               ImVec2(ImGui::GetContentRegionAvail().x, 0.0f))) {
    return;
  }

  for (const rxn::ReactionTemplate* reaction : reactions) {
    ImGui::PushID(reaction->id.c_str());
    const bool selected = state.selectedReactionId == reaction->id;
    widgets::dataRow(selected ? style::col::Accent : ImVec4(0.0f, 0.0f, 0.0f, 0.0f));
    widgets::dataCell(reaction->name.c_str());
    if (ImGui::IsItemClicked(ImGuiMouseButton_Left)) {
      state.selectedReactionId = reaction->id;
    }
    const std::string type = rxn::reactionType(*reaction);
    widgets::dataCell(type.c_str());
    widgets::dataCell(valueOrDash(reaction->substrate));
    widgets::dataCell(reaction->reagents.empty() ? "—"
                                                 : reaction->reagents.front().c_str());
    ImGui::PopID();
  }
  widgets::endDataTable();
}

void drawDetailPane(const rxn::ReactionTemplate& reaction, ImVec2 size) {
  if (!widgets::beginCard("##selected_reaction", size, style::col::BgSurface,
                          ImGuiWindowFlags_AlwaysVerticalScrollbar)) {
    return;
  }
  const std::string type = rxn::reactionType(reaction);
  widgets::cardHeader(icons::Icon::Reaction, reaction.name.c_str(), type.c_str(),
                      style::col::Accent);
  drawReactionDetails(reaction);
  widgets::endCard();
}

void drawGrid(const std::vector<const rxn::ReactionTemplate*>& reactions) {
  std::size_t groupBegin = 0;
  while (groupBegin < reactions.size()) {
    const std::string type = rxn::reactionType(*reactions[groupBegin]);
    std::size_t groupEnd = groupBegin + 1;
    while (groupEnd < reactions.size() &&
           rxn::reactionType(*reactions[groupEnd]) == type) {
      ++groupEnd;
    }

    const std::string heading = sentenceCase(type);
    widgets::sectionHeader(heading.c_str(), style::col::Data);
    const layout::Frame frame = layout::measure();
    int columnCount = layout::columnsThatFit(frame, 24.0f);
    const int densityCap = frame.density == layout::Density::Roomy ? 4 : 3;
    columnCount = std::clamp(columnCount, 1, densityCap);
    const float cardWidth = layout::columnWidth(frame, columnCount);

    for (std::size_t index = groupBegin; index < groupEnd; ++index) {
      const int column = static_cast<int>(index - groupBegin) % columnCount;
      if (column != 0) ImGui::SameLine(0.0f, frame.gap);
      ImGui::BeginGroup();
      drawReactionCard(*reactions[index], type, cardWidth);
      ImGui::EndGroup();
    }
    groupBegin = groupEnd;
  }
}

}  // namespace

void drawToolbox(AppState& st) {
  (void)st;
  static ToolboxState state;

  const std::vector<rxn::ReactionTemplate>* knowledge = nullptr;
  try {
    knowledge = &rxn::knowledgeBase();
  } catch (const std::exception& error) {
    widgets::sectionHeader("Reaction library", style::col::Danger);
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

  widgets::sectionHeader("Reaction library", style::col::Data);
  const layout::Frame pageFrame = layout::measure();
  // Two independent questions, both of which used to be answered by "is the
  // page compact": whether the toolbar run fits on one row, and whether the
  // overview charts fit at all.
  const bool compact =
      pageFrame.density == layout::Density::Compact || !combinedToolbarFits(pageFrame);

  std::vector<std::string> histogramAnnotations(state.types.size());
  std::vector<charts::BarRow> histogramRows;
  histogramRows.reserve(state.types.size());
  for (std::size_t i = 0; i < state.types.size(); ++i) {
    histogramAnnotations[i] = std::to_string(state.typeCounts[i]);
    const bool selected = state.selectedType == state.types[i];
    charts::BarRow row;
    row.label = state.types[i].c_str();
    row.value = static_cast<double>(state.typeCounts[i]);
    row.annotation = histogramAnnotations[i].c_str();
    row.accent = selected ? style::col::Accent : style::col::Data;
    row.selected = selected;
    histogramRows.push_back(row);
  }

  int unaryCount = 0;
  int binaryCount = 0;
  int otherArityCount = 0;
  for (const rxn::ReactionTemplate& reaction : *knowledge) {
    if (reaction.arity == 1) {
      ++unaryCount;
    } else if (reaction.arity == 2) {
      ++binaryCount;
    } else {
      ++otherArityCount;
    }
  }
  charts::StackSegment aritySegments[3];
  int aritySegmentCount = 0;
  if (unaryCount > 0) {
    aritySegments[aritySegmentCount++] =
        {"One reactant", static_cast<double>(unaryCount), style::col::Data};
  }
  if (binaryCount > 0) {
    aritySegments[aritySegmentCount++] =
        {"Two reactants", static_cast<double>(binaryCount), style::col::DataBright};
  }
  if (otherArityCount > 0) {
    aritySegments[aritySegmentCount++] =
        {"Other", static_cast<double>(otherArityCount), style::col::DataDim};
  }
  const std::string libraryTotal = std::to_string(knowledge->size());

  const float toolbarBand = pageFrame.control + style::metrics().gap;
  const float toolbarHeight =
      compact ? toolbarBand * 2.0f + pageFrame.gap : toolbarBand;
  const bool overviewColumns = layout::columnsThatFit(pageFrame, 18.0f) >= 2;
  const float overviewHeight =
      pageFrame.row * (static_cast<float>(state.types.size()) + 2.0f) +
      (overviewColumns ? 0.0f : pageFrame.row * 3.0f + pageFrame.gap);
  // The overview cards are supplementary; the search and the results are the
  // page. On a short page -- a small window, or a normal one at a high display
  // scale -- the cards would be squeezed until they clipped their own rows, and
  // half a histogram is worse than none, so they are dropped instead.
  const float page = layout::pageHeight();
  const bool showOverview =
      overviewHeight + pageFrame.gap + toolbarHeight <= page * 0.55f;
  const float headerMinimum =
      showOverview ? overviewHeight + pageFrame.gap + toolbarHeight : toolbarHeight;
  const float pageWeights[] = {0.0f, 1.0f};
  const float pageMinimums[] = {headerMinimum, pageFrame.row * 5.0f};
  float pageRows[2] = {};
  layout::distribute(page, pageWeights, pageMinimums, 2, pageFrame.gap, pageRows);

  const float headerY = ImGui::GetCursorPosY();
  if (ImGui::BeginChild("##toolbox_header", ImVec2(pageFrame.size.x, pageRows[0]),
                        ImGuiChildFlags_None,
                        ImGuiWindowFlags_NoScrollbar |
                            ImGuiWindowFlags_NoScrollWithMouse)) {
    const layout::Frame headerFrame = layout::measure();
    if (showOverview) {
      const float headerWeights[] = {1.0f, 0.0f};
      const float headerMinimums[] = {headerFrame.row * 4.0f, toolbarHeight};
      float headerRows[2] = {};
      layout::distribute(headerFrame.size.y, headerWeights, headerMinimums, 2,
                         headerFrame.gap, headerRows);
      const float overviewY = ImGui::GetCursorPosY();
      drawOverview(state, histogramRows, aritySegments, aritySegmentCount, libraryTotal,
                   ImVec2(headerFrame.size.x, headerRows[0]));
      layout::nextRow(overviewY + headerRows[0] + headerFrame.gap);
    }
    if (compact) {
      const float searchY = ImGui::GetCursorPosY();
      drawSearchToolbar(state, headerFrame);
      const float searchBand = headerFrame.control + style::metrics().gap;
      layout::nextRow(searchY + searchBand + headerFrame.gap);
      drawFilterToolbar(state, headerFrame, true);
    } else {
      drawCombinedToolbar(state, headerFrame);
    }
  }
  ImGui::EndChild();
  layout::nextRow(headerY + pageRows[0] + pageFrame.gap);

  auto matchesFilters = [&](const rxn::ReactionTemplate& reaction) {
    return (state.selectedType.empty() ||
            rxn::reactionType(reaction) == state.selectedType) &&
           (state.selectedSubstrate.empty() ||
            reaction.substrate == state.selectedSubstrate) &&
           matchesQuery(reaction, state.query);
  };

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

  const rxn::ReactionTemplate* selectedReaction = nullptr;
  for (const rxn::ReactionTemplate* reaction : visible) {
    if (reaction->id == state.selectedReactionId) {
      selectedReaction = reaction;
      break;
    }
  }
  if (!state.selectedReactionId.empty() && !selectedReaction) {
    state.selectedReactionId.clear();
  }

  if (ImGui::BeginChild("##toolbox_results", ImVec2(pageFrame.size.x, pageRows[1]),
                        ImGuiChildFlags_None,
                        ImGuiWindowFlags_NoScrollbar |
                            ImGuiWindowFlags_NoScrollWithMouse)) {
    const layout::Frame resultsFrame = layout::measure();
    const float resultWeights[] = {0.0f, 1.0f};
    // The summary band holds a badge inside a child window: a badge is taller
    // than a text line, and the child adds its own padding. A bare `row` here
    // clipped the badge at high display scales.
    const float summaryBand = ImGui::GetFontSize() + style::metrics().gap * 0.55f +
                              ImGui::GetStyle().WindowPadding.y * 2.0f;
    const float resultMinimums[] = {summaryBand, resultsFrame.row * 3.0f};
    float resultRows[2] = {};
    layout::distribute(resultsFrame.size.y, resultWeights, resultMinimums, 2,
                       resultsFrame.gap, resultRows);

    const float summaryY = ImGui::GetCursorPosY();
    if (ImGui::BeginChild("##result_summary", ImVec2(resultsFrame.size.x, resultRows[0]),
                          ImGuiChildFlags_None,
                          ImGuiWindowFlags_NoScrollbar |
                              ImGuiWindowFlags_NoScrollWithMouse)) {
      std::set<std::string> visibleTypes;
      for (const rxn::ReactionTemplate* reaction : visible) {
        visibleTypes.insert(rxn::reactionType(*reaction));
      }
      const std::string resultSummary =
          std::to_string(visible.size()) + " reactions · " +
          std::to_string(visibleTypes.size()) + " types";
      widgets::badge(resultSummary.c_str(), style::col::Data);
      if (!state.selectedType.empty()) {
        ImGui::SameLine(0.0f, resultsFrame.gap);
        const std::string typeSummary = "Type: " + state.selectedType;
        widgets::badge(typeSummary.c_str(), style::col::Accent);
      }
    }
    ImGui::EndChild();
    layout::nextRow(summaryY + resultRows[0] + resultsFrame.gap);

    if (ImGui::BeginChild("##results_content",
                          ImVec2(resultsFrame.size.x, resultRows[1]),
                          ImGuiChildFlags_None,
                          ImGuiWindowFlags_NoScrollbar |
                              ImGuiWindowFlags_NoScrollWithMouse)) {
      const layout::Frame contentFrame = layout::measure();
      if (visible.empty()) {
        if (ImGui::BeginChild("##empty_results_scroll", contentFrame.size,
                              ImGuiChildFlags_None)) {
          std::vector<std::string> activeNames;
          if (!state.query.empty()) activeNames.push_back("search \"" + state.query + "\"");
          if (!state.selectedType.empty()) {
            activeNames.push_back("type " + state.selectedType);
          }
          if (!state.selectedSubstrate.empty()) {
            activeNames.push_back("substrate " + state.selectedSubstrate);
          }
          const std::string emptyBody =
              activeNames.empty()
                  ? "The loaded reaction library contains no templates."
                  : "Active filters: " + join(activeNames) +
                        ". Clear them to browse the library.";
          widgets::emptyState(icons::Icon::Search, "No reactions match",
                              emptyBody.c_str());
          if (!activeNames.empty() &&
              widgets::actionButton("##clear_empty_filters", icons::Icon::Close,
                                    "Clear filters",
                                    ImVec2(widgets::actionButtonWidth(
                                               icons::Icon::Close, "Clear filters"),
                                           contentFrame.control),
                                    true)) {
            resetFilters(state);
          }
        }
        ImGui::EndChild();
      } else if (state.layout == 0) {
        if (ImGui::BeginChild("##reaction_grid_scroll", contentFrame.size,
                              ImGuiChildFlags_None)) {
          drawGrid(visible);
        }
        ImGui::EndChild();
      } else {
        const bool sideBySide = selectedReaction &&
                                layout::columnsThatFit(contentFrame, 28.0f) >= 2 &&
                                contentFrame.density != layout::Density::Compact;
        if (sideBySide) {
          const float listWidth = layout::columnWidth(contentFrame, 3, 2);
          const float detailWidth = layout::columnWidth(contentFrame, 3);
          if (ImGui::BeginChild("##reaction_table_scroll",
                                ImVec2(listWidth, contentFrame.size.y),
                                ImGuiChildFlags_Borders,
                                ImGuiWindowFlags_HorizontalScrollbar)) {
            drawReactionTable(visible, state);
          }
          ImGui::EndChild();
          ImGui::SameLine(0.0f, contentFrame.gap);
          drawDetailPane(*selectedReaction,
                         ImVec2(detailWidth, contentFrame.size.y));
        } else if (ImGui::BeginChild("##reaction_table_scroll", contentFrame.size,
                                     ImGuiChildFlags_Borders,
                                     ImGuiWindowFlags_HorizontalScrollbar)) {
          if (selectedReaction &&
              widgets::disclosure("##selected_reaction_details",
                                  selectedReaction->name.c_str(),
                                  "Selected reaction details", true,
                                  icons::Icon::Reaction, style::col::Accent)) {
            ImGui::Indent();
            drawReactionDetails(*selectedReaction);
            ImGui::Unindent();
          }
          drawReactionTable(visible, state);
          ImGui::EndChild();
        } else {
          ImGui::EndChild();
        }
      }
    }
    ImGui::EndChild();
  }
  ImGui::EndChild();
}

}  // namespace chemcad::ui

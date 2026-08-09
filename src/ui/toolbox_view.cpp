#include "imgui.h"

#include <algorithm>
#include <cctype>
#include <exception>
#include <set>
#include <string>
#include <vector>

#include "rxn/kb.hpp"
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
  bool indexed = false;
};

int resizeStringInput(ImGuiInputTextCallbackData* data) {
  if (data->EventFlag != ImGuiInputTextFlags_CallbackResize) return 0;
  auto* value = static_cast<std::string*>(data->UserData);
  value->resize(static_cast<std::size_t>(data->BufTextLen));
  data->Buf = value->data();
  return 0;
}

std::string lowerCopy(const std::string& value) {
  std::string lowered;
  lowered.reserve(value.size());
  for (unsigned char ch : value) lowered.push_back(static_cast<char>(std::tolower(ch)));
  return lowered;
}

bool containsCaseInsensitive(const std::string& text, const std::string& lowercaseNeedle) {
  return std::search(text.begin(), text.end(), lowercaseNeedle.begin(), lowercaseNeedle.end(),
                     [](unsigned char textChar, unsigned char needleChar) {
                       return std::tolower(textChar) == needleChar;
                     }) != text.end();
}

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
  if (containsCaseInsensitive(reaction.name, query) ||
      containsCaseInsensitive(reaction.substrate, query) ||
      containsCaseInsensitive(reaction.notes, query)) {
    return true;
  }
  for (const std::string& reagent : reaction.reagents) {
    if (containsCaseInsensitive(reagent, query)) return true;
  }
  return false;
}

bool toggleChip(const char* label, bool selected) {
  const style::Metrics& metrics = style::metrics();
  const ImGuiStyle& imguiStyle = ImGui::GetStyle();
  const ImVec2 textSize = ImGui::CalcTextSize(label);
  const ImVec2 size(textSize.x + imguiStyle.FramePadding.x * 2.0f, ImGui::GetFrameHeight());
  const ImVec2 min = ImGui::GetCursorScreenPos();

  ImGui::PushID(label);
  ImGui::InvisibleButton("##chip", size);
  const bool clicked = ImGui::IsItemClicked();
  const float hover = widgets::hoverT(ImGui::GetItemID(), ImGui::IsItemHovered());
  ImGui::PopID();

  const ImVec2 max(min.x + size.x, min.y + size.y);
  ImDrawList* drawList = ImGui::GetWindowDrawList();
  const ImU32 fill = selected
                         ? style::u32(style::col::Accent)
                         : style::mix(style::col::BgSurface, style::col::BgRaised, hover);
  const ImU32 border = selected
                           ? style::u32(style::col::AccentActive)
                           : style::mix(style::col::Border, style::col::BorderStrong, hover);
  const ImU32 text = selected ? style::u32(style::col::OnAccent) : style::u32(style::col::TextDim);
  drawList->AddRectFilled(min, max, fill, metrics.radiusSm);
  drawList->AddRect(min, max, border, metrics.radiusSm, 0, metrics.hairline);
  drawList->AddText(ImVec2((min.x + max.x - textSize.x) * 0.5f,
                           (min.y + max.y - textSize.y) * 0.5f),
                    text, label);
  return clicked;
}

void drawChipRow(const char* id, const std::vector<std::string>& values,
                 std::string& selectedValue) {
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
    if (toggleChip(value.c_str(), selectedValue == value)) {
      selectedValue = selectedValue == value ? std::string() : value;
    }
    usedWidth += chipWidth;
  }
  ImGui::PopID();
}

void drawField(const char* label, const std::string& value) {
  ImGui::TextColored(style::col::TextDim, "%s", label);
  ImGui::SameLine();
  ImGui::TextWrapped("%s", value.empty() ? "—" : value.c_str());
}

void drawReactionCard(const rxn::ReactionTemplate& reaction) {
  ImGui::PushID(reaction.id.c_str());
  const float width = ImGui::GetContentRegionAvail().x;
  if (widgets::beginCard("##reaction_card", ImVec2(width, 0.0f), style::col::BgSurface)) {
    const bool headingFont = style::pushFont(style::fonts::semibold());
    ImGui::TextColored(style::col::Text, "%s", reaction.name.c_str());
    style::popFont(headingFont);
    ImGui::SameLine();
    widgets::badge(reaction.substrate.c_str(), style::col::Teal);

    drawField("Reagents", join(reaction.reagents));
    drawField("Conditions", reaction.conditions);
    drawField("Outcome", reaction.outcome);
    const std::string displayedSource = reaction.source == "standard practice"
                                            ? "standard practice (no procedure citation)"
                                            : reaction.source;
    drawField("Source", displayedSource);

    const bool notesOpen = ImGui::TreeNodeEx(
        "Notes", ImGuiTreeNodeFlags_SpanAvailWidth | ImGuiTreeNodeFlags_NoTreePushOnOpen |
                     (reaction.notes.empty() ? ImGuiTreeNodeFlags_Leaf : 0));
    if (ImGui::IsItemHovered() && !reaction.notes.empty()) {
      ImGui::SetTooltip("%s", reaction.notes.c_str());
    }
    if (notesOpen && !reaction.notes.empty()) {
      ImGui::Indent();
      ImGui::TextWrapped("%s", reaction.notes.c_str());
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
    if (widgets::beginCard("##toolbox_error", ImVec2(ImGui::GetContentRegionAvail().x, 0.0f),
                           style::col::BgSurface)) {
      ImGui::TextColored(style::col::Danger, "Reaction library unavailable");
      ImGui::TextWrapped("%s", error.what());
      widgets::endCard();
    }
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
    state.indexed = true;
  }

  widgets::sectionHeader("REACTION TOOLBOX", style::col::Accent);
  ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x);
  ImGuiInputTextFlags inputFlags = ImGuiInputTextFlags_CallbackResize;
  ImGui::InputTextWithHint("##toolbox_search", "Search names, reagents, substrates or notes",
                           state.query.data(), state.query.capacity() + 1, inputFlags,
                           resizeStringInput, &state.query);

  widgets::sectionHeader("REACTION TYPE", style::col::Violet);
  drawChipRow("reaction_types", state.types, state.selectedType);
  widgets::sectionHeader("SUBSTRATE", style::col::Teal);
  // ~70 substrate classes would otherwise eat half the workspace: cap the
  // chip wall at roughly three rows and let it scroll.
  const float substrateH = ImGui::GetFrameHeight() * 3.4f;
  if (ImGui::BeginChild("##substrate_chips", ImVec2(0.0f, substrateH))) {
    drawChipRow("reaction_substrates", state.substrates, state.selectedSubstrate);
  }
  ImGui::EndChild();

  const std::string query = lowerCopy(state.query);
  std::vector<const rxn::ReactionTemplate*> visible;
  visible.reserve(knowledge->size());
  for (const rxn::ReactionTemplate& reaction : *knowledge) {
    if (!state.selectedType.empty() && rxn::reactionType(reaction) != state.selectedType) continue;
    if (!state.selectedSubstrate.empty() && reaction.substrate != state.selectedSubstrate) continue;
    if (!matchesQuery(reaction, query)) continue;
    visible.push_back(&reaction);
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
  ImGui::TextColored(style::col::TextDim, "%zu reactions  ·  %zu types", visible.size(),
                     visibleTypes.size());
  if (!state.query.empty() || !state.selectedType.empty() || !state.selectedSubstrate.empty()) {
    ImGui::SameLine();
    if (widgets::ghostButton("Clear filters")) {
      state.query.clear();
      state.selectedType.clear();
      state.selectedSubstrate.clear();
    }
  }

  if (visible.empty()) {
    if (widgets::beginCard("##toolbox_empty", ImVec2(ImGui::GetContentRegionAvail().x, 0.0f),
                           style::col::BgSurface)) {
      ImGui::TextColored(style::col::TextDim, "No reactions match these filters.");
      widgets::endCard();
    }
    return;
  }

  std::string currentType;
  for (const rxn::ReactionTemplate* reaction : visible) {
    const std::string type = rxn::reactionType(*reaction);
    if (type != currentType) {
      currentType = type;
      const std::string heading = upperCopy(type);
      widgets::sectionHeader(heading.c_str(), style::col::Accent);
    }
    drawReactionCard(*reaction);
    ImGui::Spacing();
  }
}

}  // namespace chemcad::ui

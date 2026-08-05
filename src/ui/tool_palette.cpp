#include <algorithm>
#include <array>

#include "imgui.h"

#include "ui/ui.hpp"

namespace chemcad::ui {
namespace {

struct ToolButton {
  Tool tool;
  const char* glyph;
  const char* name;
  const char* shortcut;
};

constexpr std::array kTools{
    ToolButton{Tool::Select, "->", "Select", "Esc"},
    ToolButton{Tool::Eraser, "X", "Eraser", "E"},
    ToolButton{Tool::Bond, "--", "Bond", "B"},
    ToolButton{Tool::Chain, "/\\/", "Chain", "K"},
    ToolButton{Tool::RingTemplate, "()", "Ring template", "R"},
    ToolButton{Tool::Atom, "A", "Atom", "A"},
    ToolButton{Tool::ChargePlus, "+", "Positive charge", "+"},
    ToolButton{Tool::ChargeMinus, "-", "Negative charge", "-"},
};

template <typename T, size_t N>
void drawEnumCombo(const char* id, const char* preview, T& value,
                   const std::array<std::pair<T, const char*>, N>& choices) {
  ImGui::SetNextItemWidth(-1.0f);
  if (ImGui::BeginCombo(id, preview)) {
    for (const auto& [candidate, label] : choices) {
      const bool selected = value == candidate;
      if (ImGui::Selectable(label, selected)) value = candidate;
      if (selected) ImGui::SetItemDefaultFocus();
    }
    ImGui::EndCombo();
  }
}

}  // namespace

void drawToolPalette(AppState& st) {
  const float width = std::max(32.0f, ImGui::GetContentRegionAvail().x);
  const ImVec4 active = ImGui::GetStyleColorVec4(ImGuiCol_HeaderActive);
  for (const auto& item : kTools) {
    ImGui::PushID(static_cast<int>(item.tool));
    const bool isActive = st.tool == item.tool;
    if (isActive) {
      ImGui::PushStyleColor(ImGuiCol_Button, active);
      ImGui::PushStyleColor(ImGuiCol_ButtonHovered,
                            ImGui::GetStyleColorVec4(ImGuiCol_HeaderHovered));
    }
    if (ImGui::Button(item.glyph, ImVec2(width, 0))) st.tool = item.tool;
    if (isActive) ImGui::PopStyleColor(2);
    if (ImGui::IsItemHovered()) {
      ImGui::SetTooltip("%s (%s)", item.name, item.shortcut);
    }
    ImGui::PopID();
  }

  ImGui::SeparatorText("Bond");
  constexpr std::array bondOrders{
      std::pair{core::BondOrder::Single, "Single"},
      std::pair{core::BondOrder::Double, "Double"},
      std::pair{core::BondOrder::Triple, "Triple"},
      std::pair{core::BondOrder::Aromatic, "Aromatic"},
  };
  const char* orderPreview = bondOrders.front().second;
  for (const auto& [value, label] : bondOrders) {
    if (value == st.currentOrder) orderPreview = label;
  }
  drawEnumCombo("##bond_order", orderPreview, st.currentOrder, bondOrders);

  ImGui::SeparatorText("Stereo");
  constexpr std::array stereoKinds{
      std::pair{core::BondStereo::None, "None"},
      std::pair{core::BondStereo::Wedge, "Wedge"},
      std::pair{core::BondStereo::Hash, "Hash"},
  };
  const char* stereoPreview = stereoKinds.front().second;
  for (const auto& [value, label] : stereoKinds) {
    if (value == st.currentStereo) stereoPreview = label;
  }
  drawEnumCombo("##bond_stereo", stereoPreview, st.currentStereo, stereoKinds);

  ImGui::SeparatorText("Ring");
  constexpr std::array ringKinds{
      std::pair{RingKind::Cyclopropane, "Cyclopropane"},
      std::pair{RingKind::Cyclobutane, "Cyclobutane"},
      std::pair{RingKind::Cyclopentane, "Cyclopentane"},
      std::pair{RingKind::Cyclohexane, "Cyclohexane"},
      std::pair{RingKind::Cycloheptane, "Cycloheptane"},
      std::pair{RingKind::Cyclooctane, "Cyclooctane"},
      std::pair{RingKind::Benzene, "Benzene"},
      std::pair{RingKind::Cyclopentadiene, "Cyclopentadiene"},
      std::pair{RingKind::Naphthalene, "Naphthalene"},
  };
  const char* ringPreview = ringKinds.front().second;
  for (const auto& [value, label] : ringKinds) {
    if (value == st.currentRing) ringPreview = label;
  }
  drawEnumCombo("##ring_kind", ringPreview, st.currentRing, ringKinds);
}

}  // namespace chemcad::ui

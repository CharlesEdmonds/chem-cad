#include <algorithm>
#include <array>

#include "imgui.h"

#include "ui/icons.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"
#include "ui/widgets.hpp"

namespace chemcad::ui {
namespace {

struct ToolEntry {
  Tool tool;
  icons::Icon icon;
  const char* name;
  const char* shortcut;
  const char* hint;
};

constexpr std::array kTools{
    ToolEntry{Tool::Select, icons::Icon::Select, "Select", "Esc",
              "Drag a box around atoms and bonds to select them"},
    ToolEntry{Tool::Eraser, icons::Icon::Eraser, "Eraser", "E",
              "Click an atom or bond to delete it"},
    ToolEntry{Tool::Bond, icons::Icon::Bond, "Bond", "B",
              "Click or drag to draw a bond"},
    ToolEntry{Tool::Chain, icons::Icon::Chain, "Chain", "K",
              "Drag to draw a zig-zag carbon chain"},
    ToolEntry{Tool::RingTemplate, icons::Icon::Ring, "Ring template", "R",
              "Click to stamp the ring chosen below"},
    ToolEntry{Tool::Atom, icons::Icon::Atom, "Atom", "A",
              "Click to place the element picked in the periodic table"},
    ToolEntry{Tool::ChargePlus, icons::Icon::ChargePlus, "Positive charge", "+",
              "Click an atom to raise its formal charge"},
    ToolEntry{Tool::ChargeMinus, icons::Icon::ChargeMinus, "Negative charge", "-",
              "Click an atom to lower its formal charge"},
};

struct OrderEntry {
  core::BondOrder order;
  icons::Icon icon;
  const char* name;
};
constexpr std::array kOrders{
    OrderEntry{core::BondOrder::Single, icons::Icon::BondSingle, "Single bond"},
    OrderEntry{core::BondOrder::Double, icons::Icon::BondDouble, "Double bond"},
    OrderEntry{core::BondOrder::Triple, icons::Icon::BondTriple, "Triple bond"},
    OrderEntry{core::BondOrder::Aromatic, icons::Icon::BondAromatic, "Aromatic bond"},
};

struct StereoEntry {
  core::BondStereo stereo;
  icons::Icon icon;
  const char* name;
  const char* hint;
};
constexpr std::array kStereos{
    StereoEntry{core::BondStereo::None, icons::Icon::StereoNone, "Plain bond",
                "No stereochemistry"},
    StereoEntry{core::BondStereo::Wedge, icons::Icon::StereoWedge, "Wedge",
                "Bond comes out of the plane towards you"},
    StereoEntry{core::BondStereo::Hash, icons::Icon::StereoHash, "Hashed wedge",
                "Bond goes behind the plane"},
};

// Two-column tool grid. Cells are square but capped, so an unusually wide
// dock never yields giant buttons.
void drawToolGrid(AppState& st, float avail) {
  const float spacing = ImGui::GetStyle().ItemSpacing.x;
  const float cap = style::metrics().iconSize * 3.0f;
  const float cell = std::min((avail - spacing) * 0.5f, cap);
  const float gridWidth = cell * 2.0f + spacing;
  const float indent = std::max(0.0f, (avail - gridWidth) * 0.5f);

  for (size_t i = 0; i < kTools.size(); ++i) {
    const ToolEntry& entry = kTools[i];
    if (i % 2 == 0) {
      ImGui::SetCursorPosX(ImGui::GetCursorPosX() + indent);
    } else {
      ImGui::SameLine(0.0f, spacing);
    }
    ImGui::PushID(static_cast<int>(entry.tool));
    char tooltip[192];
    std::snprintf(tooltip, sizeof(tooltip), "%s (%s)\n%s", entry.name, entry.shortcut,
                  entry.hint);
    if (widgets::iconButton("##tool", entry.icon, ImVec2(cell, cell),
                            st.tool == entry.tool, tooltip)) {
      st.tool = entry.tool;
    }
    ImGui::PopID();
  }
}

template <typename T, size_t N, typename IsCurrent, typename Apply>
void drawIconRow(const char* idPrefix, const std::array<T, N>& entries, float avail,
                 icons::Icon T::*iconMember, const char* T::*nameMember,
                 IsCurrent&& isCurrent, Apply&& apply, AppState& st,
                 const char* T::*hintMember = nullptr) {
  const float spacing = ImGui::GetStyle().ItemSpacing.x;
  const float width = (avail - spacing * static_cast<float>(N - 1)) / static_cast<float>(N);
  const float height = style::metrics().iconSize * 1.5f;
  for (size_t i = 0; i < N; ++i) {
    const T& entry = entries[i];
    if (i > 0) ImGui::SameLine(0.0f, spacing);
    ImGui::PushID(idPrefix);
    ImGui::PushID(static_cast<int>(i));
    const char* name = entry.*nameMember;
    const char* hint = hintMember ? entry.*hintMember : nullptr;
    char tooltip[192];
    if (hint)
      std::snprintf(tooltip, sizeof(tooltip), "%s\n%s", name, hint);
    else
      std::snprintf(tooltip, sizeof(tooltip), "%s", name);
    if (widgets::iconButton("##choice", entry.*iconMember, ImVec2(width, height),
                            isCurrent(entry, st), tooltip)) {
      apply(entry, st);
    }
    ImGui::PopID();
    ImGui::PopID();
  }
}

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
  const float avail = std::max(32.0f, ImGui::GetContentRegionAvail().x);

  drawToolGrid(st, avail);

  widgets::sectionHeader("Bond");
  drawIconRow(
      "##order", kOrders, avail, &OrderEntry::icon, &OrderEntry::name,
      [](const OrderEntry& e, const AppState& s) { return s.currentOrder == e.order; },
      [](const OrderEntry& e, AppState& s) { s.currentOrder = e.order; }, st);

  widgets::sectionHeader("Stereo");
  drawIconRow(
      "##stereo", kStereos, avail, &StereoEntry::icon, &StereoEntry::name,
      [](const StereoEntry& e, const AppState& s) { return s.currentStereo == e.stereo; },
      [](const StereoEntry& e, AppState& s) { s.currentStereo = e.stereo; }, st,
      &StereoEntry::hint);

  widgets::sectionHeader("Ring");
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
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "Ring stamped by the Ring template tool (R)");
}

}  // namespace chemcad::ui

#include <algorithm>
#include <array>
#include <string>
#include <utility>

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
    ToolEntry{Tool::Bond, icons::Icon::Bond, "Bond", "B / M",
              "Click an atom to attach CH3; drag to place or connect carbon. M resets a plain methyl-ready bond"},
    ToolEntry{Tool::Chain, icons::Icon::Chain, "Chain", "K",
              "Drag to draw a zig-zag carbon chain"},
    // Tool::RingTemplate has no grid button: the Ring drop-down below both
    // picks the template and arms the tool, so a button would be a second,
    // redundant way to enter the same mode.
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

// The ring picker doubles as the ring tool's activator, so it cannot use a
// plain value-binding combo.
void drawRingPicker(AppState& st) {
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
  const char* preview = ringKinds.front().second;
  for (const auto& [value, label] : ringKinds) {
    if (value == st.currentRing) preview = label;
  }

  // Armed state gets the same accent border the icon rows use, so the palette
  // still shows at a glance which mode the canvas is in.
  const bool armed = st.tool == Tool::RingTemplate;
  if (armed) ImGui::PushStyleColor(ImGuiCol_Border, style::u32(style::col::Accent));
  ImGui::SetNextItemWidth(-1.0f);
  if (ImGui::BeginCombo("##ring_kind", preview)) {
    for (const auto& [value, label] : ringKinds) {
      const bool selected = value == st.currentRing;
      if (ImGui::Selectable(label, selected)) {
        st.currentRing = value;
        st.tool = Tool::RingTemplate;
        st.statusMessage = std::string("Ring template: ") + label;
      }
      if (selected) ImGui::SetItemDefaultFocus();
    }
    ImGui::EndCombo();
  }
  if (armed) ImGui::PopStyleColor();
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "Pick a ring to arm the ring tool (R), then click the canvas");
}

}  // namespace

void drawToolPalette(AppState& st) {
  const float avail = std::max(32.0f, ImGui::GetContentRegionAvail().x);

  drawToolGrid(st, avail);

  widgets::sectionHeader("Bond  ·  M = CH3");
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
  drawRingPicker(st);
}

}  // namespace chemcad::ui

#include <algorithm>
#include <array>
#include <cstdio>
#include <string>

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
    // The ring cell is a split button drawn by drawRingCell: the body arms the
    // tool, the caret opens the ring gallery. The icon comes from kRings via
    // st.currentRing, so this entry only supplies name/shortcut/hint.
    ToolEntry{Tool::RingTemplate, icons::Icon::RingCyclohexane, "Ring", "R",
              "Stamp the selected ring onto the canvas"},
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
    StereoEntry{core::BondStereo::Wavy, icons::Icon::StereoWavy, "Wavy",
                "Stereochemistry unknown or a mixture"},
};

// ---------------------------------------------------------------- rings
// Every stampable ring with its own glyph, so the grid button and the gallery
// always show what will actually be drawn.
struct RingEntry {
  RingKind kind;
  icons::Icon icon;
  const char* name;
};
constexpr std::array kRings{
    RingEntry{RingKind::Cyclopropane, icons::Icon::RingCyclopropane, "Cyclopropane"},
    RingEntry{RingKind::Cyclobutane, icons::Icon::RingCyclobutane, "Cyclobutane"},
    RingEntry{RingKind::Cyclopentane, icons::Icon::RingCyclopentane, "Cyclopentane"},
    RingEntry{RingKind::Cyclohexane, icons::Icon::RingCyclohexane, "Cyclohexane"},
    RingEntry{RingKind::Cycloheptane, icons::Icon::RingCycloheptane, "Cycloheptane"},
    RingEntry{RingKind::Cyclooctane, icons::Icon::RingCyclooctane, "Cyclooctane"},
    RingEntry{RingKind::Benzene, icons::Icon::RingBenzene, "Benzene"},
    RingEntry{RingKind::Cyclopentadiene, icons::Icon::RingCyclopentadiene,
              "Cyclopentadiene"},
    RingEntry{RingKind::Naphthalene, icons::Icon::RingNaphthalene, "Naphthalene"},
};

const RingEntry& ringEntry(RingKind kind) {
  for (const RingEntry& entry : kRings) {
    if (entry.kind == kind) return entry;
  }
  return kRings.front();
}

// Grid of every ring, drawn inside the gallery popup. Each tile shows the
// ring's glyph over its name; the current ring carries the accent.
void drawRingGallery(AppState& st) {
  const style::Metrics& m = style::metrics();
  const float pad = m.gap * 0.9f;
  const float glyphBox = m.iconSize * 1.5f;
  const float labelH = ImGui::GetFontSize();

  float tileW = 0.0f;
  for (const RingEntry& entry : kRings) {
    tileW = std::max(tileW, ImGui::CalcTextSize(entry.name).x);
  }
  tileW = std::max(tileW + pad * 2.0f, glyphBox + pad * 1.5f);
  const float tileH = glyphBox + labelH + pad * 2.0f;
  const float spacing = ImGui::GetStyle().ItemSpacing.x;

  ImDrawList* dl = ImGui::GetWindowDrawList();
  for (size_t i = 0; i < kRings.size(); ++i) {
    const RingEntry& entry = kRings[i];
    if (i % 3 != 0) ImGui::SameLine(0.0f, spacing);
    ImGui::PushID(static_cast<int>(entry.kind));
    const ImVec2 min = ImGui::GetCursorScreenPos();
    ImGui::InvisibleButton("##ring_tile", ImVec2(tileW, tileH));
    const bool clicked = ImGui::IsItemClicked();
    const float t = widgets::hoverT(ImGui::GetItemID(), ImGui::IsItemHovered());
    ImGui::PopID();
    const ImVec2 max(min.x + tileW, min.y + tileH);
    const bool current = entry.kind == st.currentRing;

    const ImU32 fill = current ? style::mix(style::col::Accent, style::col::BgRaised, 0.82f)
                               : style::mix(style::col::BgSurface, style::col::BgRaised, t);
    dl->AddRectFilled(min, max, fill, m.radiusMd);
    dl->AddRect(min, max,
                current ? style::u32(style::col::Accent, 0.9f)
                        : style::mix(style::col::Border, style::col::BorderStrong, t),
                m.radiusMd, 0, m.hairline);

    const ImU32 glyphColor =
        current ? style::u32(style::col::Accent)
                : style::mix(style::col::TextDim, style::col::Text, 0.4f + 0.6f * t);
    icons::draw(dl, entry.icon, ImVec2((min.x + max.x) * 0.5f, min.y + pad + glyphBox * 0.5f),
                glyphBox * 0.8f, glyphColor);
    const ImVec2 labelSize = ImGui::CalcTextSize(entry.name);
    dl->AddText(ImVec2((min.x + max.x - labelSize.x) * 0.5f, max.y - pad - labelH),
                style::u32(current ? style::col::Text : style::col::TextDim), entry.name);

    if (clicked) {
      st.currentRing = entry.kind;
      st.tool = Tool::RingTemplate;
      st.statusMessage = std::string("Ring template: ") + entry.name;
      ImGui::CloseCurrentPopup();
    }
  }
}

// The ring tool's grid cell: the body arms Tool::RingTemplate with the current
// ring (and shows that ring's glyph), the caret in the corner opens the
// gallery popup. One invisible button with a manual sub-rect hit test for the
// caret -- a second overlapping item would corrupt CursorPosPrevLine and
// break the grid's SameLine stride.
void drawRingCell(AppState& st, ImVec2 cell) {
  const style::Metrics& m = style::metrics();
  const RingEntry& ring = ringEntry(st.currentRing);
  const bool armed = st.tool == Tool::RingTemplate;
  const ImVec2 min = ImGui::GetCursorScreenPos();

  ImGui::InvisibleButton("##ring", cell);
  const bool clicked = ImGui::IsItemClicked();
  const bool hovered = ImGui::IsItemHovered();
  const float t = widgets::hoverT(ImGui::GetItemID(), hovered);
  const ImVec2 max(min.x + cell.x, min.y + cell.y);

  const float caretW = std::max(14.0f, cell.x * 0.36f);
  const float caretH = std::max(12.0f, cell.y * 0.32f);
  const ImVec2 mouse = ImGui::GetMousePos();
  const bool inCaret = hovered && mouse.x >= max.x - caretW && mouse.y >= max.y - caretH;

  ImDrawList* dl = ImGui::GetWindowDrawList();
  ImU32 glyphColor;
  if (armed) {
    dl->AddRectFilled(min, max, style::mix(style::col::Accent, style::col::AccentHover, t),
                      m.radiusMd);
    dl->AddRect(min, max, style::u32(style::col::AccentActive), m.radiusMd, 0, m.hairline);
    glyphColor = style::u32(style::col::OnAccent);
  } else {
    dl->AddRectFilled(min, max, style::mix(style::col::BgSurface, style::col::BgRaised, t),
                      m.radiusMd);
    dl->AddRect(min, max, style::mix(style::col::Border, style::col::BorderStrong, t),
                m.radiusMd, 0, m.hairline);
    glyphColor = style::mix(style::col::TextDim, style::col::Text, 0.35f + 0.65f * t);
  }

  const ImVec2 centre((min.x + max.x) * 0.5f, (min.y + max.y) * 0.5f);
  const float glyph = std::min(m.iconSize * 1.25f, std::min(cell.x, cell.y) * 0.80f);
  icons::draw(dl, ring.icon, centre, glyph, glyphColor);

  // Caret: a filled triangle tucked into the corner, brightening when the
  // pointer is inside its zone so it reads as a separate click target.
  const ImU32 caretColor =
      armed ? style::u32(style::col::OnAccent, inCaret ? 1.0f : 0.65f)
            : style::mix(style::col::TextFaint, style::col::Text, inCaret ? 1.0f : 0.0f);
  const float cx = max.x - m.gap * 0.7f;
  const float cy = max.y - m.gap * 0.7f;
  const float s = std::max(3.5f, m.gap * 0.55f);
  dl->AddTriangleFilled(ImVec2(cx - s, cy - s), ImVec2(cx, cy), ImVec2(cx - s, cy), caretColor);

  if (clicked && !inCaret) st.tool = Tool::RingTemplate;
  if (hovered) {
    if (inCaret) {
      ImGui::SetTooltip("%s", "Choose a different ring");
    } else {
      char tip[192];
      std::snprintf(tip, sizeof(tip), "%s (R)\nClick to arm the ring tool", ring.name);
      ImGui::SetTooltip("%s", tip);
    }
  }

  if (clicked && inCaret) ImGui::OpenPopup("##ring_gallery");
  ImGui::SetNextWindowPos(ImVec2(min.x, max.y + m.gap * 0.5f));
  if (ImGui::BeginPopup("##ring_gallery")) {
    drawRingGallery(st);
    ImGui::EndPopup();
  }
}

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
    if (entry.tool == Tool::RingTemplate) {
      drawRingCell(st, ImVec2(cell, cell));
    } else {
      char tooltip[192];
      std::snprintf(tooltip, sizeof(tooltip), "%s (%s)\n%s", entry.name, entry.shortcut,
                    entry.hint);
      if (widgets::iconButton("##tool", entry.icon, ImVec2(cell, cell),
                              st.tool == entry.tool, tooltip)) {
        st.tool = entry.tool;
      }
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
}

}  // namespace chemcad::ui

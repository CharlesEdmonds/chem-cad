// The left dock's tool column: a two-column grid of tool buttons. The Bond,
// Ring and Atom tools are split buttons -- the body arms the tool with the
// current choice's own glyph, the caret opens a gallery popup where every
// option has an accurate icon. Everything is procedural, no image assets.

#include <algorithm>
#include <array>
#include <cstdio>
#include <string>

#include "imgui.h"

#include "chem/bridge.hpp"
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
    // Bond, Ring and Atom are split buttons drawn by drawSplitCell; their
    // icons come from the current order/ring/element, so these entries only
    // supply name/shortcut/hint.
    ToolEntry{Tool::Bond, icons::Icon::Bond, "Bond", "B / M",
              "Click an atom to attach CH3; drag to place or connect carbon. M resets a plain methyl-ready bond"},
    ToolEntry{Tool::Chain, icons::Icon::Chain, "Chain", "K",
              "Drag to draw a zig-zag carbon chain"},
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

// The elements a working chemist reaches for ninety percent of the time;
// the full grid stays in the Periodic Table panel.
struct QuickElement {
  uint8_t z;
  const char* name;
};
constexpr std::array kQuickElements{
    QuickElement{1, "Hydrogen"},  QuickElement{6, "Carbon"},   QuickElement{7, "Nitrogen"},
    QuickElement{8, "Oxygen"},    QuickElement{9, "Fluorine"}, QuickElement{15, "Phosphorus"},
    QuickElement{16, "Sulfur"},   QuickElement{17, "Chlorine"}, QuickElement{35, "Bromine"},
    QuickElement{53, "Iodine"},   QuickElement{5, "Boron"},    QuickElement{14, "Silicon"},
};

// --------------------------------------------------------------- tile grid
// One gallery glyph: an icon, or a short text (element symbols) when set.
struct TileGlyph {
  icons::Icon icon = icons::Icon::Atom;
  const char* text = nullptr;
};

// A uniform grid of icon-over-label tiles inside a gallery popup. Returns
// the clicked entry's index, or -1.
template <typename T, size_t N, typename GlyphOf, typename NameOf, typename SelectedOf>
int drawTileGrid(const char* idPrefix, const std::array<T, N>& entries, int columns,
                 GlyphOf glyphOf, NameOf nameOf, SelectedOf selectedOf) {
  const style::Metrics& m = style::metrics();
  const float pad = m.gap * 0.9f;
  const float glyphBox = m.iconSize * 1.5f;
  const float labelH = ImGui::GetFontSize();

  float tileW = 0.0f;
  for (const T& entry : entries) tileW = std::max(tileW, ImGui::CalcTextSize(nameOf(entry)).x);
  tileW = std::max(tileW + pad * 2.0f, glyphBox + pad * 1.5f);
  const float tileH = glyphBox + labelH + pad * 2.0f;
  const float spacing = ImGui::GetStyle().ItemSpacing.x;

  int clicked = -1;
  ImDrawList* dl = ImGui::GetWindowDrawList();
  ImGui::PushID(idPrefix);
  for (size_t i = 0; i < entries.size(); ++i) {
    const T& entry = entries[i];
    if (static_cast<int>(i) % columns != 0) ImGui::SameLine(0.0f, spacing);
    ImGui::PushID(static_cast<int>(i));
    const ImVec2 min = ImGui::GetCursorScreenPos();
    ImGui::InvisibleButton("##tile", ImVec2(tileW, tileH));
    const bool wasClicked = ImGui::IsItemClicked();
    const float t = widgets::hoverT(ImGui::GetItemID(), ImGui::IsItemHovered());
    ImGui::PopID();
    const ImVec2 max(min.x + tileW, min.y + tileH);
    const bool selected = selectedOf(entry);

    const ImU32 fill = selected ? style::mix(style::col::Accent, style::col::BgRaised, 0.82f)
                                : style::mix(style::col::BgSurface, style::col::BgRaised, t);
    dl->AddRectFilled(min, max, fill, m.radiusMd);
    dl->AddRect(min, max,
                selected ? style::u32(style::col::Accent, 0.9f)
                         : style::mix(style::col::Border, style::col::BorderStrong, t),
                m.radiusMd, 0, m.hairline);

    const ImU32 glyphColor =
        selected ? style::u32(style::col::Accent)
                 : style::mix(style::col::TextDim, style::col::Text, 0.4f + 0.6f * t);
    const TileGlyph glyph = glyphOf(entry);
    const ImVec2 glyphCentre((min.x + max.x) * 0.5f, min.y + pad + glyphBox * 0.5f);
    if (glyph.text) {
      const float textSize = glyphBox * 0.55f;
      const ImVec2 w = ImGui::CalcTextSize(glyph.text);
      dl->AddText(style::fonts::mono(), textSize,
                  ImVec2(glyphCentre.x - w.x * (textSize / ImGui::GetFontSize()) * 0.5f,
                         glyphCentre.y - textSize * 0.5f),
                  glyphColor, glyph.text);
    } else {
      icons::draw(dl, glyph.icon, glyphCentre, glyphBox * 0.8f, glyphColor);
    }
    const ImVec2 labelSize = ImGui::CalcTextSize(nameOf(entry));
    dl->AddText(ImVec2((min.x + max.x - labelSize.x) * 0.5f, max.y - pad - labelH),
                style::u32(selected ? style::col::Text : style::col::TextDim), nameOf(entry));

    if (wasClicked) clicked = static_cast<int>(i);
  }
  ImGui::PopID();
  return clicked;
}

void galleryCaption(const char* text) {
  ImGui::TextDisabled("%s", text);
}

// ------------------------------------------------------------- split cell
// A grid cell that is two controls in one: the body arms its tool with the
// current choice's glyph on it, the caret in the corner opens a gallery
// popup. One invisible button with a manual sub-rect hit test for the
// caret -- a second overlapping item would corrupt CursorPosPrevLine and
// break the grid's SameLine stride.
template <typename Arm, typename Popup>
void drawSplitCell(ImVec2 cell, TileGlyph bodyGlyph, bool armed, const char* bodyTip,
                   const char* caretTip, const char* popupId, Arm&& onArm,
                   Popup&& popupContent) {
  const style::Metrics& m = style::metrics();
  const ImVec2 min = ImGui::GetCursorScreenPos();

  ImGui::InvisibleButton("##split", cell);
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
  if (bodyGlyph.text) {
    const float textSize = std::min(m.iconSize * 0.85f, cell.y * 0.5f);
    const ImVec2 w = ImGui::CalcTextSize(bodyGlyph.text);
    dl->AddText(style::fonts::mono(), textSize,
                ImVec2(centre.x - w.x * (textSize / ImGui::GetFontSize()) * 0.5f,
                       centre.y - textSize * 0.5f),
                glyphColor, bodyGlyph.text);
  } else {
    const float glyph = std::min(m.iconSize * 1.25f, std::min(cell.x, cell.y) * 0.80f);
    icons::draw(dl, bodyGlyph.icon, centre, glyph, glyphColor);
  }

  const ImU32 caretColor =
      armed ? style::u32(style::col::OnAccent, inCaret ? 1.0f : 0.65f)
            : style::mix(style::col::TextFaint, style::col::Text, inCaret ? 1.0f : 0.0f);
  const float cx = max.x - m.gap * 0.7f;
  const float cy = max.y - m.gap * 0.7f;
  const float s = std::max(3.5f, m.gap * 0.55f);
  dl->AddTriangleFilled(ImVec2(cx - s, cy - s), ImVec2(cx, cy), ImVec2(cx - s, cy), caretColor);

  if (clicked && !inCaret) onArm();
  if (hovered) ImGui::SetTooltip("%s", inCaret ? caretTip : bodyTip);

  if (clicked && inCaret) ImGui::OpenPopup(popupId);
  ImGui::SetNextWindowPos(ImVec2(min.x, max.y + m.gap * 0.5f));
  if (ImGui::BeginPopup(popupId)) {
    popupContent();
    ImGui::EndPopup();
  }
}

// ---------------------------------------------------------------- galleries
void drawRingGallery(AppState& st) {
  const int hit = drawTileGrid(
      "##rings", kRings, 3, [](const RingEntry& e) { return TileGlyph{e.icon, nullptr}; },
      [](const RingEntry& e) { return e.name; },
      [&](const RingEntry& e) { return e.kind == st.currentRing; });
  if (hit >= 0) {
    st.currentRing = kRings[static_cast<size_t>(hit)].kind;
    st.tool = Tool::RingTemplate;
    st.statusMessage = std::string("Ring template: ") + kRings[static_cast<size_t>(hit)].name;
    ImGui::CloseCurrentPopup();
  }
}

void drawBondGallery(AppState& st) {
  galleryCaption("Bond order");
  const int order = drawTileGrid(
      "##orders", kOrders, 4, [](const OrderEntry& e) { return TileGlyph{e.icon, nullptr}; },
      [](const OrderEntry& e) { return e.name; },
      [&](const OrderEntry& e) {
        return st.currentStereo == core::BondStereo::None && st.currentOrder == e.order;
      });
  ImGui::Spacing();
  galleryCaption("Stereochemistry");
  const int stereo = drawTileGrid(
      "##stereos", kStereos, 4, [](const StereoEntry& e) { return TileGlyph{e.icon, nullptr}; },
      [](const StereoEntry& e) { return e.name; },
      [&](const StereoEntry& e) { return st.currentStereo == e.stereo; });
  ImGui::Spacing();
  ImGui::TextDisabled("%s", "M resets a plain single bond (CH3-ready).");

  if (order >= 0) {
    st.currentOrder = kOrders[static_cast<size_t>(order)].order;
    st.tool = Tool::Bond;
    st.statusMessage = kOrders[static_cast<size_t>(order)].name;
    ImGui::CloseCurrentPopup();
  } else if (stereo >= 0) {
    st.currentStereo = kStereos[static_cast<size_t>(stereo)].stereo;
    st.tool = Tool::Bond;
    st.statusMessage = kStereos[static_cast<size_t>(stereo)].name;
    ImGui::CloseCurrentPopup();
  }
}

void drawAtomGallery(AppState& st) {
  galleryCaption("Common elements");
  const int hit = drawTileGrid(
      "##elements", kQuickElements, 4,
      [](const QuickElement& e) { return TileGlyph{icons::Icon::Atom, chem::symbolFor(e.z)}; },
      [](const QuickElement& e) { return e.name; },
      [&](const QuickElement& e) { return st.currentElement == e.z; });
  ImGui::Spacing();
  ImGui::TextDisabled("%s", "Every element lives in the Periodic Table panel.");
  if (hit >= 0) {
    st.currentElement = kQuickElements[static_cast<size_t>(hit)].z;
    st.tool = Tool::Atom;
    st.statusMessage = std::string("Element: ") + kQuickElements[static_cast<size_t>(hit)].name;
    ImGui::CloseCurrentPopup();
  }
}

// ------------------------------------------------------------- grid cells
TileGlyph bondGlyph(const AppState& st) {
  if (st.currentStereo != core::BondStereo::None) {
    for (const StereoEntry& e : kStereos) {
      if (e.stereo == st.currentStereo) return TileGlyph{e.icon, nullptr};
    }
  }
  for (const OrderEntry& e : kOrders) {
    if (e.order == st.currentOrder) return TileGlyph{e.icon, nullptr};
  }
  return TileGlyph{icons::Icon::Bond, nullptr};
}

const char* bondStyleName(const AppState& st) {
  if (st.currentStereo != core::BondStereo::None) {
    for (const StereoEntry& e : kStereos) {
      if (e.stereo == st.currentStereo) return e.name;
    }
  }
  for (const OrderEntry& e : kOrders) {
    if (e.order == st.currentOrder) return e.name;
  }
  return "Bond";
}

void drawRingCell(AppState& st, ImVec2 cell) {
  char tip[192];
  std::snprintf(tip, sizeof(tip), "%s (R)\nClick to arm the ring tool", ringEntry(st.currentRing).name);
  drawSplitCell(
      cell, TileGlyph{ringEntry(st.currentRing).icon, nullptr}, st.tool == Tool::RingTemplate,
      tip, "Choose a different ring", "##ring_gallery", [&st] { st.tool = Tool::RingTemplate; },
      [&st] { drawRingGallery(st); });
}

void drawBondCell(AppState& st, ImVec2 cell) {
  char tip[192];
  std::snprintf(tip, sizeof(tip), "%s (B / M)\nClick to arm the bond tool", bondStyleName(st));
  drawSplitCell(
      cell, bondGlyph(st), st.tool == Tool::Bond, tip, "Choose bond type and stereo",
      "##bond_gallery", [&st] { st.tool = Tool::Bond; }, [&st] { drawBondGallery(st); });
}

void drawAtomCell(AppState& st, ImVec2 cell) {
  char tip[192];
  std::snprintf(tip, sizeof(tip), "%s (A)\nClick to arm the atom tool",
                chem::symbolFor(st.currentElement));
  drawSplitCell(
      cell, TileGlyph{icons::Icon::Atom, chem::symbolFor(st.currentElement)},
      st.tool == Tool::Atom, tip, "Choose a different element", "##atom_gallery",
      [&st] { st.tool = Tool::Atom; }, [&st] { drawAtomGallery(st); });
}

// Two-column tool grid. Cells are square but capped, so an unusually wide
// dock never yields giant buttons.
void drawToolGrid(AppState& st, float avail) {
  const style::Metrics& m = style::metrics();
  const float spacing = ImGui::GetStyle().ItemSpacing.x;
  const float cell = std::min((avail - spacing) * 0.5f, m.iconSize * 3.0f);
  const float gridWidth = cell * 2.0f + spacing;
  const float indent = std::max(0.0f, (avail - gridWidth) * 0.5f);

  for (size_t i = 0; i < kTools.size(); ++i) {
    const ToolEntry& entry = kTools[i];
    // Indent() is the sanctioned way to shift layout; raw SetCursorPosX here
    // trips ImGui's "extends parent boundaries" sanity check during the
    // docked window's hidden measurement frames.
    if (i % 2 == 0) {
      ImGui::Indent(indent);
    } else {
      ImGui::SameLine(0.0f, spacing);
    }

    ImGui::PushID(static_cast<int>(entry.tool));
    if (entry.tool == Tool::RingTemplate) {
      drawRingCell(st, ImVec2(cell, cell));
    } else if (entry.tool == Tool::Bond) {
      drawBondCell(st, ImVec2(cell, cell));
    } else if (entry.tool == Tool::Atom) {
      drawAtomCell(st, ImVec2(cell, cell));
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
    if (i % 2 == 1) ImGui::Unindent(indent);
  }
  // A trailing Dummy keeps the window boundary honest after the grid's
  // SameLine stride (SetCursorPos alone trips ImGui's boundary check).
  ImGui::Dummy(ImVec2(0.0f, ImGui::GetStyle().ItemSpacing.y));
}

}  // namespace

void drawToolPalette(AppState& st) {
  const float avail = std::max(32.0f, ImGui::GetContentRegionAvail().x);
  drawToolGrid(st, avail);
}

}  // namespace chemcad::ui

// Professional structure-tool rail: responsive labelled commands, clear active
// state, and side-opening galleries for compound controls. The rail is one
// column by design, so it remains usable in narrow dock layouts and at high DPI.

#include <algorithm>
#include <array>
#include <cstdio>
#include <string>

#include "imgui.h"

#include "chem/bridge.hpp"
#include "ui/icons.hpp"
#include "ui/element_data.hpp"
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
              "Move, inspect, and box-select structure fragments"},
    ToolEntry{Tool::Eraser, icons::Icon::Eraser, "Eraser", "E",
              "Delete atoms or bonds from the structure"},
    ToolEntry{Tool::Bond, icons::Icon::Bond, "Bond", "B / M",
              "Draw bonds; M restores a plain single bond"},
    ToolEntry{Tool::Chain, icons::Icon::Chain, "Chain", "K",
              "Drag out a carbon zig-zag chain"},
    ToolEntry{Tool::RingTemplate, icons::Icon::RingCyclohexane, "Ring", "R",
              "Stamp a ring template onto the canvas"},
    ToolEntry{Tool::Atom, icons::Icon::Atom, "Atom", "A",
              "Place the selected element"},
    ToolEntry{Tool::ChargePlus, icons::Icon::ChargePlus, "Charge +", "+",
              "Increase an atom's formal charge"},
    ToolEntry{Tool::ChargeMinus, icons::Icon::ChargeMinus, "Charge -", "-",
              "Decrease an atom's formal charge"},
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
};
constexpr std::array kStereos{
    StereoEntry{core::BondStereo::None, icons::Icon::StereoNone, "Plain"},
    StereoEntry{core::BondStereo::Wedge, icons::Icon::StereoWedge, "Wedge"},
    StereoEntry{core::BondStereo::Hash, icons::Icon::StereoHash, "Hashed"},
    StereoEntry{core::BondStereo::Wavy, icons::Icon::StereoWavy, "Wavy"},
};

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

struct QuickElement {
  uint8_t z;
  const char* name;
};
constexpr std::array kQuickElements{
    QuickElement{1, "Hydrogen"},   QuickElement{6, "Carbon"},
    QuickElement{7, "Nitrogen"},   QuickElement{8, "Oxygen"},
    QuickElement{9, "Fluorine"},   QuickElement{15, "Phosphorus"},
    QuickElement{16, "Sulfur"},    QuickElement{17, "Chlorine"},
    QuickElement{35, "Bromine"},   QuickElement{53, "Iodine"},
    QuickElement{5, "Boron"},      QuickElement{14, "Silicon"},
};

struct TileGlyph {
  icons::Icon icon = icons::Icon::Atom;
  const char* text = nullptr;
};

const RingEntry& ringEntry(RingKind kind) {
  for (const RingEntry& entry : kRings) {
    if (entry.kind == kind) return entry;
  }
  return kRings.front();
}

TileGlyph bondGlyph(const AppState& st) {
  if (st.currentStereo != core::BondStereo::None) {
    for (const StereoEntry& entry : kStereos) {
      if (entry.stereo == st.currentStereo) return {entry.icon, nullptr};
    }
  }
  for (const OrderEntry& entry : kOrders) {
    if (entry.order == st.currentOrder) return {entry.icon, nullptr};
  }
  return {icons::Icon::BondSingle, nullptr};
}

const char* bondStyleName(const AppState& st) {
  if (st.currentStereo != core::BondStereo::None) {
    for (const StereoEntry& entry : kStereos) {
      if (entry.stereo == st.currentStereo) return entry.name;
    }
  }
  for (const OrderEntry& entry : kOrders) {
    if (entry.order == st.currentOrder) return entry.name;
  }
  return "Single bond";
}

const ToolEntry& toolEntry(Tool tool) {
  for (const ToolEntry& entry : kTools) {
    if (entry.tool == tool) return entry;
  }
  return kTools.front();
}

TileGlyph toolGlyph(const AppState& st, Tool tool) {
  if (tool == Tool::Bond) return bondGlyph(st);
  if (tool == Tool::RingTemplate) return {ringEntry(st.currentRing).icon, nullptr};
  if (tool == Tool::Atom) return {icons::Icon::Atom, chem::symbolFor(st.currentElement)};
  return {toolEntry(tool).icon, nullptr};
}

std::string toolDetail(const AppState& st, Tool tool) {
  switch (tool) {
    case Tool::Select:
      return "Move and inspect";
    case Tool::Eraser:
      return "Delete structure";
    case Tool::Bond:
      return bondStyleName(st);
    case Tool::Chain:
      return "Carbon zig-zag";
    case Tool::RingTemplate:
      return ringEntry(st.currentRing).name;
    case Tool::Atom:
      return std::string(chem::symbolFor(st.currentElement)) + " element";
    case Tool::ChargePlus:
      return "Increase formal charge";
    case Tool::ChargeMinus:
      return "Decrease formal charge";
  }
  return {};
}

void drawGlyph(ImDrawList* dl, TileGlyph glyph, ImVec2 centre, float size, ImU32 color) {
  if (!glyph.text) {
    icons::draw(dl, glyph.icon, centre, size, color);
    return;
  }
  const float textSize = size * 0.68f;
  const bool mono = style::pushFont(style::fonts::mono());
  const ImVec2 measured = ImGui::CalcTextSize(glyph.text);
  style::popFont(mono);
  dl->AddText(style::fonts::mono(), textSize,
              ImVec2(centre.x - measured.x * (textSize / ImGui::GetFontSize()) * 0.5f,
                     centre.y - textSize * 0.5f),
              color, glyph.text);
}

// ---------------------------------------------------------------- gallery UI
void galleryHeader(TileGlyph glyph, const char* title, const char* subtitle) {
  const style::Metrics& m = style::metrics();
  const float iconBox = m.iconSize * 1.75f;
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::Dummy(ImVec2(std::max(300.0f, ImGui::GetContentRegionAvail().x), iconBox));
  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddRectFilled(min, ImVec2(min.x + iconBox, min.y + iconBox),
                    style::u32(style::col::BgSurface), m.radiusMd);
  dl->AddRect(min, ImVec2(min.x + iconBox, min.y + iconBox), style::u32(style::col::BorderStrong),
              m.radiusMd, 0, m.hairline);
  drawGlyph(dl, glyph, ImVec2(min.x + iconBox * 0.5f, min.y + iconBox * 0.5f),
            m.iconSize * 1.05f, style::u32(style::col::Accent));

  const float textX = min.x + iconBox + m.gap;
  const bool heading = style::pushFont(style::fonts::semibold());
  dl->AddText(ImVec2(textX, min.y + 1.0f), style::u32(style::col::Text), title);
  style::popFont(heading);
  dl->AddText(ImVec2(textX, min.y + ImGui::GetFontSize() + m.gap * 0.35f),
              style::u32(style::col::TextDim), subtitle);
  ImGui::Spacing();
}

void gallerySection(const char* label) { widgets::sectionHeader(label); }

template <typename T, size_t N, typename GlyphOf, typename NameOf, typename SelectedOf>
int drawTileGrid(const char* idPrefix, const std::array<T, N>& entries, int columns,
                 GlyphOf glyphOf, NameOf nameOf, SelectedOf selectedOf) {
  const style::Metrics& m = style::metrics();
  const float pad = m.gap * 0.85f;
  const float glyphBox = m.iconSize * 1.55f;
  const float labelH = ImGui::GetFontSize();
  const float spacing = ImGui::GetStyle().ItemSpacing.x;

  float tileW = glyphBox + pad * 2.0f;
  for (const T& entry : entries)
    tileW = std::max(tileW, ImGui::CalcTextSize(nameOf(entry)).x + pad * 2.0f);
  const float tileH = glyphBox + labelH + pad * 2.15f;

  int clicked = -1;
  ImDrawList* dl = ImGui::GetWindowDrawList();
  ImGui::PushID(idPrefix);
  for (size_t i = 0; i < entries.size(); ++i) {
    if (static_cast<int>(i) % columns != 0) ImGui::SameLine(0.0f, spacing);
    const T& entry = entries[i];
    ImGui::PushID(static_cast<int>(i));
    const ImVec2 min = ImGui::GetCursorScreenPos();
    ImGui::InvisibleButton("##tile", ImVec2(tileW, tileH));
    const bool hit = ImGui::IsItemClicked();
    const bool hovered = ImGui::IsItemHovered();
    const float hover = widgets::hoverT(ImGui::GetItemID(), hovered);
    ImGui::PopID();

    const ImVec2 max(min.x + tileW, min.y + tileH);
    const bool selected = selectedOf(entry);
    const ImU32 fill = selected
                           ? style::mix(style::col::BgSurface, style::col::Accent,
                                        0.13f + hover * 0.06f)
                           : style::mix(style::col::BgSurface, style::col::BgRaised, hover);
    dl->AddRectFilled(min, max, fill, m.radiusMd);
    dl->AddRect(min, max,
                selected ? style::u32(style::col::Accent, 0.85f)
                         : style::mix(style::col::Border, style::col::BorderStrong, hover),
                m.radiusMd, 0, m.hairline);
    if (selected) {
      dl->AddCircleFilled(ImVec2(max.x - pad * 0.75f, min.y + pad * 0.75f),
                          std::max(2.5f, m.hairline * 2.5f), style::u32(style::col::Accent));
    }

    const ImU32 glyphColor = selected
                                 ? style::u32(style::col::Accent)
                                 : style::mix(style::col::TextDim, style::col::Text, hover);
    const ImVec2 glyphCentre((min.x + max.x) * 0.5f, min.y + pad + glyphBox * 0.5f);
    drawGlyph(dl, glyphOf(entry), glyphCentre, glyphBox * 0.78f, glyphColor);
    const ImVec2 label = ImGui::CalcTextSize(nameOf(entry));
    dl->AddText(ImVec2((min.x + max.x - label.x) * 0.5f, max.y - pad - labelH),
                style::u32(selected ? style::col::Text : style::col::TextDim), nameOf(entry));

    if (hit) clicked = static_cast<int>(i);
  }
  ImGui::PopID();
  return clicked;
}

// --------------------------------------------------------------- command row
struct RowResult {
  bool body = false;
  bool menu = false;
  ImVec2 min{};
  ImVec2 max{};
};

RowResult drawCommandRow(TileGlyph glyph, const char* label, const char* detail,
                         const char* shortcut, bool selected, bool hasMenu,
                         const char* tooltip) {
  const style::Metrics& m = style::metrics();
  const float width = std::max(1.0f, ImGui::GetContentRegionAvail().x);
  const float height = std::max(52.0f, ImGui::GetFrameHeight() + m.gap * 2.7f);
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::InvisibleButton("##command", ImVec2(width, height));
  const bool clicked = ImGui::IsItemClicked();
  const bool hovered = ImGui::IsItemHovered();
  const float hover = widgets::hoverT(ImGui::GetItemID(), hovered);
  const ImVec2 max(min.x + width, min.y + height);

  const float menuW = hasMenu ? std::max(28.0f, m.gap * 3.0f) : 0.0f;
  const ImVec2 mouse = ImGui::GetMousePos();
  const bool overMenu = hasMenu && hovered && mouse.x >= max.x - menuW;
  const bool compact = width < 175.0f;
  const bool showShortcut = !compact && shortcut && shortcut[0] != '\0';

  ImDrawList* dl = ImGui::GetWindowDrawList();
  const ImVec4 selectedBase{0.16f, 0.135f, 0.085f, 1.0f};
  const ImU32 fill = selected ? style::mix(selectedBase, style::col::BgRaised, hover * 0.35f)
                              : style::mix(style::col::BgSurface, style::col::BgRaised, hover);
  dl->AddRectFilled(min, max, fill, m.radiusMd);
  dl->AddRect(min, max,
              selected ? style::u32(style::col::Accent, 0.65f)
                       : style::mix(style::col::Border, style::col::BorderStrong, hover),
              m.radiusMd, 0, m.hairline);
  if (selected) {
    dl->AddRectFilled(min, ImVec2(min.x + std::max(3.0f, m.gap * 0.34f), max.y),
                      style::u32(style::col::Accent), m.radiusMd,
                      ImDrawFlags_RoundCornersLeft);
  }

  const float pad = m.gap * 0.85f;
  const float iconBox = height - pad * 1.55f;
  const ImVec2 iconMin(min.x + pad, min.y + (height - iconBox) * 0.5f);
  const ImVec2 iconMax(iconMin.x + iconBox, iconMin.y + iconBox);
  dl->AddRectFilled(iconMin, iconMax,
                    selected ? style::u32(style::col::Accent, 0.18f)
                             : style::u32(style::col::BgRaised),
                    m.radiusSm);
  dl->AddRect(iconMin, iconMax,
              selected ? style::u32(style::col::Accent, 0.45f)
                       : style::u32(style::col::Border),
              m.radiusSm, 0, m.hairline);
  drawGlyph(dl, glyph, ImVec2((iconMin.x + iconMax.x) * 0.5f, (iconMin.y + iconMax.y) * 0.5f),
            m.iconSize, selected ? style::u32(style::col::Accent)
                                 : style::mix(style::col::TextDim, style::col::Text, hover));

  const float textX = iconMax.x + pad;
  const float rightEdge = max.x - menuW - pad;
  const float titleY = compact ? min.y + (height - ImGui::GetFontSize()) * 0.5f
                               : min.y + pad * 0.72f;
  const bool heading = style::pushFont(style::fonts::semibold());
  dl->PushClipRect(ImVec2(textX, min.y), ImVec2(rightEdge, max.y), true);
  dl->AddText(ImVec2(textX, titleY), style::u32(style::col::Text), label);
  style::popFont(heading);
  if (!compact && detail && detail[0] != '\0') {
    dl->AddText(ImVec2(textX, titleY + ImGui::GetFontSize() + m.gap * 0.28f),
                style::u32(style::col::TextDim), detail);
  }
  dl->PopClipRect();

  if (showShortcut) {
    const ImVec2 measured = ImGui::CalcTextSize(shortcut);
    const float chipW = measured.x + m.gap * 1.15f;
    const float chipH = ImGui::GetFontSize() + m.gap * 0.55f;
    const ImVec2 chipMax(rightEdge, min.y + (height + chipH) * 0.5f);
    const ImVec2 chipMin(chipMax.x - chipW, chipMax.y - chipH);
    if (chipMin.x > textX + ImGui::CalcTextSize(label).x + m.gap) {
      dl->AddRectFilled(chipMin, chipMax, style::u32(style::col::BgDeep, 0.65f), chipH * 0.5f);
      dl->AddRect(chipMin, chipMax, style::u32(style::col::Border), chipH * 0.5f, 0,
                  m.hairline);
      dl->AddText(ImVec2(chipMin.x + (chipW - measured.x) * 0.5f,
                         chipMin.y + (chipH - measured.y) * 0.5f),
                  style::u32(style::col::TextFaint), shortcut);
    }
  }

  if (hasMenu) {
    const float dividerX = max.x - menuW;
    dl->AddLine(ImVec2(dividerX, min.y + pad), ImVec2(dividerX, max.y - pad),
                style::u32(overMenu ? style::col::Accent : style::col::Border), m.hairline);
    const ImVec2 c(dividerX + menuW * 0.5f, min.y + height * 0.5f);
    const float s = std::max(3.5f, m.gap * 0.42f);
    dl->AddTriangleFilled(ImVec2(c.x - s, c.y - s * 0.45f),
                          ImVec2(c.x + s, c.y - s * 0.45f), ImVec2(c.x, c.y + s * 0.55f),
                          style::u32(overMenu ? style::col::Accent : style::col::TextDim));
  }

  if (hovered && tooltip) ImGui::SetTooltip("%s", tooltip);
  return RowResult{clicked && !overMenu, clicked && overMenu, min, max};
}

template <typename Popup>
void openRowPopup(const char* id, const RowResult& row, Popup&& content) {
  const style::Metrics& m = style::metrics();
  if (row.menu) ImGui::OpenPopup(id);
  ImGui::SetNextWindowPos(ImVec2(row.max.x + m.gap, row.min.y), ImGuiCond_Appearing);
  ImGui::SetNextWindowSizeConstraints(ImVec2(330.0f, 0.0f),
                                      ImVec2(ImGui::GetMainViewport()->WorkSize.x * 0.62f,
                                             ImGui::GetMainViewport()->WorkSize.y * 0.78f));
  ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(m.gap * 1.25f, m.gap * 1.15f));
  if (ImGui::BeginPopup(id, ImGuiWindowFlags_NoMove)) {
    content();
    ImGui::EndPopup();
  }
  ImGui::PopStyleVar();
}

// ----------------------------------------------------------------- galleries
void drawRingGallery(AppState& st) {
  const RingEntry& active = ringEntry(st.currentRing);
  galleryHeader({active.icon, nullptr}, "Ring templates", "Choose the structure to stamp");
  gallerySection("CARBOCYCLES & AROMATICS");
  const int hit = drawTileGrid(
      "##rings", kRings, 3, [](const RingEntry& entry) { return TileGlyph{entry.icon, nullptr}; },
      [](const RingEntry& entry) { return entry.name; },
      [&](const RingEntry& entry) { return entry.kind == st.currentRing; });
  if (hit >= 0) {
    st.currentRing = kRings[static_cast<size_t>(hit)].kind;
    st.tool = Tool::RingTemplate;
    st.statusMessage = std::string("Ring template: ") + kRings[static_cast<size_t>(hit)].name;
    ImGui::CloseCurrentPopup();
  }
}

void drawBondGallery(AppState& st) {
  galleryHeader(bondGlyph(st), "Bond palette", "Order and stereochemical direction");
  gallerySection("BOND ORDER");
  const int order = drawTileGrid(
      "##orders", kOrders, 4,
      [](const OrderEntry& entry) { return TileGlyph{entry.icon, nullptr}; },
      [](const OrderEntry& entry) { return entry.name; },
      [&](const OrderEntry& entry) {
        return st.currentStereo == core::BondStereo::None && st.currentOrder == entry.order;
      });
  gallerySection("STEREOCHEMISTRY");
  const int stereo = drawTileGrid(
      "##stereos", kStereos, 4,
      [](const StereoEntry& entry) { return TileGlyph{entry.icon, nullptr}; },
      [](const StereoEntry& entry) { return entry.name; },
      [&](const StereoEntry& entry) { return st.currentStereo == entry.stereo; });
  ImGui::Spacing();
  ImGui::TextDisabled("M restores a methyl-ready plain single bond.");

  if (order >= 0) {
    st.currentOrder = kOrders[static_cast<size_t>(order)].order;
    st.currentStereo = core::BondStereo::None;
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
  galleryHeader({icons::Icon::Atom, chem::symbolFor(st.currentElement)}, "Element quick pick",
                "Common elements for structure drawing");
  gallerySection("COMMON ELEMENTS");
  const int hit = drawTileGrid(
      "##elements", kQuickElements, 4,
      [](const QuickElement& entry) {
        return TileGlyph{icons::Icon::Atom, chem::symbolFor(entry.z)};
      },
      [](const QuickElement& entry) { return entry.name; },
      [&](const QuickElement& entry) { return st.currentElement == entry.z; });
  ImGui::Spacing();
  ImGui::TextDisabled("The complete element set remains in the Periodic Table.");
  if (hit >= 0) {
    st.currentElement = kQuickElements[static_cast<size_t>(hit)].z;
    st.tool = Tool::Atom;
    st.statusMessage = std::string("Element: ") + kQuickElements[static_cast<size_t>(hit)].name;
    ImGui::CloseCurrentPopup();
  }
}

// ------------------------------------------------------------- rail structure
void paletteSectionLabel(const char* label) {
  const style::Metrics& m = style::metrics();
  ImGui::Spacing();
  const ImVec2 min = ImGui::GetCursorScreenPos();
  const float height = ImGui::GetFontSize() + m.gap * 0.45f;
  ImGui::Dummy(ImVec2(ImGui::GetContentRegionAvail().x, height));
  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddText(nullptr, ImGui::GetFontSize() * 0.78f,
              ImVec2(min.x + m.gap * 0.15f, min.y + m.gap * 0.15f),
              style::u32(style::col::TextFaint), label);
}

void drawActiveToolCard(const AppState& st) {
  const style::Metrics& m = style::metrics();
  const float width = std::max(1.0f, ImGui::GetContentRegionAvail().x);
  const float height = std::max(70.0f, ImGui::GetFrameHeight() + m.gap * 4.5f);
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::Dummy(ImVec2(width, height));
  const ImVec2 max(min.x + width, min.y + height);
  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddRectFilled(min, max, style::u32(style::col::BgSurface), m.radiusLg);
  dl->AddRect(min, max, style::u32(style::col::BorderStrong), m.radiusLg, 0, m.hairline);
  dl->AddLine(ImVec2(min.x + m.gap, min.y), ImVec2(max.x - m.gap, min.y),
              style::u32(style::col::Accent), std::max(2.0f, m.hairline * 2.0f));

  const float iconBox = height - m.gap * 2.0f;
  const ImVec2 iconMin(min.x + m.gap, min.y + m.gap);
  const ImVec2 iconMax(iconMin.x + iconBox, iconMin.y + iconBox);
  dl->AddRectFilled(iconMin, iconMax, style::u32(style::col::Accent, 0.12f), m.radiusMd);
  dl->AddRect(iconMin, iconMax, style::u32(style::col::Accent, 0.42f), m.radiusMd, 0,
              m.hairline);
  drawGlyph(dl, toolGlyph(st, st.tool),
            ImVec2((iconMin.x + iconMax.x) * 0.5f, (iconMin.y + iconMax.y) * 0.5f),
            m.iconSize * 1.2f, style::u32(style::col::Accent));

  const float textX = iconMax.x + m.gap;
  dl->PushClipRect(ImVec2(textX, min.y), ImVec2(max.x - m.gap, max.y), true);
  dl->AddText(nullptr, ImGui::GetFontSize() * 0.72f, ImVec2(textX, min.y + m.gap * 0.9f),
              style::u32(style::col::TextFaint), "ACTIVE TOOL");
  const bool heading = style::pushFont(style::fonts::semibold());
  dl->AddText(ImVec2(textX, min.y + m.gap * 2.2f), style::u32(style::col::Text),
              toolEntry(st.tool).name);
  style::popFont(heading);
  const std::string detail = toolDetail(st, st.tool);
  dl->AddText(ImVec2(textX, min.y + m.gap * 2.2f + ImGui::GetFontSize() + m.gap * 0.35f),
              style::u32(style::col::TextDim), detail.c_str());
  dl->PopClipRect();
}

void drawSimpleTool(AppState& st, const ToolEntry& entry) {
  ImGui::PushID(static_cast<int>(entry.tool));
  char tooltip[256];
  std::snprintf(tooltip, sizeof(tooltip), "%s (%s)\n%s", entry.name, entry.shortcut, entry.hint);
  const std::string detail = toolDetail(st, entry.tool);
  const RowResult row = drawCommandRow({entry.icon, nullptr}, entry.name, detail.c_str(),
                                       entry.shortcut, st.tool == entry.tool, false, tooltip);
  if (row.body) st.tool = entry.tool;
  ImGui::PopID();
}

void drawBondTool(AppState& st) {
  const ToolEntry& entry = toolEntry(Tool::Bond);
  ImGui::PushID(static_cast<int>(entry.tool));
  const RowResult row = drawCommandRow(bondGlyph(st), entry.name, bondStyleName(st),
                                       entry.shortcut, st.tool == Tool::Bond, true,
                                       "Bond tool\nClick to draw; open the palette to change style");
  if (row.body) st.tool = Tool::Bond;
  openRowPopup("##bond_gallery", row, [&] { drawBondGallery(st); });
  ImGui::PopID();
}

void drawRingTool(AppState& st) {
  const ToolEntry& entry = toolEntry(Tool::RingTemplate);
  ImGui::PushID(static_cast<int>(entry.tool));
  const RingEntry& active = ringEntry(st.currentRing);
  const RowResult row = drawCommandRow({active.icon, nullptr}, entry.name, active.name,
                                       entry.shortcut, st.tool == Tool::RingTemplate, true,
                                       "Ring tool\nClick to stamp; open the palette to choose a template");
  if (row.body) st.tool = Tool::RingTemplate;
  openRowPopup("##ring_gallery", row, [&] { drawRingGallery(st); });
  ImGui::PopID();
}

void drawAtomTool(AppState& st) {
  const ToolEntry& entry = toolEntry(Tool::Atom);
  ImGui::PushID(static_cast<int>(entry.tool));
  const ElementData* element = findElement(st.currentElement);
  const std::string detail = std::string(chem::symbolFor(st.currentElement)) + " · " +
                             (element ? element->name : "Element");
  const RowResult row = drawCommandRow({icons::Icon::Atom, chem::symbolFor(st.currentElement)},
                                       entry.name, detail.c_str(), entry.shortcut,
                                       st.tool == Tool::Atom, true,
                                       "Atom tool\nClick to place; open the palette for common elements");
  if (row.body) st.tool = Tool::Atom;
  openRowPopup("##atom_gallery", row, [&] { drawAtomGallery(st); });
  ImGui::PopID();
}

}  // namespace

void drawToolPalette(AppState& st) {
  drawActiveToolCard(st);

  paletteSectionLabel("NAVIGATE");
  drawSimpleTool(st, toolEntry(Tool::Select));
  drawSimpleTool(st, toolEntry(Tool::Eraser));

  paletteSectionLabel("BUILD");
  drawBondTool(st);
  drawSimpleTool(st, toolEntry(Tool::Chain));
  drawRingTool(st);
  drawAtomTool(st);

  paletteSectionLabel("MODIFY");
  drawSimpleTool(st, toolEntry(Tool::ChargePlus));
  drawSimpleTool(st, toolEntry(Tool::ChargeMinus));

  ImGui::Spacing();
  ImGui::PushTextWrapPos();
  ImGui::TextDisabled("Keyboard shortcuts stay active while the sketch canvas is focused.");
  ImGui::PopTextWrapPos();
}

}  // namespace chemcad::ui

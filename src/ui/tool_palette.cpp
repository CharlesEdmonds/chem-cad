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
  ImGui::Dummy(
      ImVec2(std::max(ImGui::GetFontSize() * 18.0f, ImGui::GetContentRegionAvail().x), iconBox));
  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddRectFilled(min, ImVec2(min.x + iconBox, min.y + iconBox),
                    style::u32(style::col::BgSurface), m.radiusMd);
  dl->AddRect(min, ImVec2(min.x + iconBox, min.y + iconBox), style::u32(style::col::BorderStrong),
              m.radiusMd, 0, m.hairline);
  drawGlyph(dl, glyph, ImVec2(min.x + iconBox * 0.5f, min.y + iconBox * 0.5f),
            m.iconSize * 1.05f, style::u32(style::col::Accent));

  const float textX = min.x + iconBox + m.gap;
  const bool heading = style::pushFont(style::fonts::semibold());
  dl->AddText(ImVec2(textX, min.y + ImGui::GetFontSize() * 0.06f),
              style::u32(style::col::Text), title);
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
                          std::max(ImGui::GetFontSize() * 0.14f, m.hairline * 2.5f),
                          style::u32(style::col::Accent));
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

  constexpr std::array orderGlyphs{
      icons::Icon::BondSingle,
      icons::Icon::BondDouble,
      icons::Icon::BondTriple,
      icons::Icon::BondAromatic,
  };
  constexpr std::array<const char*, 4> orderTips{
      "Single bond",
      "Double bond",
      "Triple bond",
      "Aromatic bond",
  };
  int orderIndex = 0;
  for (size_t i = 0; i < kOrders.size(); ++i) {
    if (kOrders[i].order == st.currentOrder) orderIndex = static_cast<int>(i);
  }
  widgets::cardHeader(icons::Icon::Bond, "Bond order", "Mutually exclusive bond order");
  if (widgets::segmentedIcons("##bond_order", orderGlyphs.data(), orderTips.data(),
                              static_cast<int>(orderGlyphs.size()), orderIndex)) {
    st.currentOrder = kOrders[static_cast<size_t>(orderIndex)].order;
    st.currentStereo = core::BondStereo::None;
    st.tool = Tool::Bond;
    st.statusMessage = kOrders[static_cast<size_t>(orderIndex)].name;
  }

  constexpr std::array stereoGlyphs{
      icons::Icon::StereoNone,
      icons::Icon::StereoWedge,
      icons::Icon::StereoHash,
      icons::Icon::StereoWavy,
  };
  constexpr std::array<const char*, 4> stereoTips{
      "Plain bond",
      "Solid wedge",
      "Hashed wedge",
      "Wavy bond",
  };
  int stereoIndex = 0;
  for (size_t i = 0; i < kStereos.size(); ++i) {
    if (kStereos[i].stereo == st.currentStereo) stereoIndex = static_cast<int>(i);
  }
  ImGui::Spacing();
  widgets::cardHeader(icons::Icon::StereoWedge, "Stereochemistry",
                      "Directional display for the active bond");
  if (widgets::segmentedIcons("##bond_stereo", stereoGlyphs.data(), stereoTips.data(),
                              static_cast<int>(stereoGlyphs.size()), stereoIndex)) {
    st.currentStereo = kStereos[static_cast<size_t>(stereoIndex)].stereo;
    st.tool = Tool::Bond;
    st.statusMessage = kStereos[static_cast<size_t>(stereoIndex)].name;
  }

  ImGui::Spacing();
  ImGui::TextDisabled("M restores a methyl-ready plain single bond.");
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

template <size_t N>
void drawToolGrid(AppState& st, const std::array<Tool, N>& tools) {
  const style::Metrics& m = style::metrics();
  const float spacing = ImGui::GetStyle().ItemSpacing.x;
  const float available = std::max(ImGui::GetContentRegionAvail().x, ImGui::GetFrameHeight());
  const float preferred = std::max(ImGui::GetFrameHeight() * 1.85f, ImGui::GetFontSize() * 3.2f);
  const int columns =
      std::max(1, static_cast<int>((available + spacing) / (preferred + spacing)));
  const float tile = std::max(
      ImGui::GetFrameHeight(),
      (available - spacing * static_cast<float>(columns - 1)) / static_cast<float>(columns));

  for (size_t i = 0; i < tools.size(); ++i) {
    if (static_cast<int>(i) % columns != 0) ImGui::SameLine(0.0f, spacing);
    const ToolEntry& entry = toolEntry(tools[i]);
    ImGui::PushID(static_cast<int>(entry.tool));

    char tooltip[256];
    if (entry.shortcut && entry.shortcut[0] != '\0') {
      std::snprintf(tooltip, sizeof(tooltip), "%s (%s)\n%s", entry.name, entry.shortcut,
                    entry.hint);
    } else {
      std::snprintf(tooltip, sizeof(tooltip), "%s\n%s", entry.name, entry.hint);
    }

    const TileGlyph glyph = toolGlyph(st, entry.tool);
    const bool clicked = widgets::iconButton("##tool", glyph.icon, ImVec2(tile, tile),
                                             st.tool == entry.tool, tooltip);
    if (clicked) {
      st.tool = entry.tool;
      st.statusMessage = std::string(entry.name) + " tool";
      if (entry.tool == Tool::Bond) ImGui::OpenPopup("##bond_gallery");
      if (entry.tool == Tool::RingTemplate) ImGui::OpenPopup("##ring_gallery");
      if (entry.tool == Tool::Atom) ImGui::OpenPopup("##atom_gallery");
    }

    if (entry.tool == Tool::Bond || entry.tool == Tool::RingTemplate ||
        entry.tool == Tool::Atom) {
      const float font = ImGui::GetFontSize();
      ImGui::SetNextWindowSizeConstraints(
          ImVec2(font * 20.0f, 0.0f),
          ImVec2(ImGui::GetMainViewport()->WorkSize.x * 0.62f,
                 ImGui::GetMainViewport()->WorkSize.y * 0.78f));
      ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding,
                          ImVec2(m.gap * 1.25f, m.gap * 1.15f));
      if (entry.tool == Tool::Bond && ImGui::BeginPopup("##bond_gallery")) {
        drawBondGallery(st);
        ImGui::EndPopup();
      } else if (entry.tool == Tool::RingTemplate && ImGui::BeginPopup("##ring_gallery")) {
        drawRingGallery(st);
        ImGui::EndPopup();
      } else if (entry.tool == Tool::Atom && ImGui::BeginPopup("##atom_gallery")) {
        drawAtomGallery(st);
        ImGui::EndPopup();
      }
      ImGui::PopStyleVar();
    }
    ImGui::PopID();
  }
}

}  // namespace

void drawToolPalette(AppState& st) {
  const std::string activeSummary = "Active tool · " + toolDetail(st, st.tool);
  widgets::cardHeader(toolGlyph(st, st.tool).icon, toolEntry(st.tool).name,
                      activeSummary.c_str(), style::col::Accent);

  widgets::cardHeader(icons::Icon::Select, "Navigate", "Select, inspect, or remove");
  drawToolGrid(st, std::array{Tool::Select, Tool::Eraser});

  ImGui::Spacing();
  widgets::cardHeader(icons::Icon::Bond, "Build", "Bonds, chains, rings, and atoms");
  drawToolGrid(st, std::array{Tool::Bond, Tool::Chain, Tool::RingTemplate, Tool::Atom});

  ImGui::Spacing();
  if (widgets::disclosure("##charge_tools", "Formal charge", "2 tools", false,
                          icons::Icon::ChargePlus, style::col::Teal)) {
    widgets::cardHeader(icons::Icon::ChargePlus, "Charge tools",
                        "Increase or decrease formal charge", style::col::Teal);
    drawToolGrid(st, std::array{Tool::ChargePlus, Tool::ChargeMinus});
  }

  ImGui::Spacing();
  widgets::notice(icons::Icon::Info,
                  "Keyboard shortcuts stay active while the sketch canvas is focused.",
                  style::col::TextDim);
}

}  // namespace chemcad::ui

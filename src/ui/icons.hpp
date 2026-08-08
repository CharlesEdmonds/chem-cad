#pragma once
// Procedural tool glyphs drawn with ImDrawList. No icon font: every glyph is
// a handful of strokes on a normalized grid, which keeps them crisp at any UI
// scale and gives the palette its own look.

#include "imgui.h"

namespace chemcad::ui::icons {

enum class Icon {
  Select,
  Eraser,
  Bond,
  Chain,
  RingCyclopropane,
  RingCyclobutane,
  RingCyclopentane,
  RingCyclohexane,
  RingCycloheptane,
  RingCyclooctane,
  RingBenzene,
  RingCyclopentadiene,
  RingNaphthalene,
  Atom,
  ChargePlus,
  ChargeMinus,
  BondSingle,
  BondDouble,
  BondTriple,
  BondAromatic,
  StereoNone,
  StereoWedge,
  StereoHash,
  StereoWavy,
  Plus,
  Close,
  Search,
  Copy,
  ArrowRight,
  Logo,  // benzene-ring brand mark
};

// Draws `icon` centred at `centre`, inside a square of side `size`.
// `thickness` <= 0 picks a stroke proportional to `size`.
void draw(ImDrawList* dl, Icon icon, ImVec2 centre, float size, ImU32 color,
          float thickness = 0.0f);

}  // namespace chemcad::ui::icons

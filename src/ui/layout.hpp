#pragma once
// Responsive layout for panels that must fit their surface exactly.
//
// The application runs on laptop panels, 4K monitors and ultrawides, at OS
// display scales from 100% to 250%, and under its own CHEMCAD_UI_SCALE. A panel
// that hardcodes pixel thresholds is correct on exactly one of those. Two rules
// make a panel correct on all of them:
//
//   1. The only unit is the em -- ImGui::GetFontSize(). It already carries the
//      OS scale and the user scale, so an em-relative layout is resolution and
//      zoom independent by construction.
//   2. A panel that must not scroll cannot let its children choose their own
//      heights. It measures the budget once and DIVIDES it, so the parts sum to
//      the whole no matter how tall the whole is.
//
// Panels ask `measure()` what surface they were given and lay out against the
// answer. They never test raw pixel widths.

#include "imgui.h"

#include <cstdint>

namespace chemcad::ui::layout {

// How much room each element may take. Chosen from the viewport measured in
// ems, so a small window and a large window at high OS zoom land in the same
// tier -- which is the point, since they hold the same amount of text.
enum class Density : std::uint8_t {
  Compact,  // < 46 em wide, or short: single column, labels abbreviated
  Regular,  // the design target
  Roomy,    // >= 96 em wide: charts may take their preferred size
};

// The measured surface a panel has to work with.
struct Frame {
  ImVec2 size{0.0f, 0.0f};  // content region, pixels
  float em = 16.0f;         // ImGui::GetFontSize(): the one true unit
  float gap = 8.0f;         // spacing between siblings at this density
  float pad = 8.0f;         // interior padding of a card at this density
  float row = 20.0f;        // one text row including spacing
  float control = 24.0f;    // one interactive control's height
  Density density = Density::Regular;
  float aspect = 1.6f;      // size.x / size.y
  bool wide = true;         // aspect >= 1.4: side-by-side columns pay off
  bool tall = false;        // aspect < 0.9: stack, and prefer vertical charts

  // Width in ems, the number panels should branch on instead of pixels.
  float ems() const { return em > 0.0f ? size.x / em : 0.0f; }
};

// Measures the current ImGui window's content region.
Frame measure();
// Measures an explicit region, for a child whose size the caller already knows.
Frame measure(ImVec2 available);

// Columns of equal width that exactly fill `frame.size.x`, gaps included.
// `span` returns the width of several adjacent columns plus the gaps between.
float columnWidth(const Frame& frame, int columns, int span = 1);

// How many equal columns of at least `minEm` fit. Always >= 1.
int columnsThatFit(const Frame& frame, float minEm);

// Divides `budget` among `count` weighted rows so the results sum to EXACTLY
// `budget` after `count - 1` gaps -- the rounding is absorbed in the last row
// rather than accumulating into an overflow that produces a scrollbar.
// `minimums` may be null; when supplied, every row gets at least its minimum
// and the remainder is shared by weight.
void distribute(float budget, const float* weights, const float* minimums, int count,
                float gap, float* out);

// Height still available in the current window before it would scroll.
float pageHeight();

// Text below this is unreadable regardless of what the layout would prefer, so
// charts clamp their tick and legend fonts to it instead of shrinking without
// limit. Expressed in pixels at the current scale.
float minReadablePx();

// Moves the cursor to the top of the next band in a `distribute`d column.
//
// Setting the cursor alone is NOT enough: ImGui grows a window's content
// extent from submitted ITEMS, so a bare SetCursorPosY past the last item
// makes it complain that the boundary cannot grow, once per frame, forever.
// Claiming the position with a zero-size item is the sanctioned way to say
// "this space is mine", so every fixed-height band advance goes through here.
void nextRow(float yPosition);

// The size to actually draw label text at: `requested` clamped to
// [minReadablePx(), em]. Every chart axis, legend and tick goes through this.
float labelFont(float requested);

// True when `text` fits in `width` at the current font; panels use it to pick
// between a full label and its abbreviation rather than guessing from a
// character count.
bool fits(const char* text, float width);

// Longest of the candidates that fits, else the last one. Null-terminated list
// ordered longest first -- the idiom for "Interfacial area" / "Interf. area" /
// "IFA" without hardcoding breakpoints.
const char* bestLabel(float width, const char* const* candidates, int count);

}  // namespace chemcad::ui::layout

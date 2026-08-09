#pragma once
// Custom-drawn widgets built on ImDrawList + invisible buttons. These are the
// building blocks of the Bench look: outlined surfaces, animated hover lifts,
// amber selection. Behaviour matches their ImGui counterparts (IsItemClicked
// etc.) so headless interaction tests drive them exactly like stock widgets.

#include "imgui.h"

#include <string>
#include <string_view>

#include "ui/icons.hpp"
#include "ui/theme.hpp"

namespace chemcad::ui::widgets {

// Per-widget animated hover value (0..1), stored in the window's state
// storage. Call with ImGui::GetItemID() right after the behaviour item.
float hoverT(ImGuiID id, bool hovered);

// Case-insensitive substring search shared by filterable views.
bool containsCaseInsensitive(std::string_view haystack, std::string_view needle);

// ImGui text input backed by a dynamically resized std::string.
bool stringInputWithHint(const char* id, const char* hint, std::string& value);
bool stringInputWithHint(const char* id, const char* hint, std::string& value,
                         ImGuiInputTextFlags flags, bool mono = false);

// Solvent picker: a combo whose closed preview is the current solvent's name
// and whose popup leads with a search field. An empty query lists only the
// curated bench-rack solvents, because a 45-entry alphabetical wall makes the
// five solvents anyone actually reaches for hard to find; typing searches the
// whole database. Returns true when `solventId` changed.
bool solventCombo(const char* id, std::string& solventId, std::string& query);

// Section header: amber tick, small-caps-ish dim label, hairline to the edge.
void sectionHeader(const char* label, ImVec4 accent = style::col::Accent);

// Primary action: filled amber, dark label, semibold.
bool primaryButton(const char* label, ImVec2 size = ImVec2(0, 0));

// Quiet action: transparent surface, hairline border, lifts on hover.
bool ghostButton(const char* label, ImVec2 size = ImVec2(0, 0));

// Square glyph button for tool grids and icon rows. `selected` gets the amber
// fill; otherwise the surface lifts and the glyph brightens on hover.
bool iconButton(const char* id, icons::Icon icon, ImVec2 size, bool selected,
                const char* tooltip = nullptr);

// Inline pill badge (KB / AI / BEST). Advances the cursor; use between
// SameLine() calls like a text run.
void badge(const char* text, ImVec4 color);

// Dashboard stat: mono value over a dim label, on a raised card.
void statCard(const char* label, const char* value, ImVec2 size);

// Card container: BeginChild with a filled surface, border and rounding.
// `size.y == 0` auto-fits height (ImGuiChildFlags_AutoResizeY).
bool beginCard(const char* id, ImVec2 size, ImVec4 bg = style::col::BgSurface,
               ImGuiWindowFlags extraFlags = 0);
void endCard();

// --------------------------------------------------------------- structure
// Progressive disclosure is how these panels stay readable: the summary is the
// surface-level answer and the body is the derivation. Draws an animated
// chevron, an optional leading glyph, the label, and a right-aligned summary
// that survives collapse. Returns true while open; the caller indents its own
// body. State lives in the window's storage, keyed on `id`.
bool disclosure(const char* id, const char* label, const char* summary,
                bool defaultOpen = false, icons::Icon icon = icons::Icon::None,
                ImVec4 accent = style::col::Accent);

// Header for a card: leading glyph, semibold title, dim subtitle, hairline.
void cardHeader(icons::Icon icon, const char* title, const char* subtitle = nullptr,
                ImVec4 accent = style::col::Accent);

// Toolbar strip: a raised, rounded bar the caller fills with a horizontal item
// run. Pair the calls; endToolbar closes the group and draws the backing.
void beginToolbar(const char* id);
void endToolbar();

// Vertical hairline separator sized for a toolbar run.
void toolbarSeparator();

// ------------------------------------------------------------------ inputs
// One row of mutually exclusive options with a sliding amber indicator.
// Returns true when `index` changed. `width <= 0` fills the available width.
bool segmented(const char* id, const char* const* labels, int count, int& index,
               float width = 0.0f);
bool segmentedIcons(const char* id, const icons::Icon* glyphs, const char* const* tooltips,
                    int count, int& index, float width = 0.0f);

// iOS-style switch with a trailing label. Returns true when `value` changed.
bool toggle(const char* id, const char* label, bool& value, const char* tooltip = nullptr);

// Filter/tag chip. `accent` tints the selected fill.
bool chip(const char* label, bool selected, ImVec4 accent = style::col::Accent);

// Icon + label action. `primary` fills amber; otherwise it is a ghost surface.
bool actionButton(const char* id, icons::Icon icon, const char* label, ImVec2 size,
                  bool primary = false, const char* tooltip = nullptr);

// The width `actionButton` would choose for itself. Callers that must reserve
// space for a trailing action need the real number: re-deriving it by hand is
// how a label ends up spilling past its own button.
float actionButtonWidth(icons::Icon icon, const char* label, bool primary = false);

// Slider whose track carries a glyph and whose value renders inside the frame,
// so a control row does not need a separate label column.
bool glyphSlider(const char* id, icons::Icon icon, const char* label, float& value,
                 float minValue, float maxValue, const char* format,
                 const char* tooltip = nullptr);

// ---------------------------------------------------------------- readouts
// Caption over a mono value with an optional unit; `delta` is drawn as a signed
// tinted suffix when non-null.
void metric(const char* caption, const char* value, const char* unit = nullptr,
            const char* delta = nullptr, ImVec4 accent = style::col::Text);

// Key on the left, right-aligned mono value. The workhorse of the detail rows
// that live inside a disclosure body.
void keyValue(const char* key, const char* value, ImVec4 accent = style::col::Text);

// Caption, thin progress track and a right-aligned value on one line.
void progressRow(const char* caption, float fraction, const char* value, ImVec4 accent);

// Dim glyph that reveals `text` on hover. For the explanation a control needs
// but should not wear.
void helpMarker(const char* text);

// ------------------------------------------------------------- organisation
// Sub-tab strip inside a panel. This is how a workspace that holds more than
// one screen of content still fits one screen: the structural split becomes
// visible instead of becoming a scrollbar. Returns true when `index` changed.
// `glyphs` may be null.
bool subTabs(const char* id, const char* const* labels, const icons::Icon* glyphs,
             int count, int& index);

// HUD frame: a hairline rectangle with cut corners and short accent ticks, for
// instrument surfaces that should read as measurement rather than as chrome.
// Draws into the current window's draw list; reserves nothing.
void hudFrame(ImVec2 min, ImVec2 max, ImVec4 accent, float alpha = 1.0f);

// Small filled dot plus label, for live on/off state on a control surface.
void statusDot(const char* label, bool active, ImVec4 accent = style::col::Data);

// ------------------------------------------------------------------- tables
// A data table sized to its CONTENT, not stretched to the panel. Numeric
// columns are right-aligned in the mono font and measured from the widest
// value actually present, so a column of "1.33" is never as wide as its header
// plus half the panel. Row height is one text line plus the density gap -- a
// table of ten short numbers must not occupy a third of the page.
//
// `columns` describes each column once; `numeric` picks alignment and font.
struct Column {
  const char* label = nullptr;
  bool numeric = false;
  bool stretch = false;       // exactly one column may absorb slack (the name)
  const char* unit = nullptr;  // drawn dim in the header, not in every cell
  float minEm = 0.0f;          // floor for the measured width
};
bool beginDataTable(const char* id, const Column* columns, int count, ImVec2 size);
void endDataTable();
// Advances to the next cell. Text cells use the body font, numeric cells the
// mono font, right-aligned.
void dataCell(const char* text);
void dataCellf(const char* format, ...);
// Starts the next row; `accent` tints the row's left edge when non-zero alpha.
void dataRow(ImVec4 accent = ImVec4(0, 0, 0, 0));

// ------------------------------------------------------------- conditional
// Shows `body` only when `when` holds, and otherwise draws nothing at all --
// not a disabled control. Controls that cannot apply are noise; controls that
// silently vanish are confusing, so the collapsed form leaves a single dim line
// naming what would appear and why it does not.
// Returns whether the body was drawn.
bool onlyWhen(bool when, const char* absentReason);

// Empty-state / guidance block: centred glyph, headline, dim body.
void emptyState(icons::Icon icon, const char* headline, const char* body);

// Inline advisory strip: tinted surface, glyph, message.
void notice(icons::Icon icon, const char* text, ImVec4 accent = style::col::Teal);

}  // namespace chemcad::ui::widgets

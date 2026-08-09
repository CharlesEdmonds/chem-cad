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

// Empty-state / guidance block: centred glyph, headline, dim body.
void emptyState(icons::Icon icon, const char* headline, const char* body);

// Inline advisory strip: tinted surface, glyph, message.
void notice(icons::Icon icon, const char* text, ImVec4 accent = style::col::Teal);

}  // namespace chemcad::ui::widgets

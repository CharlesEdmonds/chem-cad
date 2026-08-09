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

}  // namespace chemcad::ui::widgets

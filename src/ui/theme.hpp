#pragma once
// ChemCAD "Bench" design language: palette tokens, scaled metrics and font
// handles shared by every panel and widget.
//
// Visual identity: a night-time lab bench. Deep blue-charcoal surfaces,
// hairline borders, a warm amber primary accent (the Bunsen flame) for
// selection and primary actions, teal for chemistry-positive signals, violet
// for AI-sourced results. Widgets size from these metrics or ImGui font
// metrics, never from hardcoded pixels (the app runs at 1.25 UI scale).

#include "imgui.h"

namespace chemcad::ui::style {

// ------------------------------------------------------------------ palette
namespace col {

inline constexpr ImVec4 BgDeep{0.043f, 0.055f, 0.075f, 1.00f};  // app backdrop, menu bar
inline constexpr ImVec4 BgPanel{0.067f, 0.086f, 0.114f, 1.00f};  // docked windows
inline constexpr ImVec4 BgSurface{0.094f, 0.118f, 0.153f, 1.00f};  // frames, inputs, tiles
inline constexpr ImVec4 BgRaised{0.122f, 0.153f, 0.196f, 1.00f};  // cards, hovered surfaces
inline constexpr ImVec4 Border{0.169f, 0.208f, 0.263f, 1.00f};
inline constexpr ImVec4 BorderStrong{0.243f, 0.298f, 0.369f, 1.00f};

inline constexpr ImVec4 Text{0.910f, 0.929f, 0.949f, 1.00f};
inline constexpr ImVec4 TextDim{0.604f, 0.655f, 0.706f, 1.00f};
inline constexpr ImVec4 TextFaint{0.369f, 0.420f, 0.471f, 1.00f};

inline constexpr ImVec4 Accent{0.949f, 0.663f, 0.232f, 1.00f};  // Bunsen amber
inline constexpr ImVec4 AccentHover{1.000f, 0.765f, 0.361f, 1.00f};
inline constexpr ImVec4 AccentActive{0.851f, 0.561f, 0.122f, 1.00f};
inline constexpr ImVec4 OnAccent{0.141f, 0.090f, 0.012f, 1.00f};  // text on amber

inline constexpr ImVec4 Teal{0.247f, 0.749f, 0.659f, 1.00f};  // KB / selection secondary
inline constexpr ImVec4 Violet{0.608f, 0.482f, 0.910f, 1.00f};  // AI results
inline constexpr ImVec4 Danger{0.898f, 0.337f, 0.306f, 1.00f};
inline constexpr ImVec4 Success{0.341f, 0.729f, 0.420f, 1.00f};

// Amber acts, cyan informs. Anything the user can press or that is selected
// wears the Bunsen amber above; anything that REPORTS -- traces, gauges, axes,
// numeric readouts, instrument frames -- wears cyan. Keeping the two families
// disjoint is what stops a dense instrument panel reading as noise: colour then
// tells you whether a thing is a control or a measurement before you read it.
inline constexpr ImVec4 Data{0.302f, 0.816f, 0.945f, 1.00f};       // primary data accent
inline constexpr ImVec4 DataDim{0.204f, 0.549f, 0.647f, 1.00f};    // secondary series, axes
inline constexpr ImVec4 DataBright{0.541f, 0.914f, 1.000f, 1.00f}; // hover, peak markers
inline constexpr ImVec4 GridLine{0.145f, 0.196f, 0.255f, 1.00f};   // chart grid, HUD rules
}  // namespace col

ImU32 u32(ImVec4 c, float alphaMul = 1.0f);
ImU32 mix(ImVec4 a, ImVec4 b, float t, float alphaMul = 1.0f);  // linear RGB lerp

// ------------------------------------------------------------------ metrics
// Scaled once by applyTheme(); read these instead of literals.
struct Metrics {
  float radiusSm = 4.0f;    // small controls, pills
  float radiusMd = 6.0f;    // frames, tiles, cards
  float radiusLg = 10.0f;   // popups, standalone cards
  float hairline = 1.0f;    // 1px borders at scale
  float gap = 8.0f;         // base spacing unit
  float iconSize = 22.0f;   // tool glyph bounding box
  float animSpeed = 14.0f;  // hover lerp, per second
};
const Metrics& metrics();

// ------------------------------------------------------------------- fonts
namespace fonts {
// Loads Inter (UI text) and JetBrains Mono (SMILES/code) from the assets dir.
// Safe to skip: every accessor falls back to ImGui's default font when the
// files are missing, so headless tests never need the TTFs.
void load();
ImFont* body();      // Inter Regular (or fallback)
ImFont* semibold();  // Inter SemiBold (or body fallback) — headings, emphasis
ImFont* mono();      // JetBrains Mono (or body fallback) — SMILES, formulas
}  // namespace fonts

// Pushes `font` when non-null; returns whether it pushed (pair with PopFont).
bool pushFont(ImFont* font);
void popFont(bool pushed);

}  // namespace chemcad::ui::style

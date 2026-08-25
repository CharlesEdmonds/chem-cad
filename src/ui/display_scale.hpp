#pragma once
// The one factor that turns ChemCAD's design units into physical pixels.
//
// The bench is drawn in two units and no others:
//
//   * the em -- ImGui::GetFontSize(), see ui/layout.hpp -- for anything that
//     sits next to text;
//   * the design pixel, dp(), for the handful of lengths that are not
//     text-relative: stroke widths, hit radii, the sketch canvas' bond length.
//
// Both ride on displayScale(), which is the product of
//
//   * the content scale of the monitor the window is on: Windows display zoom,
//     the GNOME text-scaling-factor, an X server's Xft.dpi -- divided by the
//     framebuffer scale, because on macOS and Wayland the window system has
//     already applied part of that factor to the framebuffer itself;
//   * CHEMCAD_UI_SCALE, the user's own preference.
//
// It is recomputed whenever the window changes monitor or the OS zoom changes,
// so a window keeps its PHYSICAL size across displays instead of its pixel
// size. Nothing in the app may read the OS scale directly; there is one factor
// and this is it.
//
// Deliberately free of ImGui and GLFW: the sketch camera depends on it, and the
// camera is pure geometry that tests compile on its own.

namespace chemcad::ui {

// What the rest of the app can actually lay out. A 4K TV at 200% zoom with a
// 2.0 user preference would ask for 4.0; past that, panels have less than the
// ~46 em they need for their Compact tier and would clip instead of reflow.
inline constexpr float kMinDisplayScale = 0.5f;
inline constexpr float kMaxDisplayScale = 4.0f;

// The window's design size, in design units. Multiplied by the scale and then
// clamped to the monitor's work area by fitWindow().
inline constexpr int kBaseWindowWidth = 1280;
inline constexpr int kBaseWindowHeight = 800;

// Below this the integrated chrome (menu bar, caption buttons, status bar) eats
// the workspace, so edge-resizing stops here.
inline constexpr int kMinWindowWidth = 600;
inline constexpr int kMinWindowHeight = 400;

// Physical pixels per design unit. 1.0 until the application sets it, so
// headless tests see the unscaled geometry they assert against.
float displayScale();
void setDisplayScale(float scale);

// `designPixels` at the current scale. The only sanctioned way to write a
// length that is not em-relative.
float dp(float designPixels);

// CHEMCAD_UI_SCALE, or `fallback` when absent, malformed or out of range.
// Accepts the raw environment string so the parse is testable.
float parseUserScale(const char* env, float fallback);

// contentScale / framebufferScale * userScale, clamped to the range above.
// Non-finite or non-positive inputs fall back to 1.0 rather than poisoning the
// whole layout: a window that is briefly zero-sized during a monitor change
// reports a framebuffer scale of 0.
float resolveDisplayScale(float contentScale, float framebufferScale, float userScale);

// Where to put a window of kBaseWindow* design units on a monitor work area.
// A window that would not fit is shrunk to the work area rather than opening
// with its own title bar and caption buttons off-screen -- which, since the app
// draws its own chrome, would leave it unclosable and unmovable.
struct WindowPlacement {
  int x = 0, y = 0, width = 0, height = 0;
};
WindowPlacement fitWindow(int workX, int workY, int workWidth, int workHeight, int baseWidth,
                          int baseHeight, float scale);

}  // namespace chemcad::ui

// The display-scale contract: ChemCAD has to be usable on whatever surface it
// is opened on -- a 1366x768 laptop panel at 100%, a 4K TV at 200% display
// zoom, a retina panel where the window system already scaled the framebuffer.
//
// The two things that can go visibly wrong are checked here: the factor itself
// (a bad one makes every panel the wrong size), and the first window placement
// (an oversized window puts the app's own caption buttons off-screen, and since
// it draws its own chrome that leaves it unclosable).
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <cmath>

#include "ui/camera.hpp"
#include "ui/display_scale.hpp"

using namespace chemcad;

namespace {

// Every test leaves the process-wide factor as it found it; the camera reads it.
struct ScopedScale {
  explicit ScopedScale(float scale) { ui::setDisplayScale(scale); }
  ~ScopedScale() { ui::setDisplayScale(1.0f); }
  ScopedScale(const ScopedScale&) = delete;
  ScopedScale& operator=(const ScopedScale&) = delete;
};

}  // namespace

TEST_CASE("the OS content scale reaches the layout, and the framebuffer scale is not double-counted") {
  // Windows display zoom / GNOME text scaling: the framebuffer is 1:1 with the
  // window, so the content scale is entirely ours to apply.
  CHECK(ui::resolveDisplayScale(1.0f, 1.0f, 1.0f) == doctest::Approx(1.0f));
  CHECK(ui::resolveDisplayScale(1.5f, 1.0f, 1.0f) == doctest::Approx(1.5f));
  CHECK(ui::resolveDisplayScale(2.0f, 1.0f, 1.25f) == doctest::Approx(2.5f));

  // macOS / Wayland: the window system already scaled the framebuffer, so the
  // app must NOT scale again -- that was the double-size bug.
  CHECK(ui::resolveDisplayScale(2.0f, 2.0f, 1.0f) == doctest::Approx(1.0f));
  CHECK(ui::resolveDisplayScale(2.0f, 2.0f, 1.25f) == doctest::Approx(1.25f));

  // A window mid-resize reports a zero framebuffer; that must not poison the
  // factor with an infinity or a zero.
  CHECK(ui::resolveDisplayScale(1.5f, 0.0f, 1.0f) == doctest::Approx(1.5f));
  CHECK(ui::resolveDisplayScale(0.0f, 1.0f, 1.0f) == doctest::Approx(1.0f));
  const float nan = std::nanf("");
  CHECK(ui::resolveDisplayScale(nan, nan, nan) == doctest::Approx(1.0f));

  // Clamped to what the panels can still lay out.
  CHECK(ui::resolveDisplayScale(8.0f, 1.0f, 2.0f) == doctest::Approx(ui::kMaxDisplayScale));
  CHECK(ui::resolveDisplayScale(0.1f, 1.0f, 0.5f) == doctest::Approx(ui::kMinDisplayScale));
}

TEST_CASE("CHEMCAD_UI_SCALE is honoured, and junk falls back") {
  CHECK(ui::parseUserScale("1.5", 1.25f) == doctest::Approx(1.5f));
  CHECK(ui::parseUserScale("0.75", 1.25f) == doctest::Approx(0.75f));
  CHECK(ui::parseUserScale(nullptr, 1.25f) == doctest::Approx(1.25f));
  CHECK(ui::parseUserScale("", 1.25f) == doctest::Approx(1.25f));
  CHECK(ui::parseUserScale("large", 1.25f) == doctest::Approx(1.25f));
  CHECK(ui::parseUserScale("0", 1.25f) == doctest::Approx(1.25f));
  CHECK(ui::parseUserScale("-2", 1.25f) == doctest::Approx(1.25f));
  CHECK(ui::parseUserScale("99", 1.25f) == doctest::Approx(1.25f));
}

TEST_CASE("dp() and the sketch camera track the display scale") {
  {
    ScopedScale unscaled(1.0f);
    CHECK(ui::dp(8.0f) == doctest::Approx(8.0f));
    ui::Camera2D camera;
    CHECK(camera.scale() == doctest::Approx(ui::kPixelsPerUnit));
  }
  {
    // At 150% zoom a bond must be 1.5x as many pixels, so it keeps its
    // PHYSICAL length instead of shrinking to two thirds.
    ScopedScale scaled(1.5f);
    CHECK(ui::dp(8.0f) == doctest::Approx(12.0f));
    ui::Camera2D camera;
    CHECK(camera.scale() == doctest::Approx(ui::kPixelsPerUnit * 1.5f));
    camera.zoom = 2.0f;
    CHECK(camera.scale() == doctest::Approx(ui::kPixelsPerUnit * 3.0f));
  }
  CHECK(ui::displayScale() == doctest::Approx(1.0f));
}

TEST_CASE("fit-to-window is scale-aware and never lands off-screen") {
  // Roomy 1440p at 100%: the design size, centred.
  {
    const ui::WindowPlacement p = ui::fitWindow(0, 0, 2560, 1440, 1280, 800, 1.0f);
    CHECK(p.width == 1280);
    CHECK(p.height == 800);
    CHECK(p.x == 640);
    CHECK(p.y == 320);
  }
  // 1080p at 150% zoom: the design size would want 1920x1200, which is taller
  // than the screen. It must be clamped, not merely centred.
  {
    const ui::WindowPlacement p = ui::fitWindow(0, 0, 1920, 1040, 1280, 800, 1.5f);
    CHECK(p.width <= 1920);
    CHECK(p.height <= 1040);
    CHECK(p.x >= 0);
    CHECK(p.y >= 0);
    CHECK(p.x + p.width <= 1920);
    CHECK(p.y + p.height <= 1040);
  }
  // A small laptop panel at 200%: entirely filled, origin respected.
  {
    const ui::WindowPlacement p = ui::fitWindow(0, 40, 1366, 728, 1280, 800, 2.0f);
    CHECK(p.width == 1366);
    CHECK(p.height == 728);
    CHECK(p.x == 0);
    CHECK(p.y == 40);
  }
  // A secondary monitor left of the primary has a negative origin; the window
  // has to follow it rather than snapping to the primary.
  {
    const ui::WindowPlacement p = ui::fitWindow(-1920, -200, 1920, 1080, 1280, 800, 1.0f);
    CHECK(p.x == -1920 + (1920 - 1280) / 2);
    CHECK(p.y == -200 + (1080 - 800) / 2);
  }
  // A work area GLFW could not report must still yield a usable window.
  {
    const ui::WindowPlacement p = ui::fitWindow(0, 0, 0, 0, 1280, 800, 1.25f);
    CHECK(p.width == 1600);
    CHECK(p.height == 1000);
  }
}

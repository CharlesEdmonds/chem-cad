#pragma once
// World <-> screen mapping for the sketch canvas.
// One world unit == one standard bond length == kPixelsPerUnit DESIGN pixels at
// zoom 1, so a benzene ring keeps its physical size at any display zoom.

#include "core/model.hpp"
#include "ui/display_scale.hpp"

namespace chemcad::ui {

inline constexpr float kPixelsPerUnit = 55.0f;
inline constexpr float kMinZoom = 0.15f;
inline constexpr float kMaxZoom = 8.0f;

struct Camera2D {
  core::Vec2 pan{0, 0};  // world point shown at the canvas origin
  float zoom = 1.0f;

  float scale() const { return dp(kPixelsPerUnit) * zoom; }

  // `origin` is the canvas top-left in screen space.
  core::Vec2 worldToScreen(core::Vec2 w, core::Vec2 origin) const {
    return {origin.x + (w.x - pan.x) * scale(), origin.y - (w.y - pan.y) * scale()};
  }
  core::Vec2 screenToWorld(core::Vec2 s, core::Vec2 origin) const {
    return {pan.x + (s.x - origin.x) / scale(), pan.y - (s.y - origin.y) / scale()};
  }

  void zoomAt(float factor, core::Vec2 screenPoint, core::Vec2 origin);
  // Fit all atoms of the document into a canvas of `sizePx`, with margin.
  void fit(const core::Document&, core::Vec2 sizePx);
};

}  // namespace chemcad::ui

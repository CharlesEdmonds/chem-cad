#include "ui/camera.hpp"

#include <algorithm>
#include <limits>

namespace chemcad::ui {

void Camera2D::zoomAt(float factor, core::Vec2 screenPoint, core::Vec2 origin) {
  const core::Vec2 before = screenToWorld(screenPoint, origin);
  zoom = std::clamp(zoom * factor, kMinZoom, kMaxZoom);
  const core::Vec2 after = screenToWorld(screenPoint, origin);
  pan.x += before.x - after.x;
  pan.y += before.y - after.y;
}

void Camera2D::fit(const core::Document& doc, core::Vec2 sizePx) {
  float minX = std::numeric_limits<float>::max(), minY = minX;
  float maxX = std::numeric_limits<float>::lowest(), maxY = maxX;
  bool any = false;
  for (const core::Molecule& m : doc.molecules) {
    for (const core::Atom& a : m.atoms()) {
      any = true;
      minX = std::min(minX, a.pos.x);
      maxX = std::max(maxX, a.pos.x);
      minY = std::min(minY, a.pos.y);
      maxY = std::max(maxY, a.pos.y);
    }
  }
  if (!any || sizePx.x <= 1 || sizePx.y <= 1) {
    pan = {-2.0f, -2.0f};
    zoom = 1.0f;
    return;
  }
  const float w = std::max(maxX - minX, 1.0f) + 2.0f;  // 1 unit margin each side
  const float h = std::max(maxY - minY, 1.0f) + 2.0f;
  zoom = std::clamp(std::min(sizePx.x / (w * kPixelsPerUnit), sizePx.y / (h * kPixelsPerUnit)),
                    kMinZoom, kMaxZoom);
  const float cx = (minX + maxX) * 0.5f, cy = (minY + maxY) * 0.5f;
  pan.x = cx - sizePx.x / (2 * scale());
  pan.y = cy + sizePx.y / (2 * scale());
}

}  // namespace chemcad::ui

#include "ui/display_scale.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace chemcad::ui {

namespace {

float gDisplayScale = 1.0f;

bool usable(const float value) { return std::isfinite(value) && value > 0.0f; }

}  // namespace

float displayScale() { return gDisplayScale; }

void setDisplayScale(const float scale) {
  gDisplayScale = usable(scale) ? std::clamp(scale, kMinDisplayScale, kMaxDisplayScale) : 1.0f;
}

float dp(const float designPixels) { return designPixels * gDisplayScale; }

float parseUserScale(const char* env, const float fallback) {
  if (env == nullptr || *env == '\0') return fallback;
  char* end = nullptr;
  const float parsed = std::strtof(env, &end);
  if (end == env || !usable(parsed)) return fallback;
  if (parsed < kMinDisplayScale || parsed > kMaxDisplayScale) return fallback;
  return parsed;
}

float resolveDisplayScale(const float contentScale, const float framebufferScale,
                          const float userScale) {
  const float content = usable(contentScale) ? contentScale : 1.0f;
  const float framebuffer = usable(framebufferScale) ? framebufferScale : 1.0f;
  const float user = usable(userScale) ? userScale : 1.0f;
  return std::clamp(content / framebuffer * user, kMinDisplayScale, kMaxDisplayScale);
}

WindowPlacement fitWindow(const int workX, const int workY, const int workWidth,
                          const int workHeight, const int baseWidth, const int baseHeight,
                          const float scale) {
  const float s = usable(scale) ? std::clamp(scale, kMinDisplayScale, kMaxDisplayScale) : 1.0f;
  WindowPlacement placement;
  placement.width = std::max(1, static_cast<int>(std::lround(static_cast<float>(baseWidth) * s)));
  placement.height = std::max(1, static_cast<int>(std::lround(static_cast<float>(baseHeight) * s)));
  if (workWidth > 0) placement.width = std::min(placement.width, workWidth);
  if (workHeight > 0) placement.height = std::min(placement.height, workHeight);
  placement.x = workX + std::max(0, (workWidth - placement.width) / 2);
  placement.y = workY + std::max(0, (workHeight - placement.height) / 2);
  return placement;
}

}  // namespace chemcad::ui

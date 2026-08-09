#include "gfx/camera3d.hpp"

#include <algorithm>
#include <cmath>

namespace chemcad::gfx {
namespace {

constexpr float kPi = 3.14159265358979323846f;

float radians(float degrees) { return degrees * (kPi / 180.0f); }

Mat4 multiply(const Mat4& a, const Mat4& b) {
  Mat4 result{};
  for (int column = 0; column < 4; ++column) {
    for (int row = 0; row < 4; ++row) {
      float value = 0.0f;
      for (int k = 0; k < 4; ++k) value += a[k * 4 + row] * b[column * 4 + k];
      result[column * 4 + row] = value;
    }
  }
  return result;
}

std::array<float, 3> subtract(const std::array<float, 3>& a,
                              const std::array<float, 3>& b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

float dot(const std::array<float, 3>& a, const std::array<float, 3>& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

std::array<float, 3> cross(const std::array<float, 3>& a,
                           const std::array<float, 3>& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

std::array<float, 3> normalise(const std::array<float, 3>& value) {
  const float length = std::sqrt(dot(value, value));
  if (length <= 1.0e-12f) return {0.0f, 0.0f, 0.0f};
  return {value[0] / length, value[1] / length, value[2] / length};
}

}  // namespace

void Camera3D::orbit(float dxPixels, float dyPixels) {
  // A quarter degree per pixel keeps a full-width drag useful without making
  // fine inspection twitchy. Pitch remains away from the lookAt pole, where
  // the camera right vector would otherwise become undefined.
  yawDeg = std::remainder(yawDeg + 0.25f * dxPixels, 360.0f);
  pitchDeg = std::clamp(pitchDeg + 0.25f * dyPixels, -85.0f, 85.0f);
}

void Camera3D::zoom(float wheelSteps) {
  // Exponential zoom makes equal wheel travel feel equal at every scale. frame()
  // keeps its fit distance as the geometric midpoint of the clip planes, so the
  // limits scale with this vessel instead of assuming an absolute metre range.
  const float frameDistance =
      std::sqrt(std::max(nearM, 1.0e-6f) * std::max(farM, 1.0e-6f));
  distanceM = std::clamp(distanceM * std::exp(-0.12f * wheelSteps),
                         0.55f * frameDistance, 2.8f * frameDistance);
}

void Camera3D::frame(double vesselHeightM, double maxRadiusM, float aspect) {
  const float height = std::max(0.001f, static_cast<float>(vesselHeightM));
  const float radius = std::max(0.0005f, static_cast<float>(maxRadiusM));
  // The orbit belongs at the vessel mid-height; targeting its base is the
  // perceptual reason a correctly sized apparatus still looks off-centre.
  targetZM = 0.5f * height;

  const float verticalHalfFov = 0.5f * radians(std::clamp(fovDeg, 5.0f, 150.0f));
  const float tangent = std::max(std::tan(verticalHalfFov), 1.0e-3f);
  const float safeAspect = std::max(aspect, 0.05f);

  // Fit height and diameter independently in perspective, then honour the
  // tighter axis. This fills wide viewports without the over-conservative
  // enclosing sphere while a 12 percent margin keeps the glass rim in frame.
  const float verticalDistance = (0.5f * height) / tangent;
  const float horizontalDistance = radius / (tangent * safeAspect);
  distanceM = std::max(0.0032f, 1.12f * std::max(verticalDistance, horizontalDistance));

  // Symmetric logarithmic clip limits both contain every allowed zoom and retain
  // the framing distance as sqrt(near*far), without adding hidden camera state.
  nearM = distanceM / 32.0f;
  farM = distanceM * 32.0f;
}

std::array<float, 3> Camera3D::eye() const {
  const float yaw = radians(yawDeg);
  const float pitch = radians(pitchDeg);
  const float horizontal = distanceM * std::cos(pitch);
  return {horizontal * std::cos(yaw), horizontal * std::sin(yaw),
          targetZM + distanceM * std::sin(pitch)};
}

Mat4 Camera3D::view() const {
  const std::array<float, 3> eyePosition = eye();
  const std::array<float, 3> target{0.0f, 0.0f, targetZM};
  const std::array<float, 3> forward = normalise(subtract(target, eyePosition));
  const std::array<float, 3> side = normalise(cross(forward, {0.0f, 0.0f, 1.0f}));
  const std::array<float, 3> up = cross(side, forward);

  return {side[0], up[0], -forward[0], 0.0f,
          side[1], up[1], -forward[1], 0.0f,
          side[2], up[2], -forward[2], 0.0f,
          -dot(side, eyePosition), -dot(up, eyePosition), dot(forward, eyePosition), 1.0f};
}

Mat4 Camera3D::projection(float aspect) const {
  const float safeAspect = std::max(aspect, 0.05f);
  const float nearPlane = std::max(nearM, 1.0e-4f);
  const float farPlane = std::max(farM, nearPlane + 1.0e-3f);
  const float cotangent = 1.0f / std::tan(0.5f * radians(std::clamp(fovDeg, 5.0f, 150.0f)));

  Mat4 result{};
  result[0] = cotangent / safeAspect;
  result[5] = cotangent;
  result[10] = (farPlane + nearPlane) / (nearPlane - farPlane);
  result[11] = -1.0f;
  result[14] = (2.0f * farPlane * nearPlane) / (nearPlane - farPlane);
  return result;
}

Mat4 Camera3D::viewProjection(float aspect) const {
  return multiply(projection(aspect), view());
}

}  // namespace chemcad::gfx

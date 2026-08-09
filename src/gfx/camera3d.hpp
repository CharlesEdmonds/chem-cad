#pragma once
// Orbit camera for the 3D fluid stage. Pure maths, no GL and no ImGui, so the
// panel can run its input handling in the headless tests where the renderer is
// absent.
//
// Convention matches the fluid solver: vessel-local metres, +z up. The camera
// orbits a target on the vessel axis, which keeps the apparatus centred while
// the user spins it.

#include <array>

namespace chemcad::gfx {

using Mat4 = std::array<float, 16>;  // column-major, as GL expects

struct Camera3D {
  float yawDeg = 35.0f;
  float pitchDeg = 14.0f;
  float distanceM = 0.55f;      // eye distance from the target
  float targetZM = 0.09f;       // orbit centre height on the vessel axis
  float fovDeg = 38.0f;
  float nearM = 0.02f;
  float farM = 3.0f;

  // Applies a drag in pixels to the orbit angles, with pitch clamped so the
  // camera can never flip over the pole.
  void orbit(float dxPixels, float dyPixels);
  // Wheel zoom, multiplicative, clamped to a sane working range.
  void zoom(float wheelSteps);
  // Frames a vessel of this height and radius so it fills the viewport.
  void frame(double vesselHeightM, double maxRadiusM, float aspect);

  std::array<float, 3> eye() const;
  Mat4 view() const;
  Mat4 projection(float aspect) const;
  Mat4 viewProjection(float aspect) const;
};

}  // namespace chemcad::gfx

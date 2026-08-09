#pragma once
// The vessel's rigid motion, and the fictitious accelerations it produces.
//
// The solver integrates in VESSEL coordinates: the boundary SDF, its density
// table and the neighbour grid then all stay in a fixed frame while the vessel
// shakes, tilts or is inverted. Moving the wall through a world-space particle
// cloud instead would rebuild boundary neighbourhoods every substep and is
// markedly less stable.
//
// With x_world = p(t) + R(t) q, the acceleration a particle feels in the
// vessel frame is
//
//   a = R^T g_world                (gravity, rotated into the vessel)
//     - R^T p''(t)                 (translational forcing: this is the shake)
//     - 2 omega x u                (Coriolis)
//     - alpha x q                  (Euler)
//     - omega x (omega x q)        (centrifugal)
//
// where u is the vessel-frame velocity, omega the body-frame angular velocity
// and alpha its derivative. p''(t) comes from the analytic shake law, never
// from finite-differencing UI positions: differencing a mouse path injects
// noise straight into the momentum equation.
//
// This is what makes VERTICAL shaking possible at all: a vertical p''(t) is
// simply another component of the same term, whereas the old 2D model could
// only fake a horizontal slosh.

#include <array>

namespace chemcad::fluid {

// Rigid pose as a translation plus a unit quaternion (w, x, y, z).
struct Pose {
  std::array<double, 3> position{0.0, 0.0, 0.0};
  std::array<double, 4> orientation{1.0, 0.0, 0.0, 0.0};
};

// What the user is doing to the vessel. `shake*` describe a driven oscillation
// along an arbitrary axis; `manual*` describe a hand-driven motion, which the
// UI supplies as smoothed acceleration rather than as raw pointer positions.
struct VesselMotion {
  bool shaking = false;
  double shakeRemainingS = 0.0;
  double shakeFrequencyHz = 3.0;
  double shakeAmplitudeM = 0.05;
  std::array<double, 3> shakeAxis{0.0, 0.0, 1.0};  // unit; z = vertical shake

  // Hand motion: acceleration in world coordinates, m/s^2, already smoothed.
  std::array<double, 3> manualAcceleration{0.0, 0.0, 0.0};

  // Vessel attitude: tilt about the horizontal axes, and the inversion the
  // user performs to vent or to drain from the neck.
  Pose pose;
  std::array<double, 3> angularVelocity{0.0, 0.0, 0.0};      // body frame, rad/s
  std::array<double, 3> angularAcceleration{0.0, 0.0, 0.0};  // body frame, rad/s^2
};

// The per-substep frame terms the solver needs, evaluated at time `t`.
struct FrameAcceleration {
  std::array<double, 3> uniform{0.0, 0.0, -9.80665};  // gravity + translational
  std::array<double, 3> omega{0.0, 0.0, 0.0};
  std::array<double, 3> alpha{0.0, 0.0, 0.0};
  bool rotating = false;  // false lets the solver skip the position-dependent terms
};

// Analytic translational acceleration of a driven shake at time t:
//   p(t)   = A sin(2 pi f t) axis
//   p''(t) = -A (2 pi f)^2 sin(2 pi f t) axis
std::array<double, 3> shakeAcceleration(const VesselMotion&, double timeS);

// Assembles the frame terms: rotates gravity into the vessel, subtracts the
// translational forcing (shake plus hand motion) and carries the rotation
// rates through for the Coriolis/Euler/centrifugal terms.
FrameAcceleration frameAcceleration(const VesselMotion&, double timeS);

// Peak speed and specific power of a driven shake, reported to the UI so the
// displayed numbers are the ones the physics used:
//   u = 2 pi f A,  epsilon = u^2 f / 2  (W/kg)
double shakePeakVelocity(const VesselMotion&);
double shakeSpecificPower(const VesselMotion&);

}  // namespace chemcad::fluid

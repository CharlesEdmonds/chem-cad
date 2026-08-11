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

  // Hand motion: acceleration in world coordinates, m/s^2, already smoothed,
  // and the world displacement of the vessel that produced it. The solver only
  // needs the acceleration; the displacement is carried so the renderer can
  // draw the glassware where the hand actually put it, instead of leaving a
  // stationary vessel with its contents mysteriously sloshing.
  std::array<double, 3> manualAcceleration{0.0, 0.0, 0.0};
  std::array<double, 3> manualOffset{0.0, 0.0, 0.0};

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

// Displacement of the same driven shake, p(t) = A sin(2 pi f t) axis. This is
// the exact integral of shakeAcceleration, so the drawn vessel and the forcing
// the fluid feels can never disagree.
std::array<double, 3> shakeDisplacement(const VesselMotion&, double timeS);

// Where a hand-held vessel actually is, given where the hand is being told to
// go. The pointer commands the hand; the vessel tracks it as a spring-damper.
// Mapping a pointer delta straight to an acceleration instead gives a hand
// stroke almost no authority and amplifies the pointer's integer jitter twice
// over, which is why the UI drives this rather than differencing positions.
//
// Two regimes, because holding a funnel and letting go of one are different
// mechanical problems:
//
//   held     6.0 Hz, zeta 0.70 -- clamped in a hand. A 3 cm wiggle is felt, a
//                                 steady drag is almost nothing.
//   released 1.2 Hz, zeta 0.85 -- nothing constrains it. The vessel keeps the
//                                 momentum of the throw and coasts back to the
//                                 bench in about half a second.
//
// `excursionLimit` is deliberately far larger than the vessel: it is allowed to
// be thrown clean off the stage and swing back.
//
// `accelerationLimit` is what the FLUID feels, and it is a physical bound, not
// a safety net. The reference bench shake this solver is written against --
// 50 mm at 3 Hz -- peaks at A(2 pi f)^2 = 17.8 m/s^2, or 1.8 g, and that is
// about as hard as a hand can drive a funnel. The limit sat at 8 g, so flinging
// the vessel across the stage handed the liquid 7.3 g of coherent forcing:
// measured, the contents left at 2 m/s and looked detonated. Above roughly 2 g
// the clamp is not protecting the solve, it is inventing forcing the user never
// applied.
//
// It binds the spring too. A 0.6 m throw at 1.2 Hz wants omega^2 x = 34 m/s^2
// coming home, so the clamp is what turns the last part of a long return from a
// snap into a glide -- and the vessel is drawn from the same integration, so
// what the user sees and what the liquid feels stay the same motion.
struct HandFollower {
  std::array<double, 3> hand{0.0, 0.0, 0.0};          // commanded hand position, m
  std::array<double, 3> position{0.0, 0.0, 0.0};      // vessel position, m
  std::array<double, 3> velocity{0.0, 0.0, 0.0};      // vessel velocity, m/s
  std::array<double, 3> acceleration{0.0, 0.0, 0.0};  // vessel acceleration, m/s^2

  double heldHz = 6.0;
  double heldDampingRatio = 0.7;
  double releasedHz = 1.2;
  double releasedDampingRatio = 0.85;
  double excursionLimit = 0.60;                  // m, per axis
  double accelerationLimit = 2.0 * 9.80665;      // m/s^2, magnitude

  // Integrates one frame. `handDelta` is this frame's commanded hand movement
  // in world metres and is ignored while `held` is false, because a released
  // vessel has no hand to follow.
  void advance(const std::array<double, 3>& handDelta, bool held, double dt);

  // True once the vessel has come back to rest at the origin.
  bool atRest() const;

  void reset();
};

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

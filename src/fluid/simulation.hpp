#pragma once
// The public face of the 3D fluid: what the UI and the chemistry talk to.
//
// One object owns the particle state, the vessel boundary, the solver, and the
// diagnostics, and exposes exactly the operations a bench chemist performs:
// charge the vessel, shake it (in any direction, including vertically), move
// it by hand, tilt or invert it, let it stand, drain it.
//
// Threading: `advance` is pure computation and may run on a worker thread.
// Rendering reads `snapshot()`, which returns an immutable copy-on-write view
// of the last completed state, so the UI never observes a half-stepped
// solve and never blocks the physics.

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "fluid/diagnostics.hpp"
#include "fluid/frame.hpp"
#include "fluid/solver.hpp"
#include "fluid/types.hpp"
#include "fluid/vessel_sdf.hpp"
#include "sol/funnel.hpp"

namespace chemcad::fluid {

// What the renderer needs, and nothing else: positions, phase, speed for
// shading, plus the geometry and diagnostics it labels the scene with.
struct Snapshot {
  std::vector<float> px, py, pz;
  std::vector<float> speed;     // |v|, m/s, for velocity-tinted shading
  std::vector<float> colour;    // smoothed phase indicator, 0..1
  std::vector<uint8_t> phase;
  std::vector<PhaseMaterial> phases;
  Diagnostics diagnostics;
  Pose pose;                    // vessel attitude, for the glass and the camera
  double vesselHeightM = 0.19;
  double maxRadiusM = 0.045;
  // Solver lattice spacing (dx). The primary quantity: everything else that
  // needs a length scale derives from it rather than back-computing it.
  double particleSpacingM = 4.0e-3;
  // Radius the surface is RENDERED at, which is deliberately not dx/2. An SPH
  // free surface is the level set of the smoothing kernel, so the surface sits
  // at roughly half the kernel support (support = 2 dx), i.e. one dx. Drawing
  // tangent dx/2 spheres instead leaves a crease at every particle and the
  // screen-space curvature flow cannot close them, so the liquid reads as a
  // cluster of orbs. At one dx the union is hole-free for any packing at least
  // as dense as the cubic lattice (which needs dx*sqrt(3)/2) and the surface
  // smooths to a continuous meniscus.
  double particleRadiusM = 4.0e-3;
  double elapsedS = 0.0;
  sol::Vessel vessel = sol::Vessel::SeparatoryFunnel;
  uint64_t revision = 0;        // bumped on every completed advance
};

class Simulation {
 public:
  Simulation();
  ~Simulation();
  Simulation(const Simulation&) = delete;
  Simulation& operator=(const Simulation&) = delete;

  // ---- setup ----------------------------------------------------------
  void setVessel(sol::Vessel, double ratedVolumeMl);
  void setResolution(double spacingM);   // rebuilds and recharges, preserving quality budget
  void setQuality(const QualityProfile&); // atomically applies spacing and pressure budget
  void setPhases(const std::vector<PhaseMaterial>&, const std::vector<double>& sigmaPairs);

  // Charges the vessel: lays a lattice of particles, dense phase at the
  // bottom, and clears all motion. Deterministic.
  void charge();
  bool charged() const;

  // ---- driving --------------------------------------------------------
  // Arms a driven shake along `axis` (unit; {0,0,1} is vertical) for
  // `durationS` at `frequencyHz` and stroke half-amplitude `amplitudeM`.
  void shake(const std::array<double, 3>& axis, double durationS, double frequencyHz,
             double amplitudeM);
  // Hand motion: world-frame acceleration, already smoothed by the caller.
  void setManualAcceleration(const std::array<double, 3>&);
  // Vessel attitude, for tilting and inverting.
  void setPose(const Pose&, const std::array<double, 3>& angularVelocity,
               const std::array<double, 3>& angularAcceleration);

  // Queues simulated time for the dedicated physics worker. Requests coalesce
  // into a bounded backlog, but accepted time is never scaled down to disguise
  // a slow solve, so this call remains non-blocking and the reported rate honest.
  void requestAdvance(double simulatedSeconds);
  bool stepping() const;
  double pendingSeconds() const;
  void waitForIdle();

  // Synchronous compatibility path for tests and batch callers. It uses the
  // same worker integration path as requestAdvance(), then waits for completion.
  void advance(double dt);

  // ---- reading --------------------------------------------------------
  // Immutable view of the last completed state. Cheap: shares storage until
  // the next advance completes.
  std::shared_ptr<const Snapshot> snapshot() const;
  const Diagnostics& diagnostics() const;
  const Solver::Stats& solverStats() const;
  double elapsedS() const;
  bool shaking() const;
  // Ratio of simulated time completed to wall time spent on the latest solve.
  // Values below one mean physics is genuinely advancing slower than real time.
  double realTimeFactor() const;

  // Total charged volume, mL, and the per-phase charge.
  double totalVolumeMl() const;

  // Human-readable description of what the solver is doing, for the UI's
  // honesty panel: resolution, particles, substeps, compression, deficit, cost.
  std::string statusLine() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace chemcad::fluid

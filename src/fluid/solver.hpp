#pragma once
// Predictive-corrective incompressible SPH (PCISPH) for two immiscible
// liquids, on the CPU, in the vessel frame.
//
// Formulation, and why:
//   * PCISPH (Solenthaler & Pajarola, ACM TOG 28(3), 2009) rather than WCSPH:
//     a Tait equation of state ties the timestep to an artificial speed of
//     sound (~70 us here), while PCISPH's prediction-correction loop reaches
//     the same incompressibility at a few ms.
//   * Density-contrast SPH (Solenthaler & Pajarola, SCA 2008, eqs. 8-14) for
//     the multiphase part: densities are computed as a NUMBER density
//     delta_i = sum_j W_ij, and each particle converts it with its own mass,
//     rho_i = m_i delta_i. Standard mass-density SPH smooths the jump at a
//     phase boundary and produces a spurious gap and interface pressure.
//     Pressure force uses their corrected symmetric form
//         F_i = -sum_j (p_i + p_j)/(2 delta_i delta_j) grad W_ij
//     which is exactly antisymmetric per pair, so momentum is conserved even
//     with unequal masses.
//   * Physical dynamic viscosity with a harmonic pair mean, so the stress at
//     the interface is limited by the less viscous liquid rather than by an
//     arithmetic average.
//   * Interfacial tension by the Continuum Surface Force (Brackbill, Kothe &
//     Zemach, JCP 100, 1992; SPH form after Morris, IJNMF 33, 2000) on the
//     colour field: f = sigma kappa grad c. |grad c| integrates to the jump in
//     c across the interface, which is 1, so this reproduces dp = 2 sigma / R
//     with sigma in N/m and nothing to fit. tests/test_fluid_solver.cpp holds
//     it to that law directly.
//
// Determinism: every reduction accumulates in double over ascending particle
// index, neighbour loops gather (never scatter), and no random number touches
// the solve. Two runs from the same charge and the same motion sequence
// produce bit-identical state.

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "fluid/frame.hpp"
#include "fluid/grid.hpp"
#include "fluid/types.hpp"
#include "fluid/vessel_sdf.hpp"

namespace chemcad::fluid {

// Tunables with their justification. Defaults are the values validated by
// tests/test_fluid_solver.cpp.
struct SolverConfig {
  Resolution resolution;
  double densityTolerance = 5.0e-3;  // maximum compression (delta/delta0 - 1)
  int minPressureIterations = 3;
  int maxPressureIterations = 32;
  double cflNumber = 0.4;            // dt <= cfl * H / v_max
  double accelerationSafety = 0.25;  // dt <= safety * sqrt(H / a_max)
  double maxSubstepS = 1.0 / 60.0;  // hard sanity ceiling; resolution sets the normal cap
  double maxSpeed = 4.0;             // m/s, hard clamp against blow-up
  double contactRadiusFactor = 0.35; // wall contact distance, in units of spacing
  double wallFriction = 0.05;        // tangential loss at the glass
  double xsphSmoothing = 0.0;        // 0 = off; display-only velocity smoothing
  bool enableSurfaceTension = true;
  bool enableCoriolis = true;
};

// A named resolution/accuracy budget. Spacing and pressure convergence move
// together so callers cannot accidentally label a coarse, loose preview as a
// quality solve (or pay quality iteration costs at preview resolution).
struct QualityProfile {
  // Release, i7-9750H, one concurrent CPU load thread, 200 mL charge:
  // Interactive: 16.14 -> 5.83 iterations/substep, 1.52x real time,
  // 3.993% compression. Balanced: 17.64 -> 10.58, 0.210x, 0.575%.
  // Quality: 21.58 -> 23.44, 0.052x, 0.488%; its 0.5% limit stays hard.
  double densityTolerance = 5.0e-3;
  int minPressureIterations = 3;
  int maxPressureIterations = 32;
  double spacing = 4.0e-3;
  bool surfaceTension = true;

  static constexpr QualityProfile interactive() {
    return {4.0e-2, 3, 12, 8.0e-3, true};
  }
  static constexpr QualityProfile balanced() {
    return {1.0e-2, 3, 20, 6.0e-3, true};
  }
  static constexpr QualityProfile quality() {
    return {5.0e-3, 3, 40, 4.0e-3, true};
  }
};

// Interfacial tension, in SI, as the caller measured it.
//
// It is a TWO-PHASE model: the colour field the CSF term differentiates is
// binary, phase 0 against everything else. setPhases rejects a tension table
// spanning more than two phases rather than silently applying it to a colour
// field that cannot represent it.
struct InterfaceModel {
  // sigma[i * phaseCount + j], N/m, symmetric, diagonal unused.
  std::vector<double> sigma;

  // The single interfacial tension this two-phase model carries, N/m.
  double interfacialTension() const { return sigma.size() >= 2 ? sigma[1] : 0.0; }
};

// Thresholds the Continuum Surface Force uses to decide where an interface is,
// when its kernel-gradient correction is trustworthy, and how sharply it is
// allowed to believe the interface curves. They live in the header because the
// GPU mirror in src/gfx/fluid_gpu.cpp compiles them into its shader prelude,
// and two copies of a threshold are two thresholds.
//
// A floor on |grad c| in units of 1/support: below it a particle has no
// interface normal worth having. It gates on the GRADIENT rather than on the
// colour value because the tension is carried by |grad c|, and a colour band
// clips its tails -- measured, a 0.05..0.95 band discarded a quarter of the
// surface delta.
inline constexpr double kInterfaceGradientFloor = 0.01;

// Smallest determinant of the Bonet-Lok correction matrix that may be inverted.
// Below it the neighbourhood is genuinely one-sided -- a wall, a free surface,
// a stray particle -- and inverting a near-singular L does not recover the
// gradient, it multiplies whatever noise is in it. The uncorrected gradient is
// the safer answer there.
inline constexpr double kInterfaceCorrectionDeterminant = 0.25;

// Largest curvature the model may act on, in units of 1/support. An SPH kernel
// cannot resolve an interface radius smaller than its own smoothing length, so
// a larger reported curvature is discretisation noise rather than geometry.
// Without this bound the divergence estimate is unbounded, and sigma kappa
// grad c turns a thin film or an isolated splash particle into an acceleration
// of tens of g: a shaken vessel visibly detonates on release.
inline constexpr double kMaxInterfaceCurvature = 2.0;

class Solver {
 public:
  void configure(const SolverConfig&);
  const SolverConfig& config() const { return config_; }

  // Materials and their pairwise interfacial tensions. Must be called before
  // the first step and after any material edit.
  void setPhases(const std::vector<PhaseMaterial>& phases, const std::vector<double>& sigmaPairs);
  const std::vector<PhaseMaterial>& phases() const { return phases_; }

  // Advances the particle state by exactly `dt` seconds of vessel time,
  // internally splitting it into CFL-limited substeps. `motion` is sampled per
  // substep so a shake keeps its phase. Returns the number of substeps taken.
  int advance(Particles& particles, const VesselBoundary& boundary,
              const VesselMotion& motion, double timeS, double dt);

  // Number density and pressure of the last completed substep, for diagnostics
  // and rendering (both read-only views into the particle arrays).
  double restNumberDensity() const { return delta0_; }

  const InterfaceModel& interfaceModel() const { return interface_; }

  // Interfacial tension the solver's own force field actually produces for a
  // free droplet of this radius, recovered from the radial force balance rather
  // than assumed. Young-Laplace says it should equal sigma; a test says so out
  // loud. Costs one static droplet evaluation, so it belongs in tests and
  // diagnostics, never in a frame.
  double measuredInterfacialTension(double dropletRadiusM) const;

  // Diagnostic counters from the last advance. Density deficit is deliberately
  // separate: pressure is non-negative and cannot fill a free-surface neighbour
  // deficit, so only compression participates in PCISPH convergence.
  struct Stats {
    int substeps = 0;
    int pressureIterations = 0;          // total across substeps
    double maxDensityError = 0.0;        // compatibility: equals compression
    double maxDensityCompression = 0.0;  // max(delta/delta0 - 1, 0)
    double maxDensityDeficit = 0.0;      // max(1 - delta/delta0, 0), diagnostic
    double maxSpeed = 0.0;               // m/s
    double substepS = 0.0;               // last substep length
    double millisecondsPerSubstep = 0.0; // wall cost over the latest advance
    double neighbourMilliseconds = 0.0;  // compatibility: grid + density
    double gridMilliseconds = 0.0;       // neighbour grid and pair-cache build
    double densityMilliseconds = 0.0;    // initial number-density gather
    double forceMilliseconds = 0.0;      // frame, viscosity, and interface
    double pressureMilliseconds = 0.0;   // prediction-correction iterations
    double integrationMilliseconds = 0.0;
    int clampedParticles = 0;            // hit the speed clamp
    int rejectedSubsteps = 0;            // halved and retried
    int stalledPressureSubsteps = 0;      // correction exited after no material progress
    unsigned workerCount = 1;             // static gather partitions used by this advance
    double pressureStiffnessSubstepS = 0.0; // quantised dt used by active stiffness
    std::uint64_t pressureStiffnessCalibrations = 0; // solver lifetime total
  };
  const Stats& stats() const { return stats_; }

 private:
  SolverConfig config_;
  std::vector<PhaseMaterial> phases_;
  std::vector<double> mass_;      // per phase: rho0 * dx^3
  InterfaceModel interface_;
  NeighbourGrid grid_;
  double delta0_ = 0.0;           // rest number density of the lattice
  std::vector<double> stiffness_; // active PCISPH scaling per phase
  struct StiffnessCacheEntry {
    double substepS = 0.0;
    std::vector<double> values;
  };
  std::vector<StiffnessCacheEntry> stiffnessCache_;
  std::uint64_t pressureStiffnessCalibrations_ = 0;
  double pressureStiffnessSubstepS_ = 0.0;
  Stats stats_;

  // Scratch arrays, kept across calls so a step allocates nothing.
  std::vector<float> predictedX_, predictedY_, predictedZ_;
  std::vector<float> predictedVX_, predictedVY_, predictedVZ_;
  std::vector<float> forceX_, forceY_, forceZ_;
  std::vector<double> densityError_;

  void ensurePressureStiffness(double substepS);
};

}  // namespace chemcad::fluid

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
//   * Interfacial tension from the colour-field curvature of the same
//     density-contrast paper (eqs. 22-24), with a cohesion term after Akinci,
//     Akinci & Teschner (ACM TOG 32(6), 2013). Its coefficient is resolution
//     dependent and is CALIBRATED against the Young-Laplace law
//     (dp = 2 sigma / R) for the working spacing -- it is never set equal to a
//     measured interfacial tension, which would be dimensionally wrong.
//
// Determinism: every reduction accumulates in double over ascending particle
// index, neighbour loops gather (never scatter), and no random number touches
// the solve. Two runs from the same charge and the same motion sequence
// produce bit-identical state.

#include <cstddef>
#include <string>
#include <vector>

#include "fluid/frame.hpp"
#include "fluid/grid.hpp"
#include "fluid/types.hpp"
#include "fluid/vessel_sdf.hpp"

namespace chemcad::fluid {

// Tunables with their justification. Defaults are the values validated by
// tests/test_fluid_solver.cpp; changing one invalidates the calibration cache.
struct SolverConfig {
  Resolution resolution;
  double densityTolerance = 5.0e-3;  // |delta/delta0 - 1| accepted by the pressure loop
  int minPressureIterations = 3;
  int maxPressureIterations = 8;
  double cflNumber = 0.4;            // dt <= cfl * H / v_max
  double accelerationSafety = 0.25;  // dt <= safety * sqrt(H / a_max)
  double maxSubstepS = 1.0 / 480.0;
  double maxSpeed = 4.0;             // m/s, hard clamp against blow-up
  double contactRadiusFactor = 0.35; // wall contact distance, in units of spacing
  double wallFriction = 0.05;        // tangential loss at the glass
  double xsphSmoothing = 0.0;        // 0 = off; display-only velocity smoothing
  bool enableSurfaceTension = true;
  bool enableCoriolis = true;
};

// Per-phase-pair interfacial tension, and the calibration that turns it into
// the solver's cohesion coefficient at this resolution.
struct InterfaceModel {
  // sigma[i * phaseCount + j], N/m, symmetric, diagonal unused.
  std::vector<double> sigma;
  // Calibration factor from the Young-Laplace test: coefficient = sigma * gain.
  double cohesionGain = 0.0;
  bool calibrated = false;
};

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

  // Young-Laplace calibration of the cohesion coefficient for the current
  // resolution and phase pair. Runs a self-contained droplet relaxation; call
  // once per resolution change. Results are cached in the interface model.
  void calibrateInterface(const VesselBoundary& boundary);
  const InterfaceModel& interfaceModel() const { return interface_; }
  const std::string& interfaceCalibrationError() const { return calibrationError_; }

  // Diagnostic counters from the last advance, surfaced so the UI can show
  // that the solver is converging rather than merely running.
  struct Stats {
    int substeps = 0;
    int pressureIterations = 0;      // total across substeps
    double maxDensityError = 0.0;    // final, dimensionless
    double maxSpeed = 0.0;           // m/s
    double substepS = 0.0;           // last substep length
    int clampedParticles = 0;        // hit the speed clamp
    int rejectedSubsteps = 0;        // halved and retried
  };
  const Stats& stats() const { return stats_; }

 private:
  SolverConfig config_;
  std::vector<PhaseMaterial> phases_;
  std::vector<double> mass_;      // per phase: rho0 * dx^3
  InterfaceModel interface_;
  NeighbourGrid grid_;
  double delta0_ = 0.0;           // rest number density of the lattice
  std::vector<double> stiffness_; // PCISPH scaling per phase, calibrated
  Stats stats_;
  std::string calibrationError_;

  // Scratch arrays, kept across calls so a step allocates nothing.
  std::vector<float> predictedX_, predictedY_, predictedZ_;
  std::vector<float> predictedVX_, predictedVY_, predictedVZ_;
  std::vector<float> forceX_, forceY_, forceZ_;
  std::vector<double> densityError_;
};

}  // namespace chemcad::fluid

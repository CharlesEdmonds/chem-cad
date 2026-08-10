#pragma once

#include <memory>
#include <string>
#include <vector>

#include "fluid/solver.hpp"

namespace chemcad::gfx {

// Optional GL 4.3 implementation of the fluid solver's substep hot path.
// Construction and every setup call are safe without a context: available()
// remains false and the owner continues to call fluid::Solver.
class FluidGpuSolver {
 public:
  FluidGpuSolver();
  ~FluidGpuSolver();
  FluidGpuSolver(const FluidGpuSolver&) = delete;
  FluidGpuSolver& operator=(const FluidGpuSolver&) = delete;
  FluidGpuSolver(FluidGpuSolver&&) noexcept;
  FluidGpuSolver& operator=(FluidGpuSolver&&) noexcept;

  bool initialise();
  bool available() const;
  const std::string& unavailableReason() const;

  void configure(const fluid::SolverConfig& config);
  const fluid::SolverConfig& config() const;
  void setPhases(const std::vector<fluid::PhaseMaterial>& phases,
                 const std::vector<double>& sigmaPairs);

  // Upload is explicit because a newly charged or reset particle set may have
  // the same size as the previous one. Once uploaded, state remains resident
  // through every pressure iteration and substep.
  bool upload(const fluid::Particles& particles,
              const fluid::VesselBoundary& boundary);
  bool resident() const;

  // Advances the resident state. Only tiny reduction counters cross to the CPU
  // for CFL and PCISPH convergence decisions.
  int advance(const fluid::VesselBoundary& boundary,
              const fluid::VesselMotion& motion, double timeS, double dt);

  // The caller invokes this only when publishing a rendering/diagnostic
  // snapshot, not after each substep.
  bool download(fluid::Particles& particles) const;

  const fluid::Solver::Stats& stats() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace chemcad::gfx

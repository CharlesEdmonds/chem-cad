#pragma once
// A device-side stand-in for Solver::advance, supplied by a layer that owns a
// graphics device.
//
// It is an interface, and it lives here, because chemcad_gfx already depends on
// chemcad_fluid: naming the concrete OpenGL implementation from this side would
// close a cycle. The fluid module stays pure CPU, headless-testable, and the
// reference against which any accelerator is judged.
//
// Threading is the whole reason this is not just a function pointer. Simulation
// owns a dedicated physics worker, and a graphics context belongs to one thread
// at a time -- the renderer keeps the main one. Every method below is therefore
// called on the worker thread, between bind() and unbind(), which is where an
// implementation makes its own context current.

#include <cstddef>
#include <string>
#include <vector>

#include "fluid/frame.hpp"
#include "fluid/solver.hpp"
#include "fluid/types.hpp"
#include "fluid/vessel_sdf.hpp"

namespace chemcad::fluid {

class Accelerator {
 public:
  virtual ~Accelerator() = default;

  // Once, on the worker thread, before anything else. A false return retires
  // the accelerator permanently and the simulation keeps the CPU solver: the
  // reasons bind() can fail -- no device, no compute support, a shader that
  // will not compile -- do not change while the program runs, so retrying would
  // only pay for the same failure every step.
  virtual bool bind() = 0;
  virtual void unbind() = 0;

  // Whether a charge of this size is worth sending to the device at all. Below
  // a few hundred particles the dispatch and convergence-readback traffic costs
  // more than the parallelism buys, and the CPU solver is simply faster;
  // measured on an RTX 2070 Max-Q the crossover sits near 600 particles.
  virtual bool worthwhile(std::size_t particleCount) const = 0;

  virtual bool configure(const SolverConfig&, const std::vector<PhaseMaterial>&,
                         const std::vector<double>& sigmaPairs) = 0;
  virtual bool upload(const Particles&, const VesselBoundary&) = 0;

  // Substeps taken, or 0 when the device could not advance -- in which case the
  // caller falls back to the CPU solver for this step and re-uploads before the
  // next one, so a transient device failure costs a frame rather than the run.
  virtual int advance(const VesselBoundary&, const VesselMotion&, double timeS, double dt) = 0;
  virtual bool download(Particles&) = 0;

  virtual const Solver::Stats& stats() const = 0;
  // Human-readable reason the device is unavailable or last failed, for the
  // status line. Empty when nothing has gone wrong.
  virtual const std::string& error() const = 0;
  // Short device identification for the same status line.
  virtual const std::string& description() const = 0;
};

}  // namespace chemcad::fluid

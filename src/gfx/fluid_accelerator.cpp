#include "gfx/fluid_accelerator.hpp"

#include <string>
#include <utility>

#include "gfx/compute.hpp"
#include "gfx/fluid_gpu.hpp"
#include "gfx/gl_api.hpp"

namespace chemcad::gfx {
namespace {

class GpuAccelerator final : public fluid::Accelerator {
 public:
  GpuAccelerator(std::function<bool()> makeCurrent, std::function<void()> clearCurrent,
                 std::size_t minimumParticles)
      : makeCurrent_(std::move(makeCurrent)),
        clearCurrent_(std::move(clearCurrent)),
        minimumParticles_(minimumParticles) {}

  bool bind() override {
    if (!makeCurrent_ || !makeCurrent_()) {
      error_ = "no second OpenGL context for the physics thread";
      return false;
    }
    current_ = true;
    const ComputeCapabilities capabilities = probeCompute();
    if (!capabilities.available) {
      error_ = capabilities.reason;
      release();
      return false;
    }
    if (!solver_.initialise()) {
      error_ = solver_.unavailableReason();
      release();
      return false;
    }
    description_ = std::string(glRendererString());
    error_.clear();
    return true;
  }

  void unbind() override { release(); }

  bool worthwhile(std::size_t particleCount) const override {
    return particleCount >= minimumParticles_;
  }

  bool configure(const fluid::SolverConfig& config,
                 const std::vector<fluid::PhaseMaterial>& phases,
                 const std::vector<double>& sigmaPairs) override {
    try {
      solver_.configure(config);
      solver_.setPhases(phases, sigmaPairs);
    } catch (const std::exception& failure) {
      // setPhases validates the material table and throws on a configuration
      // the model cannot represent. That is the CPU solver's contract too, and
      // the caller's job is to fall back, not to crash the physics worker.
      error_ = failure.what();
      return false;
    }
    error_.clear();
    return true;
  }

  bool upload(const fluid::Particles& particles,
              const fluid::VesselBoundary& boundary) override {
    if (!solver_.upload(particles, boundary)) {
      error_ = solver_.unavailableReason();
      return false;
    }
    return true;
  }

  int advance(const fluid::VesselBoundary& boundary, const fluid::VesselMotion& motion,
              double timeS, double dt) override {
    const int substeps = solver_.advance(boundary, motion, timeS, dt);
    if (substeps <= 0) error_ = solver_.unavailableReason();
    return substeps;
  }

  bool download(fluid::Particles& particles) override {
    if (!solver_.download(particles)) {
      error_ = solver_.unavailableReason();
      return false;
    }
    return true;
  }

  const fluid::Solver::Stats& stats() const override { return solver_.stats(); }
  const std::string& error() const override { return error_; }
  const std::string& description() const override { return description_; }

 private:
  void release() {
    if (!current_) return;
    current_ = false;
    if (clearCurrent_) clearCurrent_();
  }

  FluidGpuSolver solver_;
  std::function<bool()> makeCurrent_;
  std::function<void()> clearCurrent_;
  std::size_t minimumParticles_ = 600;
  bool current_ = false;
  std::string error_;
  std::string description_;
};

}  // namespace

std::shared_ptr<fluid::Accelerator> makeFluidAccelerator(std::function<bool()> makeCurrent,
                                                         std::function<void()> clearCurrent,
                                                         std::size_t minimumParticles) {
  return std::make_shared<GpuAccelerator>(std::move(makeCurrent), std::move(clearCurrent),
                                          minimumParticles);
}

}  // namespace chemcad::gfx

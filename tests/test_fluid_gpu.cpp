#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <GLFW/glfw3.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

#include "fluid/solver.hpp"
#include "fluid/vessel_sdf.hpp"
#include "gfx/compute.hpp"
#include "gfx/fluid_gpu.hpp"
#include "gfx/gl_api.hpp"
#include "sol/funnel.hpp"

// Exercises the OpenGL compute fluid backend on whatever GPU is actually
// present. It creates its own hidden context rather than borrowing the
// application's, so the whole backend -- probe, shader compilation, buffer
// residency, and every dispatch -- is covered without a window on screen.
//
// A machine with no GL 4.3 device is a legitimate configuration, not a failure:
// the runtime keeps the CPU solver in that case, and so does this suite. What
// it will not tolerate is a GPU that IS present, IS claimed, and then produces
// a different fluid from the reference implementation.

namespace fluid = chemcad::fluid;
namespace gfx = chemcad::gfx;
namespace sol = chemcad::sol;

namespace {

constexpr double kFrameS = 1.0 / 60.0;

fluid::PhaseMaterial phase(const char* label, double density, double viscosity, double volumeMl) {
  fluid::PhaseMaterial value;
  value.label = label;
  value.restDensity = density;
  value.dynamicViscosity = viscosity;
  value.volumeMl = volumeMl;
  return value;
}

// The default motion is a still, upright vessel under gravity, which is the
// state the funnel spends most of its life in.
fluid::VesselMotion atRest() { return fluid::VesselMotion{}; }

// Owns the hidden context for the whole binary. GLFW and the GL entry-point
// table are process-global, so creating one per test case would reload the
// loader against a context that had already been destroyed.
class GlFixture {
 public:
  static GlFixture& instance() {
    static GlFixture fixture;
    return fixture;
  }

  bool ready() const { return ready_; }
  const std::string& reason() const { return reason_; }

 private:
  GlFixture() {
    if (glfwInit() != GLFW_TRUE) {
      reason_ = "GLFW could not initialise (no display or driver)";
      return;
    }
    initialised_ = true;
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    window_ = glfwCreateWindow(64, 64, "chemcad gpu fluid", nullptr, nullptr);
    if (window_ == nullptr) {
      reason_ = "no OpenGL 4.3 core context available";
      return;
    }
    glfwMakeContextCurrent(window_);
    if (!gfx::loadGl(reinterpret_cast<gfx::GlProcLoader>(glfwGetProcAddress))) {
      reason_ = "required OpenGL entry points are missing";
      return;
    }
    ready_ = true;
  }

  ~GlFixture() {
    if (window_ != nullptr) glfwDestroyWindow(window_);
    if (initialised_) glfwTerminate();
  }

  GLFWwindow* window_ = nullptr;
  bool initialised_ = false;
  bool ready_ = false;
  std::string reason_;
};

struct Charge {
  fluid::SolverConfig config;
  std::vector<fluid::PhaseMaterial> phases;
  fluid::VesselBoundary boundary;
  fluid::Particles particles;
};

Charge separatoryCharge(double spacing) {
  Charge charge;
  charge.config = fluid::SolverConfig{};
  charge.config.resolution.spacing = spacing;
  charge.config.densityTolerance = 1.0e-2;
  charge.config.maxPressureIterations = 12;
  charge.phases = {phase("aqueous", 998.0, 1.0e-3, 60.0),
                   phase("organic", 850.0, 0.6e-3, 60.0)};
  charge.boundary.build(sol::Vessel::SeparatoryFunnel, 0.19,
                        charge.config.resolution.support(), charge.config.resolution.spacing);
  charge.boundary.chargeLattice(charge.phases, charge.config.resolution.spacing,
                                charge.particles);
  return charge;
}

struct Summary {
  double meanZ[2]{};
  double maxSpeed = 0.0;
  double maxCompression = 0.0;
  std::size_t count[2]{};
  bool finite = true;
  bool contained = true;
};

Summary summarise(const fluid::Particles& particles, const fluid::VesselBoundary& boundary,
                  double slack) {
  Summary summary;
  for (std::size_t i = 0; i < particles.size(); ++i) {
    const std::size_t phaseIndex = particles.phase[i] < 2 ? particles.phase[i] : 0;
    if (!std::isfinite(particles.px[i]) || !std::isfinite(particles.py[i]) ||
        !std::isfinite(particles.pz[i]) || !std::isfinite(particles.vx[i]) ||
        !std::isfinite(particles.vy[i]) || !std::isfinite(particles.vz[i])) {
      summary.finite = false;
      continue;
    }
    if (boundary.query(particles.px[i], particles.py[i], particles.pz[i]).distance > slack) {
      summary.contained = false;
    }
    summary.meanZ[phaseIndex] += particles.pz[i];
    ++summary.count[phaseIndex];
    summary.maxSpeed = std::max(
        summary.maxSpeed,
        std::sqrt(static_cast<double>(particles.vx[i]) * particles.vx[i] +
                  static_cast<double>(particles.vy[i]) * particles.vy[i] +
                  static_cast<double>(particles.vz[i]) * particles.vz[i]));
  }
  for (int index = 0; index < 2; ++index) {
    if (summary.count[index] > 0) {
      summary.meanZ[index] /= static_cast<double>(summary.count[index]);
    }
  }
  return summary;
}

}  // namespace

// One CPU-versus-GPU run of the same charge. Returned rather than asserted so
// the caller can compare across resolutions, which is the only way to say
// anything honest about a backend whose whole value is throughput at scale.
struct Comparison {
  std::size_t particles = 0;
  double cpuMs = 0.0;
  double gpuMs = 0.0;
  double speedup = 0.0;
};

Comparison compareAtSpacing(double spacing, int frames) {
  Charge charge = separatoryCharge(spacing);
  REQUIRE(charge.particles.size() > 200);
  const double slack = charge.config.contactRadiusFactor * charge.config.resolution.spacing;

  fluid::Solver reference;
  reference.configure(charge.config);
  reference.setPhases(charge.phases, {0.030});
  fluid::Particles cpuParticles = charge.particles;

  gfx::FluidGpuSolver gpu;
  gpu.configure(charge.config);
  gpu.setPhases(charge.phases, {0.030});
  INFO("gpu unavailable: ", gpu.unavailableReason());
  REQUIRE(gpu.initialise());
  REQUIRE(gpu.available());
  REQUIRE(gpu.upload(charge.particles, charge.boundary));
  REQUIRE(gpu.resident());

  const auto cpuStarted = std::chrono::steady_clock::now();
  double time = 0.0;
  for (int frame = 0; frame < frames; ++frame) {
    reference.advance(cpuParticles, charge.boundary, atRest(), time, kFrameS);
    time += kFrameS;
  }
  const double cpuMs =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - cpuStarted)
          .count();

  const auto gpuStarted = std::chrono::steady_clock::now();
  time = 0.0;
  for (int frame = 0; frame < frames; ++frame) {
    gpu.advance(charge.boundary, atRest(), time, kFrameS);
    time += kFrameS;
  }
  fluid::Particles gpuParticles;
  REQUIRE(gpu.download(gpuParticles));
  const double gpuMs =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - gpuStarted)
          .count();

  const Summary cpu = summarise(cpuParticles, charge.boundary, slack);
  const Summary device = summarise(gpuParticles, charge.boundary, slack);
  const fluid::Solver::Stats gpuStats = gpu.stats();
  std::cout << "[gpu fluid] dx=" << spacing * 1000.0 << "mm particles=" << charge.particles.size()
            << " frames=" << frames << " cpu=" << cpuMs << "ms gpu=" << gpuMs
            << "ms speedup=" << cpuMs / gpuMs << "x\n"
            << "[gpu fluid]   mean z cpu=(" << cpu.meanZ[0] << ", " << cpu.meanZ[1] << ") gpu=("
            << device.meanZ[0] << ", " << device.meanZ[1] << ")"
            << " max speed cpu=" << cpu.maxSpeed << " gpu=" << device.maxSpeed
            << " gpu compression=" << gpuStats.maxDensityCompression * 100.0 << "%\n";

  // Conservation and containment are absolute: float accumulation may move the
  // last bits, it may not lose or eject a particle.
  CHECK(gpuParticles.size() == charge.particles.size());
  CHECK(device.count[0] == cpu.count[0]);
  CHECK(device.count[1] == cpu.count[1]);
  CHECK(device.finite);
  CHECK(device.contained);
  CHECK(device.maxSpeed <= charge.config.maxSpeed);

  // The two paths are not bit-identical by design (see docs/gpu-acceleration.md):
  // GLSL reduces in float and the counting-sort scatter is atomic, so two runs
  // of the GPU path do not even agree with each other after a few hundred
  // substeps of a settling flow. What must agree is the physics -- the denser
  // phase ends up below the lighter one by the same margin.
  CHECK(device.meanZ[0] < device.meanZ[1]);
  const double separation = cpu.meanZ[1] - cpu.meanZ[0];
  REQUIRE(separation > 0.0);
  CHECK(device.meanZ[1] - device.meanZ[0] == doctest::Approx(separation).epsilon(0.30));

  return {charge.particles.size(), cpuMs, gpuMs, cpuMs / gpuMs};
}

TEST_CASE("the GPU fluid backend matches the CPU solver and scales past it") {
  GlFixture& gl = GlFixture::instance();
  if (!gl.ready()) {
    std::cout << "[gpu fluid] skipped: " << gl.reason() << '\n';
    return;
  }
  std::cout << "[gpu fluid] device: " << gfx::glRendererString() << " / "
            << gfx::glVersionString() << '\n';

  const gfx::ComputeCapabilities capabilities = gfx::probeCompute();
  if (!capabilities.available) {
    std::cout << "[gpu fluid] skipped: " << capabilities.reason << '\n';
    return;
  }
  CHECK(capabilities.majorVersion * 10 + capabilities.minorVersion >= 43);
  CHECK(capabilities.maxWorkGroupInvocations >= 128);

  // Two charges an order of magnitude apart in particle count. The small one is
  // where dispatch and convergence-readback overhead dominates and the CPU wins;
  // the large one is the regime the backend exists for. Frame counts are chosen
  // so each leg costs a few seconds, not so they are equal.
  const Comparison small = compareAtSpacing(6.0e-3, 20);
  const Comparison large = compareAtSpacing(3.0e-3, 4);
  CHECK(large.particles > 4 * small.particles);

  // The contract the backend has to earn: it must beat the reference solver on
  // the workload it was written for, and it must scale better than the CPU or
  // there is no reason to carry a second implementation.
  CHECK(large.speedup > 1.0);
  CHECK(large.speedup > small.speedup);
}

#include "fluid/simulation.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <exception>
#include <iomanip>
#include <limits>
#include <mutex>
#include <numeric>
#include <sstream>
#include <thread>
#include <utility>

namespace chemcad::fluid {
namespace {

constexpr double kMaximumAdvanceS = 0.1;
constexpr double kShakeTimerToleranceS = 1.0e-12;

std::array<double, 3> normalisedAxis(const std::array<double, 3>& axis) {
  const double lengthSquared =
      axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2];
  if (!(lengthSquared > 0.0) || !std::isfinite(lengthSquared)) {
    return {0.0, 0.0, 1.0};
  }
  const double inverseLength = 1.0 / std::sqrt(lengthSquared);
  return {axis[0] * inverseLength, axis[1] * inverseLength,
          axis[2] * inverseLength};
}

double vesselHeight(sol::Vessel vessel, double ratedVolumeMl) {
  sol::Simulation sizing;
  sizing.vessel = vessel;
  sizing.vesselVolumeMl = ratedVolumeMl;
  // columnHeightM integrates pi*r(z)^2 over the same analytic
  // sol::vesselWidthAt profile that VesselBoundary::build samples. Reusing it
  // prevents capacity calibration from drifting between the 2D and 3D models.
  return sol::columnHeightM(sizing);
}

bool hasPositiveSurfaceTension(const std::vector<double>& sigmaPairs) {
  return std::any_of(sigmaPairs.begin(), sigmaPairs.end(),
                     [](double sigma) { return sigma > 0.0; });
}

}  // namespace

struct Simulation::Impl {
  Particles particles;
  VesselBoundary boundary;
  Solver solver;
  DiagnosticsEngine diagnosticsEngine;
  NeighbourGrid diagnosticsGrid;
  VesselMotion motion;
  VesselMotion requestedMotion;
  Resolution resolution;
  std::vector<PhaseMaterial> phases;
  std::vector<double> sigmaPairs;
  std::string statusIssue;

  sol::Vessel vessel = sol::Vessel::SeparatoryFunnel;
  double ratedVolumeMl = 250.0;
  double vesselHeightM = 0.19;
  double elapsed = 0.0;
  Diagnostics latestDiagnostics;
  uint64_t revision = 0;
  uint64_t requestedShakeGeneration = 0;
  uint64_t appliedShakeGeneration = 0;

  mutable std::mutex stepMutex;
  mutable std::mutex controlMutex;
  mutable std::mutex publicationMutex;
  mutable std::mutex workMutex;
  std::condition_variable workAvailable;
  std::condition_variable idle;
  std::thread worker;
  bool stopping = false;
  bool active = false;
  bool invalidRequest = false;
  double queuedRealtimeSeconds = 0.0;
  double queuedExactSeconds = 0.0;
  double activeSeconds = 0.0;
  std::atomic<double> measuredRealTimeFactor{
      std::numeric_limits<double>::quiet_NaN()};
  std::shared_ptr<const Snapshot> published;
  std::shared_ptr<const Solver::Stats> publishedStats =
      std::make_shared<const Solver::Stats>();
  std::string publishedIssue;

  Impl() {
    SolverConfig config;
    config.resolution = resolution;
    solver.configure(config);
    vesselHeightM = vesselHeight(vessel, ratedVolumeMl);
    boundary.build(vessel, vesselHeightM, resolution.support(), resolution.spacing);
    publish();
    worker = std::thread([this] { workerLoop(); });
  }

  ~Impl() {
    {
      std::lock_guard<std::mutex> lock(workMutex);
      stopping = true;
      queuedRealtimeSeconds = 0.0;
      queuedExactSeconds = 0.0;
      invalidRequest = false;
    }
    workAvailable.notify_one();
    if (worker.joinable()) worker.join();
  }

  void rebuildBoundary() {
    boundary.build(vessel, vesselHeightM, resolution.support(), resolution.spacing);
  }

  bool calibrateIfNeeded() {
    if (phases.size() <= 1 || !hasPositiveSurfaceTension(sigmaPairs)) return true;
    try {
      // Akinci et al. (ACM TOG 32(6), 2013) require a resolution-dependent
      // cohesion coefficient. Solver calibrates that coefficient against the
      // Young-Laplace pressure jump instead of equating it to sigma in N/m.
      solver.calibrateInterface(boundary);
      return true;
    } catch (const std::exception& error) {
      statusIssue = std::string("Surface tension disabled: ") + error.what();
    } catch (...) {
      statusIssue = "Surface tension disabled: interface calibration failed";
    }
    return false;
  }

  void zeroParticleMotion() {
    std::fill(particles.vx.begin(), particles.vx.end(), 0.0f);
    std::fill(particles.vy.begin(), particles.vy.end(), 0.0f);
    std::fill(particles.vz.begin(), particles.vz.end(), 0.0f);
    std::fill(particles.ax.begin(), particles.ax.end(), 0.0f);
    std::fill(particles.ay.begin(), particles.ay.end(), 0.0f);
    std::fill(particles.az.begin(), particles.az.end(), 0.0f);
    std::fill(particles.pressure.begin(), particles.pressure.end(), 0.0f);
  }

  void recomputeDiagnostics(double intervalS) {
    diagnosticsGrid.build(particles, resolution.support());
    const Diagnostics raw = diagnosticsEngine.compute(
        particles, diagnosticsGrid, boundary, phases, resolution.spacing);
    latestDiagnostics = diagnosticsEngine.smooth(raw, intervalS);
  }

  void resetMotionForCharge() {
    std::lock_guard<std::mutex> controlLock(controlMutex);
    const Pose retainedPose = requestedMotion.pose;
    requestedMotion = {};
    requestedMotion.pose = retainedPose;
    ++requestedShakeGeneration;
    appliedShakeGeneration = requestedShakeGeneration;
    motion = {};
    motion.pose = retainedPose;
  }

  bool charge(bool preserveIssue = false) {
    if (!preserveIssue &&
        statusIssue.rfind("Surface tension disabled:", 0) != 0) {
      statusIssue.clear();
    }
    particles.clear();
    latestDiagnostics = {};
    diagnosticsEngine.reset();
    if (phases.empty()) {
      statusIssue = "Not charged: configure at least one fluid phase";
      publish();
      return false;
    }

    try {
      rebuildBoundary();
      boundary.chargeLattice(phases, resolution.spacing, particles);
      if (particles.empty()) {
        statusIssue = "Not charged: configured volumes produced no particles";
        publish();
        return false;
      }
      zeroParticleMotion();
      resetMotionForCharge();
      elapsed = 0.0;
      recomputeDiagnostics(0.0);
      publish();
      return true;
    } catch (const std::exception& error) {
      particles.clear();
      latestDiagnostics = {};
      statusIssue = std::string("Fluid charge failed: ") + error.what();
    } catch (...) {
      particles.clear();
      latestDiagnostics = {};
      statusIssue = "Fluid charge failed";
    }
    publish();
    return false;
  }

  void applyRequestedMotion(VesselMotion& motionForStep,
                            uint64_t& shakeGenerationForStep) {
    std::lock_guard<std::mutex> controlLock(controlMutex);
    motion.pose = requestedMotion.pose;
    motion.manualAcceleration = requestedMotion.manualAcceleration;
    motion.angularVelocity = requestedMotion.angularVelocity;
    motion.angularAcceleration = requestedMotion.angularAcceleration;
    if (appliedShakeGeneration != requestedShakeGeneration) {
      motion.shaking = requestedMotion.shaking;
      motion.shakeAxis = requestedMotion.shakeAxis;
      motion.shakeRemainingS = requestedMotion.shakeRemainingS;
      motion.shakeFrequencyHz = requestedMotion.shakeFrequencyHz;
      motion.shakeAmplitudeM = requestedMotion.shakeAmplitudeM;
      appliedShakeGeneration = requestedShakeGeneration;
    }
    shakeGenerationForStep = appliedShakeGeneration;
    motionForStep = motion;
  }

  void publishShakeState(uint64_t shakeGenerationForStep) {
    std::lock_guard<std::mutex> controlLock(controlMutex);
    if (requestedShakeGeneration != shakeGenerationForStep) return;
    requestedMotion.shaking = motion.shaking;
    requestedMotion.shakeRemainingS = motion.shakeRemainingS;
  }

  bool integrate(double stepS) {
    std::lock_guard<std::mutex> lock(stepMutex);
    if (particles.empty()) {
      statusIssue = phases.empty()
                        ? "Advance ignored: configure fluid phases and charge"
                        : "Advance ignored: no fluid is charged";
      publish();
      return false;
    }

    const auto started = std::chrono::steady_clock::now();
    try {
      VesselMotion motionForStep;
      uint64_t shakeGenerationForStep = 0;
      applyRequestedMotion(motionForStep, shakeGenerationForStep);
      // Keep the shake armed for the interval that begins while its timer is
      // positive. The analytic p''(t) in frame.cpp is sampled by every solver
      // substep; no pointer finite difference enters the forcing term.
      solver.advance(particles, boundary, motionForStep, elapsed, stepS);
      elapsed += stepS;
      if (motion.shaking) {
        motion.shakeRemainingS =
            std::max(0.0, motion.shakeRemainingS - stepS);
        if (motion.shakeRemainingS <= kShakeTimerToleranceS) {
          motion.shakeRemainingS = 0.0;
          motion.shaking = false;
        }
      }
      publishShakeState(shakeGenerationForStep);

      // Diagnostics are part of the immutable completed-step publication. They
      // are recomputed here so the UI never labels a new solver state with an
      // older interface estimate.
      recomputeDiagnostics(stepS);
      const double wallSeconds =
          std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
              .count();
      const double factor = wallSeconds > 0.0 ? stepS / wallSeconds : 1.0;
      measuredRealTimeFactor.store(factor, std::memory_order_release);
      publish();
      return true;
    } catch (const std::exception& error) {
      particles.clear();
      latestDiagnostics = {};
      statusIssue = std::string("Fluid advance stopped: ") + error.what();
    } catch (...) {
      particles.clear();
      latestDiagnostics = {};
      statusIssue = "Fluid advance stopped";
    }
    publish();
    return false;
  }

  void publishInvalidRequest() {
    std::lock_guard<std::mutex> lock(stepMutex);
    statusIssue = "Advance ignored: timestep must be finite and positive";
    publish();
  }

  void enqueueAdvance(double simulatedSeconds, bool realtimeRequest) {
    {
      std::lock_guard<std::mutex> lock(workMutex);
      if (!(simulatedSeconds > 0.0) || !std::isfinite(simulatedSeconds)) {
        invalidRequest = true;
      } else if (!stopping) {
        // Count the in-flight interval against the cap as well as queued
        // demand. A slow solve therefore drops excess wall-clock demand rather
        // than finishing one expensive step only to find another 0.1 s batch.
        const double queued =
            queuedRealtimeSeconds + queuedExactSeconds;
        const double room =
            std::max(0.0, kMaximumAdvanceS - activeSeconds - queued);
        double& destination =
            realtimeRequest ? queuedRealtimeSeconds : queuedExactSeconds;
        destination += std::min(simulatedSeconds, room);
      }
    }
    workAvailable.notify_one();
  }

  void workerLoop() {
    for (;;) {
      double stepS = 0.0;
      bool rejectRequest = false;
      {
        std::unique_lock<std::mutex> lock(workMutex);
        workAvailable.wait(lock, [this] {
          return stopping || invalidRequest || queuedRealtimeSeconds > 0.0 ||
                 queuedExactSeconds > 0.0;
        });
        if (stopping) break;
        rejectRequest = invalidRequest;
        invalidRequest = false;
        // Integrate every accepted queued second. The measured factor is an
        // honest diagnostic and a caller-side budget input, never a multiplier
        // that silently shortens requested simulated time.
        stepS = queuedExactSeconds + queuedRealtimeSeconds;
        queuedRealtimeSeconds = 0.0;
        queuedExactSeconds = 0.0;
        active = true;
        activeSeconds = stepS;
      }

      if (rejectRequest) publishInvalidRequest();
      if (stepS > 0.0) integrate(stepS);

      {
        std::lock_guard<std::mutex> lock(workMutex);
        active = false;
        activeSeconds = 0.0;
        if (stopping) {
          queuedRealtimeSeconds = 0.0;
          queuedExactSeconds = 0.0;
        }
      }
      idle.notify_all();
    }

    {
      std::lock_guard<std::mutex> lock(workMutex);
      active = false;
      activeSeconds = 0.0;
    }
    idle.notify_all();
  }

  void publish() {
    auto next = std::make_shared<Snapshot>();
    const std::size_t count = particles.size();
    next->px = particles.px;
    next->py = particles.py;
    next->pz = particles.pz;
    next->phase = particles.phase;
    next->speed.resize(count, 0.0f);
    next->colour.resize(count, 0.0f);
    for (std::size_t i = 0; i < count; ++i) {
      const double vx = particles.vx[i];
      const double vy = particles.vy[i];
      const double vz = particles.vz[i];
      next->speed[i] =
          static_cast<float>(std::sqrt(vx * vx + vy * vy + vz * vz));
      if (i < particles.colour.size()) {
        next->colour[i] = particles.colour[i];
      } else if (i < particles.phase.size()) {
        next->colour[i] = particles.phase[i] == 0 ? 0.0f : 1.0f;
      }
    }
    next->phases = phases;
    next->diagnostics = latestDiagnostics;
    next->pose = motion.pose;
    next->vesselHeightM = boundary.heightM();
    next->maxRadiusM = boundary.maxRadiusM();
    next->particleRadiusM = 0.5 * resolution.spacing;
    next->elapsedS = elapsed;
    next->vessel = vessel;
    next->revision = ++revision;

    auto nextStats = std::make_shared<const Solver::Stats>(solver.stats());
    std::lock_guard<std::mutex> lock(publicationMutex);
    published = std::move(next);
    publishedStats = std::move(nextStats);
    publishedIssue = statusIssue;
  }
};

Simulation::Simulation() : impl_(std::make_unique<Impl>()) {}

Simulation::~Simulation() = default;

void Simulation::setVessel(sol::Vessel vessel, double ratedVolumeMl) {
  std::lock_guard<std::mutex> lock(impl_->stepMutex);
  try {
    impl_->vessel = vessel;
    impl_->ratedVolumeMl =
        std::isfinite(ratedVolumeMl) ? std::max(ratedVolumeMl, 1.0) : 250.0;
    impl_->vesselHeightM = vesselHeight(impl_->vessel, impl_->ratedVolumeMl);
    impl_->rebuildBoundary();
    impl_->statusIssue.clear();
    impl_->calibrateIfNeeded();
    if (!impl_->phases.empty()) {
      impl_->charge(true);
    } else {
      impl_->publish();
    }
  } catch (const std::exception& error) {
    impl_->statusIssue = std::string("Vessel setup failed: ") + error.what();
    impl_->particles.clear();
    impl_->publish();
  } catch (...) {
    impl_->statusIssue = "Vessel setup failed";
    impl_->particles.clear();
    impl_->publish();
  }
}

void Simulation::setResolution(double spacingM) {
  std::lock_guard<std::mutex> lock(impl_->stepMutex);
  if (!(spacingM > 0.0) || !std::isfinite(spacingM)) {
    impl_->statusIssue = "Resolution unchanged: spacing must be finite and positive";
    impl_->publish();
    return;
  }
  try {
    impl_->resolution.spacing = spacingM;
    SolverConfig config = impl_->solver.config();
    config.resolution = impl_->resolution;
    impl_->solver.configure(config);
    if (!impl_->phases.empty()) {
      impl_->solver.setPhases(impl_->phases, impl_->sigmaPairs);
    }
    impl_->rebuildBoundary();
    impl_->statusIssue.clear();
    impl_->calibrateIfNeeded();
    if (!impl_->phases.empty()) {
      impl_->charge(true);
    } else {
      impl_->publish();
    }
  } catch (const std::exception& error) {
    impl_->statusIssue = std::string("Resolution setup failed: ") + error.what();
    impl_->particles.clear();
    impl_->publish();
  } catch (...) {
    impl_->statusIssue = "Resolution setup failed";
    impl_->particles.clear();
    impl_->publish();
  }
}

void Simulation::setPhases(const std::vector<PhaseMaterial>& phases,
                           const std::vector<double>& sigmaPairs) {
  std::lock_guard<std::mutex> lock(impl_->stepMutex);
  if (phases.empty() || phases.size() > 255) {
    impl_->phases.clear();
    impl_->sigmaPairs.clear();
    impl_->particles.clear();
    impl_->latestDiagnostics = {};
    impl_->statusIssue = "Not charged: configure between one and 255 fluid phases";
    impl_->publish();
    return;
  }
  const bool invalidMaterial = std::any_of(
      phases.begin(), phases.end(), [](const PhaseMaterial& phase) {
        return !(phase.restDensity > 0.0) || !std::isfinite(phase.restDensity) ||
               phase.dynamicViscosity < 0.0 ||
               !std::isfinite(phase.dynamicViscosity) ||
               !std::isfinite(phase.volumeMl);
      });
  const std::size_t n = phases.size();
  const bool invalidSigmaSize =
      !sigmaPairs.empty() && sigmaPairs.size() != n * n &&
      sigmaPairs.size() != n * (n - 1) / 2;
  if (invalidMaterial || invalidSigmaSize) {
    impl_->statusIssue = invalidMaterial
                             ? "Phase setup unchanged: invalid material properties"
                             : "Phase setup unchanged: invalid surface-tension table";
    impl_->publish();
    return;
  }
  try {
    impl_->solver.setPhases(phases, sigmaPairs);
    impl_->phases = phases;
    impl_->sigmaPairs = sigmaPairs;
    impl_->rebuildBoundary();
    impl_->statusIssue.clear();
    impl_->calibrateIfNeeded();
    impl_->charge(true);
  } catch (const std::exception& error) {
    impl_->statusIssue = std::string("Phase setup failed: ") + error.what();
    impl_->particles.clear();
    impl_->publish();
  } catch (...) {
    impl_->statusIssue = "Phase setup failed";
    impl_->particles.clear();
    impl_->publish();
  }
}

void Simulation::charge() {
  std::lock_guard<std::mutex> lock(impl_->stepMutex);
  // Setup order is intentionally permissive: a panel may request charge
  // before its material editor has produced a phase table. That request is a
  // reported no-op, not an exception crossing the draw call.
  impl_->charge();
}

bool Simulation::charged() const {
  std::lock_guard<std::mutex> lock(impl_->stepMutex);
  return !impl_->particles.empty();
}

void Simulation::shake(const std::array<double, 3>& axis, double durationS,
                       double frequencyHz, double amplitudeM) {
  std::lock_guard<std::mutex> lock(impl_->controlMutex);
  impl_->requestedMotion.shakeAxis = normalisedAxis(axis);
  impl_->requestedMotion.shakeRemainingS = std::max(0.0, durationS);
  impl_->requestedMotion.shakeFrequencyHz = std::max(0.0, frequencyHz);
  impl_->requestedMotion.shakeAmplitudeM = std::abs(amplitudeM);
  impl_->requestedMotion.shaking =
      impl_->requestedMotion.shakeRemainingS > 0.0 &&
      impl_->requestedMotion.shakeFrequencyHz > 0.0 &&
      impl_->requestedMotion.shakeAmplitudeM > 0.0;
  ++impl_->requestedShakeGeneration;
}

void Simulation::setManualAcceleration(
    const std::array<double, 3>& acceleration) {
  std::lock_guard<std::mutex> lock(impl_->controlMutex);
  impl_->requestedMotion.manualAcceleration = acceleration;
}

void Simulation::setPose(const Pose& pose,
                         const std::array<double, 3>& angularVelocity,
                         const std::array<double, 3>& angularAcceleration) {
  std::lock_guard<std::mutex> lock(impl_->controlMutex);
  impl_->requestedMotion.pose = pose;
  impl_->requestedMotion.angularVelocity = angularVelocity;
  impl_->requestedMotion.angularAcceleration = angularAcceleration;
}

void Simulation::requestAdvance(double simulatedSeconds) {
  impl_->enqueueAdvance(simulatedSeconds, true);
}

bool Simulation::stepping() const {
  std::lock_guard<std::mutex> lock(impl_->workMutex);
  return impl_->active;
}

double Simulation::pendingSeconds() const {
  std::lock_guard<std::mutex> lock(impl_->workMutex);
  return impl_->queuedRealtimeSeconds + impl_->queuedExactSeconds;
}

void Simulation::waitForIdle() {
  std::unique_lock<std::mutex> lock(impl_->workMutex);
  impl_->idle.wait(lock, [this] {
    return impl_->stopping ||
           (!impl_->active && !impl_->invalidRequest &&
            impl_->queuedRealtimeSeconds <= 0.0 &&
            impl_->queuedExactSeconds <= 0.0);
  });
}

void Simulation::advance(double dt) {
  // Batch callers retain the historical exact-dt contract; only live
  // non-blocking demand is reduced by the measured real-time budget.
  impl_->enqueueAdvance(dt, false);
  waitForIdle();
}

std::shared_ptr<const Snapshot> Simulation::snapshot() const {
  std::lock_guard<std::mutex> lock(impl_->publicationMutex);
  return impl_->published;
}

const Diagnostics& Simulation::diagnostics() const {
  // The thread-local owner keeps the referenced immutable snapshot alive even
  // if a worker publishes its successor immediately after this function.
  thread_local std::shared_ptr<const Snapshot> held;
  held = snapshot();
  return held->diagnostics;
}

const Solver::Stats& Simulation::solverStats() const {
  thread_local std::shared_ptr<const Solver::Stats> held;
  {
    std::lock_guard<std::mutex> lock(impl_->publicationMutex);
    held = impl_->publishedStats;
  }
  return *held;
}

double Simulation::elapsedS() const {
  const auto state = snapshot();
  return state ? state->elapsedS : 0.0;
}

bool Simulation::shaking() const {
  std::lock_guard<std::mutex> lock(impl_->controlMutex);
  return impl_->requestedMotion.shaking;
}
double Simulation::realTimeFactor() const {
  return impl_->measuredRealTimeFactor.load(std::memory_order_acquire);
}

double Simulation::totalVolumeMl() const {
  std::lock_guard<std::mutex> lock(impl_->stepMutex);
  return std::accumulate(impl_->phases.begin(), impl_->phases.end(), 0.0,
                         [](double total, const PhaseMaterial& phase) {
                           return total + phase.volumeMl;
                         });
}

std::string Simulation::statusLine() const {
  std::shared_ptr<const Snapshot> state;
  std::shared_ptr<const Solver::Stats> stats;
  std::string issue;
  {
    std::lock_guard<std::mutex> publicationLock(impl_->publicationMutex);
    state = impl_->published;
    stats = impl_->publishedStats;
    issue = impl_->publishedIssue;
  }
  if (!state || !stats) return "Fluid simulation unavailable";

  std::ostringstream text;
  if (!issue.empty()) text << issue << " | ";
  text << std::fixed << std::setprecision(1)
       << "dx " << state->particleRadiusM * 2000.0 << " mm | "
       << state->px.size() << " particles | ";
  if (state->elapsedS > 0.0) {
    text << stats->substeps << " substeps, " << stats->pressureIterations
         << " pressure iterations | worst compression "
         << stats->maxDensityCompression * 100.0 << "%, deficit "
         << stats->maxDensityDeficit * 100.0 << "% | "
         << std::setprecision(2) << stats->millisecondsPerSubstep << " ms/substep | ";
  } else {
    text << "solver step pending | ";
  }
  if (!state->diagnostics.valid) {
    text << "diagnostics pending | ";
  } else if (state->elapsedS <= 0.0) {
    text << "interface diagnostics waiting for a completed step | ";
  } else {
    text << "area " << state->diagnostics.interfacialAreaM2 * 1.0e4
         << " cm^2, dispersed " << state->diagnostics.dispersedFraction * 100.0
         << "%, d32 ";
    if (state->diagnostics.sauterDiameterM > 0.0) {
      text << state->diagnostics.sauterDiameterM * 1000.0 << " mm";
    } else {
      text << "no resolved drops";
    }
    text << " | ";
  }

  const double factor = realTimeFactor();
  if (std::isfinite(factor)) {
    text << std::setprecision(2) << "physics at " << factor << "x real time";
  } else {
    text << "physics rate pending";
  }
  return text.str();
}

}  // namespace chemcad::fluid

#include "fluid/simulation.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iomanip>
#include <mutex>
#include <numeric>
#include <sstream>
#include <utility>

namespace chemcad::fluid {
namespace {

constexpr double kDiagnosticsPeriodS = 1.0 / 30.0;
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
  Resolution resolution;
  std::vector<PhaseMaterial> phases;
  std::vector<double> sigmaPairs;
  std::string statusIssue;

  sol::Vessel vessel = sol::Vessel::SeparatoryFunnel;
  double ratedVolumeMl = 250.0;
  double vesselHeightM = 0.19;
  double elapsed = 0.0;
  double diagnosticsClock = 0.0;
  Diagnostics latestDiagnostics;
  uint64_t revision = 0;

  mutable std::mutex stepMutex;
  mutable std::mutex publicationMutex;
  std::shared_ptr<const Snapshot> published;
  std::shared_ptr<const Solver::Stats> publishedStats =
      std::make_shared<const Solver::Stats>();

  Impl() {
    SolverConfig config;
    config.resolution = resolution;
    solver.configure(config);
    vesselHeightM = vesselHeight(vessel, ratedVolumeMl);
    boundary.build(vessel, vesselHeightM, resolution.support(), resolution.spacing);
    publish();
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

      const Pose retainedPose = motion.pose;
      motion = {};
      motion.pose = retainedPose;
      elapsed = 0.0;
      diagnosticsClock = 0.0;
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
      next->speed[i] = static_cast<float>(std::sqrt(vx * vx + vy * vy + vz * vz));
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
  std::lock_guard<std::mutex> lock(impl_->stepMutex);
  impl_->motion.shakeAxis = normalisedAxis(axis);
  impl_->motion.shakeRemainingS = std::max(0.0, durationS);
  impl_->motion.shakeFrequencyHz = std::max(0.0, frequencyHz);
  impl_->motion.shakeAmplitudeM = std::abs(amplitudeM);
  impl_->motion.shaking = impl_->motion.shakeRemainingS > 0.0 &&
                          impl_->motion.shakeFrequencyHz > 0.0 &&
                          impl_->motion.shakeAmplitudeM > 0.0;
}

void Simulation::setManualAcceleration(const std::array<double, 3>& acceleration) {
  std::lock_guard<std::mutex> lock(impl_->stepMutex);
  impl_->motion.manualAcceleration = acceleration;
}

void Simulation::setPose(const Pose& pose,
                         const std::array<double, 3>& angularVelocity,
                         const std::array<double, 3>& angularAcceleration) {
  std::lock_guard<std::mutex> lock(impl_->stepMutex);
  impl_->motion.pose = pose;
  impl_->motion.angularVelocity = angularVelocity;
  impl_->motion.angularAcceleration = angularAcceleration;
}

void Simulation::advance(double dt) {
  std::lock_guard<std::mutex> lock(impl_->stepMutex);
  if (!(dt > 0.0) || !std::isfinite(dt)) {
    impl_->statusIssue = "Advance ignored: timestep must be finite and positive";
    impl_->publish();
    return;
  }
  if (impl_->particles.empty()) {
    impl_->statusIssue = impl_->phases.empty()
                             ? "Advance ignored: configure fluid phases and charge"
                             : "Advance ignored: no fluid is charged";
    impl_->publish();
    return;
  }
  const double stepS = std::min(dt, kMaximumAdvanceS);

  try {
    // Keep the shake armed for the interval that begins while its timer is
    // positive. The analytic p''(t) in frame.cpp is sampled by every solver
    // substep; no pointer finite difference enters the forcing term.
    const VesselMotion motionForStep = impl_->motion;
    impl_->solver.advance(impl_->particles, impl_->boundary, motionForStep,
                          impl_->elapsed, stepS);
    impl_->elapsed += stepS;
    if (impl_->motion.shaking) {
      impl_->motion.shakeRemainingS =
          std::max(0.0, impl_->motion.shakeRemainingS - stepS);
      if (impl_->motion.shakeRemainingS <= kShakeTimerToleranceS) {
        impl_->motion.shakeRemainingS = 0.0;
        impl_->motion.shaking = false;
      }
    }

    impl_->diagnosticsClock += stepS;
    if (!impl_->latestDiagnostics.valid ||
        impl_->diagnosticsClock >= kDiagnosticsPeriodS) {
      const double intervalS = impl_->diagnosticsClock;
      impl_->diagnosticsClock = 0.0;
      impl_->recomputeDiagnostics(intervalS);
    }
    impl_->publish();
  } catch (const std::exception& error) {
    impl_->particles.clear();
    impl_->latestDiagnostics = {};
    impl_->statusIssue = std::string("Fluid advance stopped: ") + error.what();
    impl_->publish();
  } catch (...) {
    impl_->particles.clear();
    impl_->latestDiagnostics = {};
    impl_->statusIssue = "Fluid advance stopped";
    impl_->publish();
  }
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
  std::lock_guard<std::mutex> lock(impl_->stepMutex);
  return impl_->motion.shaking;
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
  std::lock_guard<std::mutex> stepLock(impl_->stepMutex);
  issue = impl_->statusIssue;
  {
    std::lock_guard<std::mutex> publicationLock(impl_->publicationMutex);
    state = impl_->published;
    stats = impl_->publishedStats;
  }
  if (!state || !stats) return "Fluid simulation unavailable";

  std::ostringstream text;
  if (!issue.empty()) text << issue << " | ";
  text << std::fixed << std::setprecision(1)
       << "dx " << state->particleRadiusM * 2000.0 << " mm | "
       << state->px.size() << " particles | " << stats->substeps
       << " substeps | worst density error "
       << stats->maxDensityError * 100.0 << "% | ";
  if (!state->diagnostics.valid) {
    text << "diagnostics pending";
  } else {
    text << "area " << state->diagnostics.interfacialAreaM2 * 1.0e4
         << " cm^2, dispersed " << state->diagnostics.dispersedFraction * 100.0
         << "%, d32 " << state->diagnostics.sauterDiameterM * 1000.0 << " mm";
  }
  return text.str();
}

}  // namespace chemcad::fluid

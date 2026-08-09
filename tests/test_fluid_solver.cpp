#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <vector>

#include "fluid/frame.hpp"
#include "fluid/kernels.hpp"
#include "fluid/solver.hpp"
#include "fluid/simulation.hpp"
#include "fluid/vessel_sdf.hpp"
#include "sol/funnel.hpp"

namespace fluid = chemcad::fluid;
namespace sol = chemcad::sol;

namespace {

constexpr double kGravity = 9.80665;

fluid::PhaseMaterial phase(const char* label, double density, double viscosity, double volumeMl) {
  fluid::PhaseMaterial value;
  value.label = label;
  value.restDensity = density;
  value.dynamicViscosity = viscosity;
  value.volumeMl = volumeMl;
  return value;
}

fluid::SolverConfig coarseConfig() {
  fluid::SolverConfig config;
  config.resolution.spacing = 8.0e-3;
  config.maxPressureIterations = 10;
  config.maxSubstepS = 1.0 / 480.0;
  config.maxSpeed = 5.0;
  return config;
}

fluid::SolverConfig profileConfig(const fluid::QualityProfile& profile) {
  fluid::SolverConfig config;
  config.resolution.spacing = profile.spacing;
  config.densityTolerance = profile.densityTolerance;
  config.minPressureIterations = profile.minPressureIterations;
  config.maxPressureIterations = profile.maxPressureIterations;
  config.enableSurfaceTension = profile.surfaceTension;
  return config;
}

struct NamedQuality {
  const char* name;
  fluid::QualityProfile profile;
};

constexpr std::array<NamedQuality, 3> kQualityProfiles{{
    {"Interactive", fluid::QualityProfile::interactive()},
    {"Balanced", fluid::QualityProfile::balanced()},
    {"Quality", fluid::QualityProfile::quality()},
}};

fluid::VesselBoundary cylinder(const fluid::SolverConfig& config, double height = 0.12) {
  fluid::VesselBoundary boundary;
  boundary.build(sol::Vessel::GraduatedCylinder, height, config.resolution.support(),
                 config.resolution.spacing);
  return boundary;
}

fluid::VesselMotion noGravity() {
  fluid::VesselMotion motion;
  // frameAcceleration subtracts world translational acceleration. Matching
  // world gravity therefore creates a gravity-free inertial experiment.
  motion.manualAcceleration = {0.0, 0.0, -kGravity};
  return motion;
}

void runSteps(fluid::Solver& solver, fluid::Particles& particles,
              const fluid::VesselBoundary& boundary, const fluid::VesselMotion& motion,
              int count, double& time, double dt = 1.0 / 480.0) {
  for (int step = 0; step < count; ++step) {
    solver.advance(particles, boundary, motion, time, dt);
    time += dt;
  }
}

double kineticEnergy(const fluid::Particles& particles,
                     const std::vector<fluid::PhaseMaterial>& phases, double particleVolume) {
  double energy = 0.0;
  for (std::size_t i = 0; i < particles.size(); ++i) {
    const double mass = phases[particles.phase[i]].restDensity * particleVolume;
    const double speedSquared = static_cast<double>(particles.vx[i]) * particles.vx[i] +
                                static_cast<double>(particles.vy[i]) * particles.vy[i] +
                                static_cast<double>(particles.vz[i]) * particles.vz[i];
    energy += 0.5 * mass * speedSquared;
  }
  return energy;
}

std::size_t interfaceParticles(const fluid::Particles& particles) {
  return static_cast<std::size_t>(std::count_if(
      particles.colour.begin(), particles.colour.end(),
      [](float colour) { return colour > 0.2f && colour < 0.8f; }));
}

std::array<double, 3> momentum(const fluid::Particles& particles,
                               const std::vector<fluid::PhaseMaterial>& phases,
                               double particleVolume) {
  std::array<double, 3> result{};
  for (std::size_t i = 0; i < particles.size(); ++i) {
    const double mass = phases[particles.phase[i]].restDensity * particleVolume;
    result[0] += mass * particles.vx[i];
    result[1] += mass * particles.vy[i];
    result[2] += mass * particles.vz[i];
  }
  return result;
}

double norm(const std::array<double, 3>& value) {
  return std::sqrt(value[0] * value[0] + value[1] * value[1] + value[2] * value[2]);
}

void checkFiniteAndContained(const fluid::Particles& particles,
                             const fluid::VesselBoundary& boundary, double contactRadius) {
  for (std::size_t i = 0; i < particles.size(); ++i) {
    CAPTURE(i);
    CHECK(std::isfinite(particles.px[i]));
    CHECK(std::isfinite(particles.py[i]));
    CHECK(std::isfinite(particles.pz[i]));
    CHECK(std::isfinite(particles.vx[i]));
    CHECK(std::isfinite(particles.vy[i]));
    CHECK(std::isfinite(particles.vz[i]));
    CHECK(boundary.query(particles.px[i], particles.py[i], particles.pz[i]).distance <=
          contactRadius + 2.0e-6);
  }
}

void setFluidWorkers(unsigned workers) {
  const char* value = workers == 1 ? "1" : workers == 2 ? "2" : workers == 6 ? "6" : "0";
#ifdef _WIN32
  REQUIRE(_putenv_s("CHEMCAD_FLUID_WORKERS", value) == 0);
#else
  REQUIRE(setenv("CHEMCAD_FLUID_WORKERS", value, 1) == 0);
#endif
}

fluid::VesselBoundary separatoryFunnel(const fluid::SolverConfig& config) {
  sol::Simulation sizing;
  sizing.vessel = sol::Vessel::SeparatoryFunnel;
  sizing.vesselVolumeMl = 250.0;
  fluid::VesselBoundary boundary;
  boundary.build(sizing.vessel, sol::columnHeightM(sizing), config.resolution.support(),
                 config.resolution.spacing);
  return boundary;
}

}  // namespace

TEST_CASE("a coarse PCISPH column recovers hydrostatic density and pressure") {
  fluid::SolverConfig config = coarseConfig();
  config.maxPressureIterations = 12;
  config.enableSurfaceTension = false;
  fluid::Solver solver;
  solver.configure(config);
  // The 160 mm fixture has enough lattice sites to exercise wall-corrected
  // hydrostatics while retaining the production 8 mm coarse resolution.
  const std::vector<fluid::PhaseMaterial> phases{phase("water", 998.0, 0.89e-3, 70.0)};
  solver.setPhases(phases, {});
  fluid::VesselBoundary boundary = cylinder(config, 0.16);
  fluid::Particles particles;
  REQUIRE(boundary.chargeLattice(phases, config.resolution.spacing, particles) > 30);

  double time = 0.0;
  runSteps(solver, particles, boundary, {}, 480, time);

  double freeSurface = 0.0;
  for (float z : particles.pz) freeSurface = std::max(freeSurface, static_cast<double>(z));
  double maxDensityCompression = 0.0;
  std::vector<double> depths;
  std::vector<double> pressures;
  for (std::size_t i = 0; i < particles.size(); ++i) {
    if (particles.pz[i] < freeSurface - 1.5 * config.resolution.support()) {
      maxDensityCompression =
          std::max(maxDensityCompression,
                   std::max(0.0, static_cast<double>(particles.delta[i]) /
                                         solver.restNumberDensity() -
                                     1.0));
      depths.push_back(particles.pz[i]);
      pressures.push_back(particles.pressure[i]);
    }
  }
  std::cout << "[hydro diagnostic] particles=" << particles.size()
            << " interior=" << depths.size()
            << " surface-z=" << freeSurface
            << " compression=" << solver.stats().maxDensityCompression
            << " deficit=" << solver.stats().maxDensityDeficit
            << " iterations=" << solver.stats().pressureIterations << '\n';
  REQUIRE(depths.size() >= 6);
  CHECK(maxDensityCompression <= config.densityTolerance);
  CHECK(solver.stats().maxDensityCompression <= config.densityTolerance);
  CHECK(solver.stats().maxDensityDeficit > 0.1);

  double meanZ = 0.0;
  double meanP = 0.0;
  for (std::size_t i = 0; i < depths.size(); ++i) {
    meanZ += depths[i];
    meanP += pressures[i];
  }
  meanZ /= static_cast<double>(depths.size());
  meanP /= static_cast<double>(depths.size());
  double covariance = 0.0;
  double variance = 0.0;
  for (std::size_t i = 0; i < depths.size(); ++i) {
    covariance += (depths[i] - meanZ) * (pressures[i] - meanP);
    variance += (depths[i] - meanZ) * (depths[i] - meanZ);
  }
  REQUIRE(variance > 0.0);
  const double pressureGradient = covariance / variance;
  // The wall correction is locally planar at this coarse resolution. The
  // controlled quantity nevertheless meets the configured compression limit,
  // while pressure still increases monotonically with depth.
  CHECK(pressureGradient < 0.0);
  CHECK(std::abs(pressureGradient) <= 3.0 * phases[0].restDensity * kGravity);
}

TEST_CASE("corrected pressure and harmonic viscosity conserve momentum") {
  fluid::SolverConfig config = coarseConfig();
  config.enableSurfaceTension = false;
  fluid::Solver solver;
  solver.configure(config);
  const std::vector<fluid::PhaseMaterial> phases{phase("liquid", 1000.0, 4.0e-3, 0.0)};
  solver.setPhases(phases, {});
  fluid::VesselBoundary boundary = cylinder(config, 10.0);
  fluid::Particles particles;
  for (int z = -2; z <= 2; ++z) {
    for (int y = -2; y <= 2; ++y) {
      for (int x = -2; x <= 2; ++x) {
        const std::size_t i = particles.add(static_cast<float>(x * config.resolution.spacing),
                                            static_cast<float>(y * config.resolution.spacing),
                                            static_cast<float>(5.0 + z * config.resolution.spacing), 0);
        particles.vx[i] = static_cast<float>(0.75 + 0.02 * x);
        particles.vy[i] = static_cast<float>(-0.01 * y);
        particles.vz[i] = static_cast<float>(0.01 * z);
      }
    }
  }
  const std::array<double, 3> before =
      momentum(particles, phases, config.resolution.particleVolume());
  double time = 0.0;
  runSteps(solver, particles, boundary, noGravity(), 100, time);
  const std::array<double, 3> after =
      momentum(particles, phases, config.resolution.particleVolume());
  const std::array<double, 3> change{after[0] - before[0], after[1] - before[1],
                                     after[2] - before[2]};
  const double momentumScale = std::max(norm(before), 1.0e-12);
  // Pair-ordered forces leave a 1.2207e-10 kg*m/s residual after 100 Release
  // steps (2.54e-9 relative), from storing particle velocities as float.
  CHECK(norm(change) <= 5.0e-9 * momentumScale);
}

TEST_CASE("interface calibration is a safe no-op without a physical interface") {
  const fluid::SolverConfig config = coarseConfig();
  fluid::Solver solver;
  solver.configure(config);
  fluid::VesselBoundary boundary = cylinder(config);
  solver.setPhases({phase("single", 1000.0, 1.0e-3, 10.0)}, {});
  CHECK_NOTHROW(solver.calibrateInterface(boundary));
  CHECK(solver.interfaceModel().calibrated);
  CHECK(solver.interfaceModel().cohesionGain == 0.0);
  CHECK(solver.interfaceCalibrationError().empty());
  CHECK(solver.stats().interfaceCalibrations == 1);

  solver.setPhases(
      {phase("a", 1000.0, 1.0e-3, 5.0), phase("b", 900.0, 1.0e-3, 5.0)}, {0.0});
  CHECK_NOTHROW(solver.calibrateInterface(boundary));
  CHECK(solver.interfaceModel().calibrated);
  CHECK(solver.interfaceModel().cohesionGain == 0.0);
  CHECK(solver.interfaceCalibrationError().empty());
  CHECK(solver.stats().interfaceCalibrations == 2);
  CHECK_NOTHROW(solver.calibrateInterface(boundary));
  CHECK(solver.stats().interfaceCalibrations == 2);
}

TEST_CASE("pressure stiffness cache absorbs frame-clock wobble") {
  fluid::SolverConfig config = coarseConfig();
  config.enableSurfaceTension = false;
  fluid::Solver solver;
  solver.configure(config);
  const std::vector<fluid::PhaseMaterial> phases{
      phase("liquid", 1000.0, 1.0e-3, 20.0)};
  solver.setPhases(phases, {});
  fluid::VesselBoundary boundary = cylinder(config);
  fluid::Particles particles;
  REQUIRE(boundary.chargeLattice(phases, config.resolution.spacing, particles) > 0);

  const double nominal = 0.9 * std::min(
      config.maxSubstepS, config.cflNumber * config.resolution.support());
  solver.advance(particles, boundary, noGravity(), 0.0, nominal);
  const std::uint64_t calibrated = solver.stats().pressureStiffnessCalibrations;
  for (int frame = 1; frame <= 32; ++frame) {
    const double wobble = nominal + 1.0e-6 * static_cast<double>(frame % 3 - 1);
    solver.advance(particles, boundary, noGravity(), frame * nominal, wobble);
  }
  CHECK(solver.stats().pressureStiffnessCalibrations == calibrated);
}

TEST_CASE("shaking conserves particles and the analytic wall contains them") {
  fluid::SolverConfig config = coarseConfig();
  fluid::Solver solver;
  solver.configure(config);
  const std::vector<fluid::PhaseMaterial> phases{
      phase("light", 998.0, 1.0e-3, 35.0), phase("dense", 1326.0, 1.5e-3, 35.0)};
  solver.setPhases(phases, {0.030});
  fluid::VesselBoundary boundary = cylinder(config);
  fluid::Particles particles;
  boundary.chargeLattice(phases, config.resolution.spacing, particles);
  const std::size_t count = particles.size();
  const std::size_t phaseZero =
      static_cast<std::size_t>(std::count(particles.phase.begin(), particles.phase.end(), 0));
  const std::size_t phaseOne = count - phaseZero;

  fluid::VesselMotion shake;
  shake.shaking = true;
  shake.shakeRemainingS = 10.0;
  shake.shakeFrequencyHz = 3.0;
  shake.shakeAmplitudeM = 0.05;
  shake.shakeAxis = {0.0, 0.0, 1.0};
  double time = 0.0;
  runSteps(solver, particles, boundary, shake, 600, time);

  CHECK(particles.size() == count);
  CHECK(static_cast<std::size_t>(std::count(particles.phase.begin(), particles.phase.end(), 0)) ==
        phaseZero);
  CHECK(static_cast<std::size_t>(std::count(particles.phase.begin(), particles.phase.end(), 1)) ==
        phaseOne);
  checkFiniteAndContained(particles, boundary,
                          config.contactRadiusFactor * config.resolution.spacing);
}

TEST_CASE("off-axis shaking injects energy and interface, then both decay") {
  fluid::SolverConfig config = coarseConfig();
  config.wallFriction = 0.15;
  fluid::Solver solver;
  solver.configure(config);
  const std::vector<fluid::PhaseMaterial> phases{
      phase("aqueous", 998.0, 2.0e-3, 40.0), phase("organic", 850.0, 2.0e-3, 40.0)};
  solver.setPhases(phases, {0.025});
  fluid::VesselBoundary boundary = cylinder(config, 0.16);
  solver.calibrateInterface(boundary);
  INFO("calibration=", solver.interfaceCalibrationError(),
       " gain=", solver.interfaceModel().cohesionGain);
  REQUIRE(solver.interfaceModel().calibrated);
  const std::uint64_t interfaceCalibrations = solver.stats().interfaceCalibrations;
  CHECK(interfaceCalibrations == 1);
  solver.calibrateInterface(boundary);
  CHECK(solver.stats().interfaceCalibrations == interfaceCalibrations);
  fluid::Particles particles;
  boundary.chargeLattice(phases, config.resolution.spacing, particles);

  double time = 0.0;
  runSteps(solver, particles, boundary, {}, 240, time);
  CHECK(solver.stats().interfaceCalibrations == interfaceCalibrations);
  const double restingEnergy =
      kineticEnergy(particles, phases, config.resolution.particleVolume());
  const std::size_t restingInterface = interfaceParticles(particles);

  fluid::VesselMotion shake;
  shake.shaking = true;
  shake.shakeRemainingS = 2.0;
  shake.shakeFrequencyHz = 3.0;
  shake.shakeAmplitudeM = 0.05;
  // A perfectly axial lattice preserves exact rotational symmetry and cannot
  // develop the sloshing mode whose interfacial area this test measures.
  shake.shakeAxis = {1.0, 0.0, 0.25};
  double shakenEnergy = 0.0;
  std::size_t shakenInterface = 0;
  for (int step = 0; step < 480; ++step) {
    solver.advance(particles, boundary, shake, time, 1.0 / 480.0);
    time += 1.0 / 480.0;
    shakenEnergy =
        std::max(shakenEnergy,
                 kineticEnergy(particles, phases, config.resolution.particleVolume()));
    shakenInterface = std::max(shakenInterface, interfaceParticles(particles));
  }
  // A whole number of sinusoidal cycles can end at a kinetic-energy node; the
  // peak over the forcing interval is the physical injected-energy measure.
  CHECK(shakenEnergy > 2.0 * restingEnergy + 1.0e-8);
  CHECK(shakenInterface > restingInterface + particles.size() / 20);

  runSteps(solver, particles, boundary, {}, 720, time);
  const double settledEnergy =
      kineticEnergy(particles, phases, config.resolution.particleVolume());
  const std::size_t settledInterface = interfaceParticles(particles);
  CHECK(settledEnergy < 0.75 * shakenEnergy);
  CHECK(settledInterface < shakenInterface);
}

TEST_CASE("density contrast keeps the denser phase below a sharp interface") {
  fluid::SolverConfig config = coarseConfig();
  fluid::Solver solver;
  solver.configure(config);
  const std::vector<fluid::PhaseMaterial> phases{
      phase("water", 998.0, 1.0e-3, 45.0), phase("brine", 1326.0, 1.5e-3, 45.0)};
  solver.setPhases(phases, {0.040});
  fluid::VesselBoundary boundary = cylinder(config, 0.16);
  fluid::Particles particles;
  boundary.chargeLattice(phases, config.resolution.spacing, particles);
  double time = 0.0;
  runSteps(solver, particles, boundary, {}, 480, time);

  double meanZ[2]{};
  std::size_t count[2]{};
  for (std::size_t i = 0; i < particles.size(); ++i) {
    meanZ[particles.phase[i]] += particles.pz[i];
    ++count[particles.phase[i]];
  }
  REQUIRE(count[0] > 0);
  REQUIRE(count[1] > 0);
  meanZ[0] /= static_cast<double>(count[0]);
  meanZ[1] /= static_cast<double>(count[1]);
  CHECK(meanZ[1] < meanZ[0]);
  CHECK(interfaceParticles(particles) < particles.size() / 4);
}

TEST_CASE("identical fluid runs are bit deterministic") {
  fluid::SolverConfig config = coarseConfig();
  const std::vector<fluid::PhaseMaterial> phases{
      phase("light", 900.0, 1.0e-3, 30.0), phase("dense", 1100.0, 2.0e-3, 30.0)};
  fluid::VesselBoundary boundary = cylinder(config);
  fluid::Particles a;
  boundary.chargeLattice(phases, config.resolution.spacing, a);
  fluid::Particles b = a;
  fluid::Particles c = a;
  fluid::Solver first;
  fluid::Solver second;
  fluid::Solver third;
  first.configure(config);
  second.configure(config);
  third.configure(config);
  first.setPhases(phases, {0.025});
  second.setPhases(phases, {0.025});
  third.setPhases(phases, {0.025});
  fluid::VesselMotion shake;
  shake.shaking = true;
  shake.shakeRemainingS = 1.0;
  shake.shakeFrequencyHz = 3.0;
  shake.shakeAmplitudeM = 0.05;
  shake.shakeAxis = {0.25, -0.5, 1.0};
  double ta = 0.0;
  double tb = 0.0;
  setFluidWorkers(1);
  runSteps(first, a, boundary, shake, 240, ta);
  setFluidWorkers(2);
  runSteps(second, b, boundary, shake, 240, tb);
  double tc = 0.0;
  setFluidWorkers(6);
  runSteps(third, c, boundary, shake, 240, tc);
  setFluidWorkers(0);

  CHECK(a.id == b.id);
  CHECK(a.phase == b.phase);
  CHECK(a.px == b.px);
  CHECK(a.py == b.py);
  CHECK(a.pz == b.pz);
  CHECK(a.vx == b.vx);
  CHECK(a.vy == b.vy);
  CHECK(a.vz == b.vz);
  CHECK(a.id == c.id);
  CHECK(a.phase == c.phase);
  CHECK(a.px == c.px);
  CHECK(a.py == c.py);
  CHECK(a.pz == c.pz);
  CHECK(a.vx == c.vx);
  CHECK(a.vy == c.vy);
  CHECK(a.vz == c.vz);
}

TEST_CASE("each named quality profile meets its own compression limit") {
  const std::vector<fluid::PhaseMaterial> phases{
      phase("water", 998.0, 0.89e-3, 25.0),
      phase("dichloromethane", 1326.0, 0.43e-3, 25.0)};

  for (const NamedQuality& named : kQualityProfiles) {
    CAPTURE(named.name);
    fluid::SolverConfig config = profileConfig(named.profile);
    fluid::Solver solver;
    solver.configure(config);
    solver.setPhases(phases, {});
    fluid::VesselBoundary boundary = separatoryFunnel(config);
    fluid::Particles particles;
    REQUIRE(boundary.chargeLattice(phases, named.profile.spacing, particles) > 0);

    solver.advance(particles, boundary, noGravity(), 0.0, 30.0e-3);
    const fluid::Solver::Stats& stats = solver.stats();
    CHECK(stats.maxDensityCompression <= named.profile.densityTolerance);
    CHECK(stats.pressureStiffnessSubstepS ==
          doctest::Approx(stats.substepS).epsilon(1.0e-12));
    CHECK(interfaceParticles(particles) < particles.size() / 4);
    checkFiniteAndContained(
        particles, boundary,
        config.contactRadiusFactor * config.resolution.spacing);
  }
}

TEST_CASE("setQuality recharges consistently and remains deterministic") {
  const std::vector<fluid::PhaseMaterial> phases{
      phase("water", 998.0, 0.89e-3, 20.0),
      phase("dichloromethane", 1326.0, 0.43e-3, 20.0)};
  const fluid::QualityProfile quality = fluid::QualityProfile::interactive();

  fluid::Simulation changedAfterCharge;
  changedAfterCharge.setPhases(phases, {});
  const std::size_t fineCount = changedAfterCharge.snapshot()->px.size();
  changedAfterCharge.setQuality(quality);

  fluid::Simulation configuredFirst;
  configuredFirst.setQuality(quality);
  configuredFirst.setPhases(phases, {});

  const auto recharged = changedAfterCharge.snapshot();
  const auto fresh = configuredFirst.snapshot();
  REQUIRE(recharged->px.size() > 0);
  CHECK(recharged->px.size() < fineCount);
  CHECK(recharged->elapsedS == 0.0);
  CHECK(recharged->phase == fresh->phase);
  CHECK(recharged->px == fresh->px);
  CHECK(recharged->py == fresh->py);
  CHECK(recharged->pz == fresh->pz);

  changedAfterCharge.setManualAcceleration({0.0, 0.0, -kGravity});
  configuredFirst.setManualAcceleration({0.0, 0.0, -kGravity});
  changedAfterCharge.advance(30.0e-3);
  configuredFirst.advance(30.0e-3);
  const auto advancedA = changedAfterCharge.snapshot();
  const auto advancedB = configuredFirst.snapshot();
  CHECK(advancedA->phase == advancedB->phase);
  CHECK(advancedA->px == advancedB->px);
  CHECK(advancedA->py == advancedB->py);
  CHECK(advancedA->pz == advancedB->pz);
}

TEST_CASE("resolution-aware ceiling keeps an interactive frame stable") {
  const std::vector<fluid::PhaseMaterial> phases{
      phase("water", 998.0, 0.89e-3, 100.0),
      phase("dichloromethane", 1326.0, 0.43e-3, 100.0)};
  int coarseSubsteps = 0;
  int fineSubsteps = 0;

  for (double spacing : {8.0e-3, 4.0e-3}) {
    fluid::SolverConfig config;
    config.resolution.spacing = spacing;
    config.enableSurfaceTension = false;
    fluid::Solver solver;
    solver.configure(config);
    solver.setPhases(phases, {0.028});
    fluid::VesselBoundary boundary = separatoryFunnel(config);
    fluid::Particles particles;
    boundary.chargeLattice(phases, spacing, particles);
    if (spacing == 8.0e-3) REQUIRE(particles.size() > 300);
    if (spacing == 4.0e-3) REQUIRE(particles.size() == 3124);

    const int substeps =
        solver.advance(particles, boundary, noGravity(), 0.0, 30.0e-3);
    checkFiniteAndContained(
        particles, boundary, config.contactRadiusFactor * config.resolution.spacing);
    CHECK(solver.stats().maxDensityCompression <= config.densityTolerance);
    if (spacing == 8.0e-3) {
      coarseSubsteps = substeps;
      CHECK(substeps > 0);
      CHECK(substeps < 10);
    } else {
      fineSubsteps = substeps;
    }
  }
  std::cout << "[fluid substeps] dx=8mm " << coarseSubsteps
            << ", dx=4mm " << fineSubsteps << " for 30ms\n";

  CHECK(fineSubsteps > coarseSubsteps);
  CHECK(30.0e-3 / fineSubsteps < 30.0e-3 / coarseSubsteps);
}

TEST_CASE("quality-profile timing reports the loaded interactive budget") {
  const std::vector<fluid::PhaseMaterial> phases{
      phase("water", 998.0, 0.89e-3, 100.0),
      phase("dichloromethane", 1326.0, 0.43e-3, 100.0)};
  constexpr double kFrameS = 30.0e-3;

  // Keep one CPU busy while the solve runs, approximating the application
  // thread's input/layout/render submission work without launching the app.
  std::atomic<std::uint64_t> renderTicks{0};
  std::jthread renderLoad([&](std::stop_token stop) {
    while (!stop.stop_requested()) {
      renderTicks.fetch_add(1, std::memory_order_relaxed);
    }
  });
  while (renderTicks.load(std::memory_order_relaxed) < 10000) {
    std::this_thread::yield();
  }

  for (const NamedQuality& named : kQualityProfiles) {
    CAPTURE(named.name);
    fluid::SolverConfig legacy = profileConfig(named.profile);
    legacy.densityTolerance = 5.0e-3;
    legacy.minPressureIterations = 3;
    legacy.maxPressureIterations = 32;

    fluid::Solver solver;
    solver.configure(legacy);
    solver.setPhases(phases, {0.028});
    fluid::VesselBoundary boundary = separatoryFunnel(legacy);
    solver.calibrateInterface(boundary);
    INFO("calibration=", solver.interfaceCalibrationError());
    REQUIRE(solver.interfaceModel().calibrated);
    fluid::Particles charged;
    boundary.chargeLattice(phases, named.profile.spacing, charged);
    REQUIRE(charged.size() > 300);

    // Warm retained capacities, pair storage, and boundary samples. Both
    // budgets below then take the same steady-state path as the live worker.
    fluid::Particles warm = charged;
    REQUIRE(solver.advance(warm, boundary, noGravity(), 0.0, 1.0e-6) == 1);

    struct Timing {
      double realTimeFactor = 0.0;
      fluid::Solver::Stats stats;
    };
    const auto measure = [&](const char* budget,
                             const fluid::SolverConfig& config) {
      solver.configure(config);
      solver.setPhases(phases, {0.028});
      fluid::Particles particles = charged;
      const auto started = std::chrono::steady_clock::now();
      const int substeps =
          solver.advance(particles, boundary, {}, 0.0, kFrameS);
      const double elapsedMs =
          std::chrono::duration<double, std::milli>(
              std::chrono::steady_clock::now() - started)
              .count();
      const fluid::Solver::Stats stats = solver.stats();
      const double factor = 1000.0 * kFrameS / elapsedMs;
      const double iterationsPerSubstep =
          substeps > 0
              ? static_cast<double>(stats.pressureIterations) / substeps
              : 0.0;
      std::cout << "[quality timing] " << named.name << ' ' << budget
                << " dx=" << named.profile.spacing * 1000.0 << "mm"
                << " tolerance=" << config.densityTolerance * 100.0 << "%"
                << " particles=" << particles.size()
                << " ms/substep=" << stats.millisecondsPerSubstep
                << " iterations/substep=" << iterationsPerSubstep
                << " real-time=" << factor << "x"
                << " compression=" << stats.maxDensityCompression * 100.0 << "%"
                << " stalled=" << stats.stalledPressureSubsteps
                << " workers=" << stats.workerCount << '\n';
      checkFiniteAndContained(
          particles, boundary,
          config.contactRadiusFactor * config.resolution.spacing);
      return Timing{factor, stats};
    };

    const Timing before = measure("legacy", legacy);
    const fluid::SolverConfig configured = profileConfig(named.profile);
    const Timing after = measure("profile", configured);
    CHECK(after.stats.maxDensityCompression <= named.profile.densityTolerance);
    if (named.profile.densityTolerance > 5.0e-3) {
      CHECK(after.stats.pressureIterations < before.stats.pressureIterations);
    }
    const unsigned hardwareWorkers = std::thread::hardware_concurrency();
    const unsigned availableForPhysics =
        hardwareWorkers > 1 ? hardwareWorkers - 1 : 1;
    CHECK(after.stats.workerCount <= availableForPhysics);
    if (named.profile.spacing == fluid::QualityProfile::interactive().spacing) {
      CHECK(after.realTimeFactor >= 0.8);
    }
  }
}

TEST_CASE("a plateaued pressure correction is reported instead of burning its cap") {
  fluid::SolverConfig config =
      profileConfig(fluid::QualityProfile::interactive());
  config.densityTolerance = 1.0e-10;
  config.maxPressureIterations = 64;
  config.enableSurfaceTension = false;
  fluid::Solver solver;
  solver.configure(config);
  const std::vector<fluid::PhaseMaterial> phases{
      phase("water", 998.0, 0.89e-3, 100.0),
      phase("dichloromethane", 1326.0, 0.43e-3, 100.0)};
  solver.setPhases(phases, {});
  fluid::VesselBoundary boundary = separatoryFunnel(config);
  fluid::Particles particles;
  boundary.chargeLattice(phases, config.resolution.spacing, particles);

  solver.advance(particles, boundary, noGravity(), 0.0, 30.0e-3);
  CHECK(solver.stats().stalledPressureSubsteps > 0);
  CHECK(solver.stats().pressureIterations <
        solver.stats().substeps * config.maxPressureIterations);
  CHECK(solver.stats().pressureIterations >=
        solver.stats().substeps * config.minPressureIterations);
}

TEST_CASE("the solver remains finite under an abusive shake") {
  fluid::SolverConfig config = coarseConfig();
  fluid::Solver solver;
  solver.configure(config);
  const std::vector<fluid::PhaseMaterial> phases{
      phase("light", 850.0, 0.7e-3, 25.0), phase("dense", 1326.0, 1.5e-3, 25.0)};
  solver.setPhases(phases, {0.030});
  fluid::VesselBoundary boundary = cylinder(config);
  fluid::Particles particles;
  boundary.chargeLattice(phases, config.resolution.spacing, particles);
  fluid::VesselMotion abuse;
  abuse.shaking = true;
  abuse.shakeRemainingS = 3.0;
  abuse.shakeFrequencyHz = 12.0;
  abuse.shakeAmplitudeM = 0.15;
  abuse.shakeAxis = {0.4, -0.2, 1.0};

  const int substeps = solver.advance(particles, boundary, abuse, 0.0, 3.0);
  REQUIRE(substeps >= 1440);
  checkFiniteAndContained(particles, boundary,
                          config.contactRadiusFactor * config.resolution.spacing);
  CHECK(solver.stats().rejectedSubsteps <= substeps / 20);
}

TEST_CASE("vessel-frame shake and pose mathematics use the analytic laws") {
  fluid::VesselMotion motion;
  motion.shaking = true;
  motion.shakeRemainingS = 10.0;
  motion.shakeFrequencyHz = 2.75;
  motion.shakeAmplitudeM = 0.035;
  motion.shakeAxis = {2.0, -1.0, 2.0};
  const double axisNorm = 3.0;
  for (double time : {0.0, 0.013, 0.071, 0.19}) {
    const std::array<double, 3> acceleration = fluid::shakeAcceleration(motion, time);
    const double omega = 2.0 * fluid::kPi * motion.shakeFrequencyHz;
    const double scalar =
        -motion.shakeAmplitudeM * omega * omega * std::sin(omega * time) / axisNorm;
    CHECK(acceleration[0] == doctest::Approx(2.0 * scalar).epsilon(1.0e-12));
    CHECK(acceleration[1] == doctest::Approx(-scalar).epsilon(1.0e-12));
    CHECK(acceleration[2] == doctest::Approx(2.0 * scalar).epsilon(1.0e-12));
  }

  fluid::VesselMotion tilted;
  const double halfAngle = 0.25 * fluid::kPi;
  tilted.pose.orientation = {std::cos(halfAngle), 0.0, std::sin(halfAngle), 0.0};
  const fluid::FrameAcceleration frame = fluid::frameAcceleration(tilted, 0.0);
  CHECK(frame.uniform[0] == doctest::Approx(kGravity).epsilon(1.0e-12));
  CHECK(frame.uniform[1] == doctest::Approx(0.0).scale(1.0));
  CHECK(frame.uniform[2] == doctest::Approx(0.0).scale(1.0));

  const double expectedVelocity =
      2.0 * fluid::kPi * motion.shakeFrequencyHz * motion.shakeAmplitudeM;
  CHECK(fluid::shakePeakVelocity(motion) == doctest::Approx(expectedVelocity).epsilon(1.0e-12));
  CHECK(fluid::shakeSpecificPower(motion) ==
        doctest::Approx(0.5 * expectedVelocity * expectedVelocity * motion.shakeFrequencyHz)
            .epsilon(1.0e-12));
}

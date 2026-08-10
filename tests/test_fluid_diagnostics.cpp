#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <string>
#include <thread>
#include <vector>

#include "fluid/diagnostics.hpp"
#include "fluid/grid.hpp"
#include "fluid/kernels.hpp"
#include "fluid/simulation.hpp"
#include "fluid/types.hpp"
#include "fluid/vessel_sdf.hpp"
#include "sol/funnel.hpp"

namespace {
using chemcad::fluid::Diagnostics;
using chemcad::fluid::DiagnosticsEngine;
using chemcad::fluid::NeighbourGrid;
using chemcad::fluid::Particles;
using chemcad::fluid::PhaseMaterial;
using chemcad::fluid::Simulation;
using chemcad::fluid::VesselBoundary;

constexpr double kSpacing = 6.0e-3;
constexpr double kRatedVolumeMl = 250.0;

std::vector<PhaseMaterial> testPhases() {
  PhaseMaterial dense;
  dense.label = "dense";
  dense.restDensity = 1200.0;
  dense.dynamicViscosity = 2.0e-3;
  dense.volumeMl = 60.0;

  PhaseMaterial light;
  light.label = "light";
  light.restDensity = 800.0;
  light.dynamicViscosity = 1.0e-3;
  light.volumeMl = 60.0;
  return {dense, light};
}

double testVesselHeight() {
  chemcad::sol::Simulation sizing;
  sizing.vessel = chemcad::sol::Vessel::SeparatoryFunnel;
  sizing.vesselVolumeMl = kRatedVolumeMl;
  return chemcad::sol::columnHeightM(sizing);
}

VesselBoundary makeBoundary(double spacing = kSpacing) {
  VesselBoundary boundary;
  boundary.build(chemcad::sol::Vessel::SeparatoryFunnel, testVesselHeight(),
                 2.0 * spacing, spacing);
  return boundary;
}

struct DiagnosedState {
  Diagnostics diagnostics;
  std::size_t particleCount = 0;
};

DiagnosedState stratifiedState() {
  const auto phases = testPhases();
  VesselBoundary boundary = makeBoundary();
  Particles particles;
  boundary.chargeLattice(phases, kSpacing, particles);
  NeighbourGrid grid;
  grid.build(particles, 2.0 * kSpacing);
  DiagnosticsEngine engine;
  return {engine.compute(particles, grid, boundary, phases, kSpacing),
          particles.size()};
}

DiagnosedState boxState(bool dispersed) {
  const auto phases = testPhases();
  constexpr int halfCells = 9;
  constexpr int blockPeriod = 7;
  Particles particles;
  for (int iz = -halfCells; iz <= halfCells; ++iz) {
    for (int iy = -halfCells; iy <= halfCells; ++iy) {
      for (int ix = -halfCells; ix <= halfCells; ++ix) {
        uint8_t phase = iz >= 0 ? uint8_t{1} : uint8_t{0};
        if (dispersed) {
          const int bx = (ix + halfCells) % blockPeriod;
          const int by = (iy + halfCells) % blockPeriod;
          const int bz = (iz + halfCells) % blockPeriod;
          const bool inResolvedBlock =
              bx >= 1 && bx <= 3 && by >= 1 && by <= 3 &&
              bz >= 1 && bz <= 3;
          phase = inResolvedBlock ? uint8_t{1} : uint8_t{0};
        }
        particles.add(static_cast<float>(ix * kSpacing),
                      static_cast<float>(iy * kSpacing),
                      static_cast<float>(iz * kSpacing), phase);
      }
    }
  }

  NeighbourGrid grid;
  grid.build(particles, 2.0 * kSpacing);
  DiagnosticsEngine engine;
  VesselBoundary boundary = makeBoundary();
  return {engine.compute(particles, grid, boundary, phases, kSpacing),
          particles.size()};
}

float maximumSpeed(const std::shared_ptr<const chemcad::fluid::Snapshot>& snapshot) {
  return snapshot->speed.empty()
             ? 0.0f
             : *std::max_element(snapshot->speed.begin(), snapshot->speed.end());
}

void checkSnapshotShape(const chemcad::fluid::Snapshot& snapshot) {
  const std::size_t count = snapshot.px.size();
  CHECK(snapshot.py.size() == count);
  CHECK(snapshot.pz.size() == count);
  CHECK(snapshot.speed.size() == count);
  CHECK(snapshot.colour.size() == count);
  CHECK(snapshot.phase.size() == count);
}

std::array<std::size_t, 2> phaseCounts(
    const std::shared_ptr<const chemcad::fluid::Snapshot>& snapshot) {
  std::array<std::size_t, 2> counts{};
  if (!snapshot) return counts;
  for (uint8_t phase : snapshot->phase) {
    if (phase < counts.size()) ++counts[phase];
  }
  return counts;
}

void checkCharge(const Simulation& simulation, double spacing,
                 const std::array<double, 2>& volumeMl) {
  REQUIRE(simulation.charged());
  const auto state = simulation.snapshot();
  REQUIRE(state);
  const auto counts = phaseCounts(state);
  const double siteVolumeMl = spacing * spacing * spacing * 1.0e6;
  CHECK(counts[0] ==
        static_cast<std::size_t>(std::llround(volumeMl[0] / siteVolumeMl)));
  CHECK(counts[1] ==
        static_cast<std::size_t>(std::llround(volumeMl[1] / siteVolumeMl)));
}

void configureSimulation(Simulation& simulation) {
  simulation.setVessel(chemcad::sol::Vessel::SeparatoryFunnel,
                       kRatedVolumeMl);
  simulation.setResolution(10.0e-3);
  auto phases = testPhases();
  phases[0].volumeMl = 35.0;
  phases[1].volumeMl = 25.0;
  simulation.setPhases(phases, std::vector<double>(4, 0.0));
}

}  // namespace

TEST_CASE("stratified phases resolve bulk volumes and ordered layer tops") {
  const auto phases = testPhases();
  const DiagnosedState state = stratifiedState();
  const Diagnostics& diagnostics = state.diagnostics;
  REQUIRE(diagnostics.valid);
  REQUIRE(diagnostics.phases.size() == 2);

  CHECK(diagnostics.phases[0].bulkResolved);
  CHECK(diagnostics.phases[1].bulkResolved);
  CHECK(diagnostics.phases[0].bulkMl ==
        doctest::Approx(phases[0].volumeMl).epsilon(0.05));
  CHECK(diagnostics.phases[1].bulkMl ==
        doctest::Approx(phases[1].volumeMl).epsilon(0.05));
  CHECK(diagnostics.dispersedFraction < 0.02);
  CHECK(diagnostics.phases[0].layerTopM < diagnostics.phases[1].layerTopM);

  VesselBoundary boundary = makeBoundary();
  const double totalM3 = (phases[0].volumeMl + phases[1].volumeMl) * 1.0e-6;
  CHECK(std::abs(diagnostics.freeSurfaceM - boundary.heightForVolume(totalM3)) <=
        kSpacing);
}

TEST_CASE("resolved phase blocks expose more interface than a flat layer") {
  const DiagnosedState stratified = boxState(false);
  const DiagnosedState mixed = boxState(true);
  REQUIRE(mixed.diagnostics.valid);

  // Single-particle checkerboards are below H=2*dx and kernel-average back
  // toward a uniform colour, so their small reported gradient is physical
  // under-resolution rather than missing area. Resolved 3*dx blocks separated
  // by 4*dx measure 3.12 times the flat-interface result in this geometry.
  CHECK(mixed.diagnostics.dispersedFraction > 0.05);
  CHECK(mixed.diagnostics.dispersedComponents > 0);
  CHECK(mixed.diagnostics.sauterDiameterM > 0.0);
  CHECK(mixed.diagnostics.interfacialAreaM2 >
        3.0 * stratified.diagnostics.interfacialAreaM2);
}

TEST_CASE("colour-gradient area recovers a flat interface") {
  constexpr int halfCells = 9;
  const double side = (2 * halfCells + 1) * kSpacing;
  const double analyticArea = side * side;
  const Diagnostics diagnostics = boxState(false).diagnostics;
  const double measuredRatio = diagnostics.interfacialAreaM2 / analyticArea;
  CAPTURE(measuredRatio);
  // Measured A/A_flat=0.949 (5.1% low) for the 19*19-cell analytic plane.

  CHECK(diagnostics.interfacialAreaM2 ==
        doctest::Approx(analyticArea).epsilon(0.25));
}

TEST_CASE("colour-gradient area recovers the surface area of a sphere") {
  constexpr double spacing = 6.0e-3;
  constexpr double support = 2.0 * spacing;
  constexpr double radius = 3.0 * support;
  constexpr double padding = 2.0 * support;
  constexpr double centreZ = 0.09;
  const double halfWidth = radius + padding;

  Particles particles;
  for (double z = centreZ - halfWidth; z <= centreZ + halfWidth + 0.25 * spacing;
       z += spacing) {
    for (double y = -halfWidth; y <= halfWidth + 0.25 * spacing; y += spacing) {
      for (double x = -halfWidth; x <= halfWidth + 0.25 * spacing; x += spacing) {
        const double dz = z - centreZ;
        const bool inDroplet = x * x + y * y + dz * dz <= radius * radius;
        particles.add(static_cast<float>(x), static_cast<float>(y),
                      static_cast<float>(z), inDroplet ? uint8_t{1} : uint8_t{0});
      }
    }
  }

  auto phases = testPhases();
  NeighbourGrid grid;
  grid.build(particles, support);
  VesselBoundary boundary = makeBoundary(spacing);
  DiagnosticsEngine engine;
  const Diagnostics diagnostics =
      engine.compute(particles, grid, boundary, phases, spacing);
  const double analyticArea = 4.0 * chemcad::fluid::kPi * radius * radius;
  const double measuredRatio = diagnostics.interfacialAreaM2 / analyticArea;
  CAPTURE(measuredRatio);

  // R=3H is the least-resolved sphere admitted by the estimator contract.
  // Measured A/(4*pi*R^2)=0.971 (2.9% low); the 25% band also covers
  // voxelisation at this deliberately coarse test resolution.
  CHECK(diagnostics.interfacialAreaM2 ==
        doctest::Approx(analyticArea).epsilon(0.25));
}

TEST_CASE("display smoothing is gradual and leaves conserved diagnostics raw") {
  DiagnosticsEngine engine;
  Diagnostics initial;
  initial.valid = true;
  initial.phases.resize(2);
  initial.phases[0].totalMl = 40.0;
  initial.phases[0].bulkMl = 40.0;
  initial.phases[0].layerTopM = 0.02;
  initial.phases[1].totalMl = 30.0;
  initial.phases[1].bulkMl = 30.0;
  initial.phases[1].layerTopM = 0.04;
  initial.freeSurfaceM = 0.04;
  initial.interfacialAreaM2 = 0.001;

  const Diagnostics seeded = engine.smooth(initial, 0.1);
  Diagnostics step = initial;
  step.phases[0].totalMl = 41.0;
  step.phases[0].bulkMl = 39.0;
  step.phases[0].layerTopM = 0.08;
  step.phases[1].layerTopM = 0.12;
  step.freeSurfaceM = 0.12;
  step.dispersedFraction = 0.6;
  step.sauterDiameterM = 0.01;
  step.interfacialAreaM2 = 0.009;
  step.kineticEnergyJ = 2.5;

  Diagnostics filtered = engine.smooth(step, 0.03);
  CHECK(filtered.phases[0].layerTopM > seeded.phases[0].layerTopM);
  CHECK(filtered.phases[0].layerTopM < step.phases[0].layerTopM);
  CHECK(filtered.dispersedFraction > 0.0);
  CHECK(filtered.dispersedFraction < step.dispersedFraction);
  CHECK(filtered.sauterDiameterM > 0.0);
  CHECK(filtered.sauterDiameterM < step.sauterDiameterM);
  CHECK(filtered.phases[0].totalMl == step.phases[0].totalMl);
  CHECK(filtered.phases[0].bulkMl == step.phases[0].bulkMl);
  CHECK(filtered.interfacialAreaM2 == step.interfacialAreaM2);
  CHECK(filtered.kineticEnergyJ == step.kineticEnergyJ);

  for (int i = 0; i < 20; ++i) filtered = engine.smooth(step, 0.15);
  CHECK(filtered.phases[0].layerTopM ==
        doctest::Approx(step.phases[0].layerTopM).epsilon(1.0e-6));
  CHECK(filtered.dispersedFraction ==
        doctest::Approx(step.dispersedFraction).epsilon(1.0e-6));
  CHECK(filtered.sauterDiameterM ==
        doctest::Approx(step.sauterDiameterM).epsilon(1.0e-6));
}

TEST_CASE("simulation setup order is exception-free and preserves phase charges") {
  auto phases = testPhases();
  phases[0].volumeMl = 35.0;
  phases[1].volumeMl = 25.0;
  const std::array<double, 2> volumes{35.0, 25.0};
  const std::vector<double> sigma(4, 0.0);

  SUBCASE("panel order: vessel and resolution before phases") {
    Simulation simulation;
    CHECK_NOTHROW(simulation.setVessel(
        chemcad::sol::Vessel::SeparatoryFunnel, kRatedVolumeMl));
    CHECK_NOTHROW(simulation.setResolution(10.0e-3));
    CHECK_NOTHROW(simulation.charge());
    CHECK_FALSE(simulation.charged());
    const double elapsed = simulation.elapsedS();
    CHECK_NOTHROW(simulation.advance(0.02));
    CHECK(simulation.elapsedS() == elapsed);
    CHECK(simulation.statusLine().find("configure") != std::string::npos);

    CHECK_NOTHROW(simulation.setPhases(phases, sigma));
    checkCharge(simulation, 10.0e-3, volumes);
  }

  SUBCASE("phases before vessel") {
    Simulation simulation;
    CHECK_NOTHROW(simulation.setPhases(phases, sigma));
    CHECK_NOTHROW(simulation.setVessel(
        chemcad::sol::Vessel::SeparatoryFunnel, kRatedVolumeMl));
    checkCharge(simulation, 4.0e-3, volumes);
  }

  SUBCASE("resolution change after charging") {
    Simulation simulation;
    CHECK_NOTHROW(simulation.setResolution(8.0e-3));
    CHECK_NOTHROW(simulation.setPhases(phases, sigma));
    REQUIRE(simulation.charged());
    CHECK_NOTHROW(simulation.setResolution(10.0e-3));
    checkCharge(simulation, 10.0e-3, volumes);
    CHECK(simulation.totalVolumeMl() == doctest::Approx(60.0));
  }

  SUBCASE("a surface-tension table survives setup and charges") {
    Simulation simulation;
    CHECK_NOTHROW(simulation.setResolution(10.0e-3));
    CHECK_NOTHROW(simulation.setPhases(
        phases, std::vector<double>{0.0, 0.03, 0.03, 0.0}));
    checkCharge(simulation, 10.0e-3, volumes);
    CHECK_FALSE(simulation.statusLine().empty());
  }

  SUBCASE("empty phase table is reported and remains uncharged") {
    Simulation simulation;
    CHECK_NOTHROW(simulation.setPhases({}, {}));
    CHECK_NOTHROW(simulation.charge());
    CHECK_NOTHROW(simulation.advance(0.01));
    CHECK_FALSE(simulation.charged());
    CHECK(simulation.statusLine().find("configure") != std::string::npos);
  }
}

TEST_CASE("simulation charge is deterministic and reports its configured volume") {
  Simulation simulation;
  configureSimulation(simulation);
  simulation.charge();
  const auto first = simulation.snapshot();
  const Diagnostics firstDiagnostics = simulation.diagnostics();

  simulation.charge();
  const auto second = simulation.snapshot();
  const Diagnostics secondDiagnostics = simulation.diagnostics();
  REQUIRE(first);
  REQUIRE(second);
  CHECK(first->px.size() == second->px.size());
  REQUIRE(firstDiagnostics.phases.size() == secondDiagnostics.phases.size());
  for (std::size_t phase = 0; phase < firstDiagnostics.phases.size(); ++phase) {
    CHECK(firstDiagnostics.phases[phase].totalMl ==
          secondDiagnostics.phases[phase].totalMl);
    CHECK(firstDiagnostics.phases[phase].bulkMl ==
          secondDiagnostics.phases[phase].bulkMl);
    CHECK(firstDiagnostics.phases[phase].dispersedMl ==
          secondDiagnostics.phases[phase].dispersedMl);
  }
  CHECK(firstDiagnostics.interfacialAreaM2 ==
        secondDiagnostics.interfacialAreaM2);
  CHECK(simulation.totalVolumeMl() == doctest::Approx(60.0));
  const std::string status = simulation.statusLine();
  CHECK_FALSE(status.empty());
  CHECK(status.find(std::to_string(second->px.size())) != std::string::npos);
}

TEST_CASE("vertical shake publishes complete snapshots and then decays") {
  Simulation simulation;
  configureSimulation(simulation);
  const auto before = simulation.snapshot();
  REQUIRE(before);
  const uint64_t initialRevision = before->revision;

  simulation.shake({0.0, 0.0, 1.0}, 0.30, 4.0, 0.04);
  float drivenPeak = 0.0f;
  for (int i = 0; i < 6; ++i) {
    simulation.advance(0.05);
    drivenPeak = std::max(drivenPeak, maximumSpeed(simulation.snapshot()));
  }
  const auto driven = simulation.snapshot();
  REQUIRE(driven);
  CHECK(driven->revision > initialRevision);
  checkSnapshotShape(*driven);
  CHECK(drivenPeak > 0.0f);
  CHECK_FALSE(simulation.shaking());

  for (int i = 0; i < 80; ++i) simulation.advance(0.05);
  const auto settled = simulation.snapshot();
  REQUIRE(settled);
  checkSnapshotShape(*settled);
  CHECK(maximumSpeed(settled) < drivenPeak);
}

TEST_CASE("published snapshots remain internally consistent during advance") {
  Simulation simulation;
  configureSimulation(simulation);
  simulation.shake({0.0, 0.0, 1.0}, 0.5, 3.0, 0.02);
  std::atomic<bool> start{false};
  std::thread worker([&] {
    while (!start.load(std::memory_order_acquire)) std::this_thread::yield();
    for (int i = 0; i < 200; ++i) simulation.advance(0.002);
  });

  start.store(true, std::memory_order_release);
  uint64_t lastRevision = 0;
  for (int i = 0; i < 200; ++i) {
    const auto state = simulation.snapshot();
    CHECK(state != nullptr);
    if (state) {
      checkSnapshotShape(*state);
      CHECK(state->revision >= lastRevision);
      lastRevision = state->revision;
    }
    std::this_thread::yield();
  }
  worker.join();
  const auto finalState = simulation.snapshot();
  REQUIRE(finalState);
  checkSnapshotShape(*finalState);
  CHECK(finalState->revision >= lastRevision);
}

TEST_CASE("asynchronous requests use the deterministic synchronous integration path") {
  Simulation asynchronous;
  Simulation synchronous;
  configureSimulation(asynchronous);
  configureSimulation(synchronous);
  asynchronous.shake({0.0, 0.0, 1.0}, 0.10, 3.0, 0.02);
  synchronous.shake({0.0, 0.0, 1.0}, 0.10, 3.0, 0.02);

  asynchronous.requestAdvance(0.02);
  asynchronous.waitForIdle();
  synchronous.advance(0.02);

  const auto asyncState = asynchronous.snapshot();
  const auto syncState = synchronous.snapshot();
  REQUIRE(asyncState);
  REQUIRE(syncState);
  CHECK(asyncState->elapsedS == syncState->elapsedS);
  CHECK(asyncState->px == syncState->px);
  CHECK(asyncState->py == syncState->py);
  CHECK(asyncState->pz == syncState->pz);
  CHECK(asyncState->speed == syncState->speed);
  CHECK(asyncState->colour == syncState->colour);
  CHECK(asyncState->phase == syncState->phase);
  CHECK(asynchronous.solverStats().substeps ==
        synchronous.solverStats().substeps);
  CHECK(asynchronous.solverStats().maxDensityError ==
        synchronous.solverStats().maxDensityError);
}

TEST_CASE("advance bursts coalesce into a bounded backlog") {
  Simulation simulation;
  configureSimulation(simulation);
  for (int i = 0; i < 128; ++i) simulation.requestAdvance(1.0);

  CHECK(simulation.pendingSeconds() <= 0.1);
  simulation.waitForIdle();
  CHECK(simulation.pendingSeconds() == 0.0);
  CHECK_FALSE(simulation.stepping());
  CHECK(simulation.elapsedS() <= 0.100000000001);
}

TEST_CASE("destroying a simulation joins an in-flight step") {
  CHECK_NOTHROW([] {
    Simulation simulation;
    configureSimulation(simulation);
    simulation.shake({1.0, 0.0, 0.0}, 0.1, 3.0, 0.02);
    simulation.requestAdvance(0.1);
    for (int attempt = 0; attempt < 100 && !simulation.stepping(); ++attempt)
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    CHECK(simulation.stepping());
  }());
}

TEST_CASE("real-time factor is finite after a completed step") {
  Simulation simulation;
  configureSimulation(simulation);
  simulation.advance(0.01);

  CHECK(std::isfinite(simulation.realTimeFactor()));
  CHECK(simulation.realTimeFactor() > 0.0);
  // A factor above 1 is not an error: it means the solver integrated the
  // requested simulated time faster than wall clock, which is exactly what a
  // coarse test charge does. Only a non-finite or non-positive value is a bug.
  CHECK(simulation.statusLine().find("physics at") != std::string::npos);
}

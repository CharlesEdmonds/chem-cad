#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "fluid/grid.hpp"
#include "fluid/kernels.hpp"
#include "fluid/types.hpp"
#include "fluid/vessel_sdf.hpp"
#include "sol/funnel.hpp"

using namespace chemcad;

namespace {

double integrateSpherical(double support, bool wendland) {
  constexpr int kShells = 200000;
  const double dr = support / kShells;
  double integral = 0.0;
  for (int i = 0; i < kShells; ++i) {
    const double r = (static_cast<double>(i) + 0.5) * dr;
    const double w = wendland ? fluid::wendlandW(r, support)
                              : fluid::cubicW(r, support);
    integral += 4.0 * fluid::kPi * r * r * w * dr;
  }
  return integral;
}

struct FixedPositions {
  fluid::Particles particles;
  std::vector<float> xById;
  std::vector<float> yById;
  std::vector<float> zById;
};

FixedPositions makeFixedPositions(std::size_t count) {
  FixedPositions result;
  result.particles.reserve(count);
  result.xById.resize(count);
  result.yById.resize(count);
  result.zById.resize(count);
  uint32_t state = 0x6d2b79f5U;
  const auto coordinate = [&]() {
    state = state * 1664525U + 1013904223U;
    const double unit = static_cast<double>(state) /
                        static_cast<double>(std::numeric_limits<uint32_t>::max());
    return static_cast<float>((2.0 * unit - 1.0) * 0.080);
  };
  for (std::size_t i = 0; i < count; ++i) {
    const float x = coordinate();
    const float y = coordinate();
    const float z = coordinate();
    result.particles.add(x, y, z, 0);
    const uint32_t id = result.particles.id.back();
    if (id < count) {
      result.xById[id] = x;
      result.yById[id] = y;
      result.zById[id] = z;
    }
  }
  return result;
}

std::vector<uint32_t> candidateIds(const fluid::NeighbourGrid& grid,
                                   const fluid::Particles& particles,
                                   std::size_t i) {
  std::vector<uint32_t> result;
  grid.forEachCandidate(i, [&](std::size_t j) { result.push_back(particles.id[j]); });
  return result;
}

std::array<sol::Vessel, 3> allVessels() {
  return {sol::Vessel::SeparatoryFunnel, sol::Vessel::DecantingFlask,
          sol::Vessel::GraduatedCylinder};
}

double ratedHeight(sol::Vessel vessel, double capacityMl = 250.0) {
  sol::Simulation simulation;
  simulation.vessel = vessel;
  simulation.vesselVolumeMl = capacityMl;
  return sol::columnHeightM(simulation);
}

}  // namespace

TEST_CASE("SPH kernels are normalized compact supports with inward gradients") {
  constexpr double support = 8.0e-3;
  CHECK(integrateSpherical(support, true) == doctest::Approx(1.0).epsilon(0.01));
  CHECK(integrateSpherical(support, false) == doctest::Approx(1.0).epsilon(0.01));
  CHECK(fluid::wendlandW(support, support) == 0.0);
  CHECK(fluid::cubicW(support, support) == 0.0);

  CHECK(fluid::wendlandGradMagnitude(0.0, support) == 0.0);
  CHECK(fluid::wendlandGradMagnitude(support, support) == 0.0);
  for (double q : {0.10, 0.37, 0.75, 0.99}) {
    CHECK(fluid::wendlandGradMagnitude(q * support, support) < 0.0);
  }
  CHECK(fluid::cubicGradMagnitude(0.0, support) == 0.0);
  CHECK(fluid::cubicGradMagnitude(support, support) == 0.0);
  for (double q : {0.10, 0.37, 0.75, 0.99}) {
    CHECK(fluid::cubicGradMagnitude(q * support, support) < 0.0);
  }
}

TEST_CASE("rest number density follows the cubic lattice resolution") {
  constexpr double spacing = 4.0e-3;
  constexpr double support = 2.0 * spacing;
  int supportSites = 0;
  for (int z = -2; z <= 2; ++z) {
    for (int y = -2; y <= 2; ++y) {
      for (int x = -2; x <= 2; ++x) {
        const double r2 = static_cast<double>(x * x + y * y + z * z) *
                          spacing * spacing;
        if (r2 <= support * support) ++supportSites;
      }
    }
  }
  // The six sites exactly at H belong to the compact-support neighbour set and
  // carry W(H)=0 by the identity tested above; the complete interior stencil
  // therefore has the expected 30-45 candidate sites.
  CHECK(supportSites >= 30);
  CHECK(supportSites <= 45);

  const double coarse = fluid::restNumberDensity(spacing, support);
  const double fine = fluid::restNumberDensity(0.5 * spacing, 0.5 * support);
  CHECK(coarse > 0.0);
  CHECK(fine / coarse == doctest::Approx(8.0).epsilon(0.02));
}

TEST_CASE("Morton grid candidates exactly contain every true neighbour") {
  constexpr std::size_t count = 2000;
  constexpr double support = 8.0e-3;
  FixedPositions fixed = makeFixedPositions(count);
  fluid::NeighbourGrid grid;
  grid.build(fixed.particles, support);
  REQUIRE(grid.size() == count);

  for (std::size_t i = 0; i < count; ++i) {
    const uint32_t id = fixed.particles.id[i];
    REQUIRE(id < count);
    CHECK(fixed.particles.px[i] == fixed.xById[id]);
    CHECK(fixed.particles.py[i] == fixed.yById[id]);
    CHECK(fixed.particles.pz[i] == fixed.zById[id]);

    std::vector<uint32_t> found;
    grid.forEachCandidate(i, [&](std::size_t j) {
      const double dx = static_cast<double>(fixed.particles.px[i]) -
                        fixed.particles.px[j];
      const double dy = static_cast<double>(fixed.particles.py[i]) -
                        fixed.particles.py[j];
      const double dz = static_cast<double>(fixed.particles.pz[i]) -
                        fixed.particles.pz[j];
      if (dx * dx + dy * dy + dz * dz < support * support) {
        found.push_back(fixed.particles.id[j]);
      }
    });

    std::vector<uint32_t> brute;
    for (std::size_t j = 0; j < count; ++j) {
      const double dx = static_cast<double>(fixed.particles.px[i]) -
                        fixed.particles.px[j];
      const double dy = static_cast<double>(fixed.particles.py[i]) -
                        fixed.particles.py[j];
      const double dz = static_cast<double>(fixed.particles.pz[i]) -
                        fixed.particles.pz[j];
      if (dx * dx + dy * dy + dz * dz < support * support) {
        brute.push_back(fixed.particles.id[j]);
      }
    }
    std::sort(found.begin(), found.end());
    std::sort(brute.begin(), brute.end());
    CHECK(found == brute);
  }
}

TEST_CASE("Morton builds and candidate sequences are deterministic") {
  constexpr double support = 8.0e-3;
  FixedPositions fixed = makeFixedPositions(2000);
  fluid::Particles first = fixed.particles;
  fluid::Particles second = fixed.particles;
  fluid::NeighbourGrid firstGrid;
  fluid::NeighbourGrid secondGrid;
  firstGrid.build(first, support);
  secondGrid.build(second, support);

  CHECK(first.id == second.id);
  CHECK(first.px == second.px);
  CHECK(first.py == second.py);
  CHECK(first.pz == second.pz);
  for (std::size_t i = 0; i < first.size(); ++i) {
    CHECK(candidateIds(firstGrid, first, i) == candidateIds(secondGrid, second, i));
  }
}

TEST_CASE("vessel SDF is exact on its polyline and has a signed unit normal") {
  constexpr double support = 8.0e-3;
  constexpr double spacing = 4.0e-3;
  for (sol::Vessel vessel : allVessels()) {
    const double height = ratedHeight(vessel);
    fluid::VesselBoundary boundary;
    boundary.build(vessel, height, support, spacing);
    for (int sample = 0; sample < 200; ++sample) {
      const double z = height * static_cast<double>(sample) / 199.0;
      const double radius = boundary.radiusAt(z);
      const fluid::SurfaceQuery wall = boundary.query(radius, 0.0, z);
      CHECK(std::abs(wall.distance) <= 1.0e-9);
      const double normalLength = std::sqrt(wall.nx * wall.nx + wall.ny * wall.ny +
                                            wall.nz * wall.nz);
      CHECK(normalLength == doctest::Approx(1.0).epsilon(1.0e-10));
    }
    CHECK(boundary.query(0.0, 0.0, 0.5 * height).distance < 0.0);
    const double midRadius = boundary.radiusAt(0.5 * height);
    CHECK(boundary.query(midRadius + 0.02, 0.0, 0.5 * height).distance > 0.0);
  }

  fluid::VesselBoundary cylinder;
  const double height = ratedHeight(sol::Vessel::GraduatedCylinder);
  cylinder.build(sol::Vessel::GraduatedCylinder, height, support, spacing);
  const double z = 0.5 * height;
  const double radius = cylinder.radiusAt(z);
  for (double radial : {0.75 * radius, 1.25 * radius}) {
    const double expectedSignedDistance = radial - radius;
    CHECK(cylinder.query(radial, 0.0, z).distance ==
          doctest::Approx(expectedSignedDistance).epsilon(1.0e-9));
  }
}

TEST_CASE("vessel volume table reproduces capacity and inverts height") {
  constexpr double capacityMl = 250.0;
  for (sol::Vessel vessel : allVessels()) {
    const double height = ratedHeight(vessel, capacityMl);
    fluid::VesselBoundary boundary;
    boundary.build(vessel, height, 8.0e-3, 4.0e-3);
    CHECK(boundary.volumeBelow(height) ==
          doctest::Approx(capacityMl * 1.0e-6).epsilon(0.01));
    for (int sample = 0; sample < 50; ++sample) {
      const double z = height * static_cast<double>(sample) / 49.0;
      const double recovered = boundary.heightForVolume(boundary.volumeBelow(z));
      CHECK(std::abs(recovered - z) <= height / 1024.0);
    }
  }
}

TEST_CASE("lattice charging fills dense phase below light phase") {
  constexpr double spacing = 4.0e-3;
  fluid::VesselBoundary boundary;
  boundary.build(sol::Vessel::SeparatoryFunnel,
                 ratedHeight(sol::Vessel::SeparatoryFunnel, 250.0),
                 2.0 * spacing, spacing);

  std::vector<fluid::PhaseMaterial> phases(2);
  phases[0].label = "Light";
  phases[0].restDensity = 800.0;
  phases[0].volumeMl = 100.0;
  phases[1].label = "Dense";
  phases[1].restDensity = 1100.0;
  phases[1].volumeMl = 100.0;

  fluid::Particles particles;
  const std::size_t charged = boundary.chargeLattice(phases, spacing, particles);
  REQUIRE(charged == particles.size());
  double zSum[2] = {0.0, 0.0};
  std::size_t phaseCount[2] = {0, 0};
  for (std::size_t i = 0; i < particles.size(); ++i) {
    REQUIRE(particles.phase[i] < 2);
    CHECK(boundary.query(particles.px[i], particles.py[i], particles.pz[i]).distance <
          0.0);
    zSum[particles.phase[i]] += particles.pz[i];
    ++phaseCount[particles.phase[i]];
  }
  REQUIRE(phaseCount[0] > 0);
  REQUIRE(phaseCount[1] > 0);
  CHECK(zSum[1] / phaseCount[1] < zSum[0] / phaseCount[0]);

  const double expected = 200.0e-6 / (spacing * spacing * spacing);
  CHECK(static_cast<double>(charged) == doctest::Approx(expected).epsilon(0.15));
}

TEST_CASE("contact-radius charging keeps both phases in a narrow vessel") {
  constexpr double spacing = 8.0e-3;
  fluid::VesselBoundary boundary;
  boundary.build(sol::Vessel::GraduatedCylinder, 0.12, 2.0 * spacing, spacing);

  std::vector<fluid::PhaseMaterial> phases(2);
  phases[0].label = "Light";
  phases[0].restDensity = 800.0;
  phases[0].volumeMl = 100.0;
  phases[1].label = "Dense";
  phases[1].restDensity = 1200.0;
  phases[1].volumeMl = 100.0;

  fluid::Particles particles;
  const std::size_t charged = boundary.chargeLattice(phases, spacing, particles);
  std::array<std::size_t, 2> phaseCount{};
  for (uint8_t phase : particles.phase) {
    REQUIRE(phase < phaseCount.size());
    ++phaseCount[phase];
  }

  CHECK(charged > 30);
  CHECK(phaseCount[0] > 0);
  CHECK(phaseCount[1] > 0);
  const double requested =
      std::round(200.0e-6 / (spacing * spacing * spacing));
  CHECK(static_cast<double>(charged) < requested);
}

TEST_CASE("lattice charging honours each phase volume allocation") {
  constexpr double spacing = 4.0e-3;
  fluid::VesselBoundary boundary;
  boundary.build(sol::Vessel::SeparatoryFunnel,
                 ratedHeight(sol::Vessel::SeparatoryFunnel, 250.0),
                 2.0 * spacing, spacing);

  std::vector<fluid::PhaseMaterial> phases(2);
  phases[0].restDensity = 850.0;
  phases[0].volumeMl = 60.0;
  phases[1].restDensity = 1100.0;
  phases[1].volumeMl = 120.0;

  fluid::Particles particles;
  boundary.chargeLattice(phases, spacing, particles);
  std::array<std::size_t, 2> phaseCount{};
  for (uint8_t phase : particles.phase) {
    REQUIRE(phase < phaseCount.size());
    ++phaseCount[phase];
  }

  const double siteVolume = spacing * spacing * spacing;
  const double expectedLight = std::round(phases[0].volumeMl * 1.0e-6 / siteVolume);
  const double expectedDense = std::round(phases[1].volumeMl * 1.0e-6 / siteVolume);
  REQUIRE(phaseCount[0] > 0);
  REQUIRE(phaseCount[1] > 0);
  CHECK(static_cast<double>(phaseCount[0]) ==
        doctest::Approx(expectedLight).epsilon(0.10));
  CHECK(static_cast<double>(phaseCount[1]) ==
        doctest::Approx(expectedDense).epsilon(0.10));
  CHECK(static_cast<double>(phaseCount[1]) / phaseCount[0] ==
        doctest::Approx(phases[1].volumeMl / phases[0].volumeMl).epsilon(0.10));
}

TEST_CASE("charged lattice respects the solver contact-radius wall margin") {
  constexpr double spacing = 8.0e-3;
  constexpr double contactRadius = 0.35 * spacing;
  fluid::VesselBoundary boundary;
  boundary.build(sol::Vessel::GraduatedCylinder, 0.12, 2.0 * spacing, spacing);

  std::vector<fluid::PhaseMaterial> phases(1);
  phases[0].restDensity = 998.0;
  phases[0].volumeMl = 20.0;

  fluid::Particles particles;
  REQUIRE(boundary.chargeLattice(phases, spacing, particles) > 0);
  for (std::size_t i = 0; i < particles.size(); ++i) {
    CHECK(boundary.query(particles.px[i], particles.py[i], particles.pz[i]).distance <=
          -contactRadius + 2.0e-8);
  }
}

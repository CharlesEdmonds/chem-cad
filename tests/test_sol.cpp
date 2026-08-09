#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

#include "chem/bridge.hpp"
#include "core/model.hpp"
#include "sol/funnel.hpp"
#include "sol/solubility.hpp"
#include "sol/solvent.hpp"

namespace chem = chemcad::chem;
namespace core = chemcad::core;
namespace sol = chemcad::sol;

namespace {

// Parcel volume is the conserved bulk quantity. Physical radius belongs only
// to breakup and settling and must never be used for volume bookkeeping.
double dropletVolumeMl(const sol::Simulation& sim) {
  double volume = 0.0;
  for (const sol::Droplet& d : sim.droplets) volume += double(d.parcelMl);
  return volume;
}

double settledVolumeMl(const sol::Simulation& sim) {
  double volume = 0.0;
  for (double v : sim.settledMl) volume += v;
  return volume;
}

}  // namespace

TEST_CASE("solvent database loads with unique ids and known entries") {
  const std::vector<sol::Solvent>& db = sol::solvents();
  REQUIRE(db.size() >= 40);

  std::vector<std::string> ids;
  ids.reserve(db.size());
  for (const sol::Solvent& s : db) ids.push_back(s.id);
  std::sort(ids.begin(), ids.end());
  CHECK(std::adjacent_find(ids.begin(), ids.end()) == ids.end());

  CHECK(sol::findSolvent("water") != nullptr);
  CHECK(sol::findSolvent("toluene") != nullptr);
  CHECK(sol::findSolvent("definitely_not_a_solvent") == nullptr);
}

TEST_CASE("solvent molar volume agrees with molar mass over density") {
  for (const sol::Solvent& s : sol::solvents()) {
    REQUIRE(s.density > 0.0);
    const double implied = s.molarMass / s.density;
    CAPTURE(s.id);
    CHECK(s.molarVolume == doctest::Approx(implied).epsilon(0.05));
  }
}

TEST_CASE("blend combines solvent properties by volume fraction") {
  const sol::Solvent* water = sol::findSolvent("water");
  const sol::Solvent* ethanol = sol::findSolvent("ethanol");
  REQUIRE(water != nullptr);
  REQUIRE(ethanol != nullptr);

  // A single component reproduces that solvent's numbers exactly.
  const sol::Mixture pure = sol::blend({{water, 1.0}});
  CHECK(pure.hansen.dispersion == doctest::Approx(water->hansen.dispersion).epsilon(1e-9));
  CHECK(pure.hansen.polar == doctest::Approx(water->hansen.polar).epsilon(1e-9));
  CHECK(pure.hansen.hydrogenBond == doctest::Approx(water->hansen.hydrogenBond).epsilon(1e-9));
  CHECK(pure.density == doctest::Approx(water->density).epsilon(1e-9));

  // Equal fractions land exactly halfway between the two solvents.
  const sol::Mixture half = sol::blend({{water, 0.5}, {ethanol, 0.5}});
  CHECK(half.hansen.dispersion ==
        doctest::Approx((water->hansen.dispersion + ethanol->hansen.dispersion) / 2.0)
            .epsilon(1e-9));
  CHECK(half.hansen.polar ==
        doctest::Approx((water->hansen.polar + ethanol->hansen.polar) / 2.0).epsilon(1e-9));
  CHECK(half.hansen.hydrogenBond ==
        doctest::Approx((water->hansen.hydrogenBond + ethanol->hansen.hydrogenBond) / 2.0)
            .epsilon(1e-9));
  CHECK(half.density ==
        doctest::Approx((water->density + ethanol->density) / 2.0).epsilon(1e-9));

  // Unnormalised 3:1 ratios behave the same as normalised 75:25 fractions.
  const sol::Mixture skewed = sol::blend({{water, 3.0}, {ethanol, 1.0}});
  const double expectedDispersion =
      0.75 * water->hansen.dispersion + 0.25 * ethanol->hansen.dispersion;
  const double expectedDensity = 0.75 * water->density + 0.25 * ethanol->density;
  CHECK(skewed.hansen.dispersion == doctest::Approx(expectedDispersion).epsilon(1e-9));
  CHECK(skewed.density == doctest::Approx(expectedDensity).epsilon(1e-9));

  // Empty and all-null inputs must never throw and must yield a zeroed blend.
  sol::Mixture empty;
  CHECK_NOTHROW(empty = sol::blend({}));
  CHECK(empty.density == doctest::Approx(0.0).epsilon(1e-9));
  CHECK(empty.hansen.dispersion == doctest::Approx(0.0).epsilon(1e-9));
  CHECK(empty.hansen.polar == doctest::Approx(0.0).epsilon(1e-9));
  CHECK(empty.hansen.hydrogenBond == doctest::Approx(0.0).epsilon(1e-9));

  sol::Mixture allNull;
  CHECK_NOTHROW(allNull = sol::blend({{nullptr, 1.0}, {nullptr, 2.0}}));
  CHECK(allNull.density == doctest::Approx(0.0).epsilon(1e-9));
}

TEST_CASE("predict scores an exact Hansen match as fully soluble and penalises distance") {
  const sol::Solvent* water = sol::findSolvent("water");
  REQUIRE(water != nullptr);

  sol::Solute matched;
  matched.name = "matched";
  matched.molarMass = 100.0;
  matched.molarVolume = 90.0;
  matched.hansen = water->hansen;
  matched.meltingPoint = 20.0;  // below 25 C: treated as a liquid, no ideal-solubility penalty
  matched.interactionRadius = 8.0;

  const std::vector<sol::Component> pureWater{{water, 1.0}};
  const sol::Prediction hit = sol::predict(matched, pureWater);
  CHECK(hit.ra == doctest::Approx(0.0).epsilon(1e-6));
  CHECK(hit.relativeEnergyDifference == doctest::Approx(0.0).epsilon(1e-6));
  CHECK(hit.activityCoefficient == doctest::Approx(1.0).epsilon(1e-6));
  CHECK_FALSE(hit.outsideSphere);

  sol::Solute distant = matched;
  distant.hansen.dispersion += 30.0;
  distant.hansen.polar += 30.0;
  distant.hansen.hydrogenBond += 30.0;

  const sol::Prediction miss = sol::predict(distant, pureWater);
  CHECK(miss.ra > matched.interactionRadius);
  CHECK(miss.outsideSphere);
  CHECK(miss.gramsPerMillilitre < hit.gramsPerMillilitre);
}

TEST_CASE("predict stays finite for degenerate solute, mixture and temperature input") {
  const sol::Solvent* water = sol::findSolvent("water");
  REQUIRE(water != nullptr);

  sol::Solute zeroVolume;
  zeroVolume.molarMass = 50.0;
  zeroVolume.molarVolume = 0.0;
  zeroVolume.hansen = water->hansen;

  const sol::Prediction zeroVolPrediction = sol::predict(zeroVolume, {{water, 1.0}});
  CHECK(std::isfinite(zeroVolPrediction.gramsPerMillilitre));
  CHECK(std::isfinite(zeroVolPrediction.moleFraction));
  CHECK(std::isfinite(zeroVolPrediction.ra));
  CHECK(std::isfinite(zeroVolPrediction.activityCoefficient));

  sol::Solute normal;
  normal.molarMass = 120.0;
  normal.molarVolume = 100.0;
  normal.hansen = water->hansen;

  const sol::Prediction emptyMixPrediction = sol::predict(normal, {});
  CHECK(std::isfinite(emptyMixPrediction.gramsPerMillilitre));
  CHECK(std::isfinite(emptyMixPrediction.moleFraction));
  CHECK(std::isfinite(emptyMixPrediction.ra));
  CHECK(std::isfinite(emptyMixPrediction.activityCoefficient));

  const sol::Prediction hotPrediction = sol::predict(normal, {{water, 1.0}}, 5000.0);
  CHECK(std::isfinite(hotPrediction.gramsPerMillilitre));
  CHECK(std::isfinite(hotPrediction.activityCoefficient));

  const sol::Prediction coldPrediction = sol::predict(normal, {{water, 1.0}}, -270.0);
  CHECK(std::isfinite(coldPrediction.gramsPerMillilitre));
  CHECK(std::isfinite(coldPrediction.activityCoefficient));
}

TEST_CASE("melting point above 25 C strictly lowers solubility, at or below does not") {
  const sol::Solvent* water = sol::findSolvent("water");
  REQUIRE(water != nullptr);

  sol::Solute liquidAt25;
  liquidAt25.molarMass = 120.0;
  liquidAt25.molarVolume = 100.0;
  liquidAt25.hansen = water->hansen;
  liquidAt25.meltingPoint = 25.0;

  sol::Solute liquidBelow25 = liquidAt25;
  liquidBelow25.meltingPoint = -20.0;

  sol::Solute solidAbove25 = liquidAt25;
  solidAbove25.meltingPoint = 180.0;

  const std::vector<sol::Component> pureWater{{water, 1.0}};
  const sol::Prediction predAt25 = sol::predict(liquidAt25, pureWater);
  const sol::Prediction predBelow25 = sol::predict(liquidBelow25, pureWater);
  const sol::Prediction predAbove25 = sol::predict(solidAbove25, pureWater);

  CHECK(predAt25.gramsPerMillilitre ==
        doctest::Approx(predBelow25.gramsPerMillilitre).epsilon(1e-9));
  CHECK(predAbove25.gramsPerMillilitre < predAt25.gramsPerMillilitre);
}

TEST_CASE("sweep enumerates the simplex grid and validates its input") {
  const sol::Solvent* water = sol::findSolvent("water");
  const sol::Solvent* ethanol = sol::findSolvent("ethanol");
  const sol::Solvent* acetone = sol::findSolvent("acetone");
  REQUIRE(water != nullptr);
  REQUIRE(ethanol != nullptr);
  REQUIRE(acetone != nullptr);

  sol::Solute solute;
  solute.molarMass = 150.0;
  solute.molarVolume = 120.0;
  solute.hansen = water->hansen;

  const std::vector<sol::SweepPoint> one = sol::sweep(solute, {water}, 10);
  CHECK(one.size() == 1);

  const std::vector<sol::SweepPoint> two = sol::sweep(solute, {water, ethanol}, 10);
  CHECK(two.size() == 11);

  const std::vector<sol::SweepPoint> three = sol::sweep(solute, {water, ethanol, acetone}, 4);
  CHECK(three.size() == 15);

  for (const std::vector<sol::SweepPoint>* points : {&one, &two, &three}) {
    for (const sol::SweepPoint& p : *points) {
      const double total = p.fractions[0] + p.fractions[1] + p.fractions[2];
      CHECK(std::abs(total - 1.0) <= 1e-9);
    }
  }

  CHECK_THROWS_AS(sol::sweep(solute, {water, ethanol, acetone, water}, 4), sol::SolError);

  // steps below the [2, 64] floor must clamp rather than fail or misbehave.
  CHECK_NOTHROW(sol::sweep(solute, {water, ethanol}, 0));
  const std::vector<sol::SweepPoint> clampedLow = sol::sweep(solute, {water, ethanol}, 1);
  CHECK(clampedLow.size() >= 3);
}

TEST_CASE("describeSolute estimates molar mass, volume and Hansen space from a structure") {
  const core::Molecule ethanol = chem::fromSmiles("CCO");
  const sol::Solute solute = sol::describeSolute(ethanol);
  CHECK(solute.molarMass == doctest::Approx(46.0419).epsilon(0.05));
  CHECK(solute.molarVolume > 0.0);
  CHECK(std::isfinite(solute.hansen.dispersion));
  CHECK(std::isfinite(solute.hansen.polar));
  CHECK(std::isfinite(solute.hansen.hydrogenBond));

  core::Molecule empty;
  CHECK_THROWS_AS(sol::describeSolute(empty), sol::SolError);
}

TEST_CASE("describeSolute raises hydrogen bonding for hydroxyl-bearing solutes") {
  const sol::Solute benzene = sol::describeSolute(chem::fromSmiles("c1ccccc1"));
  const sol::Solute ethanol = sol::describeSolute(chem::fromSmiles("CCO"));
  CHECK(benzene.hansen.hydrogenBond < ethanol.hansen.hydrogenBond);
}

TEST_CASE("shake derives slosh velocity and power input from physical inputs") {
  sol::Phase water;
  water.label = "water";
  water.volumeMl = 100.0;
  water.density = 1.0;
  sol::Phase oil;
  oil.label = "oil";
  oil.volumeMl = 100.0;
  oil.density = 0.8;

  sol::Simulation sim;
  sim.phases = {water, oil};
  sol::reset(sim);

  const sol::ShakeParams params{5.0, 3.0, 0.05};
  sol::shake(sim, params);

  const double expectedU = 2.0 * std::numbers::pi_v<double> * 3.0 * 0.05;
  CHECK(sim.shake.peakVelocity == doctest::Approx(expectedU).epsilon(1e-9));
  CHECK(sim.shake.specificPower ==
        doctest::Approx(0.5 * expectedU * expectedU * 3.0).epsilon(1e-9));
  CHECK(sim.shake.active);
  CHECK(sim.shake.remainingS == doctest::Approx(5.0));
  CHECK(sim.shake.sauterRadiusM > 0.0);

  // A shake is a promise of motion, not instant dispersion: nothing is
  // emulsified until the clock advances.
  CHECK(sol::emulsifiedFraction(sim) == doctest::Approx(0.0));
}

TEST_CASE("harder shaking makes smaller droplets and disperses faster") {
  auto shakenAt = [](double frequencyHz, double amplitudeM, double tension) {
    sol::Phase water;
    water.label = "water";
    water.volumeMl = 100.0;
    water.density = 1.0;
    water.interfacialTension = tension;
    sol::Phase oil;
    oil.label = "oil";
    oil.volumeMl = 100.0;
    oil.density = 0.8;
    oil.interfacialTension = tension;

    sol::Simulation sim;
    sim.seed = 11u;
    sim.phases = {water, oil};
    sol::reset(sim);
    sol::shake(sim, sol::ShakeParams{8.0, frequencyHz, amplitudeM});
    return sim;
  };

  sol::Simulation gentle = shakenAt(1.0, 0.02, 30.0);
  sol::Simulation hard = shakenAt(4.0, 0.08, 30.0);

  // Hinze: d32 falls as epsilon^-0.4, so the harder shake must predict the
  // smaller mean droplet.
  CHECK(hard.shake.sauterRadiusM < gentle.shake.sauterRadiusM);

  // Same starting point: nothing dispersed yet.
  REQUIRE(sol::emulsifiedFraction(gentle) == doctest::Approx(0.0));
  REQUIRE(sol::emulsifiedFraction(hard) == doctest::Approx(0.0));

  // After the same wall time the harder shake has turned the column over
  // more often, so more volume is dispersed.
  for (int i = 0; i < 40; ++i) {  // 2 s
    sol::step(gentle, 0.05);
    sol::step(hard, 0.05);
  }
  CHECK(sol::emulsifiedFraction(hard) > sol::emulsifiedFraction(gentle));

  // Higher interfacial tension resists breakup: larger mean droplets.
  sol::Simulation slippery = shakenAt(3.0, 0.05, 5.0);
  sol::Simulation stiff = shakenAt(3.0, 0.05, 60.0);
  CHECK(stiff.shake.sauterRadiusM > slippery.shake.sauterRadiusM);
}

TEST_CASE("a shake emulsifies progressively over its duration, then stops") {
  sol::Phase water;
  water.label = "water";
  water.volumeMl = 100.0;
  water.density = 1.0;
  water.emulsionStability = 1.0;  // nothing re-settles: isolates dispersion
  sol::Phase oil;
  oil.label = "oil";
  oil.volumeMl = 100.0;
  oil.density = 0.8;
  oil.emulsionStability = 1.0;

  sol::Simulation sim;
  sim.seed = 5u;
  sim.phases = {water, oil};
  sol::reset(sim);
  sol::shake(sim, sol::ShakeParams{4.0, 3.0, 0.05});

  double previous = 0.0;
  for (int i = 0; i < 80; ++i) {  // 4 s == the shake duration
    sol::step(sim, 0.05);
    const double fraction = sol::emulsifiedFraction(sim);
    CAPTURE(i);
    CHECK(fraction >= previous - 1e-12);  // monotonic while shaking
    previous = fraction;
  }
  CHECK(previous > 0.5);  // a 4 s firm shake thoroughly emulsifies
  CHECK(!sim.shake.active);
}

TEST_CASE("funnel simulation conserves total charged volume through shaking and settling") {
  sol::Simulation sim;
  sim.vessel = sol::Vessel::SeparatoryFunnel;
  sim.vesselVolumeMl = 250.0;
  sim.seed = 42u;

  sol::Phase water;
  water.label = "water";
  water.volumeMl = 100.0;
  water.density = 1.0;
  water.viscosity = 1.0;
  water.interfacialTension = 30.0;
  water.emulsionStability = 0.3;

  sol::Phase ether;
  ether.label = "ether";
  ether.volumeMl = 80.0;
  ether.density = 0.71;
  ether.viscosity = 0.5;
  ether.interfacialTension = 15.0;
  ether.emulsionStability = 0.3;

  sim.phases = {water, ether};
  sol::reset(sim);

  const double total = sol::totalVolumeMl(sim);
  CHECK(total == doctest::Approx(180.0).epsilon(1e-9));
  CHECK(std::abs(settledVolumeMl(sim) + dropletVolumeMl(sim) - total) <= 1e-9);

  sol::shake(sim, sol::ShakeParams{4.0, 3.0, 0.06});
  CHECK(std::abs(settledVolumeMl(sim) + dropletVolumeMl(sim) - total) <= 1e-9);

  for (int i = 0; i < 500; ++i) {
    sol::step(sim, 0.05);
    CAPTURE(i);
    CHECK(std::abs(settledVolumeMl(sim) + dropletVolumeMl(sim) - total) <= 1e-9);
  }
}

TEST_CASE("reset sorts phases dense-first regardless of input order") {
  sol::Phase light;
  light.label = "light";
  light.volumeMl = 50.0;
  light.density = 0.6;

  sol::Phase mid;
  mid.label = "mid";
  mid.volumeMl = 50.0;
  mid.density = 0.9;

  sol::Phase dense;
  dense.label = "dense";
  dense.volumeMl = 50.0;
  dense.density = 1.3;

  sol::Simulation sim;
  sim.phases = {light, dense, mid};  // deliberately unsorted
  sol::reset(sim);

  REQUIRE(sim.phases.size() == 3);
  for (size_t i = 0; i + 1 < sim.phases.size(); ++i) {
    CHECK(sim.phases[i].density >= sim.phases[i + 1].density);
  }
  CHECK(sim.phases.front().label == "dense");
  CHECK(sim.phases.back().label == "light");
}

TEST_CASE("emulsion stability governs how long a shaken mixture stays dispersed") {
  auto buildShaken = [](double stability) {
    sol::Phase water;
    water.label = "water";
    water.volumeMl = 100.0;
    water.density = 1.0;
    water.interfacialTension = 20.0;
    water.emulsionStability = stability;

    sol::Phase oil;
    oil.label = "oil";
    oil.volumeMl = 100.0;
    oil.density = 0.85;
    oil.interfacialTension = 20.0;
    oil.emulsionStability = stability;

    sol::Simulation sim;
    sim.seed = 7u;
    sim.phases = {water, oil};
    sol::reset(sim);
    sol::shake(sim, sol::ShakeParams{6.0, 3.5, 0.06});
    // Run the shake to completion so both builds start settling from the
    // same fully-shaken state.
    for (int t = 0; t < 140; ++t) sol::step(sim, 0.05);
    return sim;
  };

  sol::Simulation unstable = buildShaken(0.0);
  sol::Simulation stable = buildShaken(1.0);

  constexpr double kDt = 0.05;
  constexpr int kSteps = 1200;  // 1200 * 0.05 s == 60 s
  for (int i = 0; i < kSteps; ++i) {
    sol::step(unstable, kDt);
    sol::step(stable, kDt);
  }

  CHECK(sol::emulsifiedFraction(unstable) < 0.02);
  CHECK(sol::emulsifiedFraction(stable) > 0.2);
}

TEST_CASE("identical seeds and action sequences reproduce identical outcomes") {
  auto buildCharged = []() {
    sol::Phase water;
    water.label = "water";
    water.volumeMl = 120.0;
    water.density = 1.0;

    sol::Phase oil;
    oil.label = "oil";
    oil.volumeMl = 90.0;
    oil.density = 0.8;

    sol::Simulation sim;
    sim.seed = 99u;
    sim.phases = {water, oil};
    sol::reset(sim);
    return sim;
  };

  sol::Simulation simA = buildCharged();
  sol::Simulation simB = buildCharged();

  sol::shake(simA, sol::ShakeParams{5.0, 3.0, 0.05});
  sol::shake(simB, sol::ShakeParams{5.0, 3.0, 0.05});
  for (int i = 0; i < 200; ++i) {
    sol::step(simA, 0.03);
    sol::step(simB, 0.03);
  }

  CHECK(simA.droplets.size() == simB.droplets.size());
  CHECK(sol::emulsifiedFraction(simA) ==
        doctest::Approx(sol::emulsifiedFraction(simB)).epsilon(1e-12));
}

TEST_CASE("vessel geometry stays bounded, clamps and narrows correctly at the base") {
  const std::array<sol::Vessel, 3> vessels = {
      sol::Vessel::SeparatoryFunnel, sol::Vessel::DecantingFlask, sol::Vessel::GraduatedCylinder};
  const std::array<double, 5> fractions = {0.0, 0.25, 0.5, 0.75, 1.0};

  for (sol::Vessel v : vessels) {
    for (double f : fractions) {
      const double width = sol::vesselWidthAt(v, f);
      CHECK(width > 0.0);
      CHECK(width <= 1.0);
    }

    // Out-of-range fractions must clamp rather than return garbage.
    const double below = sol::vesselWidthAt(v, -5.0);
    const double above = sol::vesselWidthAt(v, 5.0);
    CHECK(std::isfinite(below));
    CHECK(std::isfinite(above));
    CHECK(below > 0.0);
    CHECK(below <= 1.0);
    CHECK(above > 0.0);
    CHECK(above <= 1.0);

    const std::vector<core::Vec2> outline = sol::vesselOutline(v, 0.2);
    CHECK(outline.size() >= 24);
  }

  CHECK(sol::vesselWidthAt(sol::Vessel::SeparatoryFunnel, 0.0) <
        sol::vesselWidthAt(sol::Vessel::SeparatoryFunnel, 0.5));
}

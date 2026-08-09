#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <cmath>
#include <string>

#include "sol/transfer.hpp"

namespace sol = chemcad::sol;

namespace {

sol::TransferInput representativeInput(double areaM2) {
  sol::TransferInput input;
  input.soluteMassMg = 100.0;
  input.aqueousVolumeMl = 100.0;
  input.organicVolumeMl = 100.0;
  input.distributionRatio = 3.0;
  input.interfacialAreaM2 = areaM2;
  input.sauterDiameterM = 1.0e-3;
  input.dispersedFraction = 0.20;
  input.soluteMolarVolumeCm3 = 100.0;
  input.slipVelocityMPerS = 0.06;
  input.contactTimeS = 0.02;
  return input;
}

double analyticOrganicEquilibrium(const sol::TransferInput& input) {
  return input.soluteMassMg * input.distributionRatio * input.organicVolumeMl /
         (input.aqueousVolumeMl +
          input.distributionRatio * input.organicVolumeMl);
}

}  // namespace

TEST_CASE("analytic transfer reaches the exact partition equilibrium") {
  const sol::TransferInput input = representativeInput(0.06);
  const double expected = analyticOrganicEquilibrium(input);

  const sol::TransferState oneLargeStep =
      sol::advanceTransfer(input, {}, 10000.0);
  CHECK(oneLargeStep.organicMassMg ==
        doctest::Approx(expected).epsilon(0.001));
  CHECK(oneLargeStep.equilibriumOrganicMassMg ==
        doctest::Approx(expected).epsilon(1e-12));
  CHECK(oneLargeStep.aqueousMassMg + oneLargeStep.organicMassMg ==
        doctest::Approx(input.soluteMassMg).epsilon(1e-12));

  sol::TransferState manySteps;
  for (int i = 0; i < 2000; ++i) {
    manySteps = sol::advanceTransfer(input, manySteps, 0.25);
  }
  CHECK(manySteps.organicMassMg == doctest::Approx(expected).epsilon(0.001));
  CHECK(manySteps.approachToEquilibrium ==
        doctest::Approx(1.0).epsilon(0.001));
}

TEST_CASE("transfer rate scales with measured interfacial area") {
  const sol::TransferInput smallArea = representativeInput(0.02);
  const sol::TransferInput largeArea = representativeInput(0.04);

  const sol::TransferState slow = sol::advanceTransfer(smallArea, {}, 2.0);
  const sol::TransferState fast = sol::advanceTransfer(largeArea, {}, 2.0);

  REQUIRE(slow.timeConstantS > 0.0);
  REQUIRE(fast.timeConstantS > 0.0);
  CHECK(fast.timeConstantS ==
        doctest::Approx(0.5 * slow.timeConstantS).epsilon(1e-12));
  CHECK(fast.organicMassMg > slow.organicMassMg);
  CHECK(slow.areaMeasured);
  CHECK(fast.areaMeasured);
  CHECK(fast.note.find("measured by the fluid simulation") != std::string::npos);
}

TEST_CASE("zero interfacial area transfers no solute") {
  sol::TransferInput input = representativeInput(0.0);
  input.sauterDiameterM = 0.0;
  input.dispersedFraction = 0.0;

  const sol::TransferState state = sol::advanceTransfer(input, {}, 3600.0);
  CHECK(state.organicMassMg == doctest::Approx(0.0));
  CHECK(state.aqueousMassMg == doctest::Approx(input.soluteMassMg));
  CHECK(state.timeConstantS == doctest::Approx(0.0));
  CHECK_FALSE(state.areaMeasured);
  CHECK(state.note.find("transfer is disabled") != std::string::npos);
}

TEST_CASE("distribution ratio follows Henderson-Hasselbalch ionisation") {
  sol::Solute neutral;
  neutral.name = "naphthalene";
  neutral.logP = 3.30;
  neutral.ionization.ionClass = sol::IonClass::Neutral;
  CHECK(sol::distributionRatio(neutral, 7.0, 25.0) ==
        doctest::Approx(std::pow(10.0, neutral.logP)).epsilon(1e-12));

  sol::Solute base;
  base.name = "lidocaine";
  base.logP = 3.26;
  base.ionization.ionClass = sol::IonClass::Base;
  base.ionization.pKa = 7.90;
  const double baseLowPH = sol::distributionRatio(base, 2.0, 25.0);
  const double baseHighPH = sol::distributionRatio(base, 12.0, 25.0);
  CHECK(baseLowPH < baseHighPH);
  CHECK(baseLowPH < std::pow(10.0, base.logP) * 1.0e-5);

  sol::Solute acid;
  acid.name = "benzoic acid";
  acid.logP = 1.87;
  acid.ionization.ionClass = sol::IonClass::Acid;
  acid.ionization.pKa = 4.20;
  const double acidLowPH = sol::distributionRatio(acid, 1.0, 25.0);
  const double acidHighPH = sol::distributionRatio(acid, 10.0, 25.0);
  // Acid ionisation rises above pKa, so its organic/water D must fall.
  CHECK(acidHighPH < acidLowPH);
  CHECK(acidHighPH < std::pow(10.0, acid.logP) * 1.0e-5);
}

TEST_CASE("specific interfacial area uses the Sauter diameter safely") {
  CHECK(sol::specificArea(0.20, 2.0e-3) ==
        doctest::Approx(6.0 * 0.20 / 2.0e-3).epsilon(1e-12));
  CHECK(sol::specificArea(0.20, 0.0) == doctest::Approx(0.0));

  sol::TransferInput correlated = representativeInput(0.0);
  correlated.dispersedFraction = 0.20;
  correlated.sauterDiameterM = 2.0e-3;
  const sol::TransferState state = sol::advanceTransfer(correlated, {}, 1.0);
  CHECK_FALSE(state.areaMeasured);
  CHECK(state.organicMassMg > 0.0);
  CHECK(state.note.find("estimated, not measured") != std::string::npos);

  sol::TransferInput hinzeFallback = representativeInput(0.0);
  hinzeFallback.sauterDiameterM = 0.0;
  hinzeFallback.hinzeDiameterM = 3.0e-3;
  hinzeFallback.hinzeRosinRammlerShape = 3.0;
  const sol::TransferState hinzeState =
      sol::advanceTransfer(hinzeFallback, {}, 1.0);
  CHECK_FALSE(hinzeState.areaMeasured);
  CHECK(hinzeState.organicMassMg > 0.0);
  CHECK(hinzeState.note.find("converted from the Hinze d95") !=
        std::string::npos);
}

TEST_CASE("Wilke-Chang reproduces acetic acid diffusion in water") {
  // Wilke and Chang, AIChE J. 1 (1955) 264-270, report the standard
  // acetic-acid/water 25 C check near 1.21e-9 m2/s. Using water phi=2.6,
  // mu=0.89 cP and acetic-acid V_b=63.8 cm3/mol should agree within the
  // correlation's normal engineering scatter rather than be fitted exactly.
  const double predicted = sol::diffusivityWilkeChang(
      18.01528, 0.89e-3, 63.8, 298.15, 2.6);
  CHECK(predicted == doctest::Approx(1.21e-9).epsilon(0.25));
}

TEST_CASE("worked caffeine extraction approaches equilibrium through shake and settle") {
  sol::Solute caffeine;
  caffeine.name = "caffeine";
  caffeine.logP = -0.07;
  caffeine.molarVolume = 133.0;
  caffeine.ionization.ionClass = sol::IonClass::Base;
  caffeine.ionization.pKa = 0.60;

  sol::TransferInput input;
  input.soluteMassMg = 100.0;
  input.aqueousVolumeMl = 100.0;
  input.organicVolumeMl = 100.0;
  input.distributionRatio = sol::distributionRatio(caffeine, 7.0, 25.0);
  input.interfacialAreaM2 = 0.08;
  input.sauterDiameterM = 1.0e-3;
  input.dispersedFraction = 0.20;
  input.soluteMolarVolumeCm3 = caffeine.molarVolume;
  input.slipVelocityMPerS = 0.08;
  input.contactTimeS = 0.02;

  sol::TransferState state;
  double previousApproach = 0.0;
  double previousOrganic = 0.0;

  // Five seconds of measured high area while shaking.
  for (int i = 0; i < 50; ++i) {
    state = sol::advanceTransfer(input, state, 0.1);
    CHECK(state.approachToEquilibrium >= previousApproach);
    CHECK(state.organicMassMg >= previousOrganic);
    previousApproach = state.approachToEquilibrium;
    previousOrganic = state.organicMassMg;
  }

  // During a 30 s settle the measured area decays as droplets coalesce; it is
  // still the simulation's area, not the 6 phi/d32 fallback.
  for (int second = 0; second < 30; ++second) {
    input.interfacialAreaM2 = 0.02 * std::exp(-0.08 * second);
    input.sauterDiameterM = 1.0e-3 * (1.0 + 0.05 * second);
    state = sol::advanceTransfer(input, state, 1.0);
    CHECK(state.approachToEquilibrium >= previousApproach);
    CHECK(state.organicMassMg >= previousOrganic);
    previousApproach = state.approachToEquilibrium;
    previousOrganic = state.organicMassMg;
  }

  CHECK(state.organicMassMg > 0.0);
  CHECK(state.organicMassMg <= state.equilibriumOrganicMassMg);
  CHECK(state.approachToEquilibrium > 0.0);
  CHECK(state.areaMeasured);
}

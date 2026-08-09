// Ionisation path: acids, bases, salts, and the Born dielectric penalty.
//
// The neutral Flory-Huggins/Hansen/Yalkowsky chain is only half the story for
// real bench chemistry. These cases pin the half that sol/ionization.cpp adds:
// which site ionises, what a hydrochloride does to aqueous solubility, and why
// that same hydrochloride will not go into chloroform.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <cmath>
#include <string>

#include "chem/bridge.hpp"
#include "core/model.hpp"
#include "sol/ionization.hpp"
#include "sol/solubility.hpp"
#include "sol/solvent.hpp"

namespace chem = chemcad::chem;
namespace core = chemcad::core;
namespace sol = chemcad::sol;

namespace {

// A tropane alkaloid ester and the same structure drawn as its
// hydrochloride: the exact case that used to be modelled as a free base.
constexpr const char* kAlkaloid = "COC(=O)C1C(OC(=O)c2ccccc2)CC2CCC1N2C";
constexpr const char* kAlkaloidHcl = "COC(=O)C1C(OC(=O)c2ccccc2)CC2CCC1N2C.Cl";

sol::Solute soluteOf(const std::string& smiles) {
  return sol::describeSolute(chem::fromSmiles(smiles));
}

double predictIn(const sol::Solute& solute, const char* solventId,
                 double pH = sol::kAutoPH) {
  const sol::Solvent* solvent = sol::findSolvent(solventId);
  REQUIRE(solvent != nullptr);
  return sol::predict(solute, {{solvent, 1.0}}, 25.0, nullptr, 0.0, pH).gramsPerMillilitre;
}

}  // namespace

TEST_CASE("a drawn salt keeps the skeleton's descriptors and records its counter-ion") {
  const sol::Solute freeBase = soluteOf(kAlkaloid);
  const sol::Solute salt = soluteOf(kAlkaloidHcl);

  // The counter-ion must not pollute the skeleton descriptors: same molecule,
  // same molar mass and molar volume as the free base.
  CHECK(salt.molarMass == doctest::Approx(freeBase.molarMass).epsilon(1e-9));
  CHECK(salt.molarVolume == doctest::Approx(freeBase.molarVolume).epsilon(1e-9));

  CHECK(salt.ionization.fragmentCount == 2);
  CHECK(salt.ionization.saltForm);
  CHECK(salt.ionization.counterIon.find("Cl") != std::string::npos);
  CHECK(salt.ionization.ionClass == sol::IonClass::Base);
  CHECK(salt.ionization.pKa > 8.0);
  CHECK(salt.ionization.pKa < 12.0);

  // Drawn without the counter-ion it is still a base, but not a salt.
  CHECK(freeBase.ionization.ionClass == sol::IonClass::Base);
  CHECK_FALSE(freeBase.ionization.saltForm);
  CHECK(freeBase.ionization.fragmentCount == 1);
}

TEST_CASE("an amine hydrochloride is far more water soluble than its free base") {
  const sol::Solute freeBase = soluteOf(kAlkaloid);
  const sol::Solute salt = soluteOf(kAlkaloidHcl);

  const double baseWater = predictIn(freeBase, "water");
  const double saltWater = predictIn(salt, "water");
  REQUIRE(baseWater > 0.0);
  REQUIRE(saltWater > 0.0);

  // Salt formation is worth orders of magnitude in water; the measured value
  // for this compound class is a few hundred mg/mL upward.
  CHECK(std::log10(saltWater / baseWater) > 1.5);
  CHECK(saltWater > 0.1);

  // ... and it can never exceed the crystal-density mass-balance bound.
  const double ceiling = salt.molarMass / salt.molarVolume;
  CHECK(saltWater <= ceiling);

  const sol::Solvent* water = sol::findSolvent("water");
  REQUIRE(water != nullptr);
  const sol::Prediction p = sol::predict(salt, {{water, 1.0}});
  CHECK(p.ionicPath);
  CHECK(p.pHSelfBuffered);
  CHECK(p.pH < 7.0);            // BH+ hydrolyses acidic
  CHECK(p.ionisedFraction > 0.9);
  CHECK(p.bornPenaltyDecades == doctest::Approx(0.0));  // water is the reference
  CHECK_FALSE(p.ionNote.empty());
}

TEST_CASE("a hydrochloride dissolves as an ion pair in a low-dielectric solvent") {
  const sol::Solute freeBase = soluteOf(kAlkaloid);
  const sol::Solute salt = soluteOf(kAlkaloidHcl);

  const double baseChloroform = predictIn(freeBase, "chloroform");
  const double saltChloroform = predictIn(salt, "chloroform");
  const double saltHexane = predictIn(salt, "hexane");
  REQUIRE(baseChloroform > 0.0);
  REQUIRE(saltChloroform > 0.0);

  // Bench reality for an amine hydrochloride: soluble in chloroform, but
  // roughly a decade below its free base -- not the twelve decades that free
  // -ion Born electrostatics alone would predict, because below eps ~ 15 the
  // species crossing into solution is a contact ion pair.
  const double gap = std::log10(baseChloroform / saltChloroform);
  CHECK(gap > 0.4);
  CHECK(gap < 2.5);

  // A hydrocarbon solvates neither ions nor the pair.
  CHECK(saltHexane < saltChloroform * 0.05);

  const sol::Solvent* chloroform = sol::findSolvent("chloroform");
  REQUIRE(chloroform != nullptr);
  const sol::Prediction p = sol::predict(salt, {{chloroform, 1.0}});
  CHECK(p.ionicPath);
  CHECK(p.bornPenaltyDecades > 3.0);

  // Monotonic in dielectric constant: the more polar the medium, the smaller
  // the penalty.
  const double rSolute = sol::bornRadiusFromMolarVolumeNm(salt.molarVolume);
  const double inChloroform = sol::bornPenaltyDecades(4.8, rSolute, 0.181, 298.15);
  const double inEthanol = sol::bornPenaltyDecades(24.3, rSolute, 0.181, 298.15);
  const double inWater = sol::bornPenaltyDecades(78.4, rSolute, 0.181, 298.15);
  CHECK(inChloroform > inEthanol);
  CHECK(inEthanol > inWater);
  CHECK(inWater == doctest::Approx(0.0));
}

TEST_CASE("pH moves an acid and a base in opposite directions") {
  const sol::Solute benzoicAcid = soluteOf("OC(=O)c1ccccc1");
  REQUIRE(benzoicAcid.ionization.ionClass == sol::IonClass::Acid);
  CHECK(benzoicAcid.ionization.pKa < 5.0);

  const double acidLowPH = predictIn(benzoicAcid, "water", 2.0);
  const double acidHighPH = predictIn(benzoicAcid, "water", 10.0);
  CHECK(acidHighPH > acidLowPH * 10.0);  // the carboxylate is far more soluble

  const sol::Solute amine = soluteOf("CCCCCCN");
  REQUIRE(amine.ionization.ionClass == sol::IonClass::Base);
  const double baseLowPH = predictIn(amine, "water", 2.0);
  const double baseHighPH = predictIn(amine, "water", 11.0);
  CHECK(baseLowPH > baseHighPH);
}

TEST_CASE("ring carbonyls strip basicity so a xanthine is not treated as a base") {
  // Caffeine's imidazole nitrogen looks like a pKa-7 azole on ring type
  // alone; the two xanthine carbonyls pull its real pKa to 0.6. The model
  // must not hand caffeine an aqueous ionisation bonus at pH 7.
  const sol::Solute caffeine = soluteOf("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
  if (caffeine.ionization.ionClass == sol::IonClass::Base) {
    CHECK(caffeine.ionization.pKa < 5.0);
  }
  const sol::Solvent* water = sol::findSolvent("water");
  REQUIRE(water != nullptr);
  const sol::Prediction p = sol::predict(caffeine, {{water, 1.0}});
  CHECK(p.gramsPerMillilitre == doctest::Approx(0.0216).epsilon(0.02));  // still the anchor
}

TEST_CASE("neutral and zwitterionic solutes are untouched by the ionisation path") {
  const sol::Solute naphthalene = soluteOf("c1ccc2ccccc2c1");
  CHECK(naphthalene.ionization.ionClass == sol::IonClass::Neutral);
  const sol::Solvent* toluene = sol::findSolvent("toluene");
  REQUIRE(toluene != nullptr);
  const sol::Prediction neutral = sol::predict(naphthalene, {{toluene, 1.0}});
  CHECK_FALSE(neutral.ionicPath);
  // An explicit pH cannot move a species with nothing to ionise.
  const sol::Prediction forced =
      sol::predict(naphthalene, {{toluene, 1.0}}, 25.0, nullptr, 0.0, 1.0);
  CHECK(forced.gramsPerMillilitre == doctest::Approx(neutral.gramsPerMillilitre));

  // Glycine is a zwitterion: both sites strong, net neutral, no enhancement.
  const sol::Solute glycine = soluteOf("NCC(=O)O");
  CHECK(glycine.ionization.ionClass == sol::IonClass::Zwitterion);
  const sol::Solvent* water = sol::findSolvent("water");
  REQUIRE(water != nullptr);
  CHECK_FALSE(sol::predict(glycine, {{water, 1.0}}).ionicPath);
}

TEST_CASE("self-buffered pH follows the species that is dissolving") {
  sol::Ionization base;
  base.ionClass = sol::IonClass::Base;
  base.pKa = 9.8;
  CHECK(sol::selfBufferedPH(base, 0.01) > 7.0);  // free base: basic solution
  base.saltForm = true;
  CHECK(sol::selfBufferedPH(base, 1.0) < 7.0);   // its hydrochloride: acidic

  sol::Ionization acid;
  acid.ionClass = sol::IonClass::Acid;
  acid.pKa = 4.2;
  CHECK(sol::selfBufferedPH(acid, 0.03) < 7.0);  // free acid: acidic
  acid.saltForm = true;
  CHECK(sol::selfBufferedPH(acid, 1.0) > 7.0);   // its sodium salt: basic

  sol::Ionization neutral;
  CHECK(sol::selfBufferedPH(neutral, 0.1) == doctest::Approx(7.0));
  CHECK(sol::ionisedRatio(neutral, 3.0) == doctest::Approx(0.0));
}

TEST_CASE("component splitting orders fragments largest first and copies bonds") {
  const core::Molecule structure = chem::fromSmiles(kAlkaloidHcl);
  const std::vector<core::Molecule> parts = sol::splitComponents(structure);
  REQUIRE(parts.size() == 2);
  CHECK(parts[0].atomCount() > parts[1].atomCount());
  // The skeleton keeps every bond it had; the chloride has none.
  CHECK(parts[0].bondCount() >= parts[0].atomCount() - 1);
  CHECK(parts[1].bondCount() == 0);

  size_t total = 0;
  for (const core::Molecule& part : parts) total += part.atomCount();
  CHECK(total == structure.atomCount());
}

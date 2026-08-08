// Tests for the literature-anchor and aqueous-ionic solubility paths:
// measured values must win over the model where they exist, and the Ksp
// machinery must reproduce the common-ion effect.
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <cmath>
#include <string>
#include <vector>

#include "chem/bridge.hpp"
#include "sol/anchors.hpp"
#include "sol/solubility.hpp"
#include "sol/solvent.hpp"

namespace chem = chemcad::chem;
namespace sol = chemcad::sol;

namespace {

sol::Solute soluteOf(const char* smiles) {
  return sol::describeSolute(chem::fromSmiles(smiles));
}

sol::Prediction pureWater(const sol::Solute& solute, const sol::Electrolyte* bg = nullptr,
                          double bgM = 0.0) {
  return sol::predict(solute, {sol::Component{sol::findSolvent("water"), 1.0}}, 25.0, bg, bgM);
}

}  // namespace

TEST_CASE("the anchor database loads and keys match however the solute was drawn") {
  REQUIRE(sol::anchors().size() >= 14);
  // Two different valid SMILES for caffeine must hit the same anchor.
  const std::string a = chem::canonicalize("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
  const std::string b = chem::canonicalize("CN1C=NC2=C1C(=O)N(C)C(=O)N2C");
  REQUIRE(sol::findAnchor(a, "water", 25.0) != nullptr);
  CHECK(sol::findAnchor(b, "water", 25.0) == sol::findAnchor(a, "water", 25.0));
  CHECK(sol::findAnchor(a, "water", 25.0)->gramsPerMillilitre ==
        doctest::Approx(0.0216));
  // Far outside the temperature tolerance there is no anchor.
  CHECK(sol::findAnchor(a, "water", 120.0) == nullptr);
}

TEST_CASE("an anchored pure solvent returns the measured value exactly") {
  const sol::Solute caffeine = soluteOf("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
  const sol::Prediction water = pureWater(caffeine);
  CHECK(water.anchored);
  CHECK(!water.saltPath);
  CHECK(water.gramsPerMillilitre == doctest::Approx(0.0216).epsilon(0.01));

  const sol::Prediction ethanol = sol::predict(
      caffeine, {sol::Component{sol::findSolvent("ethanol"), 1.0}}, 25.0);
  CHECK(ethanol.gramsPerMillilitre == doctest::Approx(0.0075).epsilon(0.01));
}

TEST_CASE("an anchored blend interpolates between the measured endpoints") {
  const sol::Solute caffeine = soluteOf("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
  const sol::Prediction mid = sol::predict(caffeine,
                                           {sol::Component{sol::findSolvent("water"), 1.0},
                                            sol::Component{sol::findSolvent("ethanol"), 1.0}},
                                           25.0);
  CHECK(mid.anchored);
  // Geometric interpolation: the midpoint sits between the two measured
  // pure values, not above or below them.
  CHECK(mid.gramsPerMillilitre < 0.0216);
  CHECK(mid.gramsPerMillilitre > 0.0075);
}

TEST_CASE("an unanchored pair keeps the pure model and reports it") {
  const sol::Solute caffeine = soluteOf("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
  const sol::Prediction heptane = sol::predict(
      caffeine, {sol::Component{sol::findSolvent("heptane"), 1.0}}, 25.0);
  CHECK(!heptane.anchored);
}

TEST_CASE("sodium chloride uses the measured Ksp path in water") {
  const sol::Solute nacl = soluteOf("[Na+].[Cl-]");
  REQUIRE(!nacl.canonicalSmiles.empty());
  REQUIRE(sol::findSalt(nacl.canonicalSmiles) != nullptr);

  const sol::Prediction prediction = pureWater(nacl);
  CHECK(prediction.saltPath);
  CHECK(prediction.anchored);
  // 359 g/L at 25 C, exact by construction of the measured endpoint.
  CHECK(prediction.gramsPerMillilitre == doctest::Approx(0.359).epsilon(0.02));
}

TEST_CASE("the common-ion effect depresses salt solubility") {
  const sol::Solute nacl = soluteOf("[Na+].[Cl-]");
  const sol::Electrolyte* kcl = sol::findElectrolyte("kcl");
  const sol::Electrolyte* kno3 = sol::findElectrolyte("kno3");
  REQUIRE(kcl != nullptr);
  REQUIRE(kno3 != nullptr);

  const double pure = pureWater(nacl).gramsPerMillilitre;
  // KCl shares the chloride: less NaCl dissolves.
  const double commonIon = pureWater(nacl, kcl, 1.0).gramsPerMillilitre;
  CHECK(commonIon < pure);
  // KNO3 shares nothing: ionic strength alone slightly increases it.
  const double inertIon = pureWater(nacl, kno3, 1.0).gramsPerMillilitre;
  CHECK(inertIon >= pure);
}

TEST_CASE("endothermic salts dissolve better hot (van't Hoff)") {
  const sol::Solute kno3 = soluteOf("[K+].[O-][N+](=O)[O-]");
  const double hot = sol::predict(kno3, {sol::Component{sol::findSolvent("water"), 1.0}}, 60.0)
                         .gramsPerMillilitre;
  const double cold = sol::predict(kno3, {sol::Component{sol::findSolvent("water"), 1.0}}, 5.0)
                          .gramsPerMillilitre;
  CHECK(hot > cold);
}

TEST_CASE("a salt is effectively insoluble without water") {
  const sol::Solute nacl = soluteOf("[Na+].[Cl-]");
  const sol::Prediction hexane = sol::predict(
      nacl, {sol::Component{sol::findSolvent("hexane"), 1.0}}, 25.0);
  CHECK(hexane.saltPath);
  CHECK(hexane.gramsPerMillilitre == doctest::Approx(0.0).epsilon(1e-9));
}

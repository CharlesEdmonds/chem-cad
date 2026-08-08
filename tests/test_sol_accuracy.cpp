#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "chem/bridge.hpp"
#include "core/model.hpp"
#include "sol/solubility.hpp"
#include "sol/solvent.hpp"

namespace chem = chemcad::chem;
namespace core = chemcad::core;
namespace sol = chemcad::sol;

namespace {

// Builds a solute from a SMILES structure via the real group-contribution
// pipeline, then swaps in the literature melting point so accuracy cases
// isolate the Flory-Huggins/Hansen solubility model from the Joback melting
// point estimate (that estimate is validated on its own in case 4).
sol::Solute soluteFromLiterature(const std::string& smiles, double literatureMeltingPointC) {
  const core::Molecule molecule = chem::fromSmiles(smiles);
  sol::Solute solute = sol::describeSolute(molecule);
  solute.meltingPoint = literatureMeltingPointC;
  return solute;
}

struct LiteratureCase {
  const char* soluteName;
  const char* smiles;
  double meltingPointC;
  const char* solventId;
  double literatureGramsPerMl;
};

// Real solubility data points spanning polar, protic and aprotic solvents so
// the order-of-magnitude check exercises the whole Hansen space, not just
// the easy cases.
const LiteratureCase kLiteratureCases[] = {
    {"caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", 235.0, "water", 2.16e-2},
    {"caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", 235.0, "ethanol", 7.5e-3},
    {"caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", 235.0, "chloroform", 1.8e-1},
    {"benzoic acid", "OC(=O)c1ccccc1", 122.0, "water", 3.4e-3},
    {"benzoic acid", "OC(=O)c1ccccc1", 122.0, "ethanol", 4.3e-1},
    {"benzoic acid", "OC(=O)c1ccccc1", 122.0, "toluene", 1.0e-1},
    {"naphthalene", "c1ccc2ccccc2c1", 80.0, "ethanol", 1.0e-1},
    {"naphthalene", "c1ccc2ccccc2c1", 80.0, "toluene", 3.2e-1},
    {"paracetamol", "CC(=O)Nc1ccc(O)cc1", 169.0, "water", 1.4e-2},
    {"paracetamol", "CC(=O)Nc1ccc(O)cc1", 169.0, "ethanol", 1.7e-1},
};

}  // namespace

// 1. Order-of-magnitude accuracy against real literature solubilities.
TEST_CASE("predicted solubility is within one and a half decades of literature") {
  for (const LiteratureCase& c : kLiteratureCases) {
    const sol::Solvent* solvent = sol::findSolvent(c.solventId);
    REQUIRE(solvent != nullptr);

    const sol::Solute solute = soluteFromLiterature(c.smiles, c.meltingPointC);
    const sol::Prediction prediction = sol::predict(solute, {{solvent, 1.0}});
    CHECK(prediction.converged);

    const double predictedLog = std::log10(std::max(prediction.gramsPerMillilitre, 1e-300));
    const double literatureLog = std::log10(c.literatureGramsPerMl);
    const double logDiff = std::abs(predictedLog - literatureLog);

    std::cout << c.soluteName << " in " << c.solventId << ": predicted "
              << prediction.gramsPerMillilitre << " g/mL vs literature "
              << c.literatureGramsPerMl << " g/mL (log10 diff " << logDiff << ")\n";
    CAPTURE(c.soluteName);
    CAPTURE(c.solventId);
    CHECK(logDiff <= 1.5);
  }
}

// 2. Ranking within a solute is a stronger, strictly-ordered claim.
TEST_CASE("solubility ranking across solvents matches known polarity trends") {
  const sol::Solvent* water = sol::findSolvent("water");
  const sol::Solvent* ethanol = sol::findSolvent("ethanol");
  const sol::Solvent* toluene = sol::findSolvent("toluene");
  REQUIRE(water != nullptr);
  REQUIRE(ethanol != nullptr);
  REQUIRE(toluene != nullptr);

  const sol::Solute benzoicAcid = soluteFromLiterature("OC(=O)c1ccccc1", 122.0);
  const double baWater = sol::predict(benzoicAcid, {{water, 1.0}}).gramsPerMillilitre;
  const double baEthanol = sol::predict(benzoicAcid, {{ethanol, 1.0}}).gramsPerMillilitre;
  const double baToluene = sol::predict(benzoicAcid, {{toluene, 1.0}}).gramsPerMillilitre;
  CHECK(baEthanol > baToluene);
  CHECK(baToluene > baWater);

  const sol::Solute naphthalene = soluteFromLiterature("c1ccc2ccccc2c1", 80.0);
  const double nWater = sol::predict(naphthalene, {{water, 1.0}}).gramsPerMillilitre;
  const double nEthanol = sol::predict(naphthalene, {{ethanol, 1.0}}).gramsPerMillilitre;
  const double nToluene = sol::predict(naphthalene, {{toluene, 1.0}}).gramsPerMillilitre;
  CHECK(nToluene > nEthanol);
  CHECK(nEthanol > nWater);
}

// 3. Melting point must drive both the magnitude and the liquid/crystalline
// transition of the ideal solubility term.
TEST_CASE("melting point drives solubility magnitude and liquid transition") {
  const sol::Solvent* ethanol = sol::findSolvent("ethanol");
  REQUIRE(ethanol != nullptr);

  const sol::Solute hot = soluteFromLiterature("OC(=O)c1ccccc1", 250.0);
  const sol::Solute cold = soluteFromLiterature("OC(=O)c1ccccc1", 100.0);
  const sol::Prediction hotPrediction = sol::predict(hot, {{ethanol, 1.0}});
  const sol::Prediction coldPrediction = sol::predict(cold, {{ethanol, 1.0}});
  CHECK(hotPrediction.gramsPerMillilitre < coldPrediction.gramsPerMillilitre);

  // 20 C is below the default 25 C working temperature: liquid, no melting penalty.
  const sol::Solute liquid = soluteFromLiterature("OC(=O)c1ccccc1", 20.0);
  const sol::Prediction liquidPrediction = sol::predict(liquid, {{ethanol, 1.0}});
  CHECK(liquidPrediction.idealMoleFraction == doctest::Approx(1.0).epsilon(1e-9));
}

// 4. Joback group-contribution melting point is coarse; 80 C is the honest bar.
TEST_CASE("Joback melting point estimate lands within eighty degrees of literature") {
  struct JobackCase {
    const char* smiles;
    double literatureMeltingPointC;
  };
  const JobackCase cases[] = {
      {"OC(=O)c1ccccc1", 122.0},
      {"c1ccc2ccccc2c1", 80.0},
      {"Cn1cnc2c1c(=O)n(C)c(=O)n2C", 235.0},
      {"CC(=O)Nc1ccc(O)cc1", 169.0},
  };
  for (const JobackCase& c : cases) {
    const sol::Solute solute = sol::describeSolute(chem::fromSmiles(c.smiles));
    CAPTURE(c.smiles);
    CHECK(solute.meltingPointEstimated);
    CHECK(std::abs(solute.meltingPoint - c.literatureMeltingPointC) <= 80.0);
  }
}

// 5. Co-solvency maximum: solubility must rise then fall across a binary
// blend, peaking strictly inside the interval, never at an endpoint.
TEST_CASE("co-solvency maximum appears at an interior blend ratio and is unimodal") {
  const sol::Solvent* water = sol::findSolvent("water");
  const sol::Solvent* toluene = sol::findSolvent("toluene");
  REQUIRE(water != nullptr);
  REQUIRE(toluene != nullptr);

  const sol::Solute paracetamol = soluteFromLiterature("CC(=O)Nc1ccc(O)cc1", 169.0);

  // blend() interpolates Hansen parameters linearly in volume fraction, so
  // Ra^2(t) along the water-toluene line is a weighted quadratic in t (the
  // dispersion axis carries the model's 4x weight). Its unconstrained
  // minimiser tStar must fall strictly inside (0, 1) for the solute to be
  // bracketed -- otherwise Ra, and hence solubility, would be monotone
  // across the whole blend and no interior peak could exist.
  const sol::Hansen& a = water->hansen;
  const sol::Hansen& b = toluene->hansen;
  const sol::Hansen& s = paracetamol.hansen;
  const double dDelta = b.dispersion - a.dispersion;
  const double pDelta = b.polar - a.polar;
  const double hDelta = b.hydrogenBond - a.hydrogenBond;
  const double denom = 4.0 * dDelta * dDelta + pDelta * pDelta + hDelta * hDelta;
  REQUIRE(denom > 0.0);
  const double numer = 4.0 * dDelta * (s.dispersion - a.dispersion) +
                        pDelta * (s.polar - a.polar) + hDelta * (s.hydrogenBond - a.hydrogenBond);
  const double tStar = numer / denom;
  REQUIRE(tStar > 0.0);
  REQUIRE(tStar < 1.0);

  const std::vector<sol::SweepPoint> points = sol::sweep(paracetamol, {water, toluene}, 40);
  REQUIRE(points.size() == 41);

  std::size_t peak = 0;
  for (std::size_t i = 1; i < points.size(); ++i) {
    if (points[i].prediction.gramsPerMillilitre > points[peak].prediction.gramsPerMillilitre) {
      peak = i;
    }
  }
  CHECK(peak > 0);
  CHECK(peak < points.size() - 1);

  const double peakValue = points[peak].prediction.gramsPerMillilitre;
  const double endpointA = points.front().prediction.gramsPerMillilitre;
  const double endpointB = points.back().prediction.gramsPerMillilitre;
  CHECK(peakValue > endpointA * 1.05);
  CHECK(peakValue > endpointB * 1.05);

  constexpr double kSlack = 1e-12;
  for (std::size_t i = 1; i <= peak; ++i) {
    CAPTURE(i);
    CHECK(points[i].prediction.gramsPerMillilitre >=
          points[i - 1].prediction.gramsPerMillilitre - kSlack);
  }
  for (std::size_t i = peak + 1; i < points.size(); ++i) {
    CAPTURE(i);
    CHECK(points[i].prediction.gramsPerMillilitre <=
          points[i - 1].prediction.gramsPerMillilitre + kSlack);
  }
}

// 6. predict() must be internally consistent with the concentration it reports.
TEST_CASE("predict is self-consistent with its own reported mole fraction") {
  const sol::Solvent* water = sol::findSolvent("water");
  const sol::Solvent* ethanol = sol::findSolvent("ethanol");
  const sol::Solvent* chloroform = sol::findSolvent("chloroform");
  const sol::Solvent* toluene = sol::findSolvent("toluene");
  REQUIRE(water != nullptr);
  REQUIRE(ethanol != nullptr);
  REQUIRE(chloroform != nullptr);
  REQUIRE(toluene != nullptr);

  struct Case {
    sol::Solute solute;
    const sol::Solvent* solvent;
  };
  const std::vector<Case> cases = {
      {soluteFromLiterature("Cn1cnc2c1c(=O)n(C)c(=O)n2C", 235.0), water},
      {soluteFromLiterature("OC(=O)c1ccccc1", 122.0), ethanol},
      {soluteFromLiterature("c1ccc2ccccc2c1", 80.0), toluene},
      {soluteFromLiterature("CC(=O)Nc1ccc(O)cc1", 169.0), chloroform},
  };

  for (const Case& c : cases) {
    const sol::Prediction prediction = sol::predict(c.solute, {{c.solvent, 1.0}});
    CHECK(prediction.converged);

    // Pure solvent, so the mixture molar volume is exactly the solvent's own.
    const double x = prediction.moleFraction;
    const double solutionMolarVolume =
        x * c.solute.molarVolume + (1.0 - x) * c.solvent->molarVolume;
    const double reconstructed = x * c.solute.molarMass / solutionMolarVolume;

    CAPTURE(c.solvent->id);
    CHECK(reconstructed == doctest::Approx(prediction.gramsPerMillilitre).epsilon(1e-9));
  }
}

// 7. Solubility of a crystalline solute must increase monotonically with temperature.
TEST_CASE("solubility rises monotonically with temperature for a crystalline solute") {
  const sol::Solvent* ethanol = sol::findSolvent("ethanol");
  REQUIRE(ethanol != nullptr);
  const sol::Solute benzoicAcid = soluteFromLiterature("OC(=O)c1ccccc1", 122.0);

  double previous = -1.0;
  for (double tC = 5.0; tC <= 60.0 + 1e-9; tC += 5.0) {
    const sol::Prediction prediction = sol::predict(benzoicAcid, {{ethanol, 1.0}}, tC);
    CAPTURE(tC);
    CHECK(prediction.converged);
    CHECK(prediction.gramsPerMillilitre > previous);
    previous = prediction.gramsPerMillilitre;
  }
}

// 8. Degenerate inputs must stay finite and must never throw.
TEST_CASE("degenerate inputs stay finite and never throw") {
  const sol::Solvent* water = sol::findSolvent("water");
  REQUIRE(water != nullptr);
  const sol::Solute solute = soluteFromLiterature("OC(=O)c1ccccc1", 122.0);

  sol::Prediction emptyPrediction;
  CHECK_NOTHROW(emptyPrediction = sol::predict(solute, {}));
  CHECK(std::isfinite(emptyPrediction.gramsPerMillilitre));
  CHECK(std::isfinite(emptyPrediction.moleFraction));
  CHECK(std::isfinite(emptyPrediction.ra));
  CHECK(std::isfinite(emptyPrediction.activityCoefficient));

  sol::Prediction nullSolventPrediction;
  CHECK_NOTHROW(nullSolventPrediction = sol::predict(solute, {{nullptr, 1.0}}));
  CHECK(std::isfinite(nullSolventPrediction.gramsPerMillilitre));
  CHECK(std::isfinite(nullSolventPrediction.moleFraction));

  sol::Solute zeroVolumeSolute = solute;
  zeroVolumeSolute.molarVolume = 0.0;
  sol::Prediction zeroVolumePrediction;
  CHECK_NOTHROW(zeroVolumePrediction = sol::predict(zeroVolumeSolute, {{water, 1.0}}));
  CHECK(std::isfinite(zeroVolumePrediction.gramsPerMillilitre));
  CHECK(std::isfinite(zeroVolumePrediction.moleFraction));

  sol::Prediction coldPrediction;
  CHECK_NOTHROW(coldPrediction = sol::predict(solute, {{water, 1.0}}, -273.0));
  CHECK(std::isfinite(coldPrediction.gramsPerMillilitre));
  CHECK(std::isfinite(coldPrediction.moleFraction));
}

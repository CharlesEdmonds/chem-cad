// Tests for the solvent screen, miscibility rule and liquid-liquid partition
// helper added to the sol layer for the Solubility Suite overhaul.
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <cmath>
#include <string>

#include "chem/bridge.hpp"
#include "sol/solubility.hpp"
#include "sol/solvent.hpp"

namespace chem = chemcad::chem;
namespace sol = chemcad::sol;

namespace {

sol::Solute benzoicAcid() {
  return sol::describeSolute(chem::fromSmiles("OC(=O)c1ccccc1"));
}

sol::Solute caffeine() {
  return sol::describeSolute(chem::fromSmiles("Cn1c(=O)c2c(ncn2C)n(C)c1=O"));
}

const sol::ScreenRow* rowFor(const std::vector<sol::ScreenRow>& rows, const std::string& id) {
  for (const sol::ScreenRow& row : rows) {
    if (row.solvent->id == id) return &row;
  }
  return nullptr;
}

}  // namespace

TEST_CASE("screen covers the whole database, best solvent first") {
  const std::vector<sol::ScreenRow> rows = sol::screen(benzoicAcid());
  REQUIRE(rows.size() == sol::solvents().size());
  for (size_t i = 1; i < rows.size(); ++i) {
    CHECK(rows[i - 1].prediction.gramsPerMillilitre >=
          rows[i].prediction.gramsPerMillilitre);
  }
  // Every row points at a live database entry.
  for (const sol::ScreenRow& row : rows) {
    CHECK(sol::findSolvent(row.solvent->id) == row.solvent);
  }
}

TEST_CASE("screen agrees with a direct prediction for the same solvent") {
  const sol::Solute solute = benzoicAcid();
  const std::vector<sol::ScreenRow> rows = sol::screen(solute, 40.0);
  const sol::ScreenRow* water = rowFor(rows, "water");
  REQUIRE(water != nullptr);
  const sol::Prediction direct =
      sol::predict(solute, {sol::Component{sol::findSolvent("water"), 1.0}}, 40.0);
  CHECK(water->prediction.gramsPerMillilitre ==
        doctest::Approx(direct.gramsPerMillilitre).epsilon(1e-9));
}

TEST_CASE("polar solutes screen better in water than in hexane") {
  // Caffeine (H-bond acceptor) and glycine (zwitterion) are both essentially
  // insoluble in alkanes; the screen must put water ahead of hexane for
  // both. This guards the asymmetric H-bond penalty in chi -- a regression
  // to a symmetric or zero C3 term flips these rankings (glycine once came
  // out MORE soluble in hexane than in water).
  const std::vector<sol::ScreenRow> cafRows = sol::screen(caffeine());
  const sol::ScreenRow* cafWater = rowFor(cafRows, "water");
  const sol::ScreenRow* cafHexane = rowFor(cafRows, "hexane");
  REQUIRE(cafWater != nullptr);
  REQUIRE(cafHexane != nullptr);
  CHECK(cafWater->prediction.gramsPerMillilitre > cafHexane->prediction.gramsPerMillilitre);

  const std::vector<sol::ScreenRow> glyRows =
      sol::screen(sol::describeSolute(chem::fromSmiles("NCC(=O)O")));
  const sol::ScreenRow* glyWater = rowFor(glyRows, "water");
  const sol::ScreenRow* glyHexane = rowFor(glyRows, "hexane");
  REQUIRE(glyWater != nullptr);
  REQUIRE(glyHexane != nullptr);
  CHECK(glyWater->prediction.gramsPerMillilitre >
        5.0 * glyHexane->prediction.gramsPerMillilitre);

  // The hydrophobic end of the rule: naphthalene must strongly prefer
  // ethanol over water.
  const std::vector<sol::ScreenRow> naphRows =
      sol::screen(sol::describeSolute(chem::fromSmiles("c1ccc2ccccc2c1")));
  const sol::ScreenRow* naphWater = rowFor(naphRows, "water");
  const sol::ScreenRow* naphEthanol = rowFor(naphRows, "ethanol");
  REQUIRE(naphWater != nullptr);
  REQUIRE(naphEthanol != nullptr);
  CHECK(naphEthanol->prediction.gramsPerMillilitre >
        5.0 * naphWater->prediction.gramsPerMillilitre);
}

TEST_CASE("partition conserves mass across logP and volume ratios") {
  const double massMg = 250.0;
  const double logPs[] = {-1.0, 0.5, 4.0};
  const double ratios[][2] = {{50.0, 50.0}, {80.0, 20.0}, {10.0, 120.0}};
  for (double logP : logPs) {
    for (const auto& v : ratios) {
      const sol::Partition p = sol::partition(massMg, logP, v[0], v[1]);
      CHECK(p.mgAqueous + p.mgOrganic == doctest::Approx(massMg).epsilon(1e-9));
      CHECK(p.fractionOrganic ==
            doctest::Approx(p.mgOrganic / massMg).epsilon(1e-9));
    }
  }
}

TEST_CASE("partition direction follows logP") {
  CHECK(sol::partition(100.0, -3.0, 50.0, 50.0).fractionOrganic < 0.01);
  CHECK(sol::partition(100.0, 3.0, 50.0, 50.0).fractionOrganic > 0.99);
  CHECK(sol::partition(100.0, 0.0, 50.0, 50.0).fractionOrganic ==
        doctest::Approx(0.5).epsilon(0.01));
  // Doubling the organic volume pulls more solute across.
  CHECK(sol::partition(100.0, 1.0, 50.0, 100.0).fractionOrganic >
        sol::partition(100.0, 1.0, 100.0, 50.0).fractionOrganic);
}

TEST_CASE("partition handles degenerate input without inventing solute") {
  const sol::Partition zeroMass = sol::partition(0.0, 2.0, 50.0, 50.0);
  CHECK(zeroMass.mgAqueous == 0.0);
  CHECK(zeroMass.mgOrganic == 0.0);
  CHECK(zeroMass.fractionOrganic == 0.0);

  const sol::Partition noOrganic = sol::partition(100.0, 2.0, 50.0, 0.0);
  CHECK(noOrganic.mgAqueous == doctest::Approx(100.0));
  CHECK(noOrganic.mgOrganic == 0.0);
  CHECK(noOrganic.fractionOrganic == 0.0);

  const sol::Partition noPhases = sol::partition(100.0, 2.0, 0.0, 0.0);
  CHECK(noPhases.mgAqueous + noPhases.mgOrganic == 0.0);
}

TEST_CASE("miscibleWith encodes the water-miscibility rule") {
  const sol::Solvent* water = sol::findSolvent("water");
  const sol::Solvent* ethanol = sol::findSolvent("ethanol");
  const sol::Solvent* hexane = sol::findSolvent("hexane");
  REQUIRE(water != nullptr);
  REQUIRE(ethanol != nullptr);
  REQUIRE(hexane != nullptr);
  CHECK(sol::miscibleWith(*water, *ethanol));
  CHECK(sol::miscibleWith(*ethanol, *water));  // symmetric
  CHECK(!sol::miscibleWith(*water, *hexane));
  CHECK(!sol::miscibleWith(*hexane, *water));
  CHECK(sol::miscibleWith(*ethanol, *hexane));
  CHECK(sol::miscibleWith(*water, *water));
}

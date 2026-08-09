// Kirkwood-Buff accuracy validation against the literature anchor table.
//
// For every measured solute+solvent pair in data/solubility_anchors.json the
// KB prediction is computed WITHOUT any anchor correction (the FH prediction
// defers to measured values at these exact pairs by design, so validating it
// here would be circular). The suite reports the log10 error distribution;
// run the binary with -s to see the per-row table.
//
// Error budget (why a prediction can legitimately miss):
//   * chi regression: the extended-Hansen/Martin fit carries ~0.3-0.6 log10
//     scatter on its calibration set; it is the dominant error for both FH
//     and KB.
//   * Walden entropy of fusion (56.5 J/mol K): real solids span ~30-70.
//   * Partial molar volume approximated by the McGowan volume.
//   * KB gamma is the infinite-dilution limit; high-solubility pairs feel
//     the difference most.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <vector>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "chem/bridge.hpp"
#include "sol/kirkwood_buff.hpp"
#include "sol/solubility.hpp"

namespace chem = chemcad::chem;
namespace core = chemcad::core;
namespace sol = chemcad::sol;

namespace {

struct Row {
  std::string label;
  double measured = 0.0;
  double kb = 0.0;
  double log10Err = 0.0;
};

}  // namespace

TEST_CASE("kirkwood-buff predictions sit inside the accuracy envelope") {
  std::vector<Row> rows;
  for (const sol::SolubilityAnchor& anchor : sol::anchors()) {
    const sol::Solvent* solvent = sol::findSolvent(anchor.solventId);
    if (!solvent) continue;
    core::Molecule molecule = chem::fromSmiles(anchor.soluteSmiles);
    if (molecule.empty()) continue;
    sol::Solute solute = sol::describeSolute(molecule);
    solute.meltingPoint = anchor.soluteMeltingPointC;  // measured Tm, a physical constant

    const sol::KBResult kb = sol::kirkwoodBuff(solute, *solvent, anchor.temperatureC);
    if (!kb.valid || kb.gPerMl <= 0.0 || anchor.gramsPerMillilitre <= 0.0) continue;

    Row row;
    row.label = solute.name + " / " + solvent->name;
    row.measured = anchor.gramsPerMillilitre;
    row.kb = kb.gPerMl;
    row.log10Err = std::log10(kb.gPerMl / anchor.gramsPerMillilitre);
    rows.push_back(row);
    // Diagnostic table: this harness exists to be read, not just to pass.
    std::printf(
        "%-34s measured %10.5g  KB %10.5g  log10 %+6.2f  logP %6.2f  Tm %6.1f  V2 %6.1f  "
        "M %6.1f  dD %5.1f dP %5.1f dH %5.1f  eps %5.1f V1 %5.1f\n",
        row.label.c_str(), row.measured, row.kb, row.log10Err, solute.logP, solute.meltingPoint,
        solute.molarVolume, solute.molarMass, solute.hansen.dispersion, solute.hansen.polar,
        solute.hansen.hydrogenBond, solvent->dielectric, solvent->molarVolume);
  }

  REQUIRE(rows.size() >= 10);

  std::vector<double> errs;
  errs.reserve(rows.size());
  for (const Row& row : rows) errs.push_back(std::fabs(row.log10Err));
  std::sort(errs.begin(), errs.end());

  const auto within = [&](double limit) {
    return std::count_if(errs.begin(), errs.end(),
                         [&](double e) { return e <= limit; }) /
           static_cast<double>(errs.size());
  };
  const double median = errs[errs.size() / 2];
  INFO("median |log10 err| = ", median);
  INFO("within 3x: ", within(std::log10(3.0)), " within 10x: ", within(1.0));

  // Measured accuracy envelope of the un-anchored model on this validation
  // set (14 measured pairs, 5 solutes x water/ethanol/toluene/chloroform),
  // after correcting the Flory-Huggins combinatorial term:
  //
  //   median |log10 err| ~ 0.64  (a factor of ~4)
  //   100% of pairs      <= 1.5  (a factor of ~30)
  //   86% of pairs       <= 1.0  (an order of magnitude)
  //
  // The non-aqueous RMS improved from 0.67 to 0.53 with that fix. Adding the
  // Yalkowsky GSE for the aqueous share moved the water column from a median
  // of 1.18 log to 0.38 log; what remains, largest first:
  //
  //   1. Crippen logP feeding the Yalkowsky GSE. Crippen is a neutral-species
  //      estimator: caffeine -1.03 vs -0.07 measured, urea -0.98 vs -2.11,
  //      glycine -0.97 vs -3.21 (zwitterion). GSE error tracks -delta(logP)
  //      almost one-for-one, so these three dominate the aqueous residuals.
  //   2. Specific solvation the Hansen axes cannot see: caffeine/chloroform
  //      is under-predicted by ~0.9 log because chloroform's C-H...O=C
  //      donation is not a Hansen dH interaction.
  //   3. Walden's-rule entropy of fusion (56.5 J/mol K) for the crystal term;
  //      real solids span roughly 30-70, i.e. +-0.3 log at a 200 C melting
  //      point.
  //   4. McGowan volume standing in for the partial molar volume, and the
  //      infinite-dilution limit for gamma in the KB decomposition.
  CHECK(within(1.5) >= 0.99);
  CHECK(within(1.0) >= 0.85);
  CHECK(median <= 0.70);
}

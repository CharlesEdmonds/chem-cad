// Observable contracts for automatic solvent selection and its explanations.
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "chem/bridge.hpp"
#include "sol/selection.hpp"
#include "sol/solubility.hpp"
#include "sol/solvent.hpp"

namespace chem = chemcad::chem;
namespace sol = chemcad::sol;

namespace {

sol::SpeciesRole role(const char* smiles, const char* label, bool keep = true,
                      double weight = 1.0) {
  sol::SpeciesRole result;
  result.solute = sol::describeSolute(chem::fromSmiles(smiles));
  result.label = label;
  result.keep = keep;
  result.weight = weight;
  return result;
}

const sol::SolventCandidate* candidateFor(const std::vector<sol::SolventCandidate>& rows,
                                          const std::string& id) {
  for (const sol::SolventCandidate& row : rows) {
    if (row.solvent && row.solvent->id == id) return &row;
  }
  return nullptr;
}

const sol::Criterion* criterionFor(const sol::SolventCandidate& row, const std::string& name) {
  for (const sol::Criterion& criterion : row.criteria) {
    if (criterion.name == name) return &criterion;
  }
  return nullptr;
}

size_t positionOf(const std::vector<sol::SolventCandidate>& rows, const std::string& id) {
  for (size_t i = 0; i < rows.size(); ++i) {
    if (rows[i].solvent && rows[i].solvent->id == id) return i;
  }
  return rows.size();
}
bool warningContains(const sol::SolventCandidate& row, const std::string& text) {
  return std::any_of(row.warnings.begin(), row.warnings.end(),
                     [&](const std::string& warning) {
                       return warning.find(text) != std::string::npos;
                     });
}

const sol::SolventCandidate* firstDirectionFailure(
    const std::vector<sol::SolventCandidate>& rows) {
  for (const sol::SolventCandidate& row : rows) {
    if (row.selectivity < 1.0) return &row;
  }
  return nullptr;
}

void checkWorkingCandidatesPrecedeFailures(
    const std::vector<sol::SolventCandidate>& rows) {
  bool sawWorking = false;
  bool sawFailure = false;
  size_t lastWorking = 0;
  size_t firstFailure = rows.size();
  double minimumWorkingScore = 1.0;
  double maximumFailureScore = 0.0;
  for (size_t i = 0; i < rows.size(); ++i) {
    if (rows[i].selectivity >= 1.0) {
      sawWorking = true;
      lastWorking = i;
      minimumWorkingScore = std::min(minimumWorkingScore, rows[i].score);
    } else {
      sawFailure = true;
      firstFailure = std::min(firstFailure, i);
      maximumFailureScore = std::max(maximumFailureScore, rows[i].score);
    }
  }
  REQUIRE(sawWorking);
  REQUIRE(sawFailure);
  CHECK(lastWorking < firstFailure);
  CHECK(minimumWorkingScore > maximumFailureScore);
}

sol::OperationSpec extractionSpec() {
  sol::OperationSpec spec;
  spec.kind = sol::OperationKind::LiquidLiquidExtraction;
  spec.species = {role("Cn1c(=O)c2c(ncn2C)n(C)c1=O", "caffeine")};
  spec.weightGreenness = 0.0;
  spec.weightPracticality = 0.0;
  return spec;
}

}  // namespace

TEST_CASE("recrystallisation rewards hot-cold contrast and rejects boiling below hot") {
  sol::OperationSpec spec;
  spec.kind = sol::OperationKind::Recrystallisation;
  spec.species = {role("CC(=O)Nc1ccc(O)cc1", "paracetamol")};
  spec.coldTemperatureC = 20.0;
  spec.hotTemperatureC = 70.0;
  spec.weightSelectivity = 0.0;
  spec.weightRecovery = 1.0;
  spec.weightGreenness = 0.0;
  spec.weightPracticality = 0.0;

  const std::vector<sol::SolventCandidate> rows = sol::rankSolvents(spec);
  REQUIRE(rows.size() >= 2);
  CHECK(candidateFor(rows, "acetone") == nullptr);  // bp 56.1 C is below the 70 C dissolve step.

  const sol::Criterion* strongest = criterionFor(rows.front(), "Recovery");
  const sol::Criterion* weakest = criterionFor(rows.back(), "Recovery");
  REQUIRE(strongest != nullptr);
  REQUIRE(weakest != nullptr);
  CHECK(strongest->score > weakest->score);
  CHECK(strongest->detail.find("S_hot") != std::string::npos);
  CHECK(strongest->detail.find("S_cold") != std::string::npos);
  for (size_t i = 1; i < rows.size(); ++i) {
    CHECK(rows[i - 1].score >= rows[i].score);
  }
}

TEST_CASE("neutral extraction ranks high distribution first and excludes miscible liquids") {
  sol::OperationSpec spec;
  spec.kind = sol::OperationKind::LiquidLiquidExtraction;
  spec.species = {role("c1ccc2ccccc2c1", "naphthalene")};
  spec.requireWaterImmiscible = true;
  spec.weightGreenness = 0.0;
  spec.weightPracticality = 0.0;

  const std::vector<sol::SolventCandidate> rows = sol::rankSolvents(spec);
  REQUIRE(!rows.empty());
  for (const sol::SolventCandidate& row : rows) {
    REQUIRE(row.solvent != nullptr);
    CHECK(!row.solvent->waterMiscible);
    CHECK(rows.front().selectivity >= row.selectivity);
  }
  const sol::Criterion* recovery = criterionFor(rows.front(), "Recovery");
  REQUIRE(recovery != nullptr);
  CHECK(recovery->detail.find("organic phase") != std::string::npos);
  CHECK(recovery->detail.find("99%") != std::string::npos);
}

TEST_CASE("a better-extracted contaminant lowers selectivity and ranking") {
  sol::OperationSpec spec = extractionSpec();
  spec.species.push_back(role("c1ccc2ccccc2c1", "naphthalene", false));
  spec.weightSelectivity = 4.0;
  spec.weightRecovery = 0.2;

  const std::vector<sol::SolventCandidate> rows = sol::rankSolvents(spec);
  REQUIRE(rows.size() >= 2);
  const sol::SolventCandidate* worst = nullptr;
  double worstSelectivity = 1.0;
  for (const sol::SolventCandidate& row : rows) {
    const sol::Criterion* selectivity = criterionFor(row, "Selectivity");
    REQUIRE(selectivity != nullptr);
    if (selectivity->score < worstSelectivity) {
      worstSelectivity = selectivity->score;
      worst = &row;
    }
  }
  REQUIRE(worst != nullptr);
  CHECK(worstSelectivity < 0.5);
  CHECK(positionOf(rows, worst->solvent->id) > 0);
  const sol::Criterion* explanation = criterionFor(*worst, "Selectivity");
  REQUIRE(explanation != nullptr);
  CHECK(explanation->detail.find("worst contaminant") != std::string::npos);
}

TEST_CASE("extraction direction gate demotes preferential contaminant extraction") {
  sol::OperationSpec spec;
  spec.kind = sol::OperationKind::LiquidLiquidExtraction;
  spec.requireWaterImmiscible = true;
  spec.species = {
      role("CC(=O)Nc1ccc(O)cc1", "paracetamol"),
      role("Oc1ccccc1", "phenol contaminant", false),
  };

  const std::vector<sol::SolventCandidate> rows = sol::rankSolvents(spec);
  REQUIRE(!rows.empty());
  checkWorkingCandidatesPrecedeFailures(rows);

  const sol::SolventCandidate* failure = firstDirectionFailure(rows);
  REQUIRE(failure != nullptr);
  const sol::Criterion* evidence = criterionFor(*failure, "Selectivity");
  REQUIRE(evidence != nullptr);
  CHECK(evidence->detail.find("preferentially extracted") != std::string::npos);
  CHECK(evidence->detail.find("x over the target") != std::string::npos);
  CHECK(evidence->detail.find("Multiplicative separation gate") !=
        std::string::npos);
  CHECK(warningContains(*failure, "Separation direction fails"));
}

TEST_CASE("greenness cannot promote a reversed separation over a working solvent") {
  sol::OperationSpec spec;
  spec.kind = sol::OperationKind::LiquidLiquidExtraction;
  spec.requireWaterImmiscible = true;
  spec.species = {
      role("CC(=O)Nc1ccc(O)cc1", "paracetamol"),
      role("Oc1ccccc1", "phenol contaminant", false),
  };
  spec.weightSelectivity = 0.001;
  spec.weightRecovery = 0.001;
  spec.weightGreenness = 1000000.0;
  spec.weightPracticality = 1.0;

  const std::vector<sol::SolventCandidate> rows = sol::rankSolvents(spec);
  REQUIRE(!rows.empty());
  checkWorkingCandidatesPrecedeFailures(rows);
}

TEST_CASE("every separation operation explains its failed direction") {
  auto checkFailure = [](sol::OperationSpec spec,
                         const std::string& operationWording) {
    const std::vector<sol::SolventCandidate> rows = sol::rankSolvents(spec);
    REQUIRE(!rows.empty());
    const sol::SolventCandidate* failure = firstDirectionFailure(rows);
    REQUIRE(failure != nullptr);
    const sol::Criterion* evidence = criterionFor(*failure, "Selectivity");
    REQUIRE(evidence != nullptr);
    CHECK(evidence->detail.find(operationWording) != std::string::npos);
    CHECK(evidence->detail.find("Multiplicative separation gate") !=
          std::string::npos);
    CHECK(warningContains(*failure, "Separation direction fails"));
  };

  SUBCASE("liquid-liquid extraction") {
    sol::OperationSpec spec = extractionSpec();
    spec.species.push_back(
        role("c1ccc2ccccc2c1", "naphthalene contaminant", false));
    checkFailure(spec, "preferentially extracted");
  }

  SUBCASE("trituration") {
    sol::OperationSpec spec;
    spec.kind = sol::OperationKind::Trituration;
    spec.species = {
        role("Cn1c(=O)c2c(ncn2C)n(C)c1=O", "caffeine"),
        role("c1ccc2ccccc2c1", "naphthalene contaminant", false),
    };
    checkFailure(spec, "preferentially dissolved");
  }

  SUBCASE("recrystallisation") {
    sol::OperationSpec spec;
    spec.kind = sol::OperationKind::Recrystallisation;
    spec.coldTemperatureC = 5.0;
    spec.hotTemperatureC = 70.0;
    spec.species = {
        role("c1ccc2ccccc2c1", "naphthalene"),
        role("OC(=O)c1ccccc1", "benzoic acid contaminant", false),
    };
    checkFailure(spec, "preferentially crystallised");
  }

  SUBCASE("anti-solvent precipitation") {
    sol::OperationSpec spec;
    spec.kind = sol::OperationKind::AntiSolventPrecipitation;
    spec.species = {
        role("Oc1ccccc1", "phenol"),
        role("c1ccc2ccccc2c1", "naphthalene contaminant", false),
    };
    checkFailure(spec, "preferentially precipitated");
  }

  SUBCASE("chromatography") {
    sol::OperationSpec spec;
    spec.kind = sol::OperationKind::ChromatographyMobilePhase;
    spec.species = {
        role("Cn1c(=O)c2c(ncn2C)n(C)c1=O", "caffeine"),
        role("Cn1c(=O)c2c(ncn2C)n(C)c1=O", "caffeine contaminant", false),
    };
    checkFailure(spec, "co-elute");
  }
}

TEST_CASE("reaction medium explicitly does not request a separation direction") {
  sol::OperationSpec spec;
  spec.kind = sol::OperationKind::ReactionMedium;
  spec.species = {
      role("Cn1c(=O)c2c(ncn2C)n(C)c1=O", "caffeine"),
      role("OC(=O)c1ccccc1", "benzoic acid reagent", false),
  };

  const std::vector<sol::SolventCandidate> rows = sol::rankSolvents(spec);
  REQUIRE(!rows.empty());
  for (const sol::SolventCandidate& row : rows) {
    const sol::Criterion* selectivity = criterionFor(row, "Selectivity");
    REQUIRE(selectivity != nullptr);
    CHECK(selectivity->detail.find("selectivity is not desired") !=
          std::string::npos);
    CHECK(selectivity->detail.find("Multiplicative separation gate") ==
          std::string::npos);
    CHECK(!warningContains(row, "Separation direction fails"));
  }
}

TEST_CASE("all failed separations remain returned closest-first and fully warned") {
  sol::OperationSpec spec = extractionSpec();
  spec.requireWaterImmiscible = true;
  spec.species = {
      role("Cn1c(=O)c2c(ncn2C)n(C)c1=O", "caffeine"),
      role("OC(=O)c1ccccc1", "benzoic acid contaminant", false),
  };

  const std::vector<sol::SolventCandidate> rows = sol::rankSolvents(spec);
  REQUIRE(rows.size() >= 2);
  for (size_t i = 0; i < rows.size(); ++i) {
    CHECK(rows[i].selectivity < 1.0);
    CHECK(warningContains(rows[i], "Separation direction fails"));
    if (i > 0) {
      CHECK(rows[i - 1].selectivity >= rows[i].selectivity);
    }
  }
}

TEST_CASE("hard solvent constraints remove every violation") {
  sol::OperationSpec base;
  base.kind = sol::OperationKind::ReactionMedium;
  base.species = {role("OC(=O)c1ccccc1", "benzoic acid")};
  const std::vector<sol::SolventCandidate> all = sol::rankSolvents(base);
  REQUIRE(!all.empty());

  SUBCASE("boiling window") {
    sol::OperationSpec constrained = base;
    constrained.minBoilingPointC = 70.0;
    constrained.maxBoilingPointC = 120.0;
    const auto rows = sol::rankSolvents(constrained);
    CHECK(rows.size() < all.size());
    for (const auto& row : rows) {
      CHECK(row.solvent->boilingPoint >= 70.0);
      CHECK(row.solvent->boilingPoint <= 120.0);
    }
  }

  SUBCASE("chlorinated") {
    sol::OperationSpec constrained = base;
    constrained.avoidChlorinated = true;
    const auto rows = sol::rankSolvents(constrained);
    CHECK(rows.size() < all.size());
    for (const auto& row : rows) {
      CHECK(row.solvent->smiles.find("Cl") == std::string::npos);
    }
  }

  SUBCASE("peroxide formers") {
    sol::OperationSpec constrained = base;
    constrained.avoidPeroxideFormers = true;
    const auto rows = sol::rankSolvents(constrained);
    CHECK(rows.size() < all.size());
    for (const auto& row : rows) CHECK(!row.solvent->peroxideFormer);
  }

  SUBCASE("CHEM21 class ceiling") {
    sol::OperationSpec constrained = base;
    constrained.worstAcceptableClass = "recommended";
    constrained.excludeUnrated = true;
    const auto rows = sol::rankSolvents(constrained);
    CHECK(rows.size() < all.size());
    for (const auto& row : rows) CHECK(row.solvent->chem21Class == "recommended");
  }
}

TEST_CASE("greenness weight can promote a safer solvent and remains explained") {
  sol::OperationSpec performance = extractionSpec();
  const std::vector<sol::SolventCandidate> performanceRows = sol::rankSolvents(performance);
  REQUIRE(performanceRows.size() >= 2);

  const sol::SolventCandidate* hazardous = nullptr;
  const sol::SolventCandidate* greener = nullptr;
  for (size_t i = 0; i < performanceRows.size() && !greener; ++i) {
    if (performanceRows[i].solvent->chem21Class != "hazardous" &&
        performanceRows[i].solvent->chem21Class != "highly hazardous") {
      continue;
    }
    for (size_t j = i + 1; j < performanceRows.size(); ++j) {
      if (performanceRows[j].solvent->chem21Class == "recommended") {
        hazardous = &performanceRows[i];
        greener = &performanceRows[j];
        break;
      }
    }
  }
  REQUIRE(hazardous != nullptr);
  REQUIRE(greener != nullptr);
  const std::string hazardousId = hazardous->solvent->id;
  const std::string greenerId = greener->solvent->id;
  CHECK(positionOf(performanceRows, hazardousId) < positionOf(performanceRows, greenerId));

  sol::OperationSpec greenWeighted = performance;
  greenWeighted.weightSelectivity = 0.05;
  greenWeighted.weightRecovery = 0.05;
  greenWeighted.weightGreenness = 10.0;
  const std::vector<sol::SolventCandidate> greenRows = sol::rankSolvents(greenWeighted);
  CHECK(positionOf(greenRows, greenerId) < positionOf(greenRows, hazardousId));

  const sol::SolventCandidate* greenRow = candidateFor(greenRows, greenerId);
  const sol::SolventCandidate* hazardRow = candidateFor(greenRows, hazardousId);
  REQUIRE(greenRow != nullptr);
  REQUIRE(hazardRow != nullptr);
  const sol::Criterion* greenCriterion = criterionFor(*greenRow, "Greenness");
  const sol::Criterion* hazardCriterion = criterionFor(*hazardRow, "Greenness");
  REQUIRE(greenCriterion != nullptr);
  REQUIRE(hazardCriterion != nullptr);
  CHECK(greenCriterion->weight == doctest::Approx(10.0));
  CHECK(greenCriterion->score > hazardCriterion->score);
  CHECK(greenCriterion->detail.find("CHEM21") != std::string::npos);
}

TEST_CASE("ranking is deterministic down to order and score") {
  sol::OperationSpec spec = extractionSpec();
  spec.species.push_back(role("OC(=O)c1ccccc1", "benzoic acid impurity", false, 2.0));
  const std::vector<sol::SolventCandidate> first = sol::rankSolvents(spec);
  const std::vector<sol::SolventCandidate> second = sol::rankSolvents(spec);
  REQUIRE(first.size() == second.size());
  for (size_t i = 0; i < first.size(); ++i) {
    REQUIRE(first[i].solvent != nullptr);
    REQUIRE(second[i].solvent != nullptr);
    CHECK(first[i].solvent->id == second[i].solvent->id);
    CHECK(first[i].score == second[i].score);
    CHECK(first[i].selectivity == second[i].selectivity);
    CHECK(first[i].criteria.size() == second[i].criteria.size());
  }
}

TEST_CASE("degenerate selection input is safe") {
  sol::OperationSpec empty;
  CHECK_NOTHROW(sol::rankSolvents(empty));
  CHECK(sol::rankSolvents(empty).empty());

  sol::OperationSpec contaminantsOnly;
  contaminantsOnly.species = {role("O", "water contaminant", false)};
  CHECK_NOTHROW(sol::rankSolvents(contaminantsOnly));
  CHECK(sol::rankSolvents(contaminantsOnly).empty());

  sol::OperationSpec impossible;
  impossible.kind = sol::OperationKind::ReactionMedium;
  impossible.species = {role("CCO", "ethanol")};
  impossible.requireWaterMiscible = true;
  impossible.requireWaterImmiscible = true;
  CHECK_NOTHROW(sol::rankSolvents(impossible));
  CHECK(sol::rankSolvents(impossible).empty());

  sol::OperationSpec degenerate;
  degenerate.kind = sol::OperationKind::ReactionMedium;
  sol::SpeciesRole zero;
  zero.label = "zero descriptors";
  zero.keep = true;
  degenerate.species = {zero};
  CHECK_NOTHROW(sol::rankSolvents(degenerate));
  const auto rows = sol::rankSolvents(degenerate);
  for (const auto& row : rows) {
    CHECK(std::isfinite(row.score));
    CHECK(row.score >= 0.0);
    CHECK(row.score <= 1.0);
  }
}

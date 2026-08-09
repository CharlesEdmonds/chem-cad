#pragma once

#include <string>
#include <vector>

#include "sol/solubility.hpp"
#include "sol/solvent.hpp"

namespace chemcad::sol {

enum class OperationKind {
  LiquidLiquidExtraction,
  Recrystallisation,
  Trituration,
  AntiSolventPrecipitation,
  ChromatographyMobilePhase,
  ReactionMedium,
};

// One chemical species involved in the operation.
struct SpeciesRole {
  Solute solute;
  std::string label;
  bool keep = true;
  double weight = 1.0;
  double amountMg = 0.0;
};

// Everything about the operation that is not a species.
struct OperationSpec {
  OperationKind kind = OperationKind::LiquidLiquidExtraction;
  std::vector<SpeciesRole> species;
  double temperatureC = 25.0;
  double hotTemperatureC = 65.0;
  double coldTemperatureC = 0.0;
  double aqueousVolumeMl = 100.0;
  double organicVolumeMl = 100.0;
  double pH = kAutoPH;
  bool requireWaterImmiscible = false;
  bool requireWaterMiscible = false;
  double minBoilingPointC = 0.0;
  double maxBoilingPointC = 0.0;
  bool avoidPeroxideFormers = false;
  bool avoidChlorinated = false;
  bool avoidAromatics = false;
  bool excludeUnrated = false;
  std::string worstAcceptableClass;
  double weightSelectivity = 1.0;
  double weightRecovery = 1.0;
  double weightGreenness = 0.6;
  double weightPracticality = 0.5;
};

// One scored contribution, so the UI can explain the ranking.
struct Criterion {
  std::string name;
  double score = 0.0;
  double weight = 0.0;
  std::string detail;
};

struct SolventCandidate {
  const Solvent* solvent = nullptr;
  const Solvent* partner = nullptr;
  double partnerFraction = 0.0;
  double score = 0.0;
  std::vector<Criterion> criteria;
  double selectivity = 0.0;
  double recoveryFraction = 0.0;
  double targetSolubilityGPerMl = 0.0;
  double contaminantSolubilityGPerMl = 0.0;
  std::vector<std::string> warnings;
  bool estimated = true;
};

// Ranked best-first. Throws SolError only when the solvent database itself fails to load.
// An operation with no target species returns an empty vector rather than throwing.
std::vector<SolventCandidate> rankSolvents(const OperationSpec&);

// Human-readable label for a kind, and the one-line description the UI shows.
const char* operationName(OperationKind);
const char* operationDescription(OperationKind);

}  // namespace chemcad::sol

#pragma once
// Persistent UI state for the Solubility Suite's `Select` mode. Construction
// stays data-only: the solvent database is loaded lazily by rankSolvents().

#include <cstdint>
#include <string>
#include <vector>

#include "sol/selection.hpp"

namespace chemcad::ui {

struct SelectionSpeciesPresentation {
  std::string formula;
  std::string smiles;
};

struct SelectionState {
  sol::OperationSpec operation;
  std::vector<SelectionSpeciesPresentation> speciesPresentation;

  int addMode = 0;  // 0 sketch, 1 SMILES, 2 chemical name
  std::string smilesInput;
  std::string nameInput;
  std::string inputError;
  bool nameLookupRunning = false;
  uint64_t nameRequest = 0;

  std::vector<sol::SolventCandidate> candidates;
  std::string rankingSignature;
  std::string pendingSignature;
  std::string observedSignature;
  double signatureChangedAt = 0.0;
  std::string rankingError;
  bool computing = false;
  uint64_t resultsRevision = 0;

  std::string resultSearch;
  int sortMode = 0;  // score, selectivity, recovery, greenness, boiling point
  bool compareMode = false;
  std::vector<std::string> comparedCandidates;
  std::string statusMessage;
};

}  // namespace chemcad::ui

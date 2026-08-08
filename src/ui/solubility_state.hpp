#pragma once
// Panel + simulation state for the Solubility Suite (prediction and the
// separatory-funnel view). Held as a member of AppState, so it must default-
// construct without touching the filesystem -- the solvent database and
// solute group-contribution tables are only loaded lazily, on first draw.

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "sol/funnel.hpp"
#include "sol/solubility.hpp"

namespace chemcad::ui {

// Payload handed from the Solubility Suite to the Extraction Lab: the suite's
// "Send to Extraction Lab" button fills this in and sets `pending`; the
// extraction view consumes it on its next frame and clears the flag.
struct ExtractionImport {
  bool pending = false;
  std::string solventIdA;  // the aqueous (or lighter-workflow) phase
  std::string solventIdB;  // the organic partner
  double volumeMlA = 50.0;
  double volumeMlB = 50.0;
  double soluteMassMg = 100.0;
};

struct SolubilityState {
  // ------------------------------------------------------------- solute
  bool useSketch = true;        // true: pull from st.doc.molecules; false: manualSmiles
  std::string manualSmiles;

  sol::Solute solute;
  bool soluteValid = false;
  std::string soluteError;

  bool overrideSolute = false;  // user-supplied melting point / R0 vs group contribution
  float meltingPointC = 25.0f;
  float interactionRadius = 8.0f;

  uint64_t sourceRevision = 0;  // last st.docRevision the sketch solute was built from
  uint64_t soluteVersion = 0;   // bumped on every successful (re)compute or override edit

  // ------------------------------------------------------------ solvents
  int solventCount = 2;                  // 1..3 active slots
  std::array<std::string, 3> solventIds{};  // sol::Solvent::id per slot, empty = unset
  std::array<float, 3> ratios{1.0f, 1.0f, 1.0f};  // volume parts, need not be normalised

  float temperatureC = 25.0f;

  // -------------------------------------------------------------- result
  sol::Prediction prediction;

  // ------------------------------------------------------------- display
  int units = 0;        // 0 g/mL, 1 mg/mL, 2 g/100 mL, 3 mol/L
  bool logScale = false;  // ratio-plot y axis: linear or log10

  // -------------------------------------------------------- ratio sweep
  std::vector<sol::SweepPoint> sweep;
  int sweepSteps = 20;
  std::string sweepSignature;  // cache key: solvent ids + steps + temperature + soluteVersion
  int sweepPeakIndex = -1;  // index into `sweep` of the max-solubility sample, or -1

  // ------------------------------------------------------- solvent screen
  std::vector<sol::ScreenRow> screening;  // pure-solvent table, best first
  std::string screeningSignature;         // cache key: soluteVersion + temperature

  // ------------------------------------------------- common-ion effect
  // Only surfaced when the solute is a 1:1 salt (sol::findSalt match).
  bool backgroundEnabled = false;
  int backgroundElectrolyte = 0;    // index into sol::electrolytes()
  float backgroundMolarity = 0.5f;  // mol/L

  // ---------------------------------------------------------- funnel sim
  sol::Simulation funnel;
  bool funnelRunning = false;
  float funnelSpeed = 1.0f;
  // Physical shake inputs (sol::ShakeParams): what a chemist actually sets.
  float shakeDurationS = 5.0f;    // s
  float shakeFrequencyHz = 3.0f;  // Hz, 2-4 is a firm hand shake
  float shakeAmplitudeCm = 5.0f;  // cm stroke half-amplitude
  int funnelVessel = 0;  // mirrors sol::Vessel, kept as int for a plain ImGui combo

  ExtractionImport extractionImport;  // suite -> extraction lab hand-off

  std::string statusMessage;
};

}  // namespace chemcad::ui

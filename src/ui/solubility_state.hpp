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

// One fixed-C slice of the ternary composition cube. Points are ordered along
// the A:B sweep; keeping the slices here makes the expensive predictions a
// signature-keyed cache rather than per-frame drawing work.
struct TernaryLayerSweep {
  float fractionC = 0.0f;
  std::vector<sol::SweepPoint> points;
  int peakIndex = -1;
};

struct SolubilityState {
  // ------------------------------------------------------------- solute
  bool useSketch = true;        // true: pull from st.doc.molecules; false: manualSmiles
  std::string manualSmiles;

  sol::Solute solute;
  bool soluteValid = false;
  std::string soluteError;
  // Non-fatal advisory about how the solute was interpreted (e.g. a sketch
  // holding several fragments, where only the largest is the solute).
  std::string soluteNote;

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

  // pH matters only for an ionisable solute. Auto means "whatever the solute's
  // own saturated solution sets" -- a free base runs basic, its hydrochloride
  // runs acidic -- which is what happens when you drop the solid into
  // unbuffered water.
  bool pHAuto = true;
  float pH = 7.0f;

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
  std::vector<TernaryLayerSweep> ternaryLayers;  // fixed-C curves, cached with `sweep`

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
  bool funnelRender3D = false;  // false: flat 2D cross-section (default); true: shaded
  int funnelVessel = 0;  // mirrors sol::Vessel, kept as int for a plain ImGui combo

  // Grab-and-shake: the user drags the vessel itself; the drag velocity is
  // the slosh velocity. Offset is render-side, velocity feeds the physics.
  bool funnelGrabbed = false;
  float funnelGrabAnchorX = 0.0f;  // mouse x at grab, px
  float funnelDragOffsetPx = 0.0f; // vessel render offset, px
  float funnelMouseVel = 0.0f;     // smoothed |dx/dt| in vessel m/s

  ExtractionImport extractionImport;  // suite -> extraction lab hand-off

  std::string statusMessage;
};

}  // namespace chemcad::ui

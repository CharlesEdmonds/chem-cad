#pragma once
// Panel + simulation state for the Solubility Suite (prediction and the
// separatory-funnel view). Held as a member of AppState, so it must default-
// construct without touching the filesystem -- the solvent database and
// solute group-contribution tables are only loaded lazily, on first draw.

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "fluid/simulation.hpp"
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

// Cached geometry for the three-solvent composition surface. Both linear and
// logarithmic heights, colours and Lambert shades are retained so changing the
// display scale never resweeps or re-tessellates the model.
struct TernarySurfaceNode {
  sol::SweepPoint point;
  std::array<float, 3> cubeLinear{};  // ratio B, normalised solubility, fraction C
  std::array<float, 3> cubeLog{};
};

struct TernarySurfaceQuad {
  std::array<uint32_t, 4> nodes{};  // front-left, front-right, back-right, back-left
  std::array<float, 4> colourLinear{};
  std::array<float, 4> colourLog{};
  float shadeLinear = 1.0f;
  float shadeLog = 1.0f;
  float depth = 0.0f;
  int ratioIndex = 0;
  int cIndex = 0;
};

struct TernarySurfaceMesh {
  int ratioQuads = 0;
  int cQuads = 0;
  std::vector<TernarySurfaceNode> nodes;
  std::vector<TernarySurfaceQuad> quads;  // cached in back-to-front order
  int peakNode = -1;
  double minimum = 0.0;
  double maximum = 0.0;
  double logFloor = 1e-12;
  double logMinimum = 0.0;
  double logMaximum = 1.0;
};

enum class ExtractionRenderMode : uint8_t {
  Fluid3D,
  Schematic2D,
};

enum class FluidShakeAxis : uint8_t {
  Vertical,
  Horizontal,
  Diagonal,
};

enum class FluidResolution : uint8_t {
  Interactive,
  Balanced,
  Quality,
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
  bool logScale = false;  // composition graph height: linear or log10

  // -------------------------------------------------------- ratio sweep
  std::vector<sol::SweepPoint> sweep;
  int sweepSteps = 20;
  std::string sweepSignature;  // cache key: solvent ids + steps + temperature + soluteVersion
  int sweepPeakIndex = -1;  // index into `sweep` of the max-solubility sample, or -1
  TernarySurfaceMesh ternarySurface;  // full cached composition surface

  // ------------------------------------------------------- solvent screen
  std::vector<sol::ScreenRow> screening;  // pure-solvent table, best first
  std::string screeningSignature;         // cache key: soluteVersion + temperature

  // ------------------------------------------------- common-ion effect
  // Only surfaced when the solute is a 1:1 salt (sol::findSalt match).
  bool backgroundEnabled = false;
  int backgroundElectrolyte = 0;    // index into sol::electrolytes()
  float backgroundMolarity = 0.5f;  // mol/L

  // ---------------------------------------------------------- extraction
  // The legacy analytic funnel remains the editable charge model used by the
  // partition and wash calculators. The non-copyable particle simulation is
  // created only when the Extraction Lab is first drawn.
  sol::Simulation funnel;
  std::unique_ptr<fluid::Simulation> fluid;
  bool funnelRunning = false;
  float funnelSpeed = 1.0f;
  int funnelVessel = 0;  // mirrors sol::Vessel, kept as int for a plain ImGui combo

  // Physical driven-shake inputs. The axis is explicit because vertical
  // shaking must be as discoverable as horizontal shaking.
  float shakeDurationS = 5.0f;
  float shakeFrequencyHz = 3.0f;
  float shakeAmplitudeCm = 5.0f;
  FluidShakeAxis shakeAxis = FluidShakeAxis::Vertical;
  ExtractionRenderMode extractionRenderMode = ExtractionRenderMode::Fluid3D;
  FluidResolution fluidResolution = FluidResolution::Interactive;

  // Measured throughput is retained per resolution so the preset picker can
  // show an observed cost alongside the particle count before switching.
  std::array<double, 3> fluidPresetRealTimeFactor{};
  std::array<bool, 3> fluidPresetRealTimeFactorValid{};
  std::array<uint64_t, 3> fluidPresetMeasuredParticles{};

  // Shake progress is measured against Simulation::elapsedS(), never wall
  // time. A slow solver therefore advances this countdown only when it
  // actually completes simulated time.
  double fluidShakeStartElapsedS = 0.0;
  double fluidShakeEndElapsedS = 0.0;
  bool fluidShakeProgressValid = false;

  // Toolbar actions are consumed by the stage after its snapshot is available.
  bool fluidReframeRequested = false;

  // A dedicated stage toggle distinguishes grab-and-shake from camera orbit.
  // Pointer deltas are converted into world acceleration and filtered before
  // they cross the Simulation API.
  bool fluidGrabMode = false;
  bool fluidGrabActive = false;
  std::array<float, 2> fluidGrabAnchorPx{};
  std::array<double, 3> fluidManualAcceleration{};

  // Pose animation state. The angular rate is retained so setPose receives an
  // angular acceleration consistent with the visible tilt animation.
  float fluidTiltTargetDeg = 0.0f;
  float fluidTiltCurrentDeg = 0.0f;
  float fluidTiltAngularVelocityRadS = 0.0f;

  ExtractionImport extractionImport;  // suite -> extraction lab hand-off

  std::string statusMessage;
};

}  // namespace chemcad::ui

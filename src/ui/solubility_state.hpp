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

namespace chemcad::app {
class TaskRunner;
}

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

// The Solubility Suite hosts two views of the same question. `Predict` answers
// "how much of this dissolves in the blend I chose"; `Select` answers "which
// solvent should I choose for this operation". They were separate top-level
// panels, which forced the user to leave the workspace to answer the second
// half of one decision, so they are now two modes of one panel and hand their
// results to each other in place.
enum class SuiteMode : uint8_t {
  Predict,
  Select,
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
  // ------------------------------------------------------------ workspace
  SuiteMode suiteMode = SuiteMode::Predict;

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
  // created only once the user is demonstrably on the Extraction workspace.
  // ImGui submits every docked panel on the frame the dock layout is built and
  // reports it focused, so a single draw proves nothing; building the
  // simulation on that guess used to be the entire cost of the application's
  // first frame. `extractionLastDrawnFrame` turns that into evidence: from the
  // second frame onwards Begin() hides unselected dock tabs correctly, so two
  // consecutive draws mean the panel really is on top.
  sol::Simulation funnel;
  // Building the particle simulation runs the interfacial calibration, which
  // costs seconds the first time a given resolution and material pair is seen
  // on a machine. It therefore happens on a worker thread and the panel shows
  // the analytic schematic until it lands; `fluidBuildSignature` is the
  // configuration the in-flight build was started for, so a build whose inputs
  // the user has already changed is discarded instead of installed.
  // shared_ptr, not unique_ptr: TaskRunner completions live in std::function,
  // which requires a copyable payload.
  std::shared_ptr<fluid::Simulation> fluid;
  app::TaskRunner* fluidTasks = nullptr;  // borrowed from AppState each frame
  bool fluidBuildPending = false;
  std::size_t fluidBuildSignature = 0;
  int extractionLastDrawnFrame = -1;
  bool fluidConstructionAllowed = false;
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

  // Which console tab the extraction workspace is showing. Panel view state
  // belongs here beside extractionRenderMode, not in a function-local static:
  // it is user-visible, it must survive a tab switch, and a headless test has
  // to be able to select a tab without synthesising a click.
  int extractionTab = 0;

  // Dragging the stage shakes the vessel, in both the 3D and the schematic
  // view; orbiting the 3D camera is on the right button. The pointer commands
  // where the hand IS and fluid::HandFollower turns that into the vessel's
  // motion, so the vessel can be flung clean off the stage and swings back on
  // release. Only the pointer bookkeeping lives here; the mechanics are physics
  // and live in fluid/frame.hpp.
  bool fluidGrabActive = false;
  std::array<float, 2> fluidGrabAnchorPx{};

  // Grabbing the vessel promotes the stage to a full-window overlay drawn above
  // every other panel, so the funnel can be flung right across the application
  // instead of being clipped at the edge of its dock rect. The presented rect
  // lags the request by one frame because the renderer composites into an FBO
  // that ImGui shows on the NEXT frame; drawing this frame's texture into this
  // frame's rect would stretch it during the transition.
  bool fluidOverlayActive = false;
  std::array<float, 4> fluidPresentedRectPx{};  // min x, min y, max x, max y
  bool fluidPresentedValid = false;
  // Vertical field of view of the docked stage. The overlay narrows it so the
  // vessel keeps its on-screen size when the render target grows.
  float fluidStageFovDeg = 38.0f;
  // Last stage size, so a resized panel re-frames instead of cropping.
  std::array<float, 2> fluidStageSizePx{};

  fluid::HandFollower fluidHand;

  // Pose animation state. The angular rate is retained so setPose receives an
  // angular acceleration consistent with the visible tilt animation.
  float fluidTiltTargetDeg = 0.0f;
  float fluidTiltCurrentDeg = 0.0f;
  float fluidTiltAngularVelocityRadS = 0.0f;

  ExtractionImport extractionImport;  // suite -> extraction lab hand-off

  std::string statusMessage;
};

}  // namespace chemcad::ui

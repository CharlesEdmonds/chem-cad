#pragma once
// The single mutable state object threaded through every panel.
// Wave-1 UI slices may APPEND fields; renaming or repurposing existing fields
// breaks sibling slices.

#include <chrono>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "app/task_runner.hpp"
#include "core/model.hpp"
#include "rxn/engine.hpp"
#include "ui/camera.hpp"

namespace chemcad::ui {

// ---------------------------------------------------------------- addressing
// Atom/bond ids are unique only within their Molecule, so document-level
// references carry the fragment index too.
struct AtomRef {
  int mol = -1;
  core::AtomId id = core::kInvalidAtom;
  bool valid() const { return mol >= 0 && id != core::kInvalidAtom; }
  bool operator==(const AtomRef&) const = default;
};
struct BondRef {
  int mol = -1;
  core::BondId id = core::kInvalidBond;
  bool valid() const { return mol >= 0 && id != core::kInvalidBond; }
  bool operator==(const BondRef&) const = default;
};

struct Selection {
  std::vector<AtomRef> atoms;
  std::vector<BondRef> bonds;
  bool empty() const { return atoms.empty() && bonds.empty(); }
  void clear() {
    atoms.clear();
    bonds.clear();
  }
  bool contains(const AtomRef& r) const;
  bool contains(const BondRef& r) const;
};

// ---------------------------------------------------------------- tools
enum class Tool {
  Select,
  Eraser,
  Bond,
  Chain,
  RingTemplate,
  Atom,
  ChargePlus,
  ChargeMinus,
};

// Ring stamped by Tool::RingTemplate.
enum class RingKind {
  Cyclopropane,
  Cyclobutane,
  Cyclopentane,
  Cyclohexane,
  Cycloheptane,
  Cyclooctane,
  Benzene,
  Cyclopentadiene,
  Naphthalene,
};

// ---------------------------------------------------------------- async status
enum class Status { Idle, Loading, Ok, Error };

struct PropertiesCache {
  std::string formula;
  double mw = 0;
  double logP = 0;
  int rings = 0;
  std::string chemError;   // valence/sanitization message, empty when clean
  std::string smiles;      // canonical SMILES of the current document
  std::string name;        // IUPAC name (async)
  Status nameStatus = Status::Idle;
  std::string nameError;
  bool autoName = true;
  // Debounce bookkeeping: the canvas bumps docRevision on every mutation.
  uint64_t computedForRevision = ~0ull;
  uint64_t nameRequestedForRevision = ~0ull;
  std::chrono::steady_clock::time_point lastEdit{};
};

// ---------------------------------------------------------------- planner
struct MaterialBox {
  std::string smiles;
  std::string nameInput;   // text box for a name lookup
  std::string label;       // resolved display name, when known
  core::Molecule preview;  // laid-out structure for the thumbnail
  bool previewValid = false;
  std::string error;
  Status status = Status::Idle;
};

struct PlannerState {
  std::vector<MaterialBox> starts{MaterialBox{}};  // at least one
  MaterialBox target;
  std::vector<rxn::Route> routes;
  bool searching = false;
  bool searched = false;  // a search has completed at least once
  std::string error;
  bool allowLlm = true;
  int maxDepth = 3;
  int maxRoutes = 5;
};

enum class MainTab { Sketch, Planner };

// ---------------------------------------------------------------- app state
struct AppState {
  core::Document doc;
  core::UndoStack undo;
  uint64_t docRevision = 1;  // bumped by every mutation; drives caches

  Tool tool = Tool::Bond;
  uint8_t currentElement = 6;
  core::BondOrder currentOrder = core::BondOrder::Single;
  core::BondStereo currentStereo = core::BondStereo::None;
  RingKind currentRing = RingKind::Benzene;

  Selection sel;
  Camera2D cam;
  AtomRef hoverAtom;
  BondRef hoverBond;

  app::TaskRunner tasks;
  PropertiesCache props;
  PlannerState planner;

  MainTab tab = MainTab::Sketch;
  bool tabChangeRequested = false;  // set with `tab` to force the tab bar

  std::string projectPath;
  bool dirty = false;
  std::string statusMessage;
  std::string clipboardSmiles;

  // Wired up by the app layer so panels stay free of GL/file-format concerns.
  std::function<void(const std::string& path)> exportPng;
  std::function<void(const std::string& path)> exportSvg;
  std::function<void(const std::string& path)> exportMol;
  std::function<void(const std::string& path)> importMol;
  std::function<void(const std::string& path)> saveProject;
  std::function<void(const std::string& path)> openProject;

  // Canvas viewport in screen coordinates, published by drawCanvas() so the PNG
  // exporter knows which pixels to grab.
  core::Vec2 canvasOrigin{0, 0};
  core::Vec2 canvasSize{0, 0};

  // Set by exportPng(); consumed by the main loop after rendering but before
  // the buffer swap, which is the only moment the canvas pixels exist.
  std::optional<std::string> pendingPngExport;

  // Display terminal single-bonded carbons as CH3 instead of relying entirely
  // on skeletal-notation line ends. View menu can restore the compact form.
  bool showTerminalMethylLabels = true;

  void touch() {
    ++docRevision;
    dirty = true;
    props.lastEdit = std::chrono::steady_clock::now();
  }
  // Snapshot for undo. Call BEFORE mutating the document.
  void snapshot() { undo.push(doc); }

  core::Molecule* molecule(int index);
  core::Atom* atomAt(const AtomRef&);
  core::Bond* bondAt(const BondRef&);
};

}  // namespace chemcad::ui

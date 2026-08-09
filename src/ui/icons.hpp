#pragma once
// Procedural tool glyphs drawn with ImDrawList. No icon font: every glyph is
// a handful of strokes on a normalized grid, which keeps them crisp at any UI
// scale and gives the palette its own look.
//
// The set covers three families: sketch tools (the original ChemDraw palette),
// bench apparatus (glassware and instruments, so a chemistry control reads as
// chemistry rather than as a generic button), and interface controls
// (transport, disclosure, layout). A control that can be named with a glyph
// should carry one -- a label-only wall of buttons is the thing this set exists
// to replace.

#include "imgui.h"

namespace chemcad::ui::icons {

enum class Icon {
  None,  // draws nothing; lets widgets take an optional glyph without a flag

  // ---------------------------------------------------------- sketch tools
  Select,
  Eraser,
  Bond,
  Chain,
  RingCyclopropane,
  RingCyclobutane,
  RingCyclopentane,
  RingCyclohexane,
  RingCycloheptane,
  RingCyclooctane,
  RingBenzene,
  RingCyclopentadiene,
  RingNaphthalene,
  Atom,
  ChargePlus,
  ChargeMinus,
  BondSingle,
  BondDouble,
  BondTriple,
  BondAromatic,
  StereoNone,
  StereoWedge,
  StereoHash,
  StereoWavy,

  // ------------------------------------------------------ bench apparatus
  Flask,          // round-bottom flask: a solution / prediction
  Beaker,         // bulk solvent
  SepFunnel,      // separatory funnel: the extraction workspace
  Droplet,        // a solvent, a phase, a dose
  Vial,           // a sample, a screening hit
  Thermometer,    // temperature
  Balance,        // mass, partition, equilibrium
  Molecule,       // a structure, as opposed to a single atom
  Reaction,       // forward reaction arrow over a baseline
  Retro,          // retrosynthetic double-shafted arrow
  Ph,             // pH: a droplet on a scale
  Energy,         // bolt: enthalpy, driving force
  Flame,          // heating, reflux
  Snowflake,      // cooling, crystallisation
  Shake,          // hand agitation: the funnel's own verb
  Timer,          // contact time, residence time

  // ------------------------------------------------------- data & display
  ChartBars,
  ChartLine,
  ChartScatter,
  Gauge,
  Layers,        // stacked phases, stacked traces
  Grid,          // grid layout
  List,          // list layout
  Cube,          // 3D view
  Ruler,         // measurement overlay
  Table,

  // ------------------------------------------------------------- controls
  Play,
  Pause,
  Stop,
  Rewind,
  StepForward,
  Minus,
  Plus,
  Close,
  Check,
  Search,
  Filter,
  Copy,
  Trash,
  Save,
  Folder,
  Undo,
  Redo,
  Gear,
  Eye,
  EyeOff,
  Lock,
  Link,
  Star,
  Warning,
  Info,
  Sparkle,        // AI-sourced result
  Book,           // knowledge-base-sourced result
  Crosshair,      // reframe / recentre
  ZoomFit,
  ArrowRight,
  ArrowLeft,
  ChevronDown,
  ChevronRight,
  ChevronLeft,
  ChevronUp,
  DragHandle,

  Logo,  // benzene-ring brand mark
};

// Draws `icon` centred at `centre`, inside a square of side `size`.
// `thickness` <= 0 picks a stroke proportional to `size`.
void draw(ImDrawList* dl, Icon icon, ImVec2 centre, float size, ImU32 color,
          float thickness = 0.0f);

}  // namespace chemcad::ui::icons

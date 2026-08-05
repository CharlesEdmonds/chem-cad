// Headless interaction tests. These drive the REAL panel code
// (drawCanvas / drawPropertiesPanel / drawReactionPlanner) through a null ImGui
// backend, synthesizing mouse and keyboard events into ImGuiIO exactly as a
// windowing backend would. No GL context and no display are involved, so the
// sketching gestures are provable in CI and over SSH.
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <cmath>
#include <cstring>
#include <string>

#include "imgui.h"
#include "imgui_internal.h"

#include "chem/bridge.hpp"
#include "core/model.hpp"
#include "core/sprout.hpp"
#include "ui/app_state.hpp"
#include "ui/ui.hpp"

using namespace chemcad;

namespace {

constexpr float kDisplayW = 1600.0f;
constexpr float kDisplayH = 1000.0f;

// Minimal null backend: an ImGui context with a built font atlas and no
// renderer. Everything ImGui needs to run a full frame, nothing more.
struct HeadlessImGui {
  HeadlessImGui() {
    IMGUI_CHECKVERSION();
    ctx = ImGui::CreateContext();
    ImGui::SetCurrentContext(ctx);
    ImGuiIO& io = ImGui::GetIO();
    io.DisplaySize = ImVec2(kDisplayW, kDisplayH);
    io.DeltaTime = 1.0f / 60.0f;
    io.IniFilename = nullptr;  // never touch the user's layout file
    io.LogFilename = nullptr;
    unsigned char* pixels = nullptr;
    int w = 0, h = 0;
    io.Fonts->GetTexDataAsRGBA32(&pixels, &w, &h);
    io.Fonts->SetTexID(static_cast<ImTextureID>(1));
  }
  ~HeadlessImGui() { ImGui::DestroyContext(ctx); }
  ImGuiContext* ctx = nullptr;
};

// Runs one full frame with the sketch canvas docked into a full-screen window.
void canvasFrame(ui::AppState& st) {
  ImGui::NewFrame();
  ImGui::SetNextWindowPos(ImVec2(0, 0));
  ImGui::SetNextWindowSize(ImVec2(kDisplayW, kDisplayH));
  ImGui::Begin("Sketch", nullptr,
               ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                   ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings);
  ui::drawCanvas(st);
  ImGui::End();
  ImGui::Render();
}

void panelFrame(ui::AppState& st, void (*draw)(ui::AppState&), const char* title) {
  ImGui::NewFrame();
  ImGui::SetNextWindowPos(ImVec2(0, 0));
  ImGui::SetNextWindowSize(ImVec2(700, 900));
  ImGui::Begin(title, nullptr,
               ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                   ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings);
  draw(st);
  ImGui::End();
  ImGui::Render();
}

void mouseTo(float x, float y) { ImGui::GetIO().AddMousePosEvent(x, y); }

// A click needs distinct frames: ImGui latches hover on one frame, the press on
// the next and the release after that.
void clickAt(ui::AppState& st, float x, float y) {
  mouseTo(x, y);
  canvasFrame(st);
  ImGui::GetIO().AddMouseButtonEvent(ImGuiMouseButton_Left, true);
  canvasFrame(st);
  ImGui::GetIO().AddMouseButtonEvent(ImGuiMouseButton_Left, false);
  canvasFrame(st);
  canvasFrame(st);
}

// Same gesture, but driving an arbitrary panel window.
void panelClickAt(ui::AppState& st, void (*draw)(ui::AppState&), const char* title,
                  float x, float y) {
  mouseTo(x, y);
  panelFrame(st, draw, title);
  ImGui::GetIO().AddMouseButtonEvent(ImGuiMouseButton_Left, true);
  panelFrame(st, draw, title);
  ImGui::GetIO().AddMouseButtonEvent(ImGuiMouseButton_Left, false);
  panelFrame(st, draw, title);
  panelFrame(st, draw, title);
}

void dragFromTo(ui::AppState& st, float x0, float y0, float x1, float y1) {
  mouseTo(x0, y0);
  canvasFrame(st);
  ImGui::GetIO().AddMouseButtonEvent(ImGuiMouseButton_Left, true);
  canvasFrame(st);
  mouseTo((x0 + x1) * 0.5f, (y0 + y1) * 0.5f);
  canvasFrame(st);
  mouseTo(x1, y1);
  canvasFrame(st);
  ImGui::GetIO().AddMouseButtonEvent(ImGuiMouseButton_Left, false);
  canvasFrame(st);
  canvasFrame(st);
}

void pressKey(ui::AppState& st, ImGuiKey key) {
  ImGui::GetIO().AddKeyEvent(key, true);
  canvasFrame(st);
  ImGui::GetIO().AddKeyEvent(key, false);
  canvasFrame(st);
}

size_t totalAtoms(const ui::AppState& st) {
  size_t n = 0;
  for (const core::Molecule& m : st.doc.molecules) n += m.atomCount();
  return n;
}
size_t totalBonds(const ui::AppState& st) {
  size_t n = 0;
  for (const core::Molecule& m : st.doc.molecules) n += m.bondCount();
  return n;
}

// Screen position of an atom, for aiming the next gesture at it.
ImVec2 screenOf(const ui::AppState& st, int mol, size_t atomIndex) {
  const core::Atom& a = st.doc.molecules[static_cast<size_t>(mol)].atoms()[atomIndex];
  const core::Vec2 s = st.cam.worldToScreen(a.pos, st.canvasOrigin);
  return ImVec2(s.x, s.y);
}

}  // namespace

TEST_CASE("bond tool draws and extends a skeleton") {
  HeadlessImGui gui;
  ui::AppState st;
  st.tool = ui::Tool::Bond;

  // Warm-up frame so the canvas publishes its viewport rect.
  canvasFrame(st);
  REQUIRE(st.canvasSize.x > 100.0f);

  // Clicking empty canvas starts a new fragment: a C-C bond.
  clickAt(st, 700.0f, 500.0f);
  REQUIRE(st.doc.molecules.size() == 1);
  CHECK(totalAtoms(st) == 2);
  CHECK(totalBonds(st) == 1);

  // The auto bond angle must be the standard 30 degrees above horizontal, and
  // the bond exactly one bond length long.
  const core::Molecule& m = st.doc.molecules[0];
  const core::Vec2 p0 = m.atoms()[0].pos;
  const core::Vec2 p1 = m.atoms()[1].pos;
  const float dx = p1.x - p0.x, dy = p1.y - p0.y;
  CHECK(std::sqrt(dx * dx + dy * dy) == doctest::Approx(core::kBondLength).epsilon(0.01));
  CHECK(std::atan2(dy, dx) * 180.0f / 3.14159265f == doctest::Approx(30.0f).epsilon(0.02));

  // Clicking the terminal atom sprouts another bond -> a 3-carbon chain.
  const ImVec2 tip = screenOf(st, 0, 1);
  clickAt(st, tip.x, tip.y);
  CHECK(totalAtoms(st) == 3);
  CHECK(totalBonds(st) == 2);

  // Auto bond angles must zig-zag, not fold back: the chain is not collinear.
  const core::Molecule& chain = st.doc.molecules[0];
  const core::Vec2 a = chain.atoms()[0].pos, b = chain.atoms()[1].pos, c = chain.atoms()[2].pos;
  const float cross = (b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x);
  CHECK(std::fabs(cross) > 0.1f);
}

TEST_CASE("M preset attaches a visible methyl-ready terminal carbon") {
  HeadlessImGui gui;
  ui::AppState st;
  st.tool = ui::Tool::Bond;
  canvasFrame(st);
  clickAt(st, 700.0f, 500.0f);
  REQUIRE(totalAtoms(st) == 2);

  // M always restores a plain single bond, even after another bond/stereo mode.
  st.tool = ui::Tool::Chain;
  st.currentOrder = core::BondOrder::Triple;
  st.currentStereo = core::BondStereo::Wedge;
  const ImVec2 tip = screenOf(st, 0, 1);
  mouseTo(tip.x, tip.y);
  canvasFrame(st);
  REQUIRE(st.hoverAtom.valid());
  pressKey(st, ImGuiKey_M);
  CHECK(st.tool == ui::Tool::Bond);
  CHECK(st.currentOrder == core::BondOrder::Single);
  CHECK(st.currentStereo == core::BondStereo::None);

  clickAt(st, tip.x, tip.y);
  REQUIRE(totalAtoms(st) == 3);
  REQUIRE(totalBonds(st) == 2);
  const core::Molecule& propane = st.doc.molecules[0];
  const core::Atom& methyl = propane.atoms().back();
  CHECK(chem::implicitHCount(propane, methyl.id) == 3);
  CHECK(chem::canonicalize(chem::toSmiles(propane)) == chem::canonicalize("CCC"));
  CHECK(st.statusMessage == "Added methyl group (CH3)");
}

TEST_CASE("terminal methyl labels can be shown or hidden") {
  HeadlessImGui gui;
  ui::AppState st;
  st.doc.molecules.push_back(chem::fromSmiles("CC"));
  st.touch();

  st.showTerminalMethylLabels = false;
  canvasFrame(st);
  const int skeletalVertices = ImGui::GetDrawData()->TotalVtxCount;

  st.showTerminalMethylLabels = true;
  canvasFrame(st);
  const int labelledVertices = ImGui::GetDrawData()->TotalVtxCount;
  CHECK(labelledVertices > skeletalVertices);
}

TEST_CASE("hovering an atom and pressing an element key retypes it") {
  HeadlessImGui gui;
  ui::AppState st;
  st.tool = ui::Tool::Bond;
  canvasFrame(st);
  clickAt(st, 700.0f, 500.0f);
  REQUIRE(totalAtoms(st) == 2);

  // Extend to a three-carbon skeleton so retyping the tip yields ethanol.
  const ImVec2 firstTip = screenOf(st, 0, 1);
  clickAt(st, firstTip.x, firstTip.y);
  REQUIRE(totalAtoms(st) == 3);

  // Park the cursor over the terminal atom, then press O.
  const ImVec2 tip = screenOf(st, 0, 2);
  mouseTo(tip.x, tip.y);
  canvasFrame(st);
  REQUIRE(st.hoverAtom.valid());
  pressKey(st, ImGuiKey_O);

  const core::Atom* atom = st.atomAt(st.hoverAtom);
  REQUIRE(atom != nullptr);
  CHECK(atom->atomicNumber == 8);

  // C-C-O is ethanol; the chemistry layer must agree with the sketch.
  const std::string smiles = chem::toSmiles(st.doc.molecules[0]);
  CHECK(chem::canonicalize(smiles) == chem::canonicalize("CCO"));
}

TEST_CASE("undo and redo restore the document through the canvas shortcuts") {
  HeadlessImGui gui;
  ui::AppState st;
  st.tool = ui::Tool::Bond;
  canvasFrame(st);
  clickAt(st, 700.0f, 500.0f);
  const size_t afterFirst = totalAtoms(st);
  REQUIRE(afterFirst == 2);

  const ImVec2 tip = screenOf(st, 0, 1);
  clickAt(st, tip.x, tip.y);
  REQUIRE(totalAtoms(st) == 3);

  ImGui::GetIO().AddKeyEvent(ImGuiMod_Ctrl, true);
  pressKey(st, ImGuiKey_Z);
  ImGui::GetIO().AddKeyEvent(ImGuiMod_Ctrl, false);
  canvasFrame(st);
  CHECK(totalAtoms(st) == afterFirst);

  ImGui::GetIO().AddKeyEvent(ImGuiMod_Ctrl, true);
  ImGui::GetIO().AddKeyEvent(ImGuiMod_Shift, true);
  pressKey(st, ImGuiKey_Z);
  ImGui::GetIO().AddKeyEvent(ImGuiMod_Shift, false);
  ImGui::GetIO().AddKeyEvent(ImGuiMod_Ctrl, false);
  canvasFrame(st);
  CHECK(totalAtoms(st) == 3);
}

TEST_CASE("eraser removes the hovered atom and its bonds") {
  HeadlessImGui gui;
  ui::AppState st;
  st.tool = ui::Tool::Bond;
  canvasFrame(st);
  clickAt(st, 700.0f, 500.0f);
  REQUIRE(totalAtoms(st) == 2);
  REQUIRE(totalBonds(st) == 1);

  const ImVec2 tip = screenOf(st, 0, 1);
  st.tool = ui::Tool::Eraser;
  clickAt(st, tip.x, tip.y);
  CHECK(totalAtoms(st) == 1);
  CHECK(totalBonds(st) == 0);
}

TEST_CASE("ring template stamps a benzene ring") {
  HeadlessImGui gui;
  ui::AppState st;
  st.tool = ui::Tool::RingTemplate;
  st.currentRing = ui::RingKind::Benzene;
  canvasFrame(st);
  clickAt(st, 700.0f, 500.0f);

  CHECK(totalAtoms(st) == 6);
  CHECK(totalBonds(st) == 6);
  int aromatic = 0;
  for (const core::Bond& b : st.doc.molecules[0].bonds())
    if (b.order == core::BondOrder::Aromatic) ++aromatic;
  CHECK(aromatic == 6);
  CHECK(chem::canonicalize(chem::toSmiles(st.doc.molecules[0])) ==
        chem::canonicalize("c1ccccc1"));
}

TEST_CASE("atom tool places the element chosen in the periodic table") {
  HeadlessImGui gui;
  ui::AppState st;
  // Selecting tungsten in the periodic table sets these two fields.
  st.currentElement = 74;
  st.tool = ui::Tool::Atom;
  canvasFrame(st);
  clickAt(st, 700.0f, 500.0f);

  REQUIRE(totalAtoms(st) == 1);
  CHECK(st.doc.molecules[0].atoms()[0].atomicNumber == 74);
}

TEST_CASE("charge tools adjust formal charge and clamp") {
  HeadlessImGui gui;
  ui::AppState st;
  st.tool = ui::Tool::Atom;
  st.currentElement = 7;
  canvasFrame(st);
  clickAt(st, 700.0f, 500.0f);
  REQUIRE(totalAtoms(st) == 1);

  const ImVec2 at = screenOf(st, 0, 0);
  st.tool = ui::Tool::ChargePlus;
  clickAt(st, at.x, at.y);
  CHECK(st.doc.molecules[0].atoms()[0].charge == 1);

  st.tool = ui::Tool::ChargeMinus;
  clickAt(st, at.x, at.y);
  clickAt(st, at.x, at.y);
  CHECK(st.doc.molecules[0].atoms()[0].charge == -1);
}

TEST_CASE("drag closes a ring onto an existing atom instead of duplicating it") {
  HeadlessImGui gui;
  ui::AppState st;
  st.tool = ui::Tool::Bond;
  canvasFrame(st);
  clickAt(st, 700.0f, 500.0f);
  REQUIRE(totalAtoms(st) == 2);

  // Drag from one end back onto the other end: that must BOND them, not add an atom.
  const ImVec2 a = screenOf(st, 0, 0);
  const ImVec2 b = screenOf(st, 0, 1);
  const size_t before = totalAtoms(st);
  dragFromTo(st, b.x, b.y, a.x, a.y);
  CHECK(totalAtoms(st) == before);  // no new atom created
}

TEST_CASE("properties panel reports formula and mass for the sketch") {
  HeadlessImGui gui;
  ui::AppState st;
  st.doc.molecules.push_back(chem::fromSmiles("CCO"));
  st.touch();
  // The panel debounces on wall-clock time, so age the last edit past the window.
  st.props.lastEdit = std::chrono::steady_clock::now() - std::chrono::seconds(3);
  st.props.autoName = false;  // no network in tests

  for (int i = 0; i < 3; ++i) panelFrame(st, &ui::drawPropertiesPanel, "Properties");

  CHECK(st.props.formula == "C2H6O");
  CHECK(st.props.mw == doctest::Approx(46.0419).epsilon(0.001));
  CHECK(chem::canonicalize(st.props.smiles) == chem::canonicalize("CCO"));
}

TEST_CASE("tool palette icon grid switches the active tool") {
  HeadlessImGui gui;
  ui::AppState st;
  st.tool = ui::Tool::Select;
  // The 2-column grid starts at the window's content origin (8,8); cells are
  // 66px squares (cap = 3 * iconSize at default scale), rows stride 66+4.
  panelFrame(st, &ui::drawToolPalette, "Tools");

  // Eraser is row 0, column 1: grid is centred in the 700px window.
  panelClickAt(st, &ui::drawToolPalette, "Tools", 387.0f, 41.0f);
  CHECK(st.tool == ui::Tool::Eraser);

  // Bond is row 1, column 0.
  panelClickAt(st, &ui::drawToolPalette, "Tools", 313.0f, 111.0f);
  CHECK(st.tool == ui::Tool::Bond);
}

TEST_CASE("bond order icon row switches the bond order") {
  HeadlessImGui gui;
  ui::AppState st;
  panelFrame(st, &ui::drawToolPalette, "Tools");
  REQUIRE(st.currentOrder == core::BondOrder::Single);

  // Four cells share the row after the Bond section header; double is cell 1.
  panelClickAt(st, &ui::drawToolPalette, "Tools", 263.0f, 329.0f);
  CHECK(st.currentOrder == core::BondOrder::Double);
}

TEST_CASE("clicking a periodic table tile arms the atom tool with that element") {
  HeadlessImGui gui;
  ui::AppState st;
  panelFrame(st, &ui::drawPeriodicTable, "Periodic Table");

  // Hydrogen is the top-left tile; its position is pinned at the grid origin
  // regardless of how the cell-fitting math resolves.
  panelClickAt(st, &ui::drawPeriodicTable, "Periodic Table", 26.0f, 49.0f);
  CHECK(st.currentElement == 1);
  CHECK(st.tool == ui::Tool::Atom);
}

TEST_CASE("hovering a tool button fires its tooltip after the hover delay") {
  HeadlessImGui gui;
  ui::AppState st;
  // Park the cursor over the eraser cell and let the hover/stationary delays
  // elapse (0.40s + 0.15s at 1/60s per frame), then look for the tooltip.
  mouseTo(387.0f, 41.0f);
  for (int i = 0; i < 60; ++i) panelFrame(st, &ui::drawToolPalette, "Tools");

  bool tooltipActive = false;
  ImGuiContext& g = *ImGui::GetCurrentContext();
  for (ImGuiWindow* window : g.Windows) {
    if (window->WasActive && std::strstr(window->Name, "Tooltip") != nullptr) {
      tooltipActive = true;
      break;
    }
  }
  CHECK(tooltipActive);
}

TEST_CASE("periodic table search selects a unique match on Enter") {
  HeadlessImGui gui;
  ui::AppState st;
  panelFrame(st, &ui::drawPeriodicTable, "Periodic Table");

  // Focus the search field, type, hit Enter: the unique match is selected.
  panelClickAt(st, &ui::drawPeriodicTable, "Periodic Table", 350.0f, 17.0f);
  ImGui::GetIO().AddInputCharactersUTF8("tungsten");
  panelFrame(st, &ui::drawPeriodicTable, "Periodic Table");
  ImGui::GetIO().AddKeyEvent(ImGuiKey_Enter, true);
  panelFrame(st, &ui::drawPeriodicTable, "Periodic Table");
  ImGui::GetIO().AddKeyEvent(ImGuiKey_Enter, false);
  panelFrame(st, &ui::drawPeriodicTable, "Periodic Table");

  CHECK(st.currentElement == 74);
  CHECK(st.tool == ui::Tool::Atom);
}

TEST_CASE("reaction planner renders routes with their side products") {
  HeadlessImGui gui;
  ui::AppState st;
  st.planner.starts.clear();
  ui::MaterialBox acid, alcohol;
  acid.smiles = "CC(=O)O";
  alcohol.smiles = "CCO";
  st.planner.starts.push_back(acid);
  st.planner.starts.push_back(alcohol);
  st.planner.target.smiles = "CCOC(C)=O";

  // Inject a finished search result: this test covers the RESULT RENDERING path,
  // which must not crash and must survive repeated frames with cached previews.
  rxn::Step step;
  step.reactionName = "Fischer esterification";
  step.reactantSmiles = {"CC(=O)O", "CCO"};
  step.reagents = {"H2SO4 (cat.)"};
  step.conditions = "reflux";
  step.productSmiles = "CCOC(C)=O";
  step.sideProductSmiles = {"O"};
  step.source = rxn::Step::Source::KB;
  rxn::Route route;
  route.steps.push_back(step);
  st.planner.routes.push_back(route);
  st.planner.searched = true;

  for (int i = 0; i < 3; ++i) panelFrame(st, &ui::drawReactionPlanner, "Reaction Planner");

  REQUIRE(st.planner.routes.size() == 1);
  CHECK(st.planner.routes[0].steps[0].sideProductSmiles[0] == "O");
  CHECK_FALSE(st.planner.routes[0].usesLlm());
}

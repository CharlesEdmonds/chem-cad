// Headless render test for the 3D viewer: drives the real drawViewer3D()
// through a null ImGui backend, proving the turntable renderer emits geometry
// and follows the sketch, and that the Orbitals mode evaluates real hydrogenic
// wavefunctions rather than a stylised blob.
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <set>
#include <string>
#include <vector>

#include "imgui.h"

#include "chem/bridge.hpp"
#include "ui/app_state.hpp"
#include "ui/charts3d.hpp"
#include "ui/ui.hpp"
#include "ui/viewer3d_state.hpp"

using namespace chemcad;

namespace {

constexpr float kDisplayW = 1400.0f;
constexpr float kDisplayH = 1000.0f;

struct HeadlessImGui {
  HeadlessImGui() {
    IMGUI_CHECKVERSION();
    ctx = ImGui::CreateContext();
    ImGui::SetCurrentContext(ctx);
    ImGuiIO& io = ImGui::GetIO();
    io.DisplaySize = ImVec2(kDisplayW, kDisplayH);
    io.DeltaTime = 1.0f / 60.0f;
    io.IniFilename = nullptr;
    io.LogFilename = nullptr;
    unsigned char* pixels = nullptr;
    int w = 0, h = 0;
    io.Fonts->GetTexDataAsRGBA32(&pixels, &w, &h);
    io.Fonts->SetTexID(static_cast<ImTextureID>(1));
  }
  ~HeadlessImGui() { ImGui::DestroyContext(ctx); }
  ImGuiContext* ctx = nullptr;
};

// What one rendered frame tells us. `vertices` says geometry was emitted at
// all; `hash` folds in every emitted vertex position and colour, so two frames
// that differ anywhere -- a lobe at a different radius, a node ring that is
// not there -- differ here too, which a vertex count alone cannot show.
// `scrollMaxY` says whether the panel overflowed the surface it was given.
struct Geometry {
  int vertices = 0;
  std::uint64_t hash = 0;
};

Geometry captureGeometry() {
  const ImDrawData* data = ImGui::GetDrawData();
  Geometry geometry;
  geometry.vertices = data->TotalVtxCount;
  std::uint64_t hash = 14695981039346656037ull;  // FNV-1a offset basis
  for (int list = 0; list < data->CmdListsCount; ++list) {
    const ImDrawList* commands = data->CmdLists[list];
    for (int index = 0; index < commands->VtxBuffer.Size; ++index) {
      const ImDrawVert& vertex = commands->VtxBuffer[index];
      std::uint32_t words[3] = {0, 0, vertex.col};
      std::memcpy(&words[0], &vertex.pos.x, sizeof(float));
      std::memcpy(&words[1], &vertex.pos.y, sizeof(float));
      for (std::uint32_t word : words) {
        hash ^= word;
        hash *= 1099511628211ull;
      }
    }
  }
  geometry.hash = hash;
  return geometry;
}

struct FrameResult {
  Geometry geometry;
  float scrollMaxY = 0.0f;
  int vertices() const { return geometry.vertices; }
};

FrameResult panelFrameIn(ui::AppState& st, ImVec2 size) {
  ImGui::NewFrame();
  ImGui::SetNextWindowPos(ImVec2(0, 0));
  ImGui::SetNextWindowSize(size);
  ImGui::Begin("3D View", nullptr,
               ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                   ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings);
  ui::drawViewer3D(st);
  FrameResult result;
  result.scrollMaxY = ImGui::GetScrollMaxY();
  ImGui::End();
  ImGui::Render();
  result.geometry = captureGeometry();
  return result;
}

int panelFrame(ui::AppState& st) {
  return panelFrameIn(st, ImVec2(kDisplayW, kDisplayH)).vertices();
}

// ImGui reports a window's content extent from the frame before, so a size
// change needs one frame to settle before GetScrollMaxY() means anything.
FrameResult settledFrame(ui::AppState& st, ImVec2 size) {
  panelFrameIn(st, size);
  return panelFrameIn(st, size);
}

// Drives charts3d::orbital() on its own, away from the panel chrome, so a
// difference in the result can only come from the orbital itself. The Orbit is
// fresh per call and its idle spin is off, so the projection is deterministic.
Geometry orbitalFrame(int n, int l, int m, const ui::charts3d::OrbitalStyle& style,
                      ImVec2 size = ImVec2(560.0f, 560.0f)) {
  ui::charts3d::Orbit orbit;
  ImGui::NewFrame();
  ImGui::SetNextWindowPos(ImVec2(0, 0));
  ImGui::SetNextWindowSize(ImVec2(size.x + 64.0f, size.y + 64.0f));
  ImGui::Begin("Orbital", nullptr,
               ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                   ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings);
  ui::charts3d::orbital("##orbital_probe", n, l, m, size, orbit, style);
  ImGui::End();
  ImGui::Render();
  return captureGeometry();
}

// Every (n, l, m) the Orbitals toolbar can reach: n in 1..5, l in 0..min(n-1,3),
// m in -l..l.
template <typename Visit>
void forEachValidOrbital(Visit visit) {
  for (int n = 1; n <= 5; ++n) {
    for (int l = 0; l <= std::min(n - 1, 3); ++l) {
      for (int m = -l; m <= l; ++m) visit(n, l, m);
    }
  }
}

}  // namespace

TEST_CASE("the viewer renders an empty state without a sketch") {
  HeadlessImGui gui;
  ui::AppState st;
  CHECK_NOTHROW(panelFrame(st));
  CHECK(!st.viewer3d.hasModel);
}

TEST_CASE("the viewer embeds the sketch and draws real geometry") {
  HeadlessImGui gui;

  ui::AppState bare;
  const int bareVertices = panelFrame(bare);

  ui::AppState st;
  st.doc.molecules.push_back(chem::fromSmiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C"));
  st.touch();
  const int vertices = panelFrame(st);

  REQUIRE(st.viewer3d.hasModel);
  CHECK(st.viewer3d.formula == "C8H10N4O2");
  CHECK(vertices > bareVertices);
}

TEST_CASE("the model re-embeds when the sketch changes") {
  HeadlessImGui gui;
  ui::AppState st;
  st.doc.molecules.push_back(chem::fromSmiles("CCO"));
  const uint64_t firstRevision = st.docRevision;
  st.touch();
  panelFrame(st);
  REQUIRE(st.viewer3d.hasModel);
  const size_t ethanolAtoms = st.viewer3d.model.atoms.size();

  // A new molecule at a new revision must replace the cached conformer.
  st.doc.molecules.clear();
  st.doc.molecules.push_back(chem::fromSmiles("c1ccccc1"));
  st.touch();
  REQUIRE(st.docRevision != firstRevision);
  panelFrame(st);
  CHECK(st.viewer3d.model.atoms.size() != ethanolAtoms);
  CHECK(st.viewer3d.formula == "C6H6");
}

TEST_CASE("every render style draws without throwing") {
  HeadlessImGui gui;
  ui::AppState st;
  st.doc.molecules.push_back(chem::fromSmiles("c1ccccc1"));
  st.touch();
  for (int style = 0; style < 3; ++style) {
    st.viewer3d.style = style;
    CHECK_NOTHROW(panelFrame(st));
  }
}

// --------------------------------------------------------------- orbitals

TEST_CASE("the Orbitals mode is reachable and draws its own geometry") {
  HeadlessImGui gui;
  ui::AppState structure;
  const int structureVertices = panelFrame(structure);

  ui::AppState orbitals;
  orbitals.viewer3d.mode = 1;
  int orbitalVertices = 0;
  CHECK_NOTHROW(orbitalVertices = panelFrame(orbitals));

  // The Orbitals mode owes nothing to the sketch, so it must draw a real
  // surface from an empty document where Structure only shows an empty state.
  CHECK(orbitalVertices > structureVertices);
  CHECK(!orbitals.viewer3d.hasModel);
  CHECK(orbitals.viewer3d.mode == 1);
}

TEST_CASE("every valid orbital renders and its shape depends on n, l and m") {
  HeadlessImGui gui;
  ui::charts3d::OrbitalStyle style;
  style.showNodes = false;
  style.showAxes = false;

  std::set<std::uint64_t> shapes;
  std::vector<int> counts;
  forEachValidOrbital([&](int n, int l, int m) {
    Geometry drawn;
    CAPTURE(n);
    CAPTURE(l);
    CAPTURE(m);
    CHECK_NOTHROW(drawn = orbitalFrame(n, l, m, style));
    CHECK(drawn.vertices > 0);
    shapes.insert(drawn.hash);
    counts.push_back(drawn.vertices);
  });
  REQUIRE(counts.size() == 46);

  // A renderer that ignored the quantum numbers would draw one picture for all
  // 46. Every hydrogenic surface here is genuinely its own: different radial
  // extent, different shell count, or the same lobes turned to a different
  // place on screen.
  CHECK(shapes.size() == 46);

  // The case the header promises: 3d_z2 must not look like 3p_z.
  CHECK(orbitalFrame(3, 2, 0, style).vertices != orbitalFrame(3, 1, 0, style).vertices);
  // Same angular part, different radial part: 2p_z has no radial node, 3p_z
  // has one, so their surfaces cannot coincide. 3p_z against 4p_z is the same
  // argument one shell further out.
  CHECK(orbitalFrame(2, 1, 0, style).vertices != orbitalFrame(3, 1, 0, style).vertices);
  CHECK(orbitalFrame(3, 1, 0, style).vertices != orbitalFrame(4, 1, 0, style).vertices);
  // Same shell and family, different m: d_z2 against d_x2-y2.
  CHECK(orbitalFrame(3, 2, 0, style).vertices != orbitalFrame(3, 2, 2, style).vertices);

  // The s family is the strict test of the radial factor: every ns orbital is
  // spherical, so nothing but its n - 1 radial nodes can tell them apart, and
  // fitting to the view cancels even the difference in size. Each further
  // shell must add one more nested iso surface.
  std::vector<int> sCounts;
  for (int n = 1; n <= 5; ++n) sCounts.push_back(orbitalFrame(n, 0, 0, style).vertices);
  const std::set<int> distinctS(sCounts.begin(), sCounts.end());
  CHECK(distinctS.size() == 5);
  for (int n = 1; n < 5; ++n) {
    CAPTURE(n);
    CHECK(sCounts[static_cast<std::size_t>(n)] >
          sCounts[static_cast<std::size_t>(n) - 1]);
  }
}

TEST_CASE("orbitalName gives the spectroscopic name for every valid orbital") {
  HeadlessImGui gui;
  using ui::charts3d::orbitalName;

  CHECK(std::string(orbitalName(1, 0, 0)) == "1s");
  CHECK(std::string(orbitalName(2, 0, 0)) == "2s");
  CHECK(std::string(orbitalName(2, 1, 0)) == "2p z");
  CHECK(std::string(orbitalName(2, 1, 1)) == "2p x");
  CHECK(std::string(orbitalName(2, 1, -1)) == "2p y");
  CHECK(std::string(orbitalName(3, 2, 0)) == "3d z2");
  CHECK(std::string(orbitalName(3, 2, 1)) == "3d xz");
  CHECK(std::string(orbitalName(3, 2, -1)) == "3d yz");
  CHECK(std::string(orbitalName(3, 2, 2)) == "3d x2-y2");
  CHECK(std::string(orbitalName(3, 2, -2)) == "3d xy");
  CHECK(std::string(orbitalName(4, 3, 0)) == "4f z3");
  CHECK(std::string(orbitalName(4, 3, 1)) == "4f xz2");
  CHECK(std::string(orbitalName(4, 3, -1)) == "4f yz2");
  CHECK(std::string(orbitalName(4, 3, 2)) == "4f z(x2-y2)");
  CHECK(std::string(orbitalName(4, 3, -2)) == "4f xyz");
  CHECK(std::string(orbitalName(5, 3, 3)) == "5f x(x2-3y2)");
  CHECK(std::string(orbitalName(5, 3, -3)) == "5f y(3x2-y2)");

  // Every reachable orbital is named, and the family letter follows l.
  static const char kFamily[] = {'s', 'p', 'd', 'f'};
  forEachValidOrbital([&](int n, int l, int m) {
    CAPTURE(n);
    CAPTURE(l);
    CAPTURE(m);
    const std::string name = orbitalName(n, l, m);
    REQUIRE(name.size() >= 2);
    CHECK(name != "invalid");
    CHECK(name[0] == static_cast<char>('0' + n));
    CHECK(name[1] == kFamily[l]);
  });

  // Outside the valid set the contract is the literal "invalid".
  CHECK(std::string(orbitalName(0, 0, 0)) == "invalid");
  CHECK(std::string(orbitalName(-1, 0, 0)) == "invalid");
  CHECK(std::string(orbitalName(1, 1, 0)) == "invalid");   // l must be < n
  CHECK(std::string(orbitalName(2, -1, 0)) == "invalid");  // l must be >= 0
  CHECK(std::string(orbitalName(5, 4, 0)) == "invalid");   // g is not drawn
  CHECK(std::string(orbitalName(3, 1, 2)) == "invalid");   // |m| must be <= l
  CHECK(std::string(orbitalName(3, 2, -3)) == "invalid");
}

TEST_CASE("the quantum-number toolbar can never produce an invalid orbital") {
  HeadlessImGui gui;
  ui::AppState st;
  ui::Viewer3DState& vs = st.viewer3d;
  vs.mode = 1;

  // The failure case: the highest f orbital, then the principal number
  // collapsed to the ground shell. l and m must follow n down.
  vs.orbitalN = 5;
  vs.orbitalL = 3;
  vs.orbitalM = 3;
  panelFrame(st);
  CHECK(std::string(ui::charts3d::orbitalName(vs.orbitalN, vs.orbitalL,
                                              vs.orbitalM)) == "5f x(x2-3y2)");

  vs.orbitalN = 1;
  panelFrame(st);
  CHECK(vs.orbitalL == 0);
  CHECK(vs.orbitalM == 0);

  // Raising n back leaves a still-valid state rather than restoring garbage.
  vs.orbitalN = 5;
  panelFrame(st);
  CHECK(vs.orbitalL == 0);
  CHECK(vs.orbitalM == 0);
  CHECK(std::string(ui::charts3d::orbitalName(vs.orbitalN, vs.orbitalL,
                                              vs.orbitalM)) == "5s");

  // Any state a caller or a stale layout could leave behind must be repaired
  // before it reaches the renderer, for every principal number.
  const int hostile[][2] = {{3, 3}, {3, -3}, {2, 2}, {-4, 0}, {7, 6}, {0, 9}};
  for (int n = -2; n <= 9; ++n) {
    for (const auto& probe : hostile) {
      vs.orbitalN = n;
      vs.orbitalL = probe[0];
      vs.orbitalM = probe[1];
      CAPTURE(n);
      CAPTURE(probe[0]);
      CAPTURE(probe[1]);
      CHECK_NOTHROW(panelFrame(st));
      CHECK(vs.orbitalN >= 1);
      CHECK(vs.orbitalN <= 5);
      CHECK(vs.orbitalL >= 0);
      CHECK(vs.orbitalL < vs.orbitalN);
      CHECK(vs.orbitalL <= 3);
      CHECK(std::abs(vs.orbitalM) <= vs.orbitalL);
      CHECK(std::string(ui::charts3d::orbitalName(
                vs.orbitalN, vs.orbitalL, vs.orbitalM)) != "invalid");
    }
  }
}

TEST_CASE("every orbital display control changes what the renderer draws") {
  HeadlessImGui gui;
  // 3p_z has one radial node and one nodal plane, so it exercises both halves
  // of the Nodes control.
  ui::charts3d::OrbitalStyle plain;
  plain.showNodes = false;
  plain.cutaway = false;
  plain.showAxes = false;
  plain.isoLevel = 0.30f;
  const Geometry base = orbitalFrame(3, 1, 0, plain);
  REQUIRE(base.vertices > 0);

  ui::charts3d::OrbitalStyle nodes = plain;
  nodes.showNodes = true;
  CHECK(orbitalFrame(3, 1, 0, nodes).vertices != base.vertices);

  ui::charts3d::OrbitalStyle cutaway = plain;
  cutaway.cutaway = true;
  CHECK(orbitalFrame(3, 1, 0, cutaway).vertices < base.vertices);  // near half culled

  ui::charts3d::OrbitalStyle axes = plain;
  axes.showAxes = true;
  CHECK(orbitalFrame(3, 1, 0, axes).vertices != base.vertices);

  // The iso level is a continuous control: every step across its slider range
  // must move the surface, not just the two extremes. Two thresholds can happen
  // to enclose the same number of quads while placing them at different radii,
  // so this compares the drawn geometry rather than only counting it.
  const float levels[] = {0.02f, 0.10f, 0.30f, 0.55f, 0.80f, 0.95f};
  std::set<std::uint64_t> isoShapes;
  std::vector<int> isoCounts;
  for (float level : levels) {
    ui::charts3d::OrbitalStyle iso = plain;
    iso.isoLevel = level;
    CAPTURE(level);
    const Geometry drawn = orbitalFrame(3, 1, 0, iso);
    CHECK(drawn.vertices > 0);
    isoShapes.insert(drawn.hash);
    isoCounts.push_back(drawn.vertices);
  }
  CHECK(isoShapes.size() == std::size(levels));
  // A higher threshold selects a smaller region of the wavefunction.
  CHECK(isoCounts.front() > isoCounts.back());

  // The panel must hand its own style through, not a default-constructed one.
  ui::AppState st;
  st.viewer3d.mode = 1;
  st.viewer3d.orbitalStyle.showAxes = false;
  const int withoutAxes = panelFrame(st);
  st.viewer3d.orbitalStyle.showAxes = true;
  CHECK(panelFrame(st) != withoutAxes);
  st.viewer3d.orbitalStyle.isoLevel = 0.85f;
  CHECK(panelFrame(st) != withoutAxes);
}

TEST_CASE("the overlay node census matches the hydrogenic count") {
  char text[96];
  // Radial nodes = n - l - 1, angular nodes = l.
  ui::formatOrbitalNodes(1, 0, text, sizeof(text));
  CHECK(std::string(text) == "0 radial nodes  |  0 angular nodes");
  ui::formatOrbitalNodes(2, 0, text, sizeof(text));
  CHECK(std::string(text) == "1 radial node  |  0 angular nodes");
  ui::formatOrbitalNodes(2, 1, text, sizeof(text));
  CHECK(std::string(text) == "0 radial nodes  |  1 angular node");
  ui::formatOrbitalNodes(3, 1, text, sizeof(text));
  CHECK(std::string(text) == "1 radial node  |  1 angular node");
  ui::formatOrbitalNodes(3, 2, text, sizeof(text));
  CHECK(std::string(text) == "0 radial nodes  |  2 angular nodes");
  ui::formatOrbitalNodes(5, 0, text, sizeof(text));
  CHECK(std::string(text) == "4 radial nodes  |  0 angular nodes");
  ui::formatOrbitalNodes(5, 3, text, sizeof(text));
  CHECK(std::string(text) == "1 radial node  |  3 angular nodes");

  // Total nodes always come to n - 1, which is the invariant the readout is
  // there to make visible.
  forEachValidOrbital([&](int n, int l, int) {
    CAPTURE(n);
    CAPTURE(l);
    CHECK(ui::orbitalRadialNodes(n, l) + ui::orbitalAngularNodes(l) == n - 1);
  });
}

TEST_CASE("the orbital panel survives the stacked and split toolbar layouts") {
  HeadlessImGui gui;
  ui::AppState st;
  st.viewer3d.mode = 1;

  // Wide: one quantum row, one options row.
  const FrameResult wide = settledFrame(st, ImVec2(kDisplayW, kDisplayH));
  CHECK(wide.vertices() > 0);
  CHECK(wide.scrollMaxY == 0.0f);

  // Narrow enough that the quantum strips stack and the options split.
  const ImVec2 narrow(320.0f, 560.0f);
  FrameResult stacked;
  CHECK_NOTHROW(stacked = settledFrame(st, narrow));
  CHECK(stacked.vertices() > 0);
  // Every row is budgeted from the surface, so a split layout must still fit
  // without scrolling -- if it scrolls, the rows are overlapping.
  CHECK(stacked.scrollMaxY == 0.0f);

  // The controls still drive the renderer once the toolbar is split.
  st.viewer3d.orbitalStyle.showNodes = false;
  const int withoutNodes = panelFrameIn(st, narrow).vertices();
  st.viewer3d.orbitalStyle.showNodes = true;
  CHECK(panelFrameIn(st, narrow).vertices() != withoutNodes);

  st.viewer3d.orbitalN = 5;
  st.viewer3d.orbitalL = 3;
  st.viewer3d.orbitalM = -3;
  CHECK_NOTHROW(panelFrameIn(st, narrow));
  CHECK(st.viewer3d.orbitalL == 3);
  CHECK(st.viewer3d.orbitalM == -3);
  st.viewer3d.orbitalN = 2;
  CHECK_NOTHROW(panelFrameIn(st, narrow));
  CHECK(st.viewer3d.orbitalL == 1);
  CHECK(std::abs(st.viewer3d.orbitalM) <= 1);

  // Extremely small: the panel must still not throw or scroll.
  FrameResult tiny;
  CHECK_NOTHROW(tiny = settledFrame(st, ImVec2(180.0f, 360.0f)));
  CHECK(tiny.scrollMaxY == 0.0f);

  // Sweeping the width crosses both responsive thresholds -- the quantum
  // strips stack below three fitting columns, the options split below 44 ems --
  // and the panel must divide its surface exactly at every one of them.
  const float widths[] = {160.0f, 240.0f, 320.0f, 420.0f, 520.0f,
                          640.0f, 800.0f, 1000.0f, 1400.0f};
  for (float width : widths) {
    CAPTURE(width);
    FrameResult sized;
    CHECK_NOTHROW(sized = settledFrame(st, ImVec2(width, 620.0f)));
    CHECK(sized.vertices() > 0);
    CHECK(sized.scrollMaxY == 0.0f);
  }
}

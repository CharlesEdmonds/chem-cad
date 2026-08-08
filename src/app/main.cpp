// ChemCAD -- entry point, GLFW/OpenGL/ImGui bootstrap and the dock layout.

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "imgui_internal.h"

#include <GLFW/glfw3.h>

#include "app/project_io.hpp"
#include "app/screenshot.hpp"
#include "chem/bridge.hpp"
#include "core/paths.hpp"
#include "ui/app_state.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"

namespace {

using chemcad::ui::AppState;

constexpr const char* kWinSketch = "Sketch";
constexpr const char* kWinPlanner = "Reaction Planner";
constexpr const char* kWinSolubility = "Solubility Suite";
constexpr const char* kWinExtraction = "Extraction Lab";
constexpr const char* kWinViewer3D = "3D View";
constexpr const char* kWinTools = "Tools";
constexpr const char* kWinPTable = "Periodic Table";
constexpr const char* kWinProps = "Properties";

void glfwErrorCallback(int code, const char* description) {
  std::fprintf(stderr, "[glfw] error %d: %s\n", code, description);
}

// Concatenates every fragment into one molecule so exporters see the whole sketch.
chemcad::core::Molecule flatten(const chemcad::core::Document& doc) {
  chemcad::core::Molecule out;
  for (const chemcad::core::Molecule& m : doc.molecules) {
    std::unordered_map<chemcad::core::AtomId, chemcad::core::AtomId> idMap;
    for (const chemcad::core::Atom& a : m.atoms()) {
      chemcad::core::Atom copy = a;
      idMap[a.id] = out.addAtom(copy);
    }
    for (const chemcad::core::Bond& b : m.bonds()) {
      auto ia = idMap.find(b.a), ib = idMap.find(b.b);
      if (ia == idMap.end() || ib == idMap.end()) continue;
      chemcad::core::BondId nb = out.addBond(ia->second, ib->second, b.order);
      if (chemcad::core::Bond* nbp = out.bond(nb)) nbp->stereo = b.stereo;
    }
  }
  return out;
}

void writeTextFile(const std::string& path, const std::string& text) {
  std::ofstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("cannot open " + path + " for writing");
  f << text;
}

void wireCallbacks(AppState& st) {
  // Deferred: the framebuffer only holds the canvas after ImGui renders.
  st.exportPng = [&st](const std::string& path) { st.pendingPngExport = path; };
  st.exportSvg = [&st](const std::string& path) {
    try {
      writeTextFile(path, chemcad::chem::toSvg(flatten(st.doc), 900, 700));
      st.statusMessage = "Exported SVG: " + path;
    } catch (const std::exception& e) {
      st.statusMessage = std::string("SVG export failed: ") + e.what();
    }
  };
  st.exportMol = [&st](const std::string& path) {
    try {
      writeTextFile(path, chemcad::chem::toMolBlock(flatten(st.doc)));
      st.statusMessage = "Exported MOL: " + path;
    } catch (const std::exception& e) {
      st.statusMessage = std::string("MOL export failed: ") + e.what();
    }
  };
  st.importMol = [&st](const std::string& path) {
    try {
      std::ifstream f(path, std::ios::binary);
      if (!f) throw std::runtime_error("cannot open " + path);
      const std::string text((std::istreambuf_iterator<char>(f)),
                             std::istreambuf_iterator<char>());
      st.snapshot();
      st.doc.molecules.push_back(chemcad::chem::fromMolBlock(text));
      st.touch();
      st.statusMessage = "Imported " + path;
    } catch (const std::exception& e) {
      st.statusMessage = std::string("Import failed: ") + e.what();
    }
  };
  st.saveProject = [&st](const std::string& path) {
    try {
      chemcad::app::Project p;
      p.doc = st.doc;
      for (const auto& b : st.planner.starts) p.plannerStarts.push_back(b.smiles);
      p.plannerTarget = st.planner.target.smiles;
      chemcad::app::saveProject(p, path);
      st.projectPath = path;
      st.dirty = false;
      st.statusMessage = "Saved " + path;
    } catch (const std::exception& e) {
      st.statusMessage = std::string("Save failed: ") + e.what();
    }
  };
  st.openProject = [&st](const std::string& path) {
    try {
      chemcad::app::Project p = chemcad::app::loadProject(path);
      st.snapshot();
      st.doc = p.doc;
      st.planner.starts.clear();
      for (const std::string& s : p.plannerStarts) {
        chemcad::ui::MaterialBox b;
        b.smiles = s;
        st.planner.starts.push_back(b);
      }
      if (st.planner.starts.empty()) st.planner.starts.push_back({});
      st.planner.target = {};
      st.planner.target.smiles = p.plannerTarget;
      st.projectPath = path;
      st.dirty = false;
      st.touch();
      st.statusMessage = "Opened " + path;
    } catch (const std::exception& e) {
      st.statusMessage = std::string("Open failed: ") + e.what();
    }
  };
}

// Builds the default dock layout once, then leaves the user in control.
void buildDefaultLayout(ImGuiID dockspaceId) {
  ImGui::DockBuilderRemoveNode(dockspaceId);
  ImGui::DockBuilderAddNode(dockspaceId, ImGuiDockNodeFlags_DockSpace);
  ImGui::DockBuilderSetNodeSize(dockspaceId, ImGui::GetMainViewport()->WorkSize);

  ImGuiID center = dockspaceId;
  // The labelled structure-tool rail needs enough width for command names,
  // current variants, shortcuts and split-menu affordances. The versioned
  // dockspace below migrates older two-column layouts into this safer default.
  const ImGuiID left = ImGui::DockBuilderSplitNode(center, ImGuiDir_Left, 0.16f, nullptr, &center);
  const ImGuiID right = ImGui::DockBuilderSplitNode(center, ImGuiDir_Right, 0.27f, nullptr, &center);
  ImGuiID rightBottom = right;
  const ImGuiID rightTop =
      ImGui::DockBuilderSplitNode(right, ImGuiDir_Up, 0.52f, nullptr, &rightBottom);

  ImGui::DockBuilderDockWindow(kWinTools, left);
  ImGui::DockBuilderDockWindow(kWinSketch, center);
  ImGui::DockBuilderDockWindow(kWinPlanner, center);
  ImGui::DockBuilderDockWindow(kWinSolubility, center);
  ImGui::DockBuilderDockWindow(kWinExtraction, center);
  ImGui::DockBuilderDockWindow(kWinViewer3D, center);
  ImGui::DockBuilderDockWindow(kWinPTable, rightTop);
  ImGui::DockBuilderDockWindow(kWinProps, rightBottom);
  ImGui::DockBuilderFinish(dockspaceId);
}

}  // namespace

int main(int, char**) {
  glfwSetErrorCallback(glfwErrorCallback);
  if (!glfwInit()) {
    std::fprintf(stderr, "failed to initialise GLFW\n");
    return 1;
  }
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);

  GLFWwindow* window = glfwCreateWindow(1600, 1000, "ChemCAD", nullptr, nullptr);
  if (!window) {
    std::fprintf(stderr, "failed to create a window\n");
    glfwTerminate();
    return 1;
  }
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);

  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO();
  io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
  io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
  static const std::string iniPath = (chemcad::core::cacheDir() / "imgui.ini").string();
  io.IniFilename = iniPath.c_str();

  // UI scale: default 1.25, overridable for accessibility or dense displays.
  float uiScale = 1.25f;
  if (const char* env = std::getenv("CHEMCAD_UI_SCALE"); env && *env) {
    if (const float parsed = std::strtof(env, nullptr); parsed >= 0.5f && parsed <= 3.0f)
      uiScale = parsed;
  }
  chemcad::ui::applyTheme(uiScale);
  chemcad::ui::style::fonts::load();
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init("#version 330");

  AppState state;
  wireCallbacks(state);

  bool layoutBuilt = false;
  while (!glfwWindowShouldClose(window)) {
    glfwPollEvents();
    // Skip rendering only when there is genuinely no surface to draw into.
    // GLFW_ICONIFIED is unreliable without a window manager (bare X servers and
    // headless sessions report it spuriously), which would leave the window
    // permanently blank.
    int fbw = 0, fbh = 0;
    glfwGetFramebufferSize(window, &fbw, &fbh);
    if (fbw <= 0 || fbh <= 0) {
      glfwWaitEventsTimeout(0.1);
      continue;
    }

    state.tasks.pump();

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    // ---- host window holding the menu bar and the dockspace
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(vp->WorkPos);
    ImGui::SetNextWindowSize(vp->WorkSize);
    ImGui::SetNextWindowViewport(vp->ID);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    ImGui::Begin("##chemcad_host", nullptr,
                 ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking |
                     ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse |
                     ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                     ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus);
    ImGui::PopStyleVar(3);

    chemcad::ui::drawMenuBar(state);

    const ImGuiID dockspaceId = ImGui::GetID("chemcad_dockspace_v2");
    if (!layoutBuilt) {
      layoutBuilt = true;
      if (!ImGui::DockBuilderGetNode(dockspaceId)) buildDefaultLayout(dockspaceId);
    }
    ImGui::DockSpace(dockspaceId, ImVec2(0, 0), ImGuiDockNodeFlags_None);
    ImGui::End();

    // ---- panels
    if (ImGui::Begin(kWinTools)) chemcad::ui::drawToolPalette(state);
    ImGui::End();

    // A panel asking for a tab switch sets `tab` + `tabChangeRequested`;
    // focus whichever window that tab lives in.
    const auto focusTabIfRequested = [&state](chemcad::ui::MainTab t) {
      if (state.tabChangeRequested && state.tab == t) {
        ImGui::SetNextWindowFocus();
        state.tabChangeRequested = false;
      }
    };

    focusTabIfRequested(chemcad::ui::MainTab::Sketch);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    const bool sketchOpen = ImGui::Begin(kWinSketch);
    ImGui::PopStyleVar();
    if (sketchOpen) {
      if (ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows))
        state.tab = chemcad::ui::MainTab::Sketch;
      chemcad::ui::drawCanvas(state);
    }
    ImGui::End();

    focusTabIfRequested(chemcad::ui::MainTab::Planner);
    if (ImGui::Begin(kWinPlanner)) {
      if (ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows))
        state.tab = chemcad::ui::MainTab::Planner;
      chemcad::ui::drawReactionPlanner(state);
    }
    ImGui::End();

    focusTabIfRequested(chemcad::ui::MainTab::Solubility);
    if (ImGui::Begin(kWinSolubility)) {
      if (ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows))
        state.tab = chemcad::ui::MainTab::Solubility;
      chemcad::ui::drawSolubilitySuite(state);
    }
    ImGui::End();

    focusTabIfRequested(chemcad::ui::MainTab::Extraction);
    if (ImGui::Begin(kWinExtraction)) {
      if (ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows))
        state.tab = chemcad::ui::MainTab::Extraction;
      chemcad::ui::drawExtractionLab(state);
    }
    ImGui::End();

    focusTabIfRequested(chemcad::ui::MainTab::Viewer3D);
    if (ImGui::Begin(kWinViewer3D)) {
      if (ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows))
        state.tab = chemcad::ui::MainTab::Viewer3D;
      chemcad::ui::drawViewer3D(state);
    }
    ImGui::End();

    if (ImGui::Begin(kWinPTable)) chemcad::ui::drawPeriodicTable(state);
    ImGui::End();

    if (ImGui::Begin(kWinProps)) chemcad::ui::drawPropertiesPanel(state);
    ImGui::End();

    chemcad::ui::drawStatusBar(state);

    // ---- render
    ImGui::Render();
    glfwGetFramebufferSize(window, &fbw, &fbh);
    glViewport(0, 0, fbw, fbh);
    glClearColor(0.043f, 0.055f, 0.075f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    // PNG export must sample the framebuffer before the swap.
    if (state.pendingPngExport.has_value()) {
      const std::string path = *state.pendingPngExport;
      state.pendingPngExport.reset();
      std::string err;
      const float dpiX = fbw / std::max(1.0f, vp->Size.x);
      const float dpiY = fbh / std::max(1.0f, vp->Size.y);
      const bool ok = chemcad::app::capturePng(
          static_cast<int>(state.canvasOrigin.x * dpiX),
          static_cast<int>(state.canvasOrigin.y * dpiY),
          static_cast<int>(state.canvasSize.x * dpiX),
          static_cast<int>(state.canvasSize.y * dpiY), fbh, path, err);
      state.statusMessage = ok ? ("Exported PNG: " + path) : ("PNG export failed: " + err);
    }

    glfwSwapBuffers(window);
  }

  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();
  glfwDestroyWindow(window);
  glfwTerminate();
  return 0;
}

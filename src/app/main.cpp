// ChemCAD -- entry point, GLFW/OpenGL/ImGui bootstrap and the dock layout.

#include <algorithm>
#include <cmath>
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

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#define GLFW_EXPOSE_NATIVE_WIN32
#endif

// Our own GL loader (gfx/gl_api.hpp) declares the GL 3.3 enums and entry
// points this translation unit needs. GLFW would otherwise drag in the
// system's GL 1.1 header, whose enum MACROS collide with those declarations,
// so the platform header is suppressed here and the three legacy calls below
// go through the loader like everything else.
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#ifdef _WIN32
#include <GLFW/glfw3native.h>
#endif

#include "app/project_io.hpp"
#include "app/screenshot.hpp"
#include "chem/bridge.hpp"
#include "core/paths.hpp"
#include "core/profiler.hpp"
#include "gfx/fluid_accelerator.hpp"
#include "gfx/fluid_renderer.hpp"
#include "gfx/gl_api.hpp"
#include "ui/app_state.hpp"
#include "ui/display_scale.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"

namespace {

using chemcad::ui::AppState;

constexpr const char* kWinSketch = "Sketch";
constexpr const char* kWinPlanner = "Reaction Planner";
constexpr const char* kWinSolubility = "Solubility Suite";
constexpr const char* kWinExtraction = "Extraction Calculator";
constexpr const char* kWinToolbox = "Toolbox";
constexpr const char* kWinPreview3D = "Preview";
constexpr const char* kWinTools = "Tools";
constexpr const char* kWinPTable = "Periodic Table";
constexpr const char* kWinProps = "Properties";

// The bench was designed at 1.25: text and controls that read comfortably at
// arm's length on a 1080p panel at 100% OS zoom. CHEMCAD_UI_SCALE replaces it;
// the monitor's content scale multiplies it.
constexpr float kDefaultUserScale = 1.25f;


void glfwErrorCallback(int code, const char* description) {
  std::fprintf(stderr, "[glfw] error %d: %s\n", code, description);
}

// ------------------------------------------------------- integrated chrome
// The app renders its own title bar (menu bar + caption buttons), so the OS
// frame is disabled and window management is implemented here.

std::pair<double, double> screenCursorPos(GLFWwindow* window) {
#ifdef _WIN32
  (void)window;
  POINT p;
  GetCursorPos(&p);
  return {static_cast<double>(p.x), static_cast<double>(p.y)};
#else
  double cx = 0.0, cy = 0.0;
  glfwGetCursorPos(window, &cx, &cy);
  int wx = 0, wy = 0;
  glfwGetWindowPos(window, &wx, &wy);
  return {wx + cx, wy + cy};
#endif
}

struct ChromeState {
  bool dragging = false;
  double grabX = 0.0, grabY = 0.0;  // screen cursor minus window pos at grab
  int resizeEdges = 0;              // 1 left, 2 right, 4 top, 8 bottom
  double startCursorX = 0.0, startCursorY = 0.0;
  int startX = 0, startY = 0, startW = 0, startH = 0;
};

// Drives window dragging (title bar) and edge resizing for the borderless
// window. Call once per frame after ImGui content is submitted.
void pumpWindowChrome(GLFWwindow* window, const AppState& st, ChromeState& chrome) {
  const ImGuiIO& io = ImGui::GetIO();
  const auto [cursorX, cursorY] = screenCursorPos(window);
  int wx = 0, wy = 0, ww = 0, wh = 0;
  glfwGetWindowPos(window, &wx, &wy);
  glfwGetWindowSize(window, &ww, &wh);
  const bool maximized = glfwGetWindowAttrib(window, GLFW_MAXIMIZED) != 0;
  // The cursor is tracked in DESKTOP coordinates (a drag moves the window, so
  // window-relative coordinates would chase themselves), while the zones ImGui
  // reports are CLIENT coordinates. Convert once, here.
  const double localX = cursorX - wx;
  const double localY = cursorY - wy;

  // ---- title-bar drag
  const auto& zone = st.titleDragZone;
  const bool inDragZone = zone.x2 > zone.x1 && localX >= zone.x1 && localX < zone.x2 &&
                          localY >= zone.y1 && localY < zone.y2;
  if (!chrome.dragging && chrome.resizeEdges == 0 && inDragZone && !maximized &&
      ImGui::IsMouseClicked(ImGuiMouseButton_Left) && !ImGui::IsPopupOpen(nullptr, ImGuiPopupFlags_AnyPopup)) {
    chrome.dragging = true;
    chrome.grabX = localX;
    chrome.grabY = localY;
  }
  if (chrome.dragging) {
    if (!io.MouseDown[ImGuiMouseButton_Left]) {
      chrome.dragging = false;
    } else {
      glfwSetWindowPos(window, static_cast<int>(cursorX - chrome.grabX),
                       static_cast<int>(cursorY - chrome.grabY));
    }
    return;  // dragging and resizing never overlap
  }

  // ---- edge resize
  // A 6 px grab band is a 3 px band at 200% display zoom -- unhittable. Every
  // chrome length is a design pixel.
  const double edgeBand = chemcad::ui::dp(6.0f);
  int edges = 0;
  if (!maximized) {
    if (localX >= 0 && localX < ww && localY >= 0 && localY < wh) {
      if (localX < edgeBand) edges |= 1;
      if (localX >= ww - edgeBand) edges |= 2;
      if (localY < edgeBand) edges |= 4;
      if (localY >= wh - edgeBand) edges |= 8;
    }
  }
  static GLFWcursor* hCursor = nullptr;
  static GLFWcursor* vCursor = nullptr;
  if (!hCursor) hCursor = glfwCreateStandardCursor(GLFW_HRESIZE_CURSOR);
  if (!vCursor) vCursor = glfwCreateStandardCursor(GLFW_VRESIZE_CURSOR);
  if (chrome.resizeEdges == 0) {
    const bool horizontal = (edges & 3) != 0;
    const bool vertical = (edges & 12) != 0;
    glfwSetCursor(window, horizontal && !vertical ? hCursor
                              : vertical && !horizontal ? vCursor
                                                        : nullptr);
  }

  if (chrome.resizeEdges == 0 && edges != 0 && ImGui::IsMouseClicked(ImGuiMouseButton_Left) &&
      !ImGui::IsAnyItemHovered()) {
    chrome.resizeEdges = edges;
    chrome.startCursorX = cursorX;
    chrome.startCursorY = cursorY;
    chrome.startX = wx;
    chrome.startY = wy;
    chrome.startW = ww;
    chrome.startH = wh;
  }
  if (chrome.resizeEdges != 0) {
    if (!io.MouseDown[ImGuiMouseButton_Left]) {
      chrome.resizeEdges = 0;
      glfwSetCursor(window, nullptr);
      return;
    }
    const double dx = cursorX - chrome.startCursorX;
    const double dy = cursorY - chrome.startCursorY;
    int nx = chrome.startX, ny = chrome.startY, nw = chrome.startW, nh = chrome.startH;
    if (chrome.resizeEdges & 1) { nx += static_cast<int>(dx); nw -= static_cast<int>(dx); }
    if (chrome.resizeEdges & 2) { nw += static_cast<int>(dx); }
    if (chrome.resizeEdges & 4) { ny += static_cast<int>(dy); nh -= static_cast<int>(dy); }
    if (chrome.resizeEdges & 8) { nh += static_cast<int>(dy); }
    nw = std::max(nw, static_cast<int>(chemcad::ui::dp(
                          static_cast<float>(chemcad::ui::kMinWindowWidth))));
    nh = std::max(nh, static_cast<int>(chemcad::ui::dp(
                          static_cast<float>(chemcad::ui::kMinWindowHeight))));
    glfwSetWindowPos(window, nx, ny);
    glfwSetWindowSize(window, nw, nh);
  }
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

// Builds the dock layout for the active workspace tab. Each tab gets only the
// panels that are useful there -- the sketch keeps the full bench (tools,
// periodic table, preview, properties), the planner and suite keep the
// molecule panels, and the extraction calculator gets the whole bench top.
//
// The splits are budgeted in EMS, not in fractions of the viewport. A fraction
// is only correct at one display scale: 30% of a 1080p viewport is a roomy
// panel at 100% zoom and an unreadably narrow one at 200%. When the ems are not
// there, the side panels are docked into the centre node instead, where they
// become tabs beside the workspaces -- every panel stays reachable at full
// width, and nothing is squeezed until it clips.
void buildLayoutForTab(ImGuiID dockspaceId, chemcad::ui::MainTab tab) {
  using chemcad::ui::MainTab;
  // Window settings persist dock-node ids across sessions and across tab
  // layouts; without clearing them, panels keep pointing at nodes from an
  // older layout and render into nothing. The builder below IS the layout,
  // so stale settings are always safe to drop here.
  ImGui::ClearIniSettings();
  ImGui::DockBuilderRemoveNode(dockspaceId);
  ImGui::DockBuilderAddNode(dockspaceId, ImGuiDockNodeFlags_DockSpace);
  const ImVec2 work = ImGui::GetMainViewport()->WorkSize;
  ImGui::DockBuilderSetNodeSize(dockspaceId, work);

  // What each column needs to hold its contents: the tool palette is one
  // button column, the right stack has to fit the periodic table's 18 groups,
  // and the workspace itself is the point of the window.
  const float em = ImGui::GetFontSize();
  const float toolsEm = 14.0f;
  const float stackEm = 28.0f;
  const float workspaceEm = 34.0f;
  const float widthEm = work.x / std::max(em, 1.0f);
  const float heightEm = work.y / std::max(em, 1.0f);
  // A short viewport cannot stack three panels either, however wide it is.
  const bool roomForStack = widthEm >= stackEm + workspaceEm && heightEm >= 30.0f;
  const bool roomForTools = widthEm >= stackEm + workspaceEm + toolsEm;

  ImGuiID center = dockspaceId;

  if (tab == MainTab::Sketch) {
    ImGuiID tools = center;
    if (roomForTools) {
      tools = ImGui::DockBuilderSplitNode(center, ImGuiDir_Left, toolsEm * em / work.x, nullptr,
                                          &center);
    }
    if (roomForStack) {
      ImGuiID rightBottom = ImGui::DockBuilderSplitNode(
          center, ImGuiDir_Right, stackEm * em / std::max(work.x - toolsEm * em, 1.0f), nullptr,
          &center);
      const ImGuiID rightTop =
          ImGui::DockBuilderSplitNode(rightBottom, ImGuiDir_Up, 0.44f, nullptr, &rightBottom);
      ImGuiID rightProps = rightBottom;
      const ImGuiID rightMid =
          ImGui::DockBuilderSplitNode(rightBottom, ImGuiDir_Up, 0.50f, nullptr, &rightProps);
      ImGui::DockBuilderDockWindow(kWinPTable, rightTop);
      ImGui::DockBuilderDockWindow(kWinPreview3D, rightMid);
      ImGui::DockBuilderDockWindow(kWinProps, rightProps);
    } else {
      ImGui::DockBuilderDockWindow(kWinPTable, center);
      ImGui::DockBuilderDockWindow(kWinPreview3D, center);
      ImGui::DockBuilderDockWindow(kWinProps, center);
    }
    ImGui::DockBuilderDockWindow(kWinTools, roomForTools ? tools : center);
  } else if (tab == MainTab::Planner || tab == MainTab::Solubility) {
    if (roomForStack) {
      ImGuiID rightBottom = ImGui::DockBuilderSplitNode(center, ImGuiDir_Right,
                                                       stackEm * em / work.x, nullptr, &center);
      const ImGuiID rightTop =
          ImGui::DockBuilderSplitNode(rightBottom, ImGuiDir_Up, 0.52f, nullptr, &rightBottom);
      ImGui::DockBuilderDockWindow(kWinPreview3D, rightTop);
      ImGui::DockBuilderDockWindow(kWinProps, rightBottom);
    } else {
      ImGui::DockBuilderDockWindow(kWinPreview3D, center);
      ImGui::DockBuilderDockWindow(kWinProps, center);
    }
  }
  // Extraction: full-width workspace, no side panels.

  ImGui::DockBuilderDockWindow(kWinSketch, center);
  ImGui::DockBuilderDockWindow(kWinPlanner, center);
  ImGui::DockBuilderDockWindow(kWinSolubility, center);
  ImGui::DockBuilderDockWindow(kWinExtraction, center);
  ImGui::DockBuilderDockWindow(kWinToolbox, center);
  ImGui::DockBuilderFinish(dockspaceId);
}


// Renders the 3D fluid the extraction panel asked for into the renderer's own
// framebuffer and publishes the resulting texture for the NEXT frame. Called
// from the render seam with the GL context current; the panel itself never
// touches GL, which is what keeps the headless panel tests able to run it.
void renderFluidStage(AppState& state, chemcad::gfx::FluidRenderer& renderer) {
  chemcad::gfx::FluidStage& stage = state.fluidStage;
  stage.available = renderer.ready();
  if (!renderer.ready()) {
    stage.status = renderer.error().empty() ? "3D renderer unavailable" : renderer.error();
    stage.texture = 0;
    stage.requested = false;
    return;
  }
  if (!stage.requested || stage.width <= 0 || stage.height <= 0 || !stage.snapshot) {
    stage.requested = false;
    return;
  }
  const std::uint32_t texture = renderer.render(*stage.snapshot, stage.pose, stage.camera,
                                                stage.width, stage.height, stage.settings);
  stage.texture = texture;
  stage.textureWidth = renderer.width();
  stage.textureHeight = renderer.height();
  stage.lastFrameMs = renderer.lastFrameMs();
  char status[192];
  std::snprintf(status, sizeof(status), "%s | %s | %.1f ms",
                chemcad::gfx::glVersionString(), chemcad::gfx::glRendererString(),
                stage.lastFrameMs);
  stage.status = status;
  stage.requested = false;
}

// ------------------------------------------------------------ display scale
// The scale the bench should be drawn at RIGHT NOW, for the monitor this window
// is on. Two window-system facts go into it:
//
//   * the content scale -- Windows display zoom, GNOME's text-scaling-factor,
//     the macOS backing factor;
//   * the framebuffer scale, framebuffer pixels per window unit. On Windows
//     that is 1 and the content scale is ours to apply; on macOS and Wayland
//     the window system has already applied it to the framebuffer, and scaling
//     again would double-count.
//
// Queried every frame: it changes when the user drags the window to another
// monitor or changes the OS zoom, and GLFW gives no promise about which of
// those emits a callback on which platform.
float currentDisplayScale(GLFWwindow* window, float userScale) {
  float contentX = 1.0f, contentY = 1.0f;
  glfwGetWindowContentScale(window, &contentX, &contentY);
  int windowW = 0, windowH = 0, framebufferW = 0, framebufferH = 0;
  glfwGetWindowSize(window, &windowW, &windowH);
  glfwGetFramebufferSize(window, &framebufferW, &framebufferH);
  const float framebufferScale =
      windowW > 0 ? static_cast<float>(framebufferW) / static_cast<float>(windowW) : 1.0f;
  return chemcad::ui::resolveDisplayScale(std::min(contentX, contentY), framebufferScale,
                                          userScale);
}
}  // namespace

int main(int, char**) {
  glfwSetErrorCallback(glfwErrorCallback);
  if (!glfwInit()) {
    std::fprintf(stderr, "failed to initialise GLFW\n");
    return 1;
  }
  // The compute fluid backend needs OpenGL 4.3; the renderer and ImGui need
  // 3.3. Ask for 4.3 and fall back, so an older GPU still gets an application
  // instead of an error dialog -- it simply keeps the CPU solver.
  int contextVersion = 43;
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);
  glfwWindowHint(GLFW_DECORATED, GLFW_FALSE);  // the app draws its own chrome
  // Placed explicitly below, once the monitor's scale and work area are known;
  // showing it first would flash a wrongly-sized window.
  glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);

  // Provisional size only. The app draws its own title bar and caption buttons,
  // so a window larger than the screen is not merely awkward, it is unclosable;
  // fitWindow() below is what makes the first frame safe on a 1366x768 laptop
  // and on a 4K TV at 200% zoom alike.
  GLFWwindow* window = glfwCreateWindow(chemcad::ui::kBaseWindowWidth,
                                        chemcad::ui::kBaseWindowHeight, "ChemCAD", nullptr,
                                        nullptr);
  if (!window) {
    contextVersion = 33;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    window = glfwCreateWindow(chemcad::ui::kBaseWindowWidth, chemcad::ui::kBaseWindowHeight,
                              "ChemCAD", nullptr, nullptr);
  }
  glfwWindowHint(GLFW_VISIBLE, GLFW_TRUE);
  if (!window) {
    std::fprintf(stderr, "failed to create a window\n");
    glfwTerminate();
    return 1;
  }
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);

  // ---- place and size the window for the monitor it actually opened on
  const float userScale =
      chemcad::ui::parseUserScale(std::getenv("CHEMCAD_UI_SCALE"), kDefaultUserScale);
  {
    const float initialScale = currentDisplayScale(window, userScale);
    GLFWmonitor* monitor = glfwGetPrimaryMonitor();
    int workX = 0, workY = 0, workW = 0, workH = 0;
    if (monitor != nullptr) glfwGetMonitorWorkarea(monitor, &workX, &workY, &workW, &workH);
    const chemcad::ui::WindowPlacement placement =
        chemcad::ui::fitWindow(workX, workY, workW, workH, chemcad::ui::kBaseWindowWidth,
                               chemcad::ui::kBaseWindowHeight, initialScale);
    glfwSetWindowPos(window, placement.x, placement.y);
    glfwSetWindowSize(window, placement.width, placement.height);
    glfwShowWindow(window);
  }

  // The fluid renderer needs modern GL entry points; ImGui's backend keeps its
  // own private loader and its header says plainly that the rest of the app
  // must use a different one, so this resolves ours through GLFW. A failure
  // here is not fatal: the extraction panel falls back to the 2D schematic.
  chemcad::gfx::FluidRenderer fluidRenderer;
  // CHEMCAD_NO_FLUID_GL=1 keeps the loader out of the picture entirely, which is
  // the first thing to try when diagnosing a driver-specific rendering problem.
  const char* disableFluidGl = std::getenv("CHEMCAD_NO_FLUID_GL");
  const bool fluidGlAllowed = !(disableFluidGl && *disableFluidGl && *disableFluidGl != '0');
  if (fluidGlAllowed &&
      chemcad::gfx::loadGl(reinterpret_cast<chemcad::gfx::GlProcLoader>(glfwGetProcAddress))) {
    fluidRenderer.initialise();
  }

  // A second, invisible context so the physics worker can run compute shaders.
  // A GL context belongs to one thread at a time and the renderer keeps this
  // one. Nothing is shared: the compute backend's buffers are private and its
  // results come back as ordinary memory, so an independent context avoids the
  // cross-thread object-lifetime rules that sharing would impose.
  //
  // CHEMCAD_NO_FLUID_GPU=1 keeps the solve on the CPU, which is the first thing
  // to try when diagnosing a driver-specific physics problem.
  const char* disableFluidGpu = std::getenv("CHEMCAD_NO_FLUID_GPU");
  const bool fluidGpuAllowed = !(disableFluidGpu && *disableFluidGpu && *disableFluidGpu != '0');
  GLFWwindow* physicsContext = nullptr;
  if (fluidGpuAllowed && contextVersion >= 43 && chemcad::gfx::glLoaded()) {
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
    physicsContext = glfwCreateWindow(16, 16, "ChemCAD physics", nullptr, nullptr);
    glfwWindowHint(GLFW_VISIBLE, GLFW_TRUE);
  }

#ifdef _WIN32
  // Belt and braces: undecorated windows occasionally composite nothing until
  // an explicit SW_SHOW + repaint (observed on Win10 when the shell restored
  // a stale show state). Cheap and idempotent.
  if (HWND hwnd = glfwGetWin32Window(window)) {
    ShowWindow(hwnd, SW_SHOW);
    UpdateWindow(hwnd);
  }
#endif

  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO();
  io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
  io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
  static const std::string iniPath = (chemcad::core::cacheDir() / "imgui.ini").string();
  io.IniFilename = iniPath.c_str();

  // One factor for the whole bench: the monitor's content scale times the
  // user's own CHEMCAD_UI_SCALE. Re-resolved after the window is placed, since
  // placing it may have moved it onto a different monitor.
  float appliedScale = currentDisplayScale(window, userScale);
  chemcad::ui::applyTheme(appliedScale);
  chemcad::ui::style::fonts::load();
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init("#version 330");

  AppState state;
  wireCallbacks(state);
  if (physicsContext != nullptr) {
    // The callbacks run on the physics worker thread. glfwMakeContextCurrent is
    // the one GLFW entry point documented as callable from any thread, which is
    // exactly why the context handoff is a callback rather than a handle.
    state.solubility.fluidAccelerator = chemcad::gfx::makeFluidAccelerator(
        [physicsContext] {
          glfwMakeContextCurrent(physicsContext);
          return glfwGetCurrentContext() == physicsContext;
        },
        [] { glfwMakeContextCurrent(nullptr); });
  }
  // The extraction workspace needs a particle simulation whose vessel SDF and
  // first charge are more work than a frame can absorb. Starting it here puts
  // that on a worker thread while the user is still on the sketch canvas, so
  // the workspace is ready when opened.
  chemcad::ui::warmExtractionPhysics(state);

  ChromeState chrome;
  chemcad::ui::MainTab layoutTab = state.tab;
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

    // The display zoom or the monitor under the window can change at any time.
    // Rebuilding the theme is the whole response: every metric, every font size
    // and the sketch canvas' pixels-per-bond are derived from this one factor.
    // Done here, between frames -- ImGui forbids restyling mid-frame.
    if (const float wanted = currentDisplayScale(window, userScale);
        std::fabs(wanted - appliedScale) > 0.001f) {
      appliedScale = wanted;
      chemcad::ui::applyTheme(appliedScale);
      layoutBuilt = false;  // dock split sizes are in pixels; re-split at the new scale
    }

    chemcad::core::profiler().beginFrame();
    {
      CHEMCAD_PROFILE_ZONE("tasks.pump");
      state.tasks.pump();
    }

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    // ---- host window holding the menu bar and the dockspace
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(vp->WorkPos);
    ImGui::SetNextWindowSize(vp->WorkSize);
    ImGui::SetNextWindowViewport(vp->ID);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 1.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    ImGui::Begin("##chemcad_host", nullptr,
                 ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking |
                     ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse |
                     ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                     ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus);
    ImGui::PopStyleVar(3);

    chemcad::ui::drawMenuBar(state);

    const ImGuiID dockspaceId = ImGui::GetID("chemcad_dockspace_v4");
    if (!layoutBuilt || layoutTab != state.tab) {
      layoutTab = state.tab;
      layoutBuilt = true;
      buildLayoutForTab(dockspaceId, state.tab);
      // A rebuilt dock tree resets the centre node's visible tab to the first
      // docked window; pull the user's actual tab back to the front.
      state.tabChangeRequested = true;
    }
    ImGui::DockSpace(dockspaceId, ImVec2(0, 0), ImGuiDockNodeFlags_None);
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

    focusTabIfRequested(chemcad::ui::MainTab::Toolbox);
    if (ImGui::Begin(kWinToolbox)) {
      if (ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows))
        state.tab = chemcad::ui::MainTab::Toolbox;
      chemcad::ui::drawToolbox(state);
    }
    ImGui::End();

    // Side panels are workspace-specific: only the panels the active tab can
    // actually use are submitted, so nothing inapplicable clutters the bench.
    const bool moleculePanels = state.tab != chemcad::ui::MainTab::Extraction &&
                                state.tab != chemcad::ui::MainTab::Toolbox;
    if (state.tab == chemcad::ui::MainTab::Sketch) {
      if (ImGui::Begin(kWinTools)) chemcad::ui::drawToolPalette(state);
      ImGui::End();
      if (ImGui::Begin(kWinPTable)) chemcad::ui::drawPeriodicTable(state);
      ImGui::End();
    }
    if (moleculePanels) {
      if (ImGui::Begin(kWinPreview3D)) chemcad::ui::drawViewer3D(state);
      ImGui::End();
      if (ImGui::Begin(kWinProps)) chemcad::ui::drawPropertiesPanel(state);
      ImGui::End();
    }

    chemcad::ui::drawStatusBar(state);

    // A floating developer surface, deliberately outside the dockspace: the
    // profiler measures the bench, it is not part of it.
    if (state.showProfiler) {
      CHEMCAD_PROFILE_ZONE("ui.profiler");
      ImGui::SetNextWindowSize(ImVec2(vp->Size.x * 0.55f, vp->Size.y * 0.55f),
                               ImGuiCond_FirstUseEver);
      if (ImGui::Begin("Performance", &state.showProfiler)) {
        chemcad::ui::drawProfiler(state);
      }
      ImGui::End();
    }

    pumpWindowChrome(window, state, chrome);

    // ---- render
    {
      CHEMCAD_PROFILE_ZONE("imgui.render");
      ImGui::Render();
    }
    glfwGetFramebufferSize(window, &fbw, &fbh);
    if (chemcad::gfx::glLoaded()) {
      chemcad::gfx::glViewport(0, 0, fbw, fbh);
      chemcad::gfx::glClearColor(0.043f, 0.055f, 0.075f, 1.0f);
      chemcad::gfx::glClear(chemcad::gfx::GL_COLOR_BUFFER_BIT);
    }
    {
      CHEMCAD_PROFILE_ZONE("imgui.submit");
      ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    }

    // The 3D fluid is rendered here, after ImGui has recorded its draw data
    // and before the backend submits it: the panel published a request during
    // submission and drew the PREVIOUS frame's texture, so the stage is always
    // exactly one frame behind, which is invisible at 60 Hz and keeps every GL
    // call for the fluid on the main thread inside one seam.
    {
      CHEMCAD_PROFILE_ZONE("fluid.stage");
      renderFluidStage(state, fluidRenderer);
    }

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
    chemcad::core::profiler().endFrame();
  }

  // Before anything GLFW owns goes away. AppState is a stack local that
  // outlives glfwTerminate(), and its physics worker may still hold the second
  // context current on its own thread.
  chemcad::ui::shutdownExtractionPhysics(state);
  state.solubility.fluidAccelerator.reset();

  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();
  if (physicsContext != nullptr) glfwDestroyWindow(physicsContext);
  glfwDestroyWindow(window);
  glfwTerminate();
  return 0;
}

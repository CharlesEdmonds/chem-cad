#include <algorithm>
#include <filesystem>
#include <limits>
#include <string>
#include <unordered_map>

#include "imgui.h"

#include <GLFW/glfw3.h>

#include "chem/bridge.hpp"
#include "ui/file_dialog.hpp"
#include "ui/icons.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"
#include "ui/widgets.hpp"

namespace chemcad::ui {
namespace {

enum class DialogAction { None, OpenProject, SaveProject, ImportMol, ExportMol, ExportSvg, ExportPng };

core::Molecule flatten(const core::Document& document) {
  core::Molecule result;
  for (const auto& fragment : document.molecules) {
    std::unordered_map<core::AtomId, core::AtomId> ids;
    ids.reserve(fragment.atomCount());
    for (const auto& atom : fragment.atoms()) ids.emplace(atom.id, result.addAtom(atom));
    for (const auto& bond : fragment.bonds()) {
      const core::BondId id = result.addBond(ids.at(bond.a), ids.at(bond.b), bond.order);
      if (core::Bond* added = result.bond(id)) added->stereo = bond.stereo;
    }
  }
  return result;
}

void offsetForPaste(const core::Document& document, core::Molecule& molecule) {
  if (molecule.empty()) return;
  bool hasExisting = false;
  float maxX = -std::numeric_limits<float>::infinity();
  for (const auto& fragment : document.molecules) {
    for (const auto& atom : fragment.atoms()) {
      maxX = std::max(maxX, atom.pos.x);
      hasExisting = true;
    }
  }
  if (!hasExisting) return;
  float minX = std::numeric_limits<float>::infinity();
  for (const auto& atom : molecule.atoms()) minX = std::min(minX, atom.pos.x);
  for (auto& atom : molecule.mutableAtoms()) atom.pos.x += maxX - minX + 1.5f;
}

void newDocument(AppState& st) {
  if (st.doc.empty() && st.projectPath.empty()) return;
  st.snapshot();
  st.doc.clear();
  st.sel.clear();
  st.projectPath.clear();
  st.touch();
  st.statusMessage = "New document";
}

void undo(AppState& st) {
  if (st.undo.undo(st.doc)) {
    st.sel.clear();
    st.touch();
    st.statusMessage = "Undo";
  }
}

void redo(AppState& st) {
  if (st.undo.redo(st.doc)) {
    st.sel.clear();
    st.touch();
    st.statusMessage = "Redo";
  }
}

void copySmiles(AppState& st) {
  try {
    core::Molecule molecule = flatten(st.doc);
    if (molecule.empty()) {
      st.statusMessage = "Nothing to copy";
      return;
    }
    st.clipboardSmiles = chem::toSmiles(molecule);
    ImGui::SetClipboardText(st.clipboardSmiles.c_str());
    st.statusMessage = "SMILES copied";
  } catch (const std::exception& e) {
    st.statusMessage = std::string("Copy failed: ") + e.what();
  }
}

void pasteSmiles(AppState& st) {
  std::string smiles = st.clipboardSmiles;
  if (const char* system = ImGui::GetClipboardText(); system && *system) smiles = system;
  if (smiles.empty()) {
    st.statusMessage = "Clipboard does not contain SMILES";
    return;
  }
  try {
    core::Molecule molecule = chem::fromSmiles(smiles);
    offsetForPaste(st.doc, molecule);
    st.snapshot();
    st.doc.molecules.push_back(std::move(molecule));
    st.sel.clear();
    st.touch();
    st.clipboardSmiles = std::move(smiles);
    st.statusMessage = "Pasted SMILES";
  } catch (const std::exception& e) {
    st.statusMessage = std::string("Paste failed: ") + e.what();
  }
}

void cleanUp(AppState& st) {
  if (st.doc.empty()) return;
  st.snapshot();
  try {
    core::Document cleaned = st.doc;
    for (auto& molecule : cleaned.molecules) chem::layout(molecule);
    st.doc = std::move(cleaned);
    st.touch();
    st.statusMessage = "Structure cleaned up";
  } catch (const std::exception& e) {
    st.statusMessage = std::string("Clean up failed: ") + e.what();
  }
}

void deleteSelected(AppState& st) {
  if (st.sel.empty()) return;
  st.snapshot();
  for (int molIndex = static_cast<int>(st.doc.molecules.size()) - 1; molIndex >= 0;
       --molIndex) {
    core::Molecule& molecule = st.doc.molecules[static_cast<size_t>(molIndex)];
    for (const BondRef& ref : st.sel.bonds) {
      if (ref.mol == molIndex) molecule.removeBond(ref.id);
    }
    for (const AtomRef& ref : st.sel.atoms) {
      if (ref.mol == molIndex) molecule.removeAtom(ref.id);
    }
  }
  std::erase_if(st.doc.molecules, [](const core::Molecule& molecule) {
    return molecule.empty();
  });
  st.sel.clear();
  st.hoverAtom = {};
  st.hoverBond = {};
  st.touch();
  st.statusMessage = "Selection deleted";
}

void clearStructure(AppState& st) {
  if (st.doc.empty()) return;
  st.snapshot();
  st.doc.clear();
  st.sel.clear();
  st.hoverAtom = {};
  st.hoverBond = {};
  st.touch();
  st.statusMessage = "Structure cleared";
}

std::string projectFilename(const AppState& st) {
  if (st.projectPath.empty()) return {};
  std::error_code ec;
  const auto name = std::filesystem::path(st.projectPath).filename();
  return ec ? std::string{} : name.string();
}

bool iconMenuItem(icons::Icon icon, const char* label, const char* shortcut,
                  bool selected = false, bool enabled = true, bool danger = false) {
  const style::Metrics& m = style::metrics();
  std::string padded = "   ";
  padded += label;
  if (danger) {
    ImGui::PushStyleColor(ImGuiCol_HeaderHovered, style::col::Danger);
    ImGui::PushStyleColor(ImGuiCol_HeaderActive, style::col::Danger);
  }
  const bool activated = ImGui::MenuItem(padded.c_str(), shortcut, selected, enabled);
  const bool hovered = ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled);
  if (danger) ImGui::PopStyleColor(2);

  const ImVec2 min = ImGui::GetItemRectMin();
  const ImVec2 max = ImGui::GetItemRectMax();
  const float iconSize = ImGui::GetFontSize() * 0.72f;
  const ImVec4 colour =
      !enabled ? style::col::TextFaint
               : (hovered && danger ? style::col::OnAccent : style::col::TextDim);
  icons::draw(ImGui::GetWindowDrawList(), icon,
              ImVec2(min.x + m.gap + iconSize * 0.5f, (min.y + max.y) * 0.5f),
              iconSize, style::u32(colour));
  return activated;
}

bool iconToggleMenuItem(icons::Icon icon, const char* label, const char* shortcut,
                        bool* selected) {
  const style::Metrics& m = style::metrics();
  std::string padded = "   ";
  padded += label;
  const bool activated = ImGui::MenuItem(padded.c_str(), shortcut, selected);
  const ImVec2 min = ImGui::GetItemRectMin();
  const ImVec2 max = ImGui::GetItemRectMax();
  const float iconSize = ImGui::GetFontSize() * 0.72f;
  icons::draw(ImGui::GetWindowDrawList(), icon,
              ImVec2(min.x + m.gap + iconSize * 0.5f, (min.y + max.y) * 0.5f),
              iconSize, style::u32(style::col::TextDim));
  return activated;
}

void openDialog(FileDialog& dialog, DialogAction& action, DialogAction requested,
                const AppState& st) {
  action = requested;
  switch (requested) {
    case DialogAction::OpenProject:
      dialog.open("Open project", FileDialogMode::Open, {".chemcad"});
      break;
    case DialogAction::SaveProject:
      dialog.open("Save project", FileDialogMode::Save, {".chemcad"}, ".chemcad",
                  projectFilename(st));
      break;
    case DialogAction::ImportMol:
      dialog.open("Import MOL", FileDialogMode::Open, {".mol", ".sdf"});
      break;
    case DialogAction::ExportMol:
      dialog.open("Export MOL", FileDialogMode::Save, {".mol"}, ".mol");
      break;
    case DialogAction::ExportSvg:
      dialog.open("Export SVG", FileDialogMode::Save, {".svg"}, ".svg");
      break;
    case DialogAction::ExportPng:
      dialog.open("Export PNG", FileDialogMode::Save, {".png"}, ".png");
      break;
    case DialogAction::None:
      break;
  }
}

void dispatchDialogResult(AppState& st, DialogAction action, const std::string& path) {
  const auto unavailable = [&st] { st.statusMessage = "That file action is not available."; };
  switch (action) {
    case DialogAction::OpenProject:
      st.openProject ? st.openProject(path) : unavailable();
      break;
    case DialogAction::SaveProject:
      st.saveProject ? st.saveProject(path) : unavailable();
      break;
    case DialogAction::ImportMol:
      st.importMol ? st.importMol(path) : unavailable();
      break;
    case DialogAction::ExportMol:
      st.exportMol ? st.exportMol(path) : unavailable();
      break;
    case DialogAction::ExportSvg:
      st.exportSvg ? st.exportSvg(path) : unavailable();
      break;
    case DialogAction::ExportPng:
      st.exportPng ? st.exportPng(path) : unavailable();
      break;
    case DialogAction::None:
      break;
  }
}

}  // namespace

void drawMenuBar(AppState& st) {
  static FileDialog dialog;
  static DialogAction dialogAction = DialogAction::None;

  if (ImGui::BeginMenuBar()) {
    // Brand mark: benzene-ring glyph + wordmark, then the menus.
    {
      const style::Metrics& m = style::metrics();
      const float h = ImGui::GetFrameHeight();
      const ImVec2 textSize = ImGui::CalcTextSize("ChemCAD");
      ImGui::Dummy(ImVec2(h + textSize.x + m.gap * 1.75f, h));
      const ImVec2 min = ImGui::GetItemRectMin();
      ImDrawList* dl = ImGui::GetWindowDrawList();
      icons::draw(dl, icons::Icon::Logo, ImVec2(min.x + h * 0.5f, min.y + h * 0.52f),
                  h * 0.68f, style::u32(style::col::Text), m.hairline * 1.6f);
      const bool pushed = style::pushFont(style::fonts::semibold());
      dl->AddText(ImVec2(min.x + h + m.gap * 0.75f, min.y + (h - textSize.y) * 0.5f),
                  style::u32(style::col::Text), "ChemCAD");
      style::popFont(pushed);
    }
    if (ImGui::BeginMenu("File")) {
      if (iconMenuItem(icons::Icon::Plus, "New", "Ctrl+N"))
        newDocument(st);
      if (iconMenuItem(icons::Icon::Folder, "Open...", "Ctrl+O")) {
        openDialog(dialog, dialogAction, DialogAction::OpenProject, st);
      }
      if (iconMenuItem(icons::Icon::Save, "Save", "Ctrl+S")) {
        if (!st.projectPath.empty() && st.saveProject) {
          st.saveProject(st.projectPath);
        } else {
          openDialog(dialog, dialogAction, DialogAction::SaveProject, st);
        }
      }
      if (iconMenuItem(icons::Icon::Save, "Save as...", "Ctrl+Shift+S")) {
        openDialog(dialog, dialogAction, DialogAction::SaveProject, st);
      }
      ImGui::Separator();
      if (iconMenuItem(icons::Icon::ArrowLeft, "Import MOL...", "Ctrl+I")) {
        openDialog(dialog, dialogAction, DialogAction::ImportMol, st);
      }
      if (iconMenuItem(icons::Icon::ArrowRight, "Export MOL...", "Ctrl+Shift+M", false,
                       !st.doc.empty())) {
        openDialog(dialog, dialogAction, DialogAction::ExportMol, st);
      }
      if (iconMenuItem(icons::Icon::ArrowRight, "Export SVG...", "Ctrl+Shift+G", false,
                       !st.doc.empty())) {
        openDialog(dialog, dialogAction, DialogAction::ExportSvg, st);
      }
      if (iconMenuItem(icons::Icon::ArrowRight, "Export PNG...", "Ctrl+Shift+P", false,
                       !st.doc.empty())) {
        openDialog(dialog, dialogAction, DialogAction::ExportPng, st);
      }
      ImGui::Separator();
      if (iconMenuItem(icons::Icon::Close, "Quit", "Ctrl+Q")) {
        if (GLFWwindow* window = glfwGetCurrentContext()) {
          glfwSetWindowShouldClose(window, GLFW_TRUE);
        }
      }
      ImGui::EndMenu();
    }

    if (ImGui::BeginMenu("Edit")) {
      if (iconMenuItem(icons::Icon::Undo, "Undo", "Ctrl+Z", false, st.undo.canUndo()))
        undo(st);
      if (iconMenuItem(icons::Icon::Redo, "Redo", "Ctrl+Shift+Z", false,
                       st.undo.canRedo()))
        redo(st);
      ImGui::Separator();
      if (iconMenuItem(icons::Icon::Copy, "Copy SMILES", "Ctrl+C", false,
                       !st.doc.empty()))
        copySmiles(st);
      if (iconMenuItem(icons::Icon::Link, "Paste SMILES", "Ctrl+V")) pasteSmiles(st);
      ImGui::Separator();
      if (iconMenuItem(icons::Icon::Trash, "Delete selection", "Del", false,
                       !st.sel.empty(), true))
        deleteSelected(st);
      ImGui::EndMenu();
    }

    if (ImGui::BeginMenu("Structure")) {
      if (iconMenuItem(icons::Icon::Sparkle, "Clean up structure", "Ctrl+L", false,
                       !st.doc.empty()))
        cleanUp(st);
      iconToggleMenuItem(icons::Icon::Book, "Auto-name", "Ctrl+Shift+A",
                         &st.props.autoName);
      ImGui::Separator();
      if (iconMenuItem(icons::Icon::Trash, "Clear structure", "Ctrl+Shift+Del", false,
                       !st.doc.empty(), true))
        clearStructure(st);
      ImGui::EndMenu();
    }

    if (ImGui::BeginMenu("View")) {
      if (iconMenuItem(icons::Icon::ZoomFit, "Fit to window", "Ctrl+F"))
        st.cam.fit(st.doc, st.canvasSize);
      if (iconMenuItem(icons::Icon::Crosshair, "Reset zoom", "Ctrl+0"))
        st.cam.zoom = 1.0f;
      ImGui::Separator();
      if (iconMenuItem(icons::Icon::Gauge, "Performance", "F2", st.showProfiler))
        st.showProfiler = !st.showProfiler;
      ImGui::EndMenu();
    }

    if (ImGui::BeginMenu("Help")) {
      if (iconMenuItem(icons::Icon::Info, "About ChemCAD", "F1"))
        ImGui::OpenPopup("About ChemCAD");
      ImGui::EndMenu();
    }

    // ---- integrated window chrome -------------------------------------
    // The stretch between the last menu and the caption buttons is the
    // window's drag handle; the app loop reads it from AppState.
    {
      const float barH = ImGui::GetFrameHeight();
      const ImVec2 winMin = ImGui::GetWindowPos();
      const float winW = ImGui::GetWindowWidth();
      const float btnW = barH * 1.55f;
      const ImVec2 afterMenus = ImGui::GetCursorScreenPos();
      st.titleDragZone = {afterMenus.x, winMin.y, winMin.x + winW - btnW * 3.0f,
                          winMin.y + barH};

      GLFWwindow* window = glfwGetCurrentContext();
      ImDrawList* dl = ImGui::GetWindowDrawList();
      const style::Metrics& m = style::metrics();
      const ImVec2 mouse = ImGui::GetMousePos();

      const auto captionButton = [&](int index, icons::Icon glyph, bool danger,
                                     const char* tooltip) -> bool {
        const ImVec2 bMin(winMin.x + winW - btnW * static_cast<float>(3 - index), winMin.y);
        const ImVec2 bMax(bMin.x + btnW, bMin.y + barH);
        const bool hovered = mouse.x >= bMin.x && mouse.x < bMax.x && mouse.y >= bMin.y &&
                             mouse.y < bMax.y;
        ImGui::SetCursorScreenPos(bMin);
        ImGui::PushID(index);
        ImGui::InvisibleButton("##caption", ImVec2(btnW, barH));
        const bool clicked = ImGui::IsItemClicked();
        ImGui::PopID();
        if (hovered) {
          dl->AddRectFilled(bMin, bMax,
                            danger ? style::u32(style::col::Danger, 0.85f)
                                   : style::u32(style::col::BgRaised));
          ImGui::SetTooltip("%s", tooltip);
        }
        icons::draw(dl, glyph, ImVec2((bMin.x + bMax.x) * 0.5f, (bMin.y + bMax.y) * 0.5f),
                    m.iconSize * 0.62f,
                    style::u32(hovered && danger ? style::col::OnAccent
                                                 : style::col::TextDim));
        return clicked;
      };

      if (captionButton(0, icons::Icon::Minus, false, "Minimize") && window) {
        glfwIconifyWindow(window);
      }
      const bool maximized = window && glfwGetWindowAttrib(window, GLFW_MAXIMIZED);
      if (captionButton(1, maximized ? icons::Icon::ChevronDown : icons::Icon::ChevronUp,
                        false, maximized ? "Restore" : "Maximize") &&
          window) {
        if (maximized) {
          glfwRestoreWindow(window);
        } else {
          glfwMaximizeWindow(window);
        }
      }
      if (captionButton(2, icons::Icon::Close, true, "Close") && window) {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
      }
    }
    ImGui::EndMenuBar();
  }

  if (std::optional<std::string> path = dialog.draw()) {
    dispatchDialogResult(st, dialogAction, *path);
    dialogAction = DialogAction::None;
  }

  if (ImGui::BeginPopupModal("About ChemCAD", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
    const float logo = ImGui::GetFontSize() * 2.4f;
    ImGui::Dummy(ImVec2(logo, logo));
    {
      const ImVec2 min = ImGui::GetItemRectMin();
      icons::draw(ImGui::GetWindowDrawList(), icons::Icon::Logo,
                  ImVec2(min.x + logo * 0.5f, min.y + logo * 0.5f), logo * 0.9f,
                  style::u32(style::col::Text), style::metrics().hairline * 2.0f);
    }
    const bool pushed = style::pushFont(style::fonts::semibold());
    ImGui::TextUnformatted("ChemCAD");
    style::popFont(pushed);
    ImGui::TextUnformatted("Chemical structure sketching and reaction planning");
    ImGui::Separator();
    ImGui::TextDisabled("C++20 / Dear ImGui / RDKit");
    ImGui::TextDisabled("Inter and JetBrains Mono fonts, SIL OFL");
    if (widgets::primaryButton("Close")) ImGui::CloseCurrentPopup();
    ImGui::EndPopup();
  }

  const ImGuiIO& io = ImGui::GetIO();
  if (!io.WantTextInput && io.KeyCtrl) {
    if (ImGui::IsKeyPressed(ImGuiKey_N, false)) {
      newDocument(st);
    } else if (ImGui::IsKeyPressed(ImGuiKey_O, false)) {
      openDialog(dialog, dialogAction, DialogAction::OpenProject, st);
    } else if (ImGui::IsKeyPressed(ImGuiKey_I, false)) {
      openDialog(dialog, dialogAction, DialogAction::ImportMol, st);
    } else if (ImGui::IsKeyPressed(ImGuiKey_S, false)) {
      if (io.KeyShift || st.projectPath.empty() || !st.saveProject) {
        openDialog(dialog, dialogAction, DialogAction::SaveProject, st);
      } else {
        st.saveProject(st.projectPath);
      }
    } else if (io.KeyShift && ImGui::IsKeyPressed(ImGuiKey_M, false) &&
               !st.doc.empty()) {
      openDialog(dialog, dialogAction, DialogAction::ExportMol, st);
    } else if (io.KeyShift && ImGui::IsKeyPressed(ImGuiKey_G, false) &&
               !st.doc.empty()) {
      openDialog(dialog, dialogAction, DialogAction::ExportSvg, st);
    } else if (io.KeyShift && ImGui::IsKeyPressed(ImGuiKey_P, false) &&
               !st.doc.empty()) {
      openDialog(dialog, dialogAction, DialogAction::ExportPng, st);
    } else if (io.KeyShift && ImGui::IsKeyPressed(ImGuiKey_A, false)) {
      st.props.autoName = !st.props.autoName;
    } else if (io.KeyShift && ImGui::IsKeyPressed(ImGuiKey_Delete, false) &&
               !st.doc.empty()) {
      clearStructure(st);
    } else if (ImGui::IsKeyPressed(ImGuiKey_L, false) && !st.doc.empty()) {
      cleanUp(st);
    } else if (ImGui::IsKeyPressed(ImGuiKey_Q, false)) {
      if (GLFWwindow* window = glfwGetCurrentContext()) {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
      }
    } else if (ImGui::IsKeyPressed(ImGuiKey_0, false)) {
      st.cam.zoom = 1.0f;
    } else if (ImGui::IsKeyPressed(ImGuiKey_F, false)) {
      st.cam.fit(st.doc, st.canvasSize);
    }
  } else if (!io.WantTextInput && ImGui::IsKeyPressed(ImGuiKey_F1, false)) {
    ImGui::OpenPopup("About ChemCAD");
  } else if (!io.WantTextInput && ImGui::IsKeyPressed(ImGuiKey_F2, false)) {
    // The View menu has advertised F2 since the profiler landed; nothing was
    // listening for it.
    st.showProfiler = !st.showProfiler;
  }
}

}  // namespace chemcad::ui

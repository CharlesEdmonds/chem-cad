#include <algorithm>
#include <filesystem>
#include <limits>
#include <string>
#include <unordered_map>

#include "imgui.h"

#include <GLFW/glfw3.h>

#include "chem/bridge.hpp"
#include "ui/file_dialog.hpp"
#include "ui/ui.hpp"

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

std::string projectFilename(const AppState& st) {
  if (st.projectPath.empty()) return {};
  std::error_code ec;
  const auto name = std::filesystem::path(st.projectPath).filename();
  return ec ? std::string{} : name.string();
}

void openDialog(FileDialog& dialog, DialogAction& action, DialogAction requested,
                const AppState& st) {
  action = requested;
  switch (requested) {
    case DialogAction::OpenProject:
      dialog.open("Open Project", FileDialogMode::Open, {".chemcad"});
      break;
    case DialogAction::SaveProject:
      dialog.open("Save Project", FileDialogMode::Save, {".chemcad"}, ".chemcad",
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
    if (ImGui::BeginMenu("File")) {
      if (ImGui::MenuItem("New", "Ctrl+N")) newDocument(st);
      if (ImGui::MenuItem("Open...", "Ctrl+O")) {
        openDialog(dialog, dialogAction, DialogAction::OpenProject, st);
      }
      if (ImGui::MenuItem("Save", "Ctrl+S")) {
        if (!st.projectPath.empty() && st.saveProject) {
          st.saveProject(st.projectPath);
        } else {
          openDialog(dialog, dialogAction, DialogAction::SaveProject, st);
        }
      }
      if (ImGui::MenuItem("Save As...", "Ctrl+Shift+S")) {
        openDialog(dialog, dialogAction, DialogAction::SaveProject, st);
      }
      ImGui::Separator();
      if (ImGui::MenuItem("Import MOL...")) {
        openDialog(dialog, dialogAction, DialogAction::ImportMol, st);
      }
      if (ImGui::MenuItem("Export MOL...", nullptr, false, !st.doc.empty())) {
        openDialog(dialog, dialogAction, DialogAction::ExportMol, st);
      }
      if (ImGui::MenuItem("Export SVG...", nullptr, false, !st.doc.empty())) {
        openDialog(dialog, dialogAction, DialogAction::ExportSvg, st);
      }
      if (ImGui::MenuItem("Export PNG...", nullptr, false, !st.doc.empty())) {
        openDialog(dialog, dialogAction, DialogAction::ExportPng, st);
      }
      ImGui::Separator();
      if (ImGui::MenuItem("Quit", "Ctrl+Q")) {
        if (GLFWwindow* window = glfwGetCurrentContext()) {
          glfwSetWindowShouldClose(window, GLFW_TRUE);
        }
      }
      ImGui::EndMenu();
    }

    if (ImGui::BeginMenu("Edit")) {
      if (ImGui::MenuItem("Undo", "Ctrl+Z", false, st.undo.canUndo())) undo(st);
      if (ImGui::MenuItem("Redo", "Ctrl+Shift+Z", false, st.undo.canRedo())) redo(st);
      ImGui::Separator();
      if (ImGui::MenuItem("Copy SMILES", "Ctrl+C", false, !st.doc.empty())) copySmiles(st);
      if (ImGui::MenuItem("Paste SMILES", "Ctrl+V")) pasteSmiles(st);
      ImGui::EndMenu();
    }

    if (ImGui::BeginMenu("Structure")) {
      if (ImGui::MenuItem("Clean Up Structure", nullptr, false, !st.doc.empty())) cleanUp(st);
      ImGui::MenuItem("Auto-name", nullptr, &st.props.autoName);
      ImGui::EndMenu();
    }

    if (ImGui::BeginMenu("View")) {
      if (ImGui::MenuItem("Fit to window", "F")) st.cam.fit(st.doc, st.canvasSize);
      if (ImGui::MenuItem("Reset zoom", "Ctrl+0")) st.cam.zoom = 1.0f;
      ImGui::EndMenu();
    }

    if (ImGui::BeginMenu("Help")) {
      if (ImGui::MenuItem("About ChemCAD")) ImGui::OpenPopup("About ChemCAD");
      ImGui::EndMenu();
    }
    ImGui::EndMenuBar();
  }

  if (std::optional<std::string> path = dialog.draw()) {
    dispatchDialogResult(st, dialogAction, *path);
    dialogAction = DialogAction::None;
  }

  if (ImGui::BeginPopupModal("About ChemCAD", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
    ImGui::TextUnformatted("ChemCAD");
    ImGui::TextUnformatted("Chemical structure sketching and reaction planning");
    ImGui::Separator();
    ImGui::TextDisabled("C++20 / Dear ImGui / RDKit");
    if (ImGui::Button("Close")) ImGui::CloseCurrentPopup();
    ImGui::EndPopup();
  }

  const ImGuiIO& io = ImGui::GetIO();
  if (!io.WantTextInput && io.KeyCtrl) {
    if (ImGui::IsKeyPressed(ImGuiKey_N, false)) {
      newDocument(st);
    } else if (ImGui::IsKeyPressed(ImGuiKey_O, false)) {
      openDialog(dialog, dialogAction, DialogAction::OpenProject, st);
    } else if (ImGui::IsKeyPressed(ImGuiKey_S, false)) {
      if (io.KeyShift || st.projectPath.empty() || !st.saveProject) {
        openDialog(dialog, dialogAction, DialogAction::SaveProject, st);
      } else {
        st.saveProject(st.projectPath);
      }
    } else if (ImGui::IsKeyPressed(ImGuiKey_Q, false)) {
      if (GLFWwindow* window = glfwGetCurrentContext()) {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
      }
    } else if (ImGui::IsKeyPressed(ImGuiKey_0, false)) {
      st.cam.zoom = 1.0f;
    }
  }
}

}  // namespace chemcad::ui

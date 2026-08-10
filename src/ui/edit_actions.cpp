#include "ui/edit_actions.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>

#include "imgui.h"

#include "chem/bridge.hpp"

namespace chemcad::ui {
namespace {

uint64_t refKey(int mol, core::AtomId id) {
  return (static_cast<uint64_t>(static_cast<uint32_t>(mol)) << 32) | id;
}

// The selection, or the whole document when nothing is selected. A bond in the
// selection carries both its atoms, so dragging a marquee over a ring and
// copying gives the ring rather than a bag of loose atoms.
core::Molecule moleculeForClipboard(const AppState& st) {
  core::Molecule result;
  std::unordered_map<uint64_t, core::AtomId> ids;
  const bool wholeDocument = st.sel.empty();

  auto selectedAtom = [&](int mol, core::AtomId id) {
    if (wholeDocument) return true;
    if (st.sel.contains(AtomRef{mol, id})) return true;
    for (const BondRef& ref : st.sel.bonds) {
      if (ref.mol != mol) continue;
      const core::Bond* bond = st.doc.molecules[static_cast<size_t>(mol)].bond(ref.id);
      if (bond && (bond->a == id || bond->b == id)) return true;
    }
    return false;
  };

  for (int mi = 0; mi < static_cast<int>(st.doc.molecules.size()); ++mi) {
    const core::Molecule& mol = st.doc.molecules[static_cast<size_t>(mi)];
    for (const core::Atom& atom : mol.atoms()) {
      if (!selectedAtom(mi, atom.id)) continue;
      core::Atom copy = atom;
      copy.id = core::kInvalidAtom;
      ids.emplace(refKey(mi, atom.id), result.addAtom(copy));
    }
  }
  for (int mi = 0; mi < static_cast<int>(st.doc.molecules.size()); ++mi) {
    const core::Molecule& mol = st.doc.molecules[static_cast<size_t>(mi)];
    for (const core::Bond& bond : mol.bonds()) {
      const auto a = ids.find(refKey(mi, bond.a));
      const auto b = ids.find(refKey(mi, bond.b));
      if (a == ids.end() || b == ids.end()) continue;
      const core::BondId newId = result.addBond(a->second, b->second, bond.order);
      if (core::Bond* copy = result.bond(newId)) copy->stereo = bond.stereo;
    }
  }
  return result;
}

// Lands the paste clear of what is already drawn instead of on top of it.
void offsetForPaste(const core::Document& document, core::Molecule& molecule) {
  float dx = 0.6f;
  float dy = -0.6f;
  if (!document.empty() && !molecule.empty()) {
    float existingMaxX = -std::numeric_limits<float>::infinity();
    float pastedMinX = std::numeric_limits<float>::infinity();
    for (const core::Molecule& mol : document.molecules) {
      for (const core::Atom& atom : mol.atoms()) existingMaxX = std::max(existingMaxX, atom.pos.x);
    }
    for (const core::Atom& atom : molecule.atoms()) {
      pastedMinX = std::min(pastedMinX, atom.pos.x);
    }
    if (std::isfinite(existingMaxX) && std::isfinite(pastedMinX)) {
      dx = existingMaxX - pastedMinX + 1.5f;
    }
  }
  for (core::Atom& atom : molecule.mutableAtoms()) {
    atom.pos.x += dx;
    atom.pos.y += dy;
  }
}

void applyHistory(AppState& st, bool redo) {
  if (!(redo ? st.undo.redo(st.doc) : st.undo.undo(st.doc))) return;
  // Both, not either: the restored document has different atom and bond ids, so
  // a surviving selection and a surviving hover are equally stale.
  clearTransientRefs(st);
  st.touch();
  st.statusMessage = redo ? "Redo" : "Undo";
}

}  // namespace

void clearTransientRefs(AppState& st) {
  st.sel.clear();
  st.hoverAtom = {};
  st.hoverBond = {};
}

void removeEmptyFragments(AppState& st) {
  std::erase_if(st.doc.molecules, [](const core::Molecule& mol) { return mol.empty(); });
}

void undoDocument(AppState& st) { applyHistory(st, false); }
void redoDocument(AppState& st) { applyHistory(st, true); }

void copySelectionAsSmiles(AppState& st) {
  try {
    const core::Molecule copy = moleculeForClipboard(st);
    if (copy.empty()) {
      st.statusMessage = "Nothing to copy";
      return;
    }
    st.clipboardSmiles = chem::toSmiles(copy);
    ImGui::SetClipboardText(st.clipboardSmiles.c_str());
    st.statusMessage = st.sel.empty() ? "Copied SMILES" : "Copied selection as SMILES";
  } catch (const std::exception& error) {
    st.statusMessage = std::string("Copy failed: ") + error.what();
  }
}

void pasteSmilesFromClipboard(AppState& st) {
  const char* system = ImGui::GetClipboardText();
  std::string smiles = system && *system ? system : st.clipboardSmiles;
  if (smiles.empty()) {
    st.statusMessage = "Clipboard does not contain SMILES";
    return;
  }
  try {
    core::Molecule pasted = chem::fromSmiles(smiles);
    offsetForPaste(st.doc, pasted);
    st.snapshot();
    st.doc.molecules.push_back(std::move(pasted));
    clearTransientRefs(st);
    st.touch();
    st.clipboardSmiles = std::move(smiles);
    st.statusMessage = "Pasted SMILES";
  } catch (const std::exception& error) {
    st.statusMessage = std::string("Paste failed: ") + error.what();
  }
}

void deleteSelection(AppState& st) {
  if (st.sel.empty()) return;
  st.snapshot();
  for (int molIndex = static_cast<int>(st.doc.molecules.size()) - 1; molIndex >= 0; --molIndex) {
    core::Molecule& mol = st.doc.molecules[static_cast<size_t>(molIndex)];
    for (const BondRef& ref : st.sel.bonds) {
      if (ref.mol == molIndex) mol.removeBond(ref.id);
    }
    for (const AtomRef& ref : st.sel.atoms) {
      if (ref.mol == molIndex) mol.removeAtom(ref.id);
    }
  }
  removeEmptyFragments(st);
  clearTransientRefs(st);
  st.touch();
}

}  // namespace chemcad::ui

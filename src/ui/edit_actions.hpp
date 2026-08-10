#pragma once
// Document-level edit actions: undo, redo, clipboard, delete.
//
// One implementation each, because there used to be two. The Edit menu and the
// sketch canvas had grown separate copies of all four, and the copies had
// drifted: the menu's undo cleared the selection but left `hoverAtom` and
// `hoverBond` pointing into a document that had just changed underneath them,
// while the canvas's undo cleared those and left the selection dangling. The
// menu's copy ignored the selection and always took the whole document; the
// canvas's copy took the selection and fell back to the whole document only
// when nothing was selected. Both were reached by the same advertised shortcut,
// so which one you got depended on where the focus happened to be.
//
// Anything that mutates the document from a menu item, a shortcut or a gesture
// belongs here, so that divergence cannot happen again.

#include "ui/app_state.hpp"

namespace chemcad::ui {

// Drops every reference into the document that a mutation can invalidate. Call
// after any edit that can remove an atom or a bond.
void clearTransientRefs(AppState&);
void removeEmptyFragments(AppState&);

void undoDocument(AppState&);
void redoDocument(AppState&);

// Copies the selection as SMILES, or the whole document when nothing is
// selected -- which is what makes one shortcut sensible on both an empty and a
// partial selection.
void copySelectionAsSmiles(AppState&);
void pasteSmilesFromClipboard(AppState&);
void deleteSelection(AppState&);

}  // namespace chemcad::ui

#pragma once
// Small read-only structure renderer used by the reaction planner cards.
// Shares no state with the interactive canvas.

#include "imgui.h"

#include "core/model.hpp"

namespace chemcad::ui {

// Draws `mol` scaled to fit the screen-space box [min,max], centred.
// A molecule with no atoms draws a dashed placeholder box.
void drawMoleculeThumb(ImDrawList* dl, const core::Molecule& mol, ImVec2 min, ImVec2 max);

// Consumes layout space at the cursor, draws a framed thumbnail and returns
// true when it was clicked.
bool moleculeThumbButton(const char* id, const core::Molecule& mol, ImVec2 size);

}  // namespace chemcad::ui

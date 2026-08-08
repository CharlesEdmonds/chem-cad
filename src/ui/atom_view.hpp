#pragma once
// Animated atomic-structure renderer used by the periodic-table hover card.

#include "imgui.h"

#include <string>

#include "ui/element_data.hpp"

namespace chemcad::ui {

// Draws the element's atom into the given rect: a rotating pseudo-3D point
// cloud of the occupied atomic orbitals (s/p/d/f) around a proton+neutron
// nucleus. bohrMode selects the backup Bohr shell model instead.
void drawAtomModel(const ElementData& element, ImVec2 min, ImVec2 max, bool bohrMode);

// Electron configuration in noble-gas-core shorthand, e.g. "[Ar] 3d10 4s2 4p4",
// with the known anomalous ground states (Cr, Cu, Pd, ...) applied.
std::string elementConfigString(const ElementData& element);

}  // namespace chemcad::ui

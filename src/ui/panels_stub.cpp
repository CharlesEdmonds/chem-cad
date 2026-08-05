// Wave-0 placeholders. Each function below is claimed by a Wave-1 slice; as a
// slice lands it deletes its stub from this file and adds its own translation
// unit (ui/*.cpp is globbed, so no CMake edit is needed).
#include "imgui.h"

#include "ui/mol_thumb.hpp"
#include "ui/ui.hpp"

namespace chemcad::ui {

void drawMenuBar(AppState&) {}
void drawToolPalette(AppState&) { ImGui::TextUnformatted("tools"); }
void drawCanvas(AppState&) { ImGui::TextUnformatted("canvas"); }
void drawPeriodicTable(AppState&) { ImGui::TextUnformatted("periodic table"); }
void drawPropertiesPanel(AppState&) { ImGui::TextUnformatted("properties"); }
void drawReactionPlanner(AppState&) { ImGui::TextUnformatted("reaction planner"); }
void drawStatusBar(AppState&) {}
void drawMoleculeThumb(ImDrawList*, const core::Molecule&, ImVec2, ImVec2) {}
bool moleculeThumbButton(const char*, const core::Molecule&, ImVec2) { return false; }

}  // namespace chemcad::ui

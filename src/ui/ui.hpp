#pragma once
// Panel entry points. One panel per translation unit so parallel work on the
// UI never collides in a shared file.

#include "ui/app_state.hpp"

namespace chemcad::ui {

void drawMenuBar(AppState&);        // ui/menu_bar.cpp
void drawToolPalette(AppState&);    // ui/tool_palette.cpp
void drawCanvas(AppState&);         // ui/canvas.cpp
void drawPeriodicTable(AppState&);  // ui/periodic_table.cpp
void drawPropertiesPanel(AppState&);// ui/properties_panel.cpp
void drawReactionPlanner(AppState&);// ui/reaction_planner.cpp
void drawToolbox(AppState&);        // ui/toolbox_view.cpp
void drawSolventSelector(AppState&);// ui/solvent_selector.cpp
void drawSolubilitySuite(AppState&);// ui/solubility_suite.cpp
void drawExtractionLab(AppState&);  // ui/extraction_view.cpp
// Starts the extraction physics build on a worker thread at startup so the
// workspace is ready the first time it is opened. ui/extraction_view.cpp
void warmExtractionPhysics(AppState&);
void drawViewer3D(AppState&);       // ui/viewer3d_view.cpp
void drawStatusBar(AppState&);      // ui/status_bar.cpp

// Applies the dark theme, sizing and fonts. Called once at startup.
void applyTheme(float uiScale);

}  // namespace chemcad::ui

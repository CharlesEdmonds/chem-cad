#pragma once
// UI state for the 3D molecule viewer (src/ui/viewer3d_view.cpp). Kept
// separate from AppState like the other panel states so the panel stays
// self-contained.

#include <cstdint>
#include <string>

#include "chem/embed3d.hpp"

namespace chemcad::ui {

struct Viewer3DState {
  float yawDeg = 35.0f;
  float pitchDeg = -18.0f;
  float zoom = 1.0f;      // multiplies the fit-to-view scale
  bool autoRotate = false;
  int style = 0;          // 0 ball-and-stick, 1 licorice, 2 space-filling

  // Cached model, re-embedded whenever the sketch changes.
  uint64_t sourceRevision = 0;
  bool hasModel = false;
  chem::Embedded3D model;
  std::string formula;    // display caption, computed alongside the embed
  double molWeight = 0.0;
  std::string errorMessage;
};

}  // namespace chemcad::ui

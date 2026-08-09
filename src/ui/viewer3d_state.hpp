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

  // Flat 2D depiction of the same molecule (sketch coordinates), used by the
  // Skeleton style: it is the sketch, and edge-on it collapses like a sheet.
  struct SketchAtom {
    uint8_t z = 6;
    float x = 0.0f, y = 0.0f;
  };
  struct SketchBond {
    int a = 0, b = 0, order = 1;
  };
  std::vector<SketchAtom> sketchAtoms;
  std::vector<SketchBond> sketchBonds;
  float sketchRadius = 1.0f;
  bool hasSketch = false;
};

}  // namespace chemcad::ui

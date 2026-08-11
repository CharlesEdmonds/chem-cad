#pragma once
// UI state for the molecule and atomic-orbital views. Kept separate from
// AppState so both cameras and their controls survive mode changes.

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "chem/embed3d.hpp"
#include "ui/charts3d.hpp"

namespace chemcad::ui {

struct Viewer3DState {
  int mode = 0;

  int orbitalN = 3;
  int orbitalL = 2;
  int orbitalM = 0;
  charts3d::OrbitalStyle orbitalStyle{
      style::col::DataBright, style::col::DataDim, true, true, false, 0.30f};
  charts3d::Orbit orbitalOrbit;

  float yawDeg = 35.0f;
  float pitchDeg = -18.0f;
  float zoom = 1.0f;      // multiplies the fit-to-view scale
  bool autoRotate = false;
  int style = 0;          // 0 ball-and-stick, 1 licorice, 2 space-filling, 3 skeleton

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

// Node census for a hydrogenic orbital, the pair the Orbitals overlay reports.
// The radial factor has n - l - 1 spherical nodes and the angular factor has l
// nodal surfaces, so the two always sum to n - 1 (Levine, *Quantum Chemistry*,
// 7th ed., section 6.6). Split out of the draw call so the arithmetic is
// checkable without scraping pixels.
int orbitalRadialNodes(int n, int l);
int orbitalAngularNodes(int l);

// Formats that census the way the overlay shows it, e.g.
// "1 radial node  |  2 angular nodes".
void formatOrbitalNodes(int n, int l, char* out, std::size_t size);

}  // namespace chemcad::ui

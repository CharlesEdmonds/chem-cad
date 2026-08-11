#pragma once
// Dimensional and multivariate charts: the pictures chemists actually reason
// with, drawn with ImDrawList like the rest of the charting idiom.
//
// These live apart from charts.hpp because they share machinery charts.hpp has
// no use for -- a common projection, depth sorting and an orbit controller --
// and because a 2D instrument tile should not have to pay for it.
//
// Everything here is painter's-algorithm software rendering into the current
// window's draw list. That is deliberate: these views sit inside docked ImGui
// panels, must obey their clip rectangles, and must keep working in the
// headless panel tests that have no GL context at all.

#include "imgui.h"

#include <cstdint>

#include "ui/theme.hpp"

namespace chemcad::ui::charts3d {

// Orbit state for any of the views below. Panels own one per view and persist
// it, so a rotation survives a rebuild of the data.
struct Orbit {
  float yawDeg = 38.0f;
  float pitchDeg = 26.0f;
  float zoom = 1.0f;
  bool spin = false;        // slow idle rotation when the user is not dragging
  float spinRateDegS = 9.0f;
};

// Applies drag-to-orbit and wheel-to-zoom for the item just submitted, and
// advances the idle spin. Call immediately after the view's InvisibleButton.
void handleOrbitInput(Orbit& orbit, bool hovered, bool active);

// ------------------------------------------------------------------ surfaces
// A scalar field z = f(u, v) on a regular grid, drawn as a shaded wireframe
// surface with height-mapped colour. `values` is row-major, `rows` x `columns`.
struct SurfaceStyle {
  const char* uLabel = nullptr;
  const char* vLabel = nullptr;
  const char* wLabel = nullptr;    // the height axis
  ImVec4 low = style::col::DataDim;
  ImVec4 high = style::col::Data;
  ImVec4 peak = style::col::Accent;
  bool logHeight = false;
  bool showGrid = true;            // wireframe over the shaded quads
  bool showFloor = true;           // projected contour map on the base plane
  bool markPeak = true;
};
// Returns the linear index of the hovered cell, else -1.
int surface(const char* id, const double* values, int columns, int rows, ImVec2 size,
            Orbit& orbit, const SurfaceStyle& style);

// Composition triangle for three components. `points` holds `count` samples of
// three barycentric weights each (they need not be normalised), `values` the
// scalar to colour by. Optionally marks one working point.
struct TernaryStyle {
  const char* aLabel = "A";
  const char* bLabel = "B";
  const char* cLabel = "C";
  ImVec4 low = style::col::DataDim;
  ImVec4 high = style::col::Data;
  bool isolines = true;
  bool hasMarker = false;
  double markerA = 0.0, markerB = 0.0, markerC = 0.0;
};
// Returns the sample index under the pointer, else -1. When the user clicks
// inside the triangle, `outA/outB/outC` receive that composition and the
// function returns -2, which is how a panel makes the chart an INPUT.
int ternary(const char* id, const double* points, const double* values, int count,
            ImVec2 size, const TernaryStyle& style, double* outA, double* outB,
            double* outC);

// ------------------------------------------------------------- point clouds
// Labelled points in a 3D property space, with an optional sphere around one of
// them -- the Hansen picture: solvents scattered in (dD, dP, dH) and the
// solute's interaction radius drawn as the sphere that should enclose the good
// ones.
struct CloudPoint {
  double x = 0.0, y = 0.0, z = 0.0;
  const char* label = nullptr;
  ImVec4 colour = style::col::Data;
  float radius = 1.0f;   // relative marker size
  bool highlighted = false;
};
struct CloudStyle {
  const char* xLabel = nullptr;
  const char* yLabel = nullptr;
  const char* zLabel = nullptr;
  bool showAxes = true;
  bool showLabels = true;   // dropped automatically when they would collide
  bool hasSphere = false;
  double sphereX = 0.0, sphereY = 0.0, sphereZ = 0.0, sphereRadius = 0.0;
  ImVec4 sphereColour = style::col::Accent;
};
// Returns the hovered point index, else -1.
int cloud(const char* id, const CloudPoint* points, int count, ImVec2 size, Orbit& orbit,
          const CloudStyle& style);

// ------------------------------------------------------- atomic structure
// Hydrogenic atomic orbitals. `n`, `l`, `m` are the usual quantum numbers; the
// view draws the angular probability lobes with their sign phase, and the
// radial nodes of the chosen shell, so 3d_z2 looks like 3d_z2 and not like a
// generic dumbbell. This is a real evaluation of the spherical harmonic, not a
// stylised blob.
struct OrbitalStyle {
  ImVec4 positive = style::col::Data;    // positive phase lobe
  ImVec4 negative = style::col::Accent;  // negative phase lobe
  bool showAxes = true;
  bool showNodes = true;      // nodal planes / radial nodes
  bool cutaway = false;       // slice the near half to expose radial structure
  // Fraction of the OUTERMOST radial lobe's peak amplitude. Anchoring on that
  // lobe rather than on max |psi| keeps every shell of a high-n orbital inside
  // the slider's range; see orbital() for why the global maximum cannot work.
  float isoLevel = 0.30f;
};
void orbital(const char* id, int n, int l, int m, ImVec2 size, Orbit& orbit,
             const OrbitalStyle& style);

// Human name for a set of quantum numbers, e.g. "3d z2". Returns a static
// string; invalid combinations yield "invalid".
const char* orbitalName(int n, int l, int m);

// ------------------------------------------------------------- multivariate
// One polyline per candidate across shared normalised axes. The right tool for
// ranking many solvents on many properties at once, where a table forces the
// reader to hold ten numbers in their head.
struct ParallelAxis {
  const char* label = nullptr;
  double min = 0.0;
  double max = 1.0;
  bool higherIsBetter = true;
  const char* unit = nullptr;
};
struct ParallelSeries {
  const char* label = nullptr;
  const double* values = nullptr;   // one per axis, in DATA units
  ImVec4 colour = style::col::DataDim;
  bool highlighted = false;
};
// Returns the hovered series index, else -1.
int parallelCoordinates(const char* id, const ParallelAxis* axes, int axisCount,
                        const ParallelSeries* series, int seriesCount, ImVec2 size);

// Stepwise contributions from a starting value to a final one: each step is a
// floating bar, gains and losses tinted apart, with running totals.
struct WaterfallStep {
  const char* label = nullptr;
  double delta = 0.0;
  bool total = false;   // draw as an absolute column from the baseline
};
void waterfall(const char* id, double start, const WaterfallStep* steps, int count,
               const char* unit, ImVec2 size);

}  // namespace chemcad::ui::charts3d

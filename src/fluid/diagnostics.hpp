#pragma once
// Turning particle state into the numbers chemistry needs.
//
// The chemistry layer must not read particles directly: it needs volumes,
// interfacial area and droplet size, and those are integrals over the state,
// not fields of it. Each quantity here has a definition that survives being
// quoted in a report:
//
//   bulk vs dispersed: same-phase connected components over the neighbour
//     graph (r <= H). The largest component of a phase that also holds at
//     least 5% of that phase's volume is its bulk; everything else is
//     dispersed. Hysteresis (4%/6%) stops the classification flickering.
//
//   Sauter mean diameter: each dispersed component c has an equivalent-sphere
//     diameter d_c = (6 V_c / pi)^(1/3), and
//         d32 = sum_c d_c^3 / sum_c d_c^2
//     which is the diameter of a monodisperse dispersion with the same
//     volume-to-surface ratio. Reported as 0 only when nothing is dispersed,
//     never as a stand-in for "infinitely fine".
//
//   interfacial area: with a phase indicator chi (0 for phase A, 1 for B), the
//     smoothed colour field c_i and its SPH gradient give
//         A = sum_i V_i |grad c_i|
//     which is the quadrature of integral |grad chi| dV, the area of the
//     transition surface. Validated against a sphere: within 10% of 4 pi R^2
//     for R >= 3 H.
//
//   layer heights: from the bulk volumes through the vessel's own revolution
//     integral, so a layer's drawn height and its volume can never disagree.
//     Output is exponentially filtered for display only; the filtered value is
//     never fed back into the physics.

#include <vector>

#include "fluid/grid.hpp"
#include "fluid/types.hpp"
#include "fluid/vessel_sdf.hpp"

namespace chemcad::fluid {

class DiagnosticsEngine {
 public:
  // Recomputes everything from the current particle state. `grid` must have
  // been built for these positions. Deterministic for a given state.
  Diagnostics compute(const Particles&, const NeighbourGrid&, const VesselBoundary&,
                      const std::vector<PhaseMaterial>& phases, double spacing);

  // Applies the display-only smoothing to layer heights and dispersed
  // fraction, using the previous result held internally. `dt` is the interval
  // since the last call.
  Diagnostics smooth(const Diagnostics& raw, double dt);

  void reset();

 private:
  Diagnostics previous_;
  bool havePrevious_ = false;
  std::vector<int> component_;      // per particle, -1 until assigned
  std::vector<double> componentVol_;
  std::vector<uint8_t> componentPhase_;
  std::vector<double> componentCentroidZ_;
};

}  // namespace chemcad::fluid

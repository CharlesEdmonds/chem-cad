#pragma once
// Shared value types for the 3D multiphase fluid solver.
//
// Frame convention for EVERY type in this module: vessel-local coordinates in
// metres, origin at the centre of the vessel's bottom, +z up. This differs
// deliberately from the 2D cross-section renderer (which uses +y up in a
// vessel-local plane); the 3D solver is the physics, the cross-section is one
// projection of it.
//
// Storage is structure-of-arrays with float state and double reductions: the
// solver's neighbour loops are memory-bound, while its convergence tests and
// diagnostics must be reproducible, so every accumulation over particles is
// performed in double in a fixed index order (see solver.hpp).

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace chemcad::fluid {

// One liquid charged into the vessel. Values are SI; the UI converts.
struct PhaseMaterial {
  std::string label;
  double restDensity = 998.0;         // kg/m^3
  double dynamicViscosity = 0.89e-3;  // Pa.s
  double volumeMl = 0.0;              // charged volume
  float colour[4] = {0.30f, 0.60f, 0.90f, 0.85f};
};

// Particle state, structure of arrays. Index i is one particle throughout.
struct Particles {
  std::vector<float> px, py, pz;  // position, m
  std::vector<float> vx, vy, vz;  // velocity, m/s
  std::vector<float> ax, ay, az;  // non-pressure acceleration, m/s^2
  std::vector<float> delta;       // number density, 1/m^3
  std::vector<float> pressure;    // Pa, >= 0
  std::vector<float> colour;      // smoothed phase indicator, 0..1
  std::vector<float> nx, ny, nz;  // interface normal (unnormalised gradient)
  std::vector<uint8_t> phase;     // index into the material table
  std::vector<uint32_t> id;       // immutable, for deterministic tie-breaking

  std::size_t size() const { return px.size(); }
  bool empty() const { return px.empty(); }

  // Every array is kept exactly parallel, so these three operations touch all
  // of them together; a mismatch would silently corrupt the neighbour search.
  void clear() { forEach([](auto& v) { v.clear(); }); }
  void reserve(std::size_t n) { forEach([n](auto& v) { v.reserve(n); }); }
  void resize(std::size_t n) { forEach([n](auto& v) { v.resize(n); }); }

  // Appends one particle at rest. Returns its index.
  std::size_t add(float x, float y, float z, uint8_t phaseIndex) {
    const std::size_t index = px.size();
    px.push_back(x);
    py.push_back(y);
    pz.push_back(z);
    vx.push_back(0.0f);
    vy.push_back(0.0f);
    vz.push_back(0.0f);
    ax.push_back(0.0f);
    ay.push_back(0.0f);
    az.push_back(0.0f);
    delta.push_back(0.0f);
    pressure.push_back(0.0f);
    colour.push_back(0.0f);
    nx.push_back(0.0f);
    ny.push_back(0.0f);
    nz.push_back(0.0f);
    phase.push_back(phaseIndex);
    id.push_back(static_cast<uint32_t>(index));
    return index;
  }

 private:
  template <typename Fn>
  void forEach(Fn&& fn) {
    fn(px); fn(py); fn(pz);
    fn(vx); fn(vy); fn(vz);
    fn(ax); fn(ay); fn(az);
    fn(delta); fn(pressure); fn(colour);
    fn(nx); fn(ny); fn(nz);
    fn(phase); fn(id);
  }
};

// Solver resolution and the constants derived from it. Everything else in the
// solver reads its lengths from here so a resolution change stays consistent.
struct Resolution {
  double spacing = 4.0e-3;  // dx, m: nominal particle lattice spacing
  // Kernel SUPPORT radius, m. H = 2*dx is the standard support-to-spacing
  // ratio (Akinci et al. 2013 use spacing = 2r with support = 4r), giving
  // roughly 30-45 interior neighbours in 3D.
  double support() const { return 2.0 * spacing; }
  double particleVolume() const { return spacing * spacing * spacing; }
};

// Diagnostics the chemistry layer consumes. Computed at a lower rate than the
// physics step (see diagnostics.hpp) because it needs connected components.
struct PhaseDiagnostics {
  double totalMl = 0.0;      // charged volume of this phase
  double bulkMl = 0.0;       // volume in the settled bulk component
  double dispersedMl = 0.0;  // volume in droplets
  double layerTopM = 0.0;    // top of this phase's bulk layer, vessel z
  bool bulkResolved = false; // false when no component qualifies as bulk
};

struct Diagnostics {
  std::vector<PhaseDiagnostics> phases;
  double interfacialAreaM2 = 0.0;   // total A-B interfacial area
  double sauterDiameterM = 0.0;     // d32 of the dispersed components, 0 when none
  int dispersedComponents = 0;      // droplet count actually resolved
  double dispersedFraction = 0.0;   // dispersed volume / total volume, 0..1
  double freeSurfaceM = 0.0;        // top of the liquid column, vessel z
  double kineticEnergyJ = 0.0;      // total, for settling detection
  bool valid = false;
};

}  // namespace chemcad::fluid

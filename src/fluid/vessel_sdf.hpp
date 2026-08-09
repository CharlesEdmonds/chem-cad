#pragma once
// Analytic signed-distance boundary for the vessel.
//
// The vessel is the surface of revolution of the analytic half-width profile
// that already defines the 2D cross-section (`sol::vesselWidthAt`), so the 3D
// fluid and the 2D schematic are guaranteed to be the same vessel. Boundary
// particles were rejected in favour of an SDF: an analytic distance has no
// sampling holes at the neck, costs one table lookup per particle, and stays
// fixed in vessel coordinates while the vessel shakes (see frame.hpp).
//
// Two things the SDF must supply beyond distance:
//   * the kernel-integrated number density of the glass side, delta_b(phi),
//     added to a particle's own number density so a wall-adjacent particle is
//     not treated as being at a free surface (which would give it a spurious
//     negative pressure and let the liquid stick to or leak through glass);
//   * the outward normal, for the contact projection and friction.

#include <vector>

#include "fluid/types.hpp"
#include "sol/funnel.hpp"

namespace chemcad::fluid {

struct SurfaceQuery {
  double distance = 0.0;  // signed, m: negative inside the liquid space
  double nx = 0.0;        // unit normal pointing from glass into the fluid
  double ny = 0.0;
  double nz = 1.0;
};

class VesselBoundary {
 public:
  // Builds the boundary for one vessel shape at one physical height, sampling
  // the analytic profile into a polyline. `support` is the kernel support
  // radius, needed for the boundary density table.
  void build(sol::Vessel vessel, double heightM, double support, double spacing);

  double heightM() const { return heightM_; }
  double maxRadiusM() const { return maxRadiusM_; }

  // Interior radius at a height, m. Clamped outside [0, height].
  double radiusAt(double z) const;

  // Nearest-surface distance and normal for a point in vessel coordinates.
  SurfaceQuery query(double x, double y, double z) const;

  // Kernel-integrated glass number density at signed distance `phi`, and its
  // derivative. Zero for phi <= -support (fully immersed in liquid space).
  double boundaryDensity(double phi) const;
  double boundaryDensitySlope(double phi) const;

  // Volume of the vessel below a height, m^3, from the same revolution
  // integral the capacity calibration uses.
  double volumeBelow(double z) const;
  // Inverse of volumeBelow: the height that holds `volume` m^3.
  double heightForVolume(double volume) const;

  // Fills `particles` with a lattice of `spacing` covering the volume needed
  // for the given phase materials, dense phase at the bottom. Returns the
  // number of particles created. Deterministic for a given input.
  std::size_t chargeLattice(const std::vector<PhaseMaterial>& phases, double spacing,
                            Particles& particles) const;

 private:
  sol::Vessel vessel_ = sol::Vessel::SeparatoryFunnel;
  double heightM_ = 0.19;
  double maxRadiusM_ = 0.045;
  double support_ = 8.0e-3;
  double spacing_ = 4.0e-3;
  std::vector<double> profileZ_;       // sampled polyline, ascending z
  std::vector<double> profileR_;       // interior radius at profileZ_
  std::vector<double> cumulativeVol_;  // revolution volume below profileZ_
  std::vector<double> boundaryDelta_;  // delta_b table over phi in [-H, +H]
  std::vector<double> boundarySlope_;  // d(delta_b)/d(phi)
};

}  // namespace chemcad::fluid

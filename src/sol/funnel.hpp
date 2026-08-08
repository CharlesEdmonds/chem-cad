#pragma once
// 2D cross-section liquid-liquid separation simulation.

#include <string>
#include <vector>

#include "core/model.hpp"

namespace chemcad::sol {

enum class Vessel { SeparatoryFunnel, DecantingFlask, GraduatedCylinder };

// One bulk liquid charged into the vessel.
struct Phase {
  std::string label;
  double volumeMl = 0.0;
  double density = 1.0;              // g/mL -- decides the stacking order
  double viscosity = 1.0;            // mPa.s
  double interfacialTension = 30.0;  // mN/m
  double emulsionStability = 0.3;    // 0 breaks instantly, 1 never settles
  float colour[4] = {0.30f, 0.60f, 0.90f, 0.70f};
};

// A dispersed droplet of one phase suspended in the column.
struct Droplet {
  core::Vec2 position{};  // vessel-local metres, origin at the bottom centre
  core::Vec2 velocity{};
  float radius = 0.0f;    // metres
  int phase = 0;          // index into Simulation::phases
};

struct Simulation {
  Vessel vessel = Vessel::SeparatoryFunnel;
  double vesselVolumeMl = 250.0;
  std::vector<Phase> phases;
  std::vector<Droplet> droplets;
  // Settled bulk volume still in the layer stack, mL, parallel to `phases`.
  std::vector<double> settledMl;
  double shakeEnergy = 0.0;  // decays each step; drives dispersion
  double elapsed = 0.0;      // seconds since the last reset
  unsigned seed = 1u;
};

// Charges the vessel: sorts phases dense-first, puts every mL in the settled
// stack and clears all droplets.
void reset(Simulation&);

// Injects mixing energy. `vigour` in [0, 1] -- converts settled volume into
// droplets whose size falls with vigour and with interfacial tension.
void shake(Simulation&, double vigour);

// Advances the simulation by dt seconds: Stokes rise/fall, drag, coalescence
// back into the settled stack, emulsion decay.
void step(Simulation&, double dt);

// Total charged volume, mL.
double totalVolumeMl(const Simulation&);

// Fraction of the charged volume currently dispersed as droplets, 0..1.
double emulsifiedFraction(const Simulation&);

// Closed polygon outline of the vessel cross-section, in vessel-local metres,
// counter-clockwise from the bottom.
std::vector<core::Vec2> vesselOutline(Vessel, double heightMetres);

// Interior half-width of the vessel at a height fraction (0 bottom, 1 top),
// as a fraction of the widest half-width.
double vesselWidthAt(Vessel, double heightFraction);

}  // namespace chemcad::sol

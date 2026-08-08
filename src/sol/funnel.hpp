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

// Physical description of one hand shake: how long, how fast, and how far the
// vessel is oscillated. These are the quantities a chemist actually controls;
// everything else (power input, droplet size, dispersion rate) is DERIVED from
// them -- there is deliberately no dimensionless "vigour" knob.
struct ShakeParams {
  double durationS = 5.0;     // seconds of shaking
  double frequencyHz = 3.0;   // oscillation frequency; 2-4 Hz is a firm hand shake
  double amplitudeM = 0.05;   // stroke half-amplitude of the vessel motion, m
};

// Derived shake state, filled by `shake()` and advanced by `step()`. All
// derived quantities are computed once at shake start so the UI can display
// the same numbers the physics uses.
struct ShakeState {
  bool active = false;
  double remainingS = 0.0;    // counts down inside step(); pauses with the sim
  double durationS = 0.0;
  double frequencyHz = 0.0;
  double amplitudeM = 0.0;
  double peakVelocity = 0.0;   // u = 2*pi*f*A, m/s -- peak slosh velocity
  double specificPower = 0.0;  // epsilon, W/kg -- energy input per unit liquid mass
  double sauterRadiusM = 0.0;  // Hinze mean droplet radius target, m
};

struct Simulation {
  Vessel vessel = Vessel::SeparatoryFunnel;
  double vesselVolumeMl = 250.0;
  std::vector<Phase> phases;
  std::vector<Droplet> droplets;
  // Settled bulk volume still in the layer stack, mL, parallel to `phases`.
  std::vector<double> settledMl;
  ShakeState shake;
  double shakeEnergy = 0.0;  // 1 while shaking, decays after; drives churn only
  double elapsed = 0.0;      // seconds since the last reset
  unsigned seed = 1u;
};

// Charges the vessel: sorts phases dense-first, puts every mL in the settled
// stack and clears all droplets.
void reset(Simulation&);

// Starts a shake: derives the slosh velocity u = 2*pi*f*A, the specific power
// input epsilon = u^2*f/2, and the Hinze Sauter mean droplet radius from the
// phase interfacial tensions, then disperses settled volume progressively
// over the shake duration inside `step`. Calling again re-arms the shake.
void shake(Simulation&, const ShakeParams&);

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

#pragma once
// Flory-Huggins / Hansen solubility model.

#include <array>
#include <string>
#include <vector>

#include "core/model.hpp"
#include "sol/solvent.hpp"

namespace chemcad::sol {

// Everything the model needs about the dissolved compound.
struct Solute {
  std::string name;
  double molarMass = 0.0;      // g/mol
  double molarVolume = 0.0;    // cm3/mol (McGowan characteristic volume)
  Hansen hansen;
  double meltingPoint = 25.0;      // deg C; <= the working temperature means liquid
  double entropyOfFusion = 56.5;   // J/(mol K); Walden's rule
  double logP = 0.0;
  double interactionRadius = 8.0;  // Hansen sphere radius R0, MPa^0.5
  bool meltingPointEstimated = false;  // true when Joback supplied it
};

// Group-contribution estimate of the solute from a sketched structure.
// Throws SolError when the structure is empty or cannot be interpreted.
Solute describeSolute(const core::Molecule&);

// One solvent component of a blend. volumeFraction need not be normalised.
struct Component {
  const Solvent* solvent = nullptr;
  double volumeFraction = 0.0;
};

// Volume-fraction weighted properties of a blend.
struct Mixture {
  Hansen hansen;
  double density = 0.0;      // g/mL
  double molarVolume = 0.0;  // cm3/mol
};

// Empty or all-zero input yields a zeroed Mixture rather than throwing.
Mixture blend(const std::vector<Component>&);

struct Prediction {
  double gramsPerMillilitre = 0.0;
  double molesPerLitre = 0.0;
  double moleFraction = 0.0;
  double idealMoleFraction = 0.0;
  double ra = 0.0;                        // Hansen distance, MPa^0.5
  double relativeEnergyDifference = 0.0;  // Ra / R0
  double activityCoefficient = 1.0;
  double chi = 0.0;                       // Flory-Huggins interaction parameter
  bool outsideSphere = false;             // RED > 1
  bool converged = true;
};

Prediction predict(const Solute&, const std::vector<Component>&, double temperatureC = 25.0);

struct SweepPoint {
  std::array<double, 3> fractions{};  // volume fractions, sums to 1
  Prediction prediction;
};

// Sweeps 1, 2 or 3 solvents over a simplex grid. `steps` is the number of
// subdivisions per axis and is clamped to [2, 64]. More than 3 solvents throws
// SolError.
std::vector<SweepPoint> sweep(const Solute&, const std::vector<const Solvent*>&, int steps,
                              double temperatureC = 25.0);

}  // namespace chemcad::sol

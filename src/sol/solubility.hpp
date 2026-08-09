#pragma once
// Flory-Huggins / Hansen solubility model.

#include <array>
#include <string>
#include <vector>

#include "core/model.hpp"
#include "sol/anchors.hpp"
#include "sol/ionization.hpp"
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
  std::string canonicalSmiles;         // identity key for literature anchors
  // Ionisable sites, counter-ions and salt form, from analyseIonization().
  Ionization ionization;
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
  bool anchored = false;    // a measured literature value contributed
  bool saltPath = false;    // Ksp equilibrium, not the organic FH model
  std::string anchorNote;   // e.g. "measured value (caffeine, CRC Handbook)"

  // ---- ionisation (see sol/ionization.hpp) ----------------------------
  bool ionicPath = false;         // an acid/base/salt correction was applied
  double pH = 7.0;                // pH the correction used
  bool pHSelfBuffered = true;     // true when the pH came from the solute itself
  double pKa = 0.0;               // dominant site, estimated unless stated
  double ionisedFraction = 0.0;   // ionised share of the dissolved solute, 0..1
  double bornPenaltyDecades = 0.0;  // decades the ion loses to a low dielectric
  bool ceilingLimited = false;    // saturated against the crystal-density bound
  std::string ionNote;            // human-readable summary of the ionic path
};

// Neutral organic solutes use the Flory-Huggins + extended-Hansen model,
// corrected toward any literature anchors for the involved pure solvents
// (log-linear in the blend fraction, so an anchored pure endpoint returns
// its measured value exactly). 1:1 salts in a water-containing blend use
// the Ksp equilibrium instead, honouring `background` for the common-ion
// effect and ionic strength.
// Sentinel for "no pH supplied": the solute's own saturated solution sets it
// (a free base runs basic, its hydrochloride runs acidic), which is what a
// bench chemist observes when dissolving the solid in unbuffered water.
constexpr double kAutoPH = -1.0;

Prediction predict(const Solute&, const std::vector<Component>&, double temperatureC = 25.0,
                   const Electrolyte* background = nullptr, double backgroundM = 0.0,
                   double pH = kAutoPH);

struct SweepPoint {
  std::array<double, 3> fractions{};  // volume fractions, sums to 1
  Prediction prediction;
};

// Sweeps 1, 2 or 3 solvents over a simplex grid. `steps` is the number of
// subdivisions per axis and is clamped to [2, 64]. More than 3 solvents throws
// SolError.
std::vector<SweepPoint> sweep(const Solute&, const std::vector<const Solvent*>&, int steps,
                              double temperatureC = 25.0,
                              const Electrolyte* background = nullptr, double backgroundM = 0.0,
                              double pH = kAutoPH);

// One row of the solvent screen: predicted solubility in a pure solvent.
struct ScreenRow {
  const Solvent* solvent = nullptr;
  Prediction prediction;
};

// Predicts solubility in every pure solvent in the database, sorted best
// (highest g/mL) first. This backs the solvent-selection table: tens of
// predictions, cheap enough to recompute whenever the solute or temperature
// changes. Throws SolError only when the database itself fails to load.
std::vector<ScreenRow> screen(const Solute&, double temperatureC = 25.0,
                              const Electrolyte* background = nullptr, double backgroundM = 0.0,
                              double pH = kAutoPH);

// Distribution of a neutral solute between water and a water-immiscible
// organic phase, using logP as the octanol/water partition proxy.
struct Partition {
  double mgAqueous = 0.0;
  double mgOrganic = 0.0;
  double fractionOrganic = 0.0;  // mgOrganic / total, 0..1
};

// Splits `massMg` between the two phases at equilibrium:
//   D = 10^logP,  organic share = D*Vorg / (D*Vorg + Vaq).
// The textbook neutral-species approximation -- no pH or ionisation
// correction, and logP stands in for the actual solvent pair's log D.
Partition partition(double massMg, double logP, double volumeAqueousMl,
                    double volumeOrganicMl);

}  // namespace chemcad::sol

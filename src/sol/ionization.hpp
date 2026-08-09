#pragma once
// Acid/base ionisation of the solute, and what it does to solubility.
//
// The Flory-Huggins/Hansen/Yalkowsky chain in solubility.cpp describes a
// NEUTRAL species. That is the wrong species for a large share of real
// chemistry: alkaloid hydrochlorides, amine salts, sodium carboxylates and
// every ionisable drug substance dissolve as ions, and their solubility
// differs from the free base or free acid by two to four orders of magnitude
// in water while collapsing in a low-dielectric organic. This module supplies
// the three pieces needed to model that honestly:
//
//   1. which site ionises and at what pKa (group-contribution estimate),
//   2. the Henderson-Hasselbalch enhancement at a given pH, including the
//      self-buffered pH a saturated solution of the species itself sets,
//   3. the Born electrostatic penalty for carrying an ion into a solvent of
//      lower dielectric constant.
//
// Everything here is an ESTIMATE and says so: `Ionization::pKaEstimated` is
// true unless a measured value was supplied, and the group table below is a
// textbook one, not a fitted regression.

#include <string>
#include <vector>

#include "core/model.hpp"

namespace chemcad::sol {

// How the solute behaves on the pH axis.
enum class IonClass {
  Neutral,     // no site ionises inside a usable pH window
  Base,        // protonates to a cation (amines, N-heterocycles)
  Acid,        // deprotonates to an anion (carboxylic acids, phenols)
  Zwitterion,  // both, and both strong: net neutral over a wide pH range
};

// Result of inspecting a drawn structure for ionisable sites and counter-ions.
struct Ionization {
  IonClass ionClass = IonClass::Neutral;
  // pKa of the site that dominates: for a base this is the pKa of the
  // conjugate acid BH+, for an acid the pKa of AH.
  double pKa = 0.0;
  bool pKaEstimated = true;
  std::string siteLabel;  // e.g. "tertiary amine", "carboxylic acid"

  // True when the sketch carried a counter-ion fragment or a net formal
  // charge, i.e. the drawn species is already a salt rather than a free
  // base/acid. A salt's saturated solution self-buffers to the ionised side,
  // which is what makes it soluble.
  bool saltForm = false;
  std::string counterIon;        // formula of the recognised counter-ion, if any
  double counterIonRadiusNm = 0.20;  // Born radius of the counter-ion
  int netCharge = 0;             // formal charge on the solute skeleton as drawn
  int fragmentCount = 1;         // connected components in the source structure
};

// Splits a structure into its connected components (a sketched salt, or a
// sketch holding several molecules, arrives as one Molecule with several).
// Components come back largest-first by heavy-atom count.
std::vector<core::Molecule> splitComponents(const core::Molecule&);

// Inspects a whole structure: the largest component is the solute skeleton,
// everything else is examined as a possible counter-ion.
Ionization analyseIonization(const core::Molecule& structure);

// Fraction-of-ionised-to-neutral ratio at a given pH:
//   base:  [BH+]/[B]  = 10^(pKa - pH)
//   acid:  [A-]/[AH]  = 10^(pH - pKa)
// Neutral species and zwitterions return 0 (no enhancement). The result is
// capped at 1e12 so a wildly out-of-range pH cannot overflow the caller.
double ionisedRatio(const Ionization&, double pH);

// pH a saturated solution of this species sets on its own, with no buffer:
//   free acid:      pH = 0.5 (pKa - log C)
//   free base:      pH = 14 - 0.5 ((14 - pKa) - log C)
//   salt of a base: pH = 0.5 (pKa - log C)        (hydrolysis of BH+)
//   salt of an acid pH = 7 + 0.5 (pKa + log C)    (hydrolysis of A-)
// `molarity` is the saturated concentration; pass a positive estimate.
// Returns 7.0 for a neutral species (pure water).
double selfBufferedPH(const Ionization&, double molarity);

// Born electrostatic penalty, in log10 solubility units, for moving a 1:1
// ion pair from water into a solvent of dielectric constant `dielectric`:
//
//   dG_Born(r, eps) = -N_A z^2 e^2 / (8 pi eps_0 r) * (1 - 1/eps)
//   ddG = sum_ions [ N_A e^2 / (8 pi eps_0 r_i) ] * (1/eps_s - 1/eps_water)
//
// returned as ddG / (2.303 R T) >= 0, i.e. the number of decades of
// solubility the ionised branch LOSES. Radii are Latimer-corrected
// (r_eff = r_ion + 0.08 nm) because the bare-ion Born radius overestimates
// transfer free energies; even so, Born ignores ion pairing, so this is an
// upper bound on the penalty -- documented, not hidden.
double bornPenaltyDecades(double dielectric, double soluteRadiusNm,
                          double counterIonRadiusNm, double temperatureK,
                          double waterDielectric = 78.4);

// Born radius of the solute skeleton from its McGowan characteristic volume
// (a sphere of the same molecular volume), in nanometres.
double bornRadiusFromMolarVolumeNm(double molarVolumeCm3PerMol);

}  // namespace chemcad::sol

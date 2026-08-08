#pragma once
// Literature anchors and aqueous ionic chemistry for the solubility model.
//
// The Flory-Huggins/Hansen model estimates solubility from structure alone;
// where a measured value exists for an exact solute+solvent pair, the
// prediction should defer to it. Anchors pin pure-solvent endpoints to
// literature data and the model interpolates blends between them
// (log-linear in the blend fraction). Salts bypass the organic model
// entirely: their aqueous solubility is a saturation-product (Ksp)
// equilibrium, with a Davies activity correction and the common-ion effect
// when a background electrolyte shares an ion.

#include <string>
#include <string_view>
#include <vector>

#include "sol/solvent.hpp"

namespace chemcad::sol {

// One measured solute+solvent pair. `soluteSmiles` is canonicalised at load,
// so lookups match however the structure was drawn.
struct SolubilityAnchor {
  std::string soluteSmiles;  // canonical
  std::string solventId;
  double temperatureC = 25.0;
  double gramsPerMillilitre = 0.0;
  double soluteMeltingPointC = 25.0;  // literature Tm of the measured solute
  std::string note;
};

// A 1:1 ionic compound: its measured pure-water solubility (the anchor) and
// its saturation product at 25 C (used for common-ion/temperature RATIOS,
// so the measured endpoint stays exact).
struct Salt {
  std::string soluteSmiles;  // canonical
  std::string name;
  std::string cation;
  std::string anion;
  double molarMass = 0.0;
  double solubilityGPerMl25 = 0.0;  // measured, in pure water at 25 C
  double ksp25 = 0.0;
  double dissolutionEnthalpyKj = 0.0;  // for van't Hoff temperature scaling
  std::string note;
};

// A background electrolyte the user can add to the aqueous phase (for the
// common-ion effect and ionic strength).
struct Electrolyte {
  std::string id;
  std::string name;
  std::string cation;
  std::string anion;
  double molarMass = 0.0;
};

// All three tables load from data/solubility_anchors.json and follow the
// solvents() contract: sticky, sorted, and a persistent load failure throws
// SolError on every call.
const std::vector<SolubilityAnchor>& anchors();
const std::vector<Salt>& salts();
const std::vector<Electrolyte>& electrolytes();

// Exact-pair lookup; temperature must be within `toleranceC` of the anchor's
// measurement temperature. nullptr when no measured value exists.
const SolubilityAnchor* findAnchor(std::string_view canonicalSmiles, std::string_view solventId,
                                   double temperatureC, double toleranceC = 5.0);
const Salt* findSalt(std::string_view canonicalSmiles);
const Electrolyte* findElectrolyte(std::string_view id);

// Davies-equation mean ionic activity coefficient for a 1:1 electrolyte at
// the given ionic strength (clamped to the equation's valid range I <= 0.5).
double daviesGamma(double ionicStrength);

// Saturated molar solubility (mol per litre of water) of a 1:1 salt with an
// optional background electrolyte at concentration backgroundM.
//
//   gamma(I)^2 * (s + c_common) * s = Ksp(T),   I = I_background + s
//
// Ksp(T) follows the van't Hoff equation from dissolutionEnthalpyKj; gamma
// is the Davies activity coefficient. A background sharing an ion depresses
// solubility (common-ion effect); an inert one raises it slightly through
// the activity coefficient (salting-in).
double saltSolubilityMolar(const Salt&, const Electrolyte* background, double backgroundM,
                           double temperatureC);

}  // namespace chemcad::sol

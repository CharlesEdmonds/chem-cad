#pragma once
// Internal declaration of the shared extended-Hansen/Martin chi implemented
// in solubility.cpp. Used by the Kirkwood-Buff module; not public API.

#include <vector>

#include "sol/solubility.hpp"
#include "sol/solvent.hpp"

namespace chemcad::sol {

// Extended-Hansen/Martin interaction parameter for a solute against a blend
// described by its mixed Hansen parameters, at temperatureK. Identical to the
// regression driving the FH prediction (single definition, exported here so
// the Kirkwood-Buff module can report the same number).
double extendedHansenChi(const Solute& solute, const Hansen& mixtureHansen,
                         double temperatureK);

// The Flory-Huggins/extended-Hansen solve WITHOUT the literature-anchor
// correction. predict() layers anchors and the salt path on top of this; the
// Kirkwood-Buff module needs the raw model so its activity coefficient is the
// calibrated one rather than a second, divergent implementation.
Prediction floryHugginsPrediction(const Solute&, const std::vector<Component>&,
                                  double temperatureC);

}  // namespace chemcad::sol

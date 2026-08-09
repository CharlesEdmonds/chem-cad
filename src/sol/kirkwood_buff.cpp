// Kirkwood-Buff solubility module. See kirkwood_buff.hpp for the theory and
// the exact identities used here.

#include "sol/kirkwood_buff.hpp"

#include <algorithm>
#include <cmath>

namespace chemcad::sol {
namespace {

constexpr double kR = 8.314;  // J/(mol K)

// Ideal (Hildebrand/Yalkowsky) ln x_2 from the fusion term: Walden's-rule
// entropy of fusion over melting-point depression. A solute liquid at the
// working temperature is ideally miscible (ln x = 0). Same expression the FH
// path uses for its ideal ceiling, restated here for the readout.
double lnIdealSolubility(const Solute& solute, double temperatureC) {
  if (solute.meltingPoint <= temperatureC) return 0.0;
  const double meltingK = solute.meltingPoint + 273.15;
  const double temperatureK = temperatureC + 273.15;
  return -(solute.entropyOfFusion / kR) * (meltingK - temperatureK) / temperatureK;
}

}  // namespace

KBResult kirkwoodBuff(const Solute& solute, const std::vector<Component>& components,
                      double temperatureC) {
  KBResult result;
  const Mixture mixture = blend(components);

  // Volume-weighted mixture compressibility over the components that carry
  // data; a blend where fewer than half the volume has kappa_T is flagged.
  double kappaWeighted = 0.0;
  double fractionWithData = 0.0;
  double totalFraction = 0.0;
  for (const Component& component : components) {
    if (!component.solvent || component.volumeFraction <= 0.0) continue;
    totalFraction += component.volumeFraction;
    if (component.solvent->kappaT > 0.0) {
      kappaWeighted += component.solvent->kappaT * component.volumeFraction;
      fractionWithData += component.volumeFraction;
      if (result.kappaSource.empty()) result.kappaSource = component.solvent->kappaTSource;
    }
  }
  result.kappaKnown = totalFraction > 0.0 && fractionWithData / totalFraction > 0.5;
  result.kappaT = fractionWithData > 0.0 ? kappaWeighted / fractionWithData : 0.0;

  const double temperatureK = temperatureC + 273.15;
  const double v1 = std::max(mixture.molarVolume, 1e-6);      // cm3/mol
  const double v2 = std::max(solute.molarVolume, 1e-6);       // cm3/mol
  const double kappaSI = std::max(result.kappaT, 0.0) * 1e-9;  // GPa^-1 -> Pa^-1

  // Solvent-solvent KBI: the exact fluctuation-theory identity for a pure
  // liquid, G_11 = R T kappa_T - V_1 (m3/mol -> cm3/mol). For water this
  // returns about -17 cm3/mol, the accepted literature value. Without
  // kappa_T the incompressible limit (-V_1) is used and flagged.
  result.g11 = (kR * temperatureK * kappaSI - v1 * 1e-6) * 1e6;

  // The activity coefficient comes from the calibrated Flory-Huggins /
  // extended-Hansen solve rather than a second, divergent implementation:
  // KB theory is exact, so its usefulness here is the DECOMPOSITION it gives
  // (G_11 from measured compressibility, G_12 from the activity model), not
  // a rival correlation. Reimplementing the activity from KBIs alone would
  // require solute-solvent radial distribution functions that no desktop
  // input provides.
  const Prediction fh = floryHugginsPrediction(solute, components, temperatureC);
  result.chi = fh.chi;
  result.lnGammaInf = std::log(std::max(fh.activityCoefficient, 1e-12));

  // Ben-Naim inversion: ln gamma = (G_11 - G_12)/V_1  ->  G_12. A G_12 above
  // G_11 means the solute is preferentially solvated (favourable mixing).
  result.g12 = result.g11 - v1 * result.lnGammaInf;

  // Yalkowsky ideal solubility and the resulting saturation.
  result.lnIdeal = lnIdealSolubility(solute, temperatureC);
  result.moleFraction = fh.moleFraction;
  result.gPerMl = fh.gramsPerMillilitre;
  result.valid = std::isfinite(result.gPerMl) && result.gPerMl >= 0.0;
  return result;
}

KBResult kirkwoodBuff(const Solute& solute, const Solvent& solvent, double temperatureC) {
  const std::vector<Component> components{Component{&solvent, 1.0}};
  return kirkwoodBuff(solute, components, temperatureC);
}

}  // namespace chemcad::sol

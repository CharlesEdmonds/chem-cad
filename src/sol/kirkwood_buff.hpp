#pragma once
// Kirkwood-Buff theory of solubility.
//
// Kirkwood-Buff (KB) theory relates the thermodynamics of a solution to the
// Kirkwood-Buff integrals (KBIs) -- the volume integrals of the pair
// correlation functions between species:
//
//   G_ij = integral of (g_ij(r) - 1) * 4 pi r^2 dr   [volume per molecule]
//
// For a solid solute (2) dissolving in a solvent blend (1), Ben-Naim's
// inversion of KB theory gives the infinite-dilution activity coefficient:
//
//   ln gamma_2^inf = rho_1 * (G_11 - G_12)      rho_1 = 1 / V_1
//
// The KBIs in this module come from two measurable/modelled routes:
//
//   * Solvent-solvent: the exact fluctuation-theory identity for a pure
//     liquid, G_11 = R T kappa_T - V_1, from the isothermal compressibility
//     kappa_T (NIST Web Thermo Tables where available).
//   * Solute-solvent: the Flory-Huggins / extended-Hansen activity model the
//     suite already uses, rearranged into KB form. FH gives
//       ln gamma_2^inf = chi + ln(V_2/V_1) + (1 - V_2/V_1),
//     so G_12 = G_11 - V_1 * ln gamma_2^inf.
//
// Solubility then follows the ideal-solubility equation corrected by the KB
// activity coefficient:
//
//   ln x_2 = -DeltaS_fus/R * (T_m/T - 1)  -  ln gamma_2^inf
//
// (fusion term via Walden's rule entropy of fusion, shared with the FH path).
//
// References:
//   Ben-Naim, A. "Molecular Theory of Solutions", Oxford 2006, ch. 4-5.
//   Kirkwood, J.G.; Buff, F.P. J. Chem. Phys. 1951, 19, 774.
//   Matteoli, E.; Lepori, L. J. Chem. Phys. 1984, 80, 2856 (compressibility
//     route to G_11).

#include "sol/chi.hpp"
#include "sol/solubility.hpp"
#include "sol/solvent.hpp"

namespace chemcad::sol {

struct KBResult {
  bool valid = false;          // false when inputs are physically degenerate
  bool kappaKnown = false;     // false when the solvent carries no kappa_T
  std::string kappaSource;     // where kappa_T came from ("NIST WTT" / ...)
  double kappaT = 0.0;         // GPa^-1, as used
  double g11 = 0.0;            // solvent-solvent KBI, cm3/mol
  double g12 = 0.0;            // solute-solvent KBI, cm3/mol
  double chi = 0.0;            // extended-Hansen chi behind G_12
  double lnGammaInf = 0.0;     // KB infinite-dilution activity coefficient
  double lnIdeal = 0.0;        // fusion-term ideal solubility, ln x_2^ideal
  double moleFraction = 0.0;   // KB solubility, x_2
  double gPerMl = 0.0;         // same as a mass concentration
};

// KB solubility of `solute` in one pure solvent at `temperatureC`.
KBResult kirkwoodBuff(const Solute&, const Solvent&, double temperatureC = 25.0);

// KB solubility in a blend: G_11 uses the volume-weighted mixture kappa_T
// and molar volume; chi comes from the same blend mixture as the FH path.
KBResult kirkwoodBuff(const Solute&, const std::vector<Component>&,
                      double temperatureC = 25.0);

}  // namespace chemcad::sol

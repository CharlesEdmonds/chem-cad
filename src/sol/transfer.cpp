#include "sol/transfer.hpp"

#include <algorithm>
#include <cmath>
#include <string>

namespace chemcad::sol {
namespace {

constexpr double kMlToM3 = 1.0e-6;
constexpr double kWilkeChangCm2PerS = 7.4e-8;
constexpr double kCm2ToM2 = 1.0e-4;
constexpr double kPaSToCentipoise = 1.0e3;
constexpr double kPi = 3.14159265358979323846;

bool positiveFinite(double value) {
  return std::isfinite(value) && value > 0.0;
}

double nonnegativeFinite(double value) {
  return std::isfinite(value) && value > 0.0 ? value : 0.0;
}

double organicEquilibriumFraction(double distribution, double aqueousVolumeM3,
                                  double organicVolumeM3) {
  if (!positiveFinite(distribution) || !positiveFinite(organicVolumeM3)) return 0.0;
  if (!positiveFinite(aqueousVolumeM3)) return 1.0;

  const double weightedOrganic = distribution * organicVolumeM3;
  if (std::isinf(weightedOrganic)) return 1.0;
  return weightedOrganic / (aqueousVolumeM3 + weightedOrganic);
}

double approachFromMass(double organicMass, double equilibriumOrganicMass,
                        double totalMass) {
  if (!(totalMass > 0.0)) return 1.0;
  if (organicMass <= equilibriumOrganicMass) {
    if (!(equilibriumOrganicMass > 0.0)) return 1.0;
    return std::clamp(organicMass / equilibriumOrganicMass, 0.0, 1.0);
  }

  const double organicSideSpan = totalMass - equilibriumOrganicMass;
  if (!(organicSideSpan > 0.0)) return 1.0;
  return std::clamp((totalMass - organicMass) / organicSideSpan, 0.0, 1.0);
}

// Hinze's C_H = 0.725 is a maximum/high-percentile breakup scale, not a
// Sauter mean (Hinze, AIChE J. 1 (1955) 289-295, eq. 21). If a measured
// Rosin-Rammler number-distribution shape q is available, its exact moments
// give d32/d95 below. Rosin and Rammler, J. Inst. Fuel 7 (1933) 29-36.
double hinzeToSauter(double hinzeDiameter, double rosinRammlerShape) {
  if (!positiveFinite(hinzeDiameter)) return 0.0;
  if (!positiveFinite(rosinRammlerShape)) return hinzeDiameter;

  const double q = rosinRammlerShape;
  const double percentileScale = std::pow(-std::log(0.05), 1.0 / q);
  const double momentRatio = std::tgamma(1.0 + 3.0 / q) /
                             std::tgamma(1.0 + 2.0 / q);
  const double d32 = hinzeDiameter * momentRatio / percentileScale;
  return positiveFinite(d32) ? d32 : hinzeDiameter;
}

}  // namespace

double distributionRatio(const Solute& solute, double pH, double temperatureC) {
  // No transfer enthalpy is available, so applying an invented van't Hoff
  // correction would be less physical than leaving logP temperature-neutral.
  static_cast<void>(temperatureC);
  const double neutral = std::pow(10.0, std::clamp(solute.logP, -300.0, 300.0));

  // kAutoPH requires a concentration to solve the self-buffered pH. This API
  // has none, so the only reproducible interpretation is the neutral logP
  // partition used by partition(); callers wanting ionisation must supply pH.
  if (pH == kAutoPH || !std::isfinite(pH)) return neutral;

  // Henderson-Hasselbalch: bases are retained in water below pKa and acids
  // above pKa. ionisedRatio implements 10^(pKa-pH) and 10^(pH-pKa),
  // respectively, and caps the ratio before it can overflow.
  const double ratio = ionisedRatio(solute.ionization, pH);
  return neutral / (1.0 + ratio);
}

double diffusivityWilkeChang(double solventMolarMassGPerMol,
                             double solventViscosityPaS,
                             double soluteMolarVolumeCm3, double temperatureK,
                             double associationFactor) {
  if (!positiveFinite(solventMolarMassGPerMol) ||
      !positiveFinite(solventViscosityPaS) ||
      !positiveFinite(soluteMolarVolumeCm3) || !positiveFinite(temperatureK) ||
      !positiveFinite(associationFactor)) {
    return 0.0;
  }

  // Wilke and Chang, AIChE J. 1 (1955) 264-270, eq. 12:
  // D[cm2/s] = 7.4e-8 sqrt(phi M) T / (mu[cP] V_b[cm3/mol]^0.6).
  // Multiplying by 1e-4 converts cm2/s to m2/s. The caller-visible contract
  // labels McGowan volume as a proxy because it is not the fitted V_b.
  const double viscosityCp = solventViscosityPaS * kPaSToCentipoise;
  const double diffusivityCm2PerS =
      kWilkeChangCm2PerS * std::sqrt(associationFactor * solventMolarMassGPerMol) *
      temperatureK / (viscosityCp * std::pow(soluteMolarVolumeCm3, 0.6));
  return diffusivityCm2PerS * kCm2ToM2;
}

double filmCoefficientHigbie(double diffusivityM2PerS, double contactTimeS) {
  if (!positiveFinite(diffusivityM2PerS) || !positiveFinite(contactTimeS)) return 0.0;

  // Higbie penetration theory: k = 2 sqrt(D/(pi t_c)); Higbie,
  // Trans. AIChE 31 (1935) 365-389. t_c is a surface-renewal estimate, not a
  // measurement of a molecular boundary-layer thickness.
  return 2.0 * std::sqrt(diffusivityM2PerS / (kPi * contactTimeS));
}

double filmCoefficientDroplet(double diffusivityM2PerS, double dropletDiameterM,
                              double slipVelocityMPerS,
                              double kinematicViscosityM2PerS) {
  if (!positiveFinite(diffusivityM2PerS) || !positiveFinite(dropletDiameterM) ||
      !positiveFinite(kinematicViscosityM2PerS)) {
    return 0.0;
  }

  // Ranz-Marshall external-film estimate for a translating sphere:
  // Sh = 2 + 0.60 Re^1/2 Sc^1/3, k = Sh D/d, with
  // Re = |u_rel|d/nu and Sc = nu/D. Ranz and Marshall,
  // Chem. Eng. Prog. 48 (1952) 173-180. A deformable liquid droplet is not a
  // rigid sphere, so advanceTransfer labels this explicitly as an estimate.
  const double reynolds = std::abs(slipVelocityMPerS) * dropletDiameterM /
                          kinematicViscosityM2PerS;
  const double schmidt = kinematicViscosityM2PerS / diffusivityM2PerS;
  const double sherwood = 2.0 + 0.60 * std::sqrt(reynolds) * std::cbrt(schmidt);
  return sherwood * diffusivityM2PerS / dropletDiameterM;
}

double specificArea(double dispersedFraction, double sauterDiameterM) {
  if (!positiveFinite(sauterDiameterM) || !positiveFinite(dispersedFraction)) return 0.0;
  return 6.0 * std::clamp(dispersedFraction, 0.0, 1.0) / sauterDiameterM;
}

TransferState advanceTransfer(const TransferInput& input,
                              const TransferState& previous, double dt) {
  TransferState state;
  const double totalMass = nonnegativeFinite(input.soluteMassMg);
  const double aqueousVolumeM3 = nonnegativeFinite(input.aqueousVolumeMl) * kMlToM3;
  const double organicVolumeM3 = nonnegativeFinite(input.organicVolumeMl) * kMlToM3;
  const double distribution = nonnegativeFinite(input.distributionRatio);

  const double equilibriumFraction =
      organicEquilibriumFraction(distribution, aqueousVolumeM3, organicVolumeM3);
  state.equilibriumOrganicMassMg = totalMass * equilibriumFraction;

  const double previousAqueous = nonnegativeFinite(previous.aqueousMassMg);
  const double previousOrganic = nonnegativeFinite(previous.organicMassMg);
  const double previousTotal = previousAqueous + previousOrganic;
  double organicMass = 0.0;
  if (previousTotal > 0.0) {
    organicMass = totalMass * previousOrganic / previousTotal;
  }
  organicMass = std::clamp(organicMass, 0.0, totalMass);

  double diameter = nonnegativeFinite(input.sauterDiameterM);
  bool diameterFromHinze = false;
  bool diameterFromArea = false;
  if (!(diameter > 0.0) && positiveFinite(input.hinzeDiameterM)) {
    diameter = hinzeToSauter(input.hinzeDiameterM,
                             input.hinzeRosinRammlerShape);
    diameterFromHinze = true;
  }

  double area = nonnegativeFinite(input.interfacialAreaM2);
  state.areaMeasured = area > 0.0;
  const double totalLiquidVolumeM3 = aqueousVolumeM3 + organicVolumeM3;
  const double dispersedFraction =
      std::clamp(nonnegativeFinite(input.dispersedFraction), 0.0, 1.0);

  if (!(area > 0.0) && diameter > 0.0 && dispersedFraction > 0.0 &&
      totalLiquidVolumeM3 > 0.0) {
    area = specificArea(dispersedFraction, diameter) * totalLiquidVolumeM3;
  }

  // A measured area can still arrive before connected-component d32 is valid.
  // Recovering 6 V_d/A is then an equivalent Sauter diameter, not a new area
  // correlation; it only supplies the film-correlation length scale.
  if (area > 0.0 && !(diameter > 0.0)) {
    double dispersedVolumeM3 = dispersedFraction * totalLiquidVolumeM3;
    if (!(dispersedVolumeM3 > 0.0)) dispersedVolumeM3 = organicVolumeM3;
    if (dispersedVolumeM3 > 0.0) {
      diameter = 6.0 * dispersedVolumeM3 / area;
      diameterFromArea = true;
    }
  }

  const double temperatureK = input.temperatureC + 273.15;
  const double aqueousDiffusivity = diffusivityWilkeChang(
      input.aqueousMolarMassGPerMol, input.aqueousViscosityPaS,
      input.soluteMolarVolumeCm3, temperatureK,
      input.aqueousAssociationFactor);
  const double organicDiffusivity = diffusivityWilkeChang(
      input.organicMolarMassGPerMol, input.organicViscosityPaS,
      input.soluteMolarVolumeCm3, temperatureK,
      input.organicAssociationFactor);
  const double kinematicViscosity =
      positiveFinite(input.aqueousDensityKgPerM3)
          ? input.aqueousViscosityPaS / input.aqueousDensityKgPerM3
          : 0.0;
  const double aqueousFilm = filmCoefficientDroplet(
      aqueousDiffusivity, diameter, input.slipVelocityMPerS,
      kinematicViscosity);
  const double organicFilm = filmCoefficientHigbie(
      organicDiffusivity, input.contactTimeS);

  // AQUEOUS-BASIS two-film convention:
  //   J_aq->org = K_aq (c_aq - c_org/D)
  //   1/K_aq = 1/k_aq + 1/(D k_org).
  // This basis is dimensionally consistent with D = c_org*/c_aq* and gives
  // zero flux exactly at c_org = D c_aq (Lewis and Whitman, Ind. Eng. Chem.
  // 16 (1924) 1215-1220). The algebraic form avoids reciprocal overflow.
  double coefficientOverDistribution = 0.0;
  if (distribution > 0.0 && aqueousFilm > 0.0 && organicFilm > 0.0) {
    const double scaledOrganicFilm = distribution * organicFilm;
    if (std::isinf(scaledOrganicFilm)) {
      state.overallCoefficientMPerS = aqueousFilm;
    } else {
      const double denominator = scaledOrganicFilm + aqueousFilm;
      state.overallCoefficientMPerS =
          aqueousFilm * scaledOrganicFilm / denominator;
      coefficientOverDistribution =
          aqueousFilm * organicFilm / denominator;
    }
  }

  if (aqueousVolumeM3 > 0.0 && organicVolumeM3 > 0.0 && area > 0.0 &&
      state.overallCoefficientMPerS > 0.0) {
    // dm_org/dt = K A (m_aq/V_aq - m_org/(D V_org)). With conserved
    // total mass this is a scalar linear ODE. Its exact solution is
    // m_o(t+dt) = m_o,eq + (m_o(t)-m_o,eq) exp(-lambda dt), where
    // lambda = K A (1/V_aq + 1/(D V_org)). No Euler stability limit exists.
    const double lambda =
        area * (state.overallCoefficientMPerS / aqueousVolumeM3 +
                coefficientOverDistribution / organicVolumeM3);
    if (positiveFinite(lambda)) {
      state.timeConstantS = 1.0 / lambda;
      if (totalMass > 0.0 && dt > 0.0 && std::isfinite(dt)) {
        const double fractionAdvanced = -std::expm1(-lambda * dt);
        organicMass += (state.equilibriumOrganicMassMg - organicMass) *
                       fractionAdvanced;
      }
    }
  }

  state.organicMassMg = std::clamp(organicMass, 0.0, totalMass);
  state.aqueousMassMg = totalMass - state.organicMassMg;
  state.approachToEquilibrium = approachFromMass(
      state.organicMassMg, state.equilibriumOrganicMassMg, totalMass);

  const std::string basisNote =
      " Aqueous-basis two-film estimate: J=K_aq(c_aq-c_org/D), "
      "1/K_aq=1/k_aq+1/(D*k_org). The aqueous film uses the "
      "Ranz-Marshall droplet Sherwood estimate; the organic film uses the "
      "Higbie penetration estimate. Wilke-Chang uses McGowan volume as a proxy.";
  if (state.areaMeasured) {
    state.note = "Interfacial area measured by the fluid simulation.";
    if (diameterFromArea) {
      state.note += " Equivalent d32 was inferred from 6*V_d/A for film transport only.";
    }
    state.note += basisNote;
  } else if (area > 0.0) {
    if (diameterFromHinze && positiveFinite(input.hinzeRosinRammlerShape)) {
      state.note =
          "Interfacial area estimated, not measured: a=6*phi/d32, with d32 "
          "converted from the Hinze d95 scale using a calibrated "
          "Rosin-Rammler number-distribution shape.";
    } else if (diameterFromHinze) {
      state.note =
          "Interfacial area estimated, not measured: a=6*phi/d_Hinze. "
          "Without a measured distribution width, d_Hinze is retained as a "
          "conservative area lower bound and is not called a measured d32.";
    } else {
      state.note =
          "Interfacial area estimated, not measured: a=6*phi/d32 using the "
          "simulation-supplied Sauter diameter.";
    }
    state.note += basisNote;
  } else {
    state.note =
        "No interfacial area was measured and no positive dispersion "
        "correlation was available; transfer is disabled.";
  }

  return state;
}

}  // namespace chemcad::sol

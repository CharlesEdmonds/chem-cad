#pragma once
// Rate-based liquid-liquid solute transfer.

#include <string>

#include "sol/solubility.hpp"

namespace chemcad::sol {

// Distribution ratio D = c_organic / c_aqueous for a possibly ionisable
// solute: log D = log P - log10(1 + 10^(pKa-pH)) for a base, and
// log D = log P - log10(1 + 10^(pH-pKa)) for an acid. Neutral logP is the
// solvent-pair proxy because no solvent-pair-specific distribution data exist.
double distributionRatio(const Solute&, double pH, double temperatureC);

struct TransferInput {
  double soluteMassMg = 0.0;
  double aqueousVolumeMl = 0.0;
  double organicVolumeMl = 0.0;
  double distributionRatio = 0.0;     // D = c_organic / c_aqueous
  double interfacialAreaM2 = 0.0;     // measured by the fluid simulation
  double sauterDiameterM = 0.0;       // measured d32; 0 when nothing is dispersed
  double dispersedFraction = 0.0;     // measured dispersed volume / total volume, 0..1
  double aqueousViscosityPaS = 0.89e-3;
  double organicViscosityPaS = 0.41e-3;
  double soluteMolarVolumeCm3 = 0.0;  // McGowan volume, used as a documented WC proxy
  double temperatureC = 25.0;

  // Wilke-Chang and the droplet Sherwood correlation need properties that
  // viscosity and volume alone cannot determine. The defaults are water and
  // dichloromethane; callers selecting other phases should replace them.
  double aqueousMolarMassGPerMol = 18.01528;
  double organicMolarMassGPerMol = 84.93;
  double aqueousAssociationFactor = 2.6;
  double organicAssociationFactor = 1.0;
  double aqueousDensityKgPerM3 = 997.0;
  double slipVelocityMPerS = 0.0;
  double contactTimeS = 1.0;  // surface-renewal estimate for the Higbie film

  // Hinze's 0.725 correlation is a limiting (approximately d95/dmax) diameter,
  // not d32. A calibrated Rosin-Rammler number-distribution shape q permits
  // an explicit moment conversion. With q == 0, the limiting diameter itself
  // is used only as a conservative lower-bound estimate of area.
  double hinzeDiameterM = 0.0;
  double hinzeRosinRammlerShape = 0.0;
};

struct TransferState {
  double aqueousMassMg = 0.0;
  double organicMassMg = 0.0;
  double equilibriumOrganicMassMg = 0.0;
  double approachToEquilibrium = 0.0;    // 0..1
  double overallCoefficientMPerS = 0.0;  // aqueous-basis K_L
  double timeConstantS = 0.0;            // 0 when there is no active transfer path
  bool areaMeasured = false;  // false whenever a correlation supplied the area
  std::string note;           // distinguishes every measurement and estimate
};

// Advances the split by dt seconds. The two-compartment ODE is integrated in
// closed form, so the function is deterministic and remains stable for dt
// much larger than the mass-transfer time constant.
TransferState advanceTransfer(const TransferInput&, const TransferState& previous, double dt);

// Wilke-Chang diffusivity in m2/s. Its empirical solute-volume input is the
// normal-boiling molar volume; the app currently supplies McGowan volume as an
// explicitly labelled proxy.
double diffusivityWilkeChang(double solventMolarMassGPerMol, double solventViscosityPaS,
                             double soluteMolarVolumeCm3, double temperatureK,
                             double associationFactor = 1.0);

// Film coefficients in m/s.
double filmCoefficientHigbie(double diffusivityM2PerS, double contactTimeS);
double filmCoefficientDroplet(double diffusivityM2PerS, double dropletDiameterM,
                              double slipVelocityMPerS, double kinematicViscosityM2PerS);

// Specific interfacial area, 1/m, for spherical droplets: a = 6 phi / d32.
double specificArea(double dispersedFraction, double sauterDiameterM);

}  // namespace chemcad::sol

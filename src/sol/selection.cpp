#include "sol/selection.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace chemcad::sol {
namespace {

constexpr double kTiny = 1e-12;
constexpr double kReferenceSolubility = 0.01;  // g/mL; 10 mg/mL is the UI's practical-solubility reference.
constexpr double kWorkableRecovery = 0.05;
constexpr double kUsefulChromatographySpacing = 0.10;

bool isSeparationOperation(OperationKind kind) {
  return kind != OperationKind::ReactionMedium;
}

double separationGate(double selectivity) {
  // A ratio is naturally judged by orders of magnitude, not by linear distance.
  // This logistic in log10(ratio) is 0.13 at the failure boundary (1x), 0.98
  // at 10x, and approaches a 0.02 floor for strongly reversed separations.
  if (std::isinf(selectivity) && selectivity > 0.0) return 1.0;
  const double logRatio =
      std::log10(std::max(std::isfinite(selectivity) ? selectivity : 0.0, kTiny));
  return std::clamp(0.02 + 0.98 / (1.0 + std::exp(-6.0 * (logRatio - 0.35))),
                    0.0, 1.0);
}

int separationTier(const SolventCandidate& candidate) {
  if (!(candidate.selectivity >= 1.0)) return 0;
  return candidate.recoveryFraction >= kWorkableRecovery ? 2 : 1;
}

struct Roles {
  std::vector<const SpeciesRole*> targets;
  std::vector<const SpeciesRole*> contaminants;
};

struct CommonScores {
  double greenness = 0.5;
  double practicality = 0.5;
  std::string greenDetail;
  std::string practicalDetail;
};

struct EligibleSolvent {
  const Solvent* solvent = nullptr;
  std::vector<std::string> warnings;
};

double unit(double value) {
  if (!std::isfinite(value)) return 0.0;
  return std::clamp(value, 0.0, 1.0);
}

double positiveWeight(double value) {
  return std::isfinite(value) && value > 0.0 ? value : 0.0;
}

std::string number(double value, int decimals = 2) {
  char buffer[64];
  std::snprintf(buffer, sizeof(buffer), "%.*f", decimals,
                std::isfinite(value) ? value : 0.0);
  return buffer;
}

std::string speciesName(const SpeciesRole& role) {
  if (!role.label.empty()) return role.label;
  if (!role.solute.name.empty()) return role.solute.name;
  return "unnamed species";
}

void addWarning(std::vector<std::string>& warnings, std::string warning) {
  if (warning.empty()) return;
  if (std::find(warnings.begin(), warnings.end(), warning) == warnings.end()) {
    warnings.push_back(std::move(warning));
  }
}

Roles splitRoles(const OperationSpec& spec) {
  Roles roles;
  for (const SpeciesRole& species : spec.species) {
    if (species.keep) {
      roles.targets.push_back(&species);
    } else {
      roles.contaminants.push_back(&species);
    }
  }
  return roles;
}

int classRank(const std::string& value) {
  if (value == "recommended") return 0;
  if (value == "problematic") return 1;
  if (value == "hazardous") return 2;
  if (value == "highly hazardous") return 3;
  return -1;
}

bool isChlorinated(const Solvent& solvent) {
  // A literal chlorine atom in the recorded SMILES is narrower and more honest
  // than rejecting the whole "halogenated" family, which also contains non-chlorinated solvents.
  return solvent.smiles.find("Cl") != std::string::npos;
}

bool isProtic(const Solvent& solvent) {
  // Family is used because the solvent table is curated at the family level;
  // water, alcohols and carboxylic acids necessarily provide an O-H proton.
  return solvent.family == "water" || solvent.family == "alcohol" || solvent.family == "acid";
}

bool passesHardConstraints(const Solvent& solvent, const OperationSpec& spec,
                           std::vector<std::string>& warnings) {
  if (spec.requireWaterImmiscible && solvent.waterMiscible) return false;
  if (spec.requireWaterMiscible && !solvent.waterMiscible) return false;
  if (spec.minBoilingPointC > 0.0 && solvent.boilingPoint < spec.minBoilingPointC) return false;
  if (spec.maxBoilingPointC > 0.0 && solvent.boilingPoint > spec.maxBoilingPointC) return false;
  if (spec.avoidPeroxideFormers && solvent.peroxideFormer) return false;
  if (spec.avoidChlorinated && isChlorinated(solvent)) return false;
  if (spec.avoidAromatics && solvent.family == "aromatic") return false;
  if (spec.excludeUnrated && solvent.chem21Class.empty()) return false;

  const int allowed = classRank(spec.worstAcceptableClass);
  const int actual = classRank(solvent.chem21Class);
  if (allowed >= 0 && actual > allowed) return false;

  if (solvent.chem21Class.empty()) {
    addWarning(warnings, "CHEM21-unrated solvent; greenness is held at neutral rather than favoured.");
  } else if (allowed >= 0 && actual == allowed) {
    addWarning(warnings, "CHEM21 class is exactly the selected acceptance limit (" +
                             solvent.chem21Class + ").");
  }
  if (spec.minBoilingPointC > 0.0 &&
      solvent.boilingPoint - spec.minBoilingPointC <= 10.0) {
    addWarning(warnings, "Boiling point is within 10 C of the requested lower limit.");
  }
  if (spec.maxBoilingPointC > 0.0 &&
      spec.maxBoilingPointC - solvent.boilingPoint <= 10.0) {
    addWarning(warnings, "Boiling point is within 10 C of the requested upper limit.");
  }
  if (solvent.meltingPointC <= spec.temperatureC &&
      spec.temperatureC - solvent.meltingPointC <= 5.0) {
    addWarning(warnings, "Melting point is within 5 C of the working temperature.");
  }
  if (solvent.peroxideFormer) {
    addWarning(warnings, "Peroxide former: storage age and inhibitor status require review.");
  }
  if (!solvent.hazardNote.empty()) addWarning(warnings, solvent.hazardNote + ".");
  return true;
}

std::vector<EligibleSolvent> eligibleSolvents(const std::vector<Solvent>& database,
                                               const OperationSpec& spec) {
  std::vector<EligibleSolvent> eligible;
  if (spec.requireWaterImmiscible && spec.requireWaterMiscible) return eligible;
  eligible.reserve(database.size());
  for (const Solvent& solvent : database) {
    std::vector<std::string> warnings;
    if (passesHardConstraints(solvent, spec, warnings)) {
      eligible.push_back(EligibleSolvent{&solvent, std::move(warnings)});
    }
  }
  return eligible;
}

Prediction safePredict(const SpeciesRole& role, const std::vector<Component>& components,
                       double temperatureC, double pH,
                       std::vector<std::string>* warnings = nullptr) {
  try {
    return predict(role.solute, components, temperatureC, nullptr, 0.0, pH);
  } catch (const SolError& error) {
    // solvents() is loaded before this helper is entered, so a later SolError is
    // a missing optional model/anchor datum. Keep ranking usable and expose the loss.
    if (warnings) {
      addWarning(*warnings, "Could not model " + speciesName(role) + ": " + error.what() +
                                "; treated as zero solubility.");
    }
    return {};
  }
}

struct Evidence {
  bool allAnchored = true;
  std::vector<std::string> warnings;

  void observe(const SpeciesRole& role, const Prediction& prediction) {
    allAnchored = allAnchored && prediction.anchored;
    if (prediction.ionicPath && !prediction.ionNote.empty()) {
      addWarning(warnings, speciesName(role) + ": " + prediction.ionNote);
    }
  }
};

double predictionValue(const Prediction& prediction) {
  return std::isfinite(prediction.gramsPerMillilitre)
             ? std::max(0.0, prediction.gramsPerMillilitre)
             : 0.0;
}

double roleWeight(const SpeciesRole& role) {
  return positiveWeight(role.weight);
}

double weightedMean(const std::vector<const SpeciesRole*>& roles,
                    const std::vector<double>& values) {
  double sum = 0.0;
  double weights = 0.0;
  for (size_t i = 0; i < roles.size() && i < values.size(); ++i) {
    const double weight = roleWeight(*roles[i]);
    sum += weight * values[i];
    weights += weight;
  }
  return weights > 0.0 ? sum / weights : 0.0;
}

double weightedGeometricMean(const std::vector<const SpeciesRole*>& roles,
                             const std::vector<double>& values) {
  double logSum = 0.0;
  double weights = 0.0;
  for (size_t i = 0; i < roles.size() && i < values.size(); ++i) {
    const double weight = roleWeight(*roles[i]);
    logSum += weight * std::log(std::max(values[i], kTiny));
    weights += weight;
  }
  return weights > 0.0 ? std::exp(logSum / weights) : 0.0;
}

CommonScores commonScores(const Solvent& solvent, double workingTemperatureC,
                          std::vector<std::string>& warnings) {
  CommonScores scores;

  if (solvent.chem21Class.empty()) {
    scores.greenness = 0.5;
    scores.greenDetail = "CHEM21 unrated; assigned the neutral score 0.50.";
    addWarning(warnings, "CHEM21-unrated solvent; greenness is held at neutral rather than favoured.");
  } else {
    // CHEM21 reports Safety/Health/Environment on 1..10 scales with larger values
    // worse (Prat et al., Green Chem. 2016, 18, 288, DOI 10.1039/C5GC01008J).
    // Reversing the affine 1..10 scale gives 1 at score 1 and 0 at score 10.
    const double meanShe = (solvent.chem21Safety + solvent.chem21Health +
                            solvent.chem21Environment) /
                           3.0;
    const double sheScore = unit(1.0 - (meanShe - 1.0) / 9.0);
    const int rank = classRank(solvent.chem21Class);
    // The four ordered CHEM21 classes are placed uniformly on [1, 0]; averaging
    // this categorical result with the three measured guide scores avoids hiding either.
    const double classScore = rank >= 0 ? 1.0 - static_cast<double>(rank) / 3.0 : 0.5;
    scores.greenness = unit(0.75 * sheScore + 0.25 * classScore);
    scores.greenDetail = "CHEM21 S/H/E " + std::to_string(solvent.chem21Safety) + "/" +
                         std::to_string(solvent.chem21Health) + "/" +
                         std::to_string(solvent.chem21Environment) + ", class " +
                         solvent.chem21Class + "; combined score " +
                         number(scores.greenness) + ".";
  }

  // Practicality is an explicit average of five observable handling terms.
  // The boiling term is a triangular bench heuristic: 80 C is convenient for
  // rotary evaporation, and the score falls linearly to zero 120 C away.
  const double boiling = unit(1.0 - std::abs(solvent.boilingPoint - 80.0) / 120.0);
  // A liquid margin of 20 C is treated as fully robust; linear interpolation
  // down to the melting point makes near-freezing operation visibly marginal.
  const double liquid = unit((workingTemperatureC - solvent.meltingPointC) / 20.0);
  // Cost tiers are ordinal 1..4, so the unique affine reversal maps cheap to 1
  // and expensive to 0 without pretending the tiers are dollar ratios.
  const double cost = unit(1.0 - (static_cast<double>(solvent.costTier) - 1.0) / 3.0);
  const double peroxide = solvent.peroxideFormer ? 0.0 : 1.0;
  const double hazard = solvent.hazardNote.empty() ? 1.0 : 0.5;
  scores.practicality = unit((boiling + liquid + cost + peroxide + hazard) / 5.0);
  scores.practicalDetail = "bp " + number(solvent.boilingPoint, 1) + " C, mp " +
                           number(solvent.meltingPointC, 1) + " C, cost tier " +
                           std::to_string(solvent.costTier) +
                           (solvent.peroxideFormer ? ", peroxide former" : ", no peroxide flag") +
                           "; score " + number(scores.practicality) + ".";
  return scores;
}

CommonScores pairCommonScores(const Solvent& primary, const Solvent& partner,
                              double workingTemperatureC,
                              std::vector<std::string>& warnings) {
  std::vector<std::string> primaryWarnings;
  std::vector<std::string> partnerWarnings;
  const CommonScores a = commonScores(primary, workingTemperatureC, primaryWarnings);
  const CommonScores b = commonScores(partner, workingTemperatureC, partnerWarnings);
  for (const std::string& warning : primaryWarnings) addWarning(warnings, warning);
  for (const std::string& warning : partnerWarnings) {
    addWarning(warnings, "Anti-solvent: " + warning);
  }
  CommonScores result;
  // Equal volume-independent averaging is deliberate: both liquids must be
  // stocked and handled even when the anti-solvent fraction is modest.
  result.greenness = 0.5 * (a.greenness + b.greenness);
  result.practicality = 0.5 * (a.practicality + b.practicality);
  result.greenDetail = primary.name + " / " + partner.name + " mean: " +
                       number(result.greenness) + ".";
  result.practicalDetail = "Both-liquid handling mean: " + number(result.practicality) +
                           "; primary " + a.practicalDetail + " partner " + b.practicalDetail;
  return result;
}

void finishCandidate(SolventCandidate& candidate, const OperationSpec& spec,
                     double selectivityScore, std::string selectivityDetail,
                     double recoveryScore, std::string recoveryDetail,
                     const CommonScores& common, bool applySeparationGate) {
  const double selectivityWeight = positiveWeight(spec.weightSelectivity);
  const double recoveryWeight = positiveWeight(spec.weightRecovery);
  const double greenWeight = positiveWeight(spec.weightGreenness);
  const double practicalWeight = positiveWeight(spec.weightPracticality);
  if (applySeparationGate && isSeparationOperation(spec.kind)) {
    selectivityDetail += " Multiplicative separation gate " +
                         number(separationGate(candidate.selectivity), 3) +
                         " from the selectivity ratio.";
  }

  candidate.criteria = {
      Criterion{"Selectivity", unit(selectivityScore), selectivityWeight,
                std::move(selectivityDetail)},
      Criterion{"Recovery", unit(recoveryScore), recoveryWeight, std::move(recoveryDetail)},
      Criterion{"Greenness", unit(common.greenness), greenWeight, common.greenDetail},
      Criterion{"Practicality", unit(common.practicality), practicalWeight,
                common.practicalDetail},
  };

  // The weighted mean remains the user's preference score. Separation operations
  // then multiply it by a log-ratio gate, so safety or convenience cannot erase
  // a chemically reversed separation.
  const double totalWeight = selectivityWeight + recoveryWeight + greenWeight + practicalWeight;
  if (totalWeight <= 0.0) {
    candidate.score = 0.0;
    return;
  }
  const double weightedScore =
      unit((selectivityWeight * unit(selectivityScore) +
            recoveryWeight * unit(recoveryScore) +
            greenWeight * unit(common.greenness) +
            practicalWeight * unit(common.practicality)) /
           totalWeight);
  candidate.score =
      applySeparationGate && isSeparationOperation(spec.kind)
          ? unit(weightedScore * separationGate(candidate.selectivity))
          : weightedScore;
}

void appendEvidence(SolventCandidate& candidate, const Evidence& evidence) {
  for (const std::string& warning : evidence.warnings) addWarning(candidate.warnings, warning);
  candidate.estimated = !evidence.allAnchored;
  if (candidate.estimated) {
    addWarning(candidate.warnings,
               "Estimated: at least one solubility in this score is model-derived, not measured.");
  }
}

std::vector<SolventCandidate> rankExtraction(const OperationSpec& spec, const Roles& roles,
                                             const std::vector<EligibleSolvent>& eligible) {
  std::vector<SolventCandidate> result;
  const Solvent* water = findSolvent("water");
  if (!water) return result;

  std::vector<Prediction> targetWater;
  std::vector<Prediction> contaminantWater;
  targetWater.reserve(roles.targets.size());
  contaminantWater.reserve(roles.contaminants.size());
  for (const SpeciesRole* role : roles.targets) {
    targetWater.push_back(safePredict(*role, {{water, 1.0}}, spec.temperatureC, spec.pH));
  }
  for (const SpeciesRole* role : roles.contaminants) {
    contaminantWater.push_back(safePredict(*role, {{water, 1.0}}, spec.temperatureC, spec.pH));
  }

  for (const EligibleSolvent& entry : eligible) {
    const Solvent& solvent = *entry.solvent;
    // Liquid-liquid extraction requires a distinct phase, irrespective of whether
    // the user also enabled the general water-immiscibility filter.
    if (&solvent == water || solvent.waterMiscible || miscibleWith(solvent, *water)) continue;

    SolventCandidate candidate;
    candidate.solvent = &solvent;
    candidate.warnings = entry.warnings;
    Evidence evidence;
    std::vector<double> targetD;
    std::vector<double> targetRecovery;
    std::vector<double> targetOrganic;
    targetD.reserve(roles.targets.size());
    targetRecovery.reserve(roles.targets.size());
    targetOrganic.reserve(roles.targets.size());

    const double aqueousVolume = std::max(0.0, spec.aqueousVolumeMl);
    const double organicVolume = std::max(0.0, spec.organicVolumeMl);
    for (size_t i = 0; i < roles.targets.size(); ++i) {
      const SpeciesRole& role = *roles.targets[i];
      const Prediction organic =
          safePredict(role, {{&solvent, 1.0}}, spec.temperatureC, spec.pH, &evidence.warnings);
      evidence.observe(role, targetWater[i]);
      evidence.observe(role, organic);
      const double waterS = predictionValue(targetWater[i]);
      const double organicS = predictionValue(organic);
      // At dilute equilibrium, equal chemical potential makes the concentration
      // ratio D = C_org/C_aq; the model's saturation concentrations supply that ratio.
      const double distribution = organicS / std::max(waterS, kTiny);
      // One-contact extraction mass balance: f = D*Vorg/(D*Vorg + Vaq)
      // (standard liquid-liquid extraction derivation from C_org=D*C_aq).
      const double fraction = distribution * organicVolume /
                              std::max(distribution * organicVolume + aqueousVolume, kTiny);
      targetD.push_back(distribution);
      targetRecovery.push_back(unit(fraction));
      targetOrganic.push_back(organicS);
    }

    double worstContaminantD = 0.0;
    double worstContaminantS = 0.0;
    std::string worstName = "none specified";
    double maxContaminantWeight = 0.0;
    for (const SpeciesRole* role : roles.contaminants) {
      maxContaminantWeight = std::max(maxContaminantWeight, roleWeight(*role));
    }
    for (size_t i = 0; i < roles.contaminants.size(); ++i) {
      const SpeciesRole& role = *roles.contaminants[i];
      const Prediction organic =
          safePredict(role, {{&solvent, 1.0}}, spec.temperatureC, spec.pH, &evidence.warnings);
      evidence.observe(role, contaminantWater[i]);
      evidence.observe(role, organic);
      const double waterS = predictionValue(contaminantWater[i]);
      const double organicS = predictionValue(organic);
      const double distribution = organicS / std::max(waterS, kTiny);
      // Relative importance scales the worst-case contaminant continuously while
      // leaving the highest-weight contaminant at its physical D value.
      const double weightedD = distribution * roleWeight(role) /
                               std::max(maxContaminantWeight, kTiny);
      if (weightedD > worstContaminantD) {
        worstContaminantD = weightedD;
        worstContaminantS = organicS;
        worstName = speciesName(role);
      }
    }

    const double representativeD = weightedGeometricMean(roles.targets, targetD);
    const double recovery = weightedMean(roles.targets, targetRecovery);
    const double ratio = roles.contaminants.empty()
                             ? representativeD
                             : representativeD / std::max(worstContaminantD, kTiny);
    // r/(1+r) is the odds-to-probability transform: equal target/contaminant
    // distribution gives 0.5, and a worse-than-contaminant target falls below 0.5.
    const double selectivityScore = ratio / (1.0 + ratio);

    const double remaining = 1.0 /
                             std::max(1.0 + representativeD * organicVolume /
                                                std::max(aqueousVolume, kTiny),
                                      1.0);
    int contacts = 0;
    if (remaining > 0.0 && remaining < 1.0) {
      // After n identical fresh-solvent contacts the aqueous fraction is q^n;
      // q^n <= 0.01 gives n >= log(0.01)/log(q).
      contacts = static_cast<int>(std::ceil(std::log(0.01) / std::log(remaining)));
    }

    candidate.selectivity = ratio;
    candidate.recoveryFraction = recovery;
    candidate.targetSolubilityGPerMl = weightedMean(roles.targets, targetOrganic);
    candidate.contaminantSolubilityGPerMl = worstContaminantS;
    std::string selectDetail = "Target D " + number(representativeD, 2) +
                               "; worst contaminant " + worstName + " D " +
                               number(worstContaminantD, 2) + ", selectivity " +
                               number(ratio, 2) + "x.";
    if (!roles.contaminants.empty() && ratio < 1.0) {
      const std::string failure =
          "Separation direction fails: " + worstName +
          " is preferentially extracted into the solvent by " +
          number(1.0 / std::max(ratio, kTiny), 2) + "x over the target.";
      selectDetail += " " + failure;
      addWarning(candidate.warnings, failure);
    }
    const std::string recoveryDetail =
        number(100.0 * recovery, 1) + "% target in the organic phase after one contact; " +
        (contacts > 0 ? std::to_string(contacts) + " contact(s) predicted for 99%."
                      : "99% is not reachable with the entered phase volumes.");
    CommonScores common = commonScores(solvent, spec.temperatureC, candidate.warnings);
    finishCandidate(candidate, spec, selectivityScore, selectDetail, recovery,
                    recoveryDetail, common, !roles.contaminants.empty());
    appendEvidence(candidate, evidence);
    result.push_back(std::move(candidate));
  }
  return result;
}

std::vector<SolventCandidate> rankRecrystallisation(
    const OperationSpec& spec, const Roles& roles,
    const std::vector<EligibleSolvent>& eligible) {
  std::vector<SolventCandidate> result;
  const double cold = std::min(spec.coldTemperatureC, spec.hotTemperatureC);
  const double hot = std::max(spec.coldTemperatureC, spec.hotTemperatureC);
  for (const EligibleSolvent& entry : eligible) {
    const Solvent& solvent = *entry.solvent;
    // A recrystallisation solvent must stay liquid from the cold endpoint through
    // the dissolve temperature; equality at either phase boundary is not operationally safe.
    if (solvent.meltingPointC >= cold || solvent.boilingPoint <= hot) continue;

    SolventCandidate candidate;
    candidate.solvent = &solvent;
    candidate.warnings = entry.warnings;
    if (solvent.meltingPointC > cold - 5.0) {
      addWarning(candidate.warnings, "Melting point is within 5 C of the cold endpoint.");
    }
    if (solvent.boilingPoint < hot + 10.0) {
      addWarning(candidate.warnings, "Boiling point is within 10 C of the hot endpoint.");
    }
    Evidence evidence;
    std::vector<double> targetHot;
    std::vector<double> targetCold;
    std::vector<double> targetYield;
    for (const SpeciesRole* role : roles.targets) {
      const Prediction hotPrediction =
          safePredict(*role, {{&solvent, 1.0}}, hot, spec.pH, &evidence.warnings);
      const Prediction coldPrediction =
          safePredict(*role, {{&solvent, 1.0}}, cold, spec.pH, &evidence.warnings);
      evidence.observe(*role, hotPrediction);
      evidence.observe(*role, coldPrediction);
      const double sHot = predictionValue(hotPrediction);
      const double sCold = predictionValue(coldPrediction);
      targetHot.push_back(sHot);
      targetCold.push_back(sCold);
      // Cooling a saturated hot solution leaves S_cold/S_hot in the mother liquor,
      // hence the ideal crystal yield is 1 - S_cold/S_hot (fractional mass balance).
      targetYield.push_back(unit(1.0 - sCold / std::max(sHot, kTiny)));
    }
    const double meanTargetHot = weightedMean(roles.targets, targetHot);
    const double meanTargetCold = weightedMean(roles.targets, targetCold);
    const double yield = weightedMean(roles.targets, targetYield);

    std::vector<double> impurityPurity;
    std::vector<double> contaminantCold;
    for (const SpeciesRole* role : roles.contaminants) {
      const Prediction hotPrediction =
          safePredict(*role, {{&solvent, 1.0}}, hot, spec.pH, &evidence.warnings);
      const Prediction coldPrediction =
          safePredict(*role, {{&solvent, 1.0}}, cold, spec.pH, &evidence.warnings);
      evidence.observe(*role, hotPrediction);
      evidence.observe(*role, coldPrediction);
      const double sHot = predictionValue(hotPrediction);
      const double sCold = predictionValue(coldPrediction);
      contaminantCold.push_back(sCold);
      // An impurity is rejected by either limiting case: high S_cold keeps it in
      // mother liquor, while low S_hot relative to target means it never dissolved.
      // a/(a+b) maps each physical solubility ratio onto [0,1] with equality at 0.5.
      const double staysDissolved = sCold / std::max(sCold + meanTargetCold, kTiny);
      const double neverDissolved = meanTargetHot / std::max(meanTargetHot + sHot, kTiny);
      impurityPurity.push_back(std::max(staysDissolved, neverDissolved));
    }
    const double purity = roles.contaminants.empty()
                              ? 1.0
                              : weightedMean(roles.contaminants, impurityPurity);
    const double purityOdds = purity / std::max(1.0 - purity, kTiny);

    candidate.selectivity = purityOdds;
    candidate.recoveryFraction = yield;
    candidate.targetSolubilityGPerMl = meanTargetHot;
    candidate.contaminantSolubilityGPerMl =
        weightedMean(roles.contaminants, contaminantCold);
    const double targetHotCold =
        meanTargetHot / std::max(meanTargetCold, kTiny);
    std::string selectDetail =
        roles.contaminants.empty()
            ? "No contaminants specified; target hot/cold solubility ratio " +
                  number(targetHotCold, 2) + "x."
            : "Target hot/cold solubility ratio " + number(targetHotCold, 2) +
                  "x; contaminant cold solubility " +
                  number(candidate.contaminantSolubilityGPerMl, 4) +
                  " g/mL gives impurity rejection " + number(100.0 * purity, 1) + "%.";
    if (!roles.contaminants.empty() && purityOdds < 1.0) {
      const std::string failure =
          "Separation direction fails: the contaminant is preferentially crystallised by " +
          number(1.0 / std::max(purityOdds, kTiny), 2) +
          "x instead of remaining in the mother liquor.";
      selectDetail += " " + failure;
      addWarning(candidate.warnings, failure);
    }
    const std::string recoveryDetail =
        "Target S_hot " + number(meanTargetHot, 4) + " g/mL, S_cold " +
        number(meanTargetCold, 4) + " g/mL; ideal cooling yield " +
        number(100.0 * yield, 1) + "%.";
    CommonScores common = commonScores(solvent, spec.temperatureC, candidate.warnings);
    finishCandidate(candidate, spec, purity, selectDetail, yield, recoveryDetail, common,
                    !roles.contaminants.empty());
    appendEvidence(candidate, evidence);
    result.push_back(std::move(candidate));
  }
  return result;
}

std::vector<SolventCandidate> rankTrituration(const OperationSpec& spec, const Roles& roles,
                                              const std::vector<EligibleSolvent>& eligible) {
  std::vector<SolventCandidate> result;
  const double washVolume = std::max(0.0, spec.organicVolumeMl);
  for (const EligibleSolvent& entry : eligible) {
    const Solvent& solvent = *entry.solvent;
    SolventCandidate candidate;
    candidate.solvent = &solvent;
    candidate.warnings = entry.warnings;
    Evidence evidence;
    std::vector<double> targetSolubility;
    std::vector<double> targetLoss;
    for (const SpeciesRole* role : roles.targets) {
      const Prediction prediction = safePredict(*role, {{&solvent, 1.0}}, spec.temperatureC,
                                                spec.pH, &evidence.warnings);
      evidence.observe(*role, prediction);
      const double solubility = predictionValue(prediction);
      targetSolubility.push_back(solubility);
      if (role->amountMg > 0.0) {
        // Dissolved capacity is S[g/mL]*V[mL]*1000[mg/g]; dividing by charged
        // mass gives the maximum target-loss fraction from one saturated wash.
        targetLoss.push_back(unit(solubility * washVolume * 1000.0 / role->amountMg));
      } else {
        // With no charge entered, S/(S+0.01 g/mL) is a bounded loss proxy whose
        // midpoint is the stated 10 mg/mL practical-solubility reference.
        targetLoss.push_back(solubility / (solubility + kReferenceSolubility));
      }
    }
    std::vector<double> contaminantSolubility;
    for (const SpeciesRole* role : roles.contaminants) {
      const Prediction prediction = safePredict(*role, {{&solvent, 1.0}}, spec.temperatureC,
                                                spec.pH, &evidence.warnings);
      evidence.observe(*role, prediction);
      contaminantSolubility.push_back(predictionValue(prediction));
    }
    const double targetS = weightedMean(roles.targets, targetSolubility);
    const double contaminantS = weightedMean(roles.contaminants, contaminantSolubility);
    const double loss = weightedMean(roles.targets, targetLoss);
    const double ratio = roles.contaminants.empty()
                             ? kReferenceSolubility / std::max(targetS, kTiny)
                             : contaminantS / std::max(targetS, kTiny);
    // Odds mapping makes equal target/contaminant solubility score 0.5 and
    // directly enforces the inverse selectivity required for a solid wash.
    const double selectivityScore = ratio / (1.0 + ratio);
    const double recovery = unit(1.0 - loss);

    candidate.selectivity = ratio;
    candidate.recoveryFraction = recovery;
    candidate.targetSolubilityGPerMl = targetS;
    candidate.contaminantSolubilityGPerMl = contaminantS;
    std::string selectDetail =
        "Contaminant/target solubility ratio " + number(ratio, 2) + "x (target " +
        number(targetS, 4) + ", contaminants " + number(contaminantS, 4) + " g/mL).";
    if (!roles.contaminants.empty() && ratio < 1.0) {
      const std::string failure =
          "Separation direction fails: the target is preferentially dissolved by " +
          number(1.0 / std::max(ratio, kTiny), 2) +
          "x over the contaminant, so the wash would remove product.";
      selectDetail += " " + failure;
      addWarning(candidate.warnings, failure);
    }
    const std::string recoveryDetail =
        "Predicted target loss in one " + number(washVolume, 1) + " mL wash: " +
        number(100.0 * loss, 1) +
        (std::any_of(roles.targets.begin(), roles.targets.end(),
                     [](const SpeciesRole* role) { return role->amountMg <= 0.0; })
             ? "% (10 mg/mL reference used where charge was unspecified)."
             : "% from entered charges.");
    CommonScores common = commonScores(solvent, spec.temperatureC, candidate.warnings);
    finishCandidate(candidate, spec, selectivityScore, selectDetail, recovery,
                    recoveryDetail, common, !roles.contaminants.empty());
    appendEvidence(candidate, evidence);
    result.push_back(std::move(candidate));
  }
  return result;
}
bool candidateRanksAhead(const SolventCandidate& a, const SolventCandidate& b,
                         OperationKind kind, bool applySeparationGate) {
  if (applySeparationGate && isSeparationOperation(kind)) {
    const int aTier = separationTier(a);
    const int bTier = separationTier(b);
    if (aTier != bTier) return aTier > bTier;
    // If every candidate reverses the separation, closeness to the 1x boundary
    // is more useful than a falsely reassuring green or handling score.
    if (aTier == 0 && a.selectivity != b.selectivity) {
      return a.selectivity > b.selectivity;
    }
  }
  if (a.score != b.score) return a.score > b.score;
  const std::string aName = a.solvent ? a.solvent->name : std::string();
  const std::string bName = b.solvent ? b.solvent->name : std::string();
  if (aName != bName) return aName < bName;
  const std::string aPartner = a.partner ? a.partner->name : std::string();
  const std::string bPartner = b.partner ? b.partner->name : std::string();
  if (aPartner != bPartner) return aPartner < bPartner;
  return a.partnerFraction < b.partnerFraction;
}

std::vector<SolventCandidate> rankAntiSolvent(const OperationSpec& spec, const Roles& roles,
                                              const std::vector<EligibleSolvent>& eligible) {
  std::vector<SolventCandidate> result;
  if (eligible.size() < 2) return result;

  struct PureRow {
    size_t eligibleIndex = 0;
    std::vector<Prediction> targetPredictions;
    std::vector<Prediction> contaminantPredictions;
    double targetSolubility = 0.0;
    double contaminantSolubility = 0.0;
  };
  std::vector<PureRow> pure;
  pure.reserve(eligible.size());
  for (size_t i = 0; i < eligible.size(); ++i) {
    PureRow row;
    row.eligibleIndex = i;
    std::vector<double> values;
    for (const SpeciesRole* role : roles.targets) {
      row.targetPredictions.push_back(
          safePredict(*role, {{eligible[i].solvent, 1.0}}, spec.temperatureC, spec.pH));
      values.push_back(predictionValue(row.targetPredictions.back()));
    }
    row.targetSolubility = weightedMean(roles.targets, values);
    std::vector<double> contaminantValues;
    for (const SpeciesRole* role : roles.contaminants) {
      row.contaminantPredictions.push_back(
          safePredict(*role, {{eligible[i].solvent, 1.0}}, spec.temperatureC, spec.pH));
      contaminantValues.push_back(predictionValue(row.contaminantPredictions.back()));
    }
    row.contaminantSolubility =
        weightedMean(roles.contaminants, contaminantValues);
    pure.push_back(std::move(row));
  }

  std::vector<size_t> primaryOrder(pure.size());
  std::iota(primaryOrder.begin(), primaryOrder.end(), 0);
  std::stable_sort(primaryOrder.begin(), primaryOrder.end(), [&](size_t a, size_t b) {
    if (pure[a].targetSolubility != pure[b].targetSolubility) {
      return pure[a].targetSolubility > pure[b].targetSolubility;
    }
    return eligible[pure[a].eligibleIndex].solvent->name <
           eligible[pure[b].eligibleIndex].solvent->name;
  });
  // Pure-solvent pre-ranking caps the UI path at 8*12*4 = 384 blend points;
  // the four fractions span onset through a strong anti-solvent excess without a fine sweep.
  constexpr size_t kPrimaryCap = 8;
  constexpr size_t kPartnerCap = 12;
  constexpr double kFractions[] = {0.25, 0.50, 0.75, 0.90};
  if (primaryOrder.size() > kPrimaryCap) primaryOrder.resize(kPrimaryCap);

  for (size_t primaryRowIndex : primaryOrder) {
    const PureRow& primaryRow = pure[primaryRowIndex];
    const EligibleSolvent& primaryEntry = eligible[primaryRow.eligibleIndex];
    const Solvent& primary = *primaryEntry.solvent;

    std::vector<size_t> partnerOrder;
    for (size_t i = 0; i < pure.size(); ++i) {
      const Solvent& partner = *eligible[pure[i].eligibleIndex].solvent;
      if (&partner != &primary && miscibleWith(primary, partner)) partnerOrder.push_back(i);
    }
    std::stable_sort(partnerOrder.begin(), partnerOrder.end(), [&](size_t a, size_t b) {
      if (pure[a].targetSolubility != pure[b].targetSolubility) {
        return pure[a].targetSolubility < pure[b].targetSolubility;
      }
      return eligible[pure[a].eligibleIndex].solvent->name <
             eligible[pure[b].eligibleIndex].solvent->name;
    });
    if (partnerOrder.size() > kPartnerCap) partnerOrder.resize(kPartnerCap);

    bool found = false;
    SolventCandidate best;
    for (size_t partnerRowIndex : partnerOrder) {
      const EligibleSolvent& partnerEntry = eligible[pure[partnerRowIndex].eligibleIndex];
      const Solvent& partner = *partnerEntry.solvent;
      for (double fraction : kFractions) {
        Evidence evidence;
        std::vector<double> targetBlend;
        std::vector<double> targetRecovery;
        for (size_t i = 0; i < roles.targets.size(); ++i) {
          const SpeciesRole& role = *roles.targets[i];
          evidence.observe(role, primaryRow.targetPredictions[i]);
          const Prediction blendPrediction =
              safePredict(role, {{&primary, 1.0 - fraction}, {&partner, fraction}},
                          spec.temperatureC, spec.pH, &evidence.warnings);
          evidence.observe(role, blendPrediction);
          const double good = predictionValue(primaryRow.targetPredictions[i]);
          const double mixed = predictionValue(blendPrediction);
          targetBlend.push_back(mixed);
          // On adding anti-solvent at fixed dissolved charge, the fraction beyond
          // the new saturation capacity is 1 - S_blend/S_good (ideal volume basis).
          targetRecovery.push_back(unit(1.0 - mixed / std::max(good, kTiny)));
        }
        std::vector<double> contaminantBlend;
        for (size_t i = 0; i < roles.contaminants.size(); ++i) {
          const SpeciesRole& role = *roles.contaminants[i];
          evidence.observe(role, primaryRow.contaminantPredictions[i]);
          const Prediction prediction =
              safePredict(role, {{&primary, 1.0 - fraction}, {&partner, fraction}},
                          spec.temperatureC, spec.pH, &evidence.warnings);
          evidence.observe(role, prediction);
          contaminantBlend.push_back(predictionValue(prediction));
        }

        const double mixedTarget = weightedMean(roles.targets, targetBlend);
        const double mixedContaminant =
            weightedMean(roles.contaminants, contaminantBlend);
        const double recovery = weightedMean(roles.targets, targetRecovery);
        const double supersaturation =
            primaryRow.targetSolubility / std::max(mixedTarget, kTiny);
        // Selectivity compares the desired target collapse with unwanted
        // contaminant collapse. A ratio below one therefore precipitates impurity
        // preferentially, while complete contaminant carry-over gives a large ratio.
        const double contaminantCarryover =
            roles.contaminants.empty()
                ? 1.0
                : unit(mixedContaminant /
                       std::max(primaryRow.contaminantSolubility, kTiny));
        const double contaminantCollapse = unit(1.0 - contaminantCarryover);
        const double directionRatio =
            roles.contaminants.empty()
                ? supersaturation
                : recovery / std::max(contaminantCollapse, kTiny);
        const double selectivityScore = directionRatio / (1.0 + directionRatio);

        SolventCandidate candidate;
        candidate.solvent = &primary;
        candidate.partner = &partner;
        candidate.partnerFraction = fraction;
        candidate.warnings = primaryEntry.warnings;
        for (const std::string& warning : partnerEntry.warnings) {
          addWarning(candidate.warnings, "Anti-solvent: " + warning);
        }
        candidate.selectivity = directionRatio;
        candidate.recoveryFraction = recovery;
        candidate.targetSolubilityGPerMl = mixedTarget;
        candidate.contaminantSolubilityGPerMl = mixedContaminant;
        std::string selectDetail =
            "At " + number(100.0 * fraction, 0) + "% " + partner.name +
            ", target collapse is " + number(100.0 * recovery, 1) +
            "% and contaminant carry-over is " +
            number(100.0 * contaminantCarryover, 1) +
            "%; precipitation selectivity " + number(directionRatio, 2) + "x.";
        if (!roles.contaminants.empty() && directionRatio < 1.0) {
          const std::string failure =
              "Separation direction fails: the contaminant is preferentially precipitated by " +
              number(1.0 / std::max(directionRatio, kTiny), 2) +
              "x over the target.";
          selectDetail += " " + failure;
          addWarning(candidate.warnings, failure);
        }
        const std::string recoveryDetail =
            "Solubility falls from " + number(primaryRow.targetSolubility, 4) + " to " +
            number(mixedTarget, 4) + " g/mL; ideal precipitated fraction " +
            number(100.0 * recovery, 1) + "%.";
        CommonScores common =
            pairCommonScores(primary, partner, spec.temperatureC, candidate.warnings);
        finishCandidate(candidate, spec, selectivityScore, selectDetail, recovery,
                        recoveryDetail, common, !roles.contaminants.empty());
        appendEvidence(candidate, evidence);
        if (!found ||
            candidateRanksAhead(candidate, best,
                                OperationKind::AntiSolventPrecipitation,
                                !roles.contaminants.empty())) {
          best = std::move(candidate);
          found = true;
        }
      }
    }
    if (found) result.push_back(std::move(best));
  }
  return result;
}

std::vector<SolventCandidate> rankChromatography(
    const OperationSpec& spec, const Roles& roles,
    const std::vector<EligibleSolvent>& eligible) {
  std::vector<SolventCandidate> result;
  for (const EligibleSolvent& entry : eligible) {
    const Solvent& solvent = *entry.solvent;
    SolventCandidate candidate;
    candidate.solvent = &solvent;
    candidate.warnings = entry.warnings;
    Evidence evidence;
    std::vector<double> targetIndex;
    std::vector<double> targetSolubility;
    std::vector<double> contaminantIndex;
    std::vector<double> contaminantSolubility;

    auto indexFor = [&](const SpeciesRole& role, std::vector<double>& indices,
                        std::vector<double>& concentrations) {
      const Prediction prediction = safePredict(role, {{&solvent, 1.0}}, spec.temperatureC,
                                                spec.pH, &evidence.warnings);
      evidence.observe(role, prediction);
      concentrations.push_back(predictionValue(prediction));
      const double red = std::max(0.0, prediction.relativeEnergyDifference);
      // exp(-RED^2/2) is a unit Gaussian similarity in Hansen-radius units;
      // eps/(eps+20) is a bounded dielectric polarity term with midpoint eps=20.
      // Their 3:1 blend is explicitly a polarity heuristic, not a fitted retention model.
      const double hansenAffinity = std::exp(-0.5 * red * red);
      const double dielectricPolarity =
          std::max(0.0, solvent.dielectric) / (std::max(0.0, solvent.dielectric) + 20.0);
      indices.push_back(unit(0.75 * hansenAffinity + 0.25 * dielectricPolarity));
    };
    for (const SpeciesRole* role : roles.targets) {
      indexFor(*role, targetIndex, targetSolubility);
    }
    for (const SpeciesRole* role : roles.contaminants) {
      indexFor(*role, contaminantIndex, contaminantSolubility);
    }

    const double target = weightedMean(roles.targets, targetIndex);
    const double contaminant = weightedMean(roles.contaminants, contaminantIndex);
    const double separation = std::abs(target - contaminant);
    const bool opposite = !roles.contaminants.empty() && (target - 0.5) * (contaminant - 0.5) < 0.0;
    // The useful target window is centred at 0.5. Linear distance to either end
    // gives an elution-window score of 1 at centre and 0 at affinity 0 or 1.
    const double targetWindow = unit(1.0 - std::abs(target - 0.5) / 0.5);
    // A 0.10 index gap is the minimum useful heuristic window. Expressing the
    // observed gap as a ratio makes <1 an explicit chromatography failure,
    // instead of allowing a same-window pair to win on solvent convenience.
    const double retentionRatio =
        roles.contaminants.empty()
            ? targetWindow / 0.5
            : separation / kUsefulChromatographySpacing;
    const double selectivityScore = retentionRatio / (1.0 + retentionRatio);

    candidate.selectivity = retentionRatio;
    candidate.recoveryFraction = targetWindow;
    candidate.targetSolubilityGPerMl = weightedMean(roles.targets, targetSolubility);
    candidate.contaminantSolubilityGPerMl =
        weightedMean(roles.contaminants, contaminantSolubility);
    std::string selectDetail =
        "Polarity heuristic only, not a retention-factor prediction: target index " +
        number(target) + ", contaminants " + number(contaminant) +
        ", index gap " + number(separation) + " (" + number(retentionRatio, 2) +
        "x the 0.10 minimum useful retention window)" +
        (opposite ? ", on opposite sides of the 0.50 window." : ", on the same side.");
    if (!roles.contaminants.empty() && retentionRatio < 1.0) {
      const std::string failure =
          "Separation direction fails: the target and contaminant co-elute within the same "
          "retention window; the gap is only " +
          number(retentionRatio, 2) + "x the minimum useful spacing.";
      selectDetail += " " + failure;
      addWarning(candidate.warnings, failure);
    }
    const std::string recoveryDetail =
        "Target mobile-phase affinity window score " + number(targetWindow) +
        " from Hansen RED and dielectric polarity.";
    CommonScores common = commonScores(solvent, spec.temperatureC, candidate.warnings);
    finishCandidate(candidate, spec, selectivityScore, selectDetail, targetWindow,
                    recoveryDetail, common, !roles.contaminants.empty());
    appendEvidence(candidate, evidence);
    result.push_back(std::move(candidate));
  }
  return result;
}

std::vector<SolventCandidate> rankReactionMedium(
    const OperationSpec& spec, const Roles& roles,
    const std::vector<EligibleSolvent>& eligible) {
  std::vector<SolventCandidate> result;
  for (const EligibleSolvent& entry : eligible) {
    const Solvent& solvent = *entry.solvent;
    SolventCandidate candidate;
    candidate.solvent = &solvent;
    candidate.warnings = entry.warnings;
    Evidence evidence;
    std::vector<double> allSolubility;
    std::vector<double> dissolveFractions;
    std::vector<double> targetSolubility;
    std::vector<double> contaminantSolubility;
    double totalWeight = 0.0;
    double weightedLog = 0.0;
    double weightedLogSquared = 0.0;

    for (const SpeciesRole& role : spec.species) {
      const Prediction prediction = safePredict(role, {{&solvent, 1.0}}, spec.temperatureC,
                                                spec.pH, &evidence.warnings);
      evidence.observe(role, prediction);
      const double solubility = predictionValue(prediction);
      allSolubility.push_back(solubility);
      (role.keep ? targetSolubility : contaminantSolubility).push_back(solubility);
      if (role.amountMg > 0.0) {
        // Capacity S*V*1000 divided by charge is the maximum dissolved fraction;
        // organicVolumeMl is the entered working-solvent volume for this operation.
        dissolveFractions.push_back(unit(solubility * std::max(0.0, spec.organicVolumeMl) *
                                         1000.0 / role.amountMg));
      } else {
        dissolveFractions.push_back(solubility / (solubility + kReferenceSolubility));
      }
      const double weight = roleWeight(role);
      const double logS = std::log10(std::max(solubility, kTiny));
      totalWeight += weight;
      weightedLog += weight * logS;
      weightedLogSquared += weight * logS * logS;
    }

    const double minimumDissolved = dissolveFractions.empty()
                                        ? 0.0
                                        : *std::min_element(dissolveFractions.begin(),
                                                            dissolveFractions.end());
    const double meanLog = totalWeight > 0.0 ? weightedLog / totalWeight : 0.0;
    const double variance = totalWeight > 0.0
                                ? std::max(0.0, weightedLogSquared / totalWeight - meanLog * meanLog)
                                : 0.0;
    // A one-decade standard deviation consumes one sixth of the consistency
    // score; six decades spans the model's practical display range.
    const double consistency = unit(1.0 - std::sqrt(variance) / 6.0);

    if (isProtic(solvent)) {
      addWarning(candidate.warnings,
                 "Protic solvent; compatibility with strongly basic or organometallic chemistry requires chemist review.");
    } else {
      addWarning(candidate.warnings,
                 "Aprotic solvent; chemical inertness toward the drawn species is not predicted.");
    }
    candidate.selectivity = consistency;
    candidate.recoveryFraction = minimumDissolved;
    candidate.targetSolubilityGPerMl = weightedMean(roles.targets, targetSolubility);
    candidate.contaminantSolubilityGPerMl =
        weightedMean(roles.contaminants, contaminantSolubility);
    const std::string selectDetail =
        "All-species dissolution consistency " + number(consistency) +
        " (selectivity is not desired in a reaction medium).";
    const std::string recoveryDetail =
        "Worst species dissolved fraction " + number(100.0 * minimumDissolved, 1) +
        "% at the entered charge/volume; unspecified charges use the 10 mg/mL reference.";
    CommonScores common = commonScores(solvent, spec.temperatureC, candidate.warnings);
    finishCandidate(candidate, spec, consistency, selectDetail, minimumDissolved,
                    recoveryDetail, common, false);
    appendEvidence(candidate, evidence);
    result.push_back(std::move(candidate));
  }
  return result;
}

void rankDeterministically(std::vector<SolventCandidate>& candidates,
                           OperationKind kind, bool applySeparationGate) {
  std::stable_sort(
      candidates.begin(), candidates.end(),
      [kind, applySeparationGate](const SolventCandidate& a,
                                  const SolventCandidate& b) {
        return candidateRanksAhead(a, b, kind, applySeparationGate);
      });
}

}  // namespace

std::vector<SolventCandidate> rankSolvents(const OperationSpec& spec) {
  const Roles roles = splitRoles(spec);
  if (roles.targets.empty()) return {};

  // Load once before model calls. This is the only uncaught SolError path: a
  // malformed/missing solvent database is fatal, while optional prediction data degrades visibly.
  const std::vector<Solvent>& database = solvents();
  const std::vector<EligibleSolvent> eligible = eligibleSolvents(database, spec);
  std::vector<SolventCandidate> result;
  switch (spec.kind) {
    case OperationKind::LiquidLiquidExtraction:
      result = rankExtraction(spec, roles, eligible);
      break;
    case OperationKind::Recrystallisation:
      result = rankRecrystallisation(spec, roles, eligible);
      break;
    case OperationKind::Trituration:
      result = rankTrituration(spec, roles, eligible);
      break;
    case OperationKind::AntiSolventPrecipitation:
      result = rankAntiSolvent(spec, roles, eligible);
      break;
    case OperationKind::ChromatographyMobilePhase:
      result = rankChromatography(spec, roles, eligible);
      break;
    case OperationKind::ReactionMedium:
      result = rankReactionMedium(spec, roles, eligible);
      break;
  }
  rankDeterministically(result, spec.kind, !roles.contaminants.empty());
  return result;
}

const char* operationName(OperationKind kind) {
  switch (kind) {
    case OperationKind::LiquidLiquidExtraction:
      return "Liquid-liquid extraction";
    case OperationKind::Recrystallisation:
      return "Recrystallisation";
    case OperationKind::Trituration:
      return "Trituration";
    case OperationKind::AntiSolventPrecipitation:
      return "Anti-solvent precipitation";
    case OperationKind::ChromatographyMobilePhase:
      return "Chromatography mobile phase";
    case OperationKind::ReactionMedium:
      return "Reaction medium";
  }
  return "Unknown operation";
}

const char* operationDescription(OperationKind kind) {
  switch (kind) {
    case OperationKind::LiquidLiquidExtraction:
      return "Pull targets into a water-immiscible phase while leaving contaminants behind.";
    case OperationKind::Recrystallisation:
      return "Dissolve hot, crystallise targets cold, and leave impurities in the mother liquor.";
    case OperationKind::Trituration:
      return "Wash a solid so contaminants dissolve while the desired material remains.";
    case OperationKind::AntiSolventPrecipitation:
      return "Add a miscible poor solvent to collapse target solubility and precipitate it.";
    case OperationKind::ChromatographyMobilePhase:
      return "Rank polarity heuristics for useful target/contaminant elution separation.";
    case OperationKind::ReactionMedium:
      return "Dissolve every species within the requested working and boiling range.";
  }
  return "";
}

}  // namespace chemcad::sol

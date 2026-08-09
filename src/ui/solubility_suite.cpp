// Solubility Suite: the workspace for every "which solvent, and how much
// dissolves" question. Two modes share it. Predict holds solute + solvent
// inputs on the left rail with the prediction and the solubility-vs-composition
// graph dominant on the right. Select hosts the operation-first solvent ranking
// (ui/solvent_selector.cpp), which hands its winning blend straight back into
// Predict rather than across a panel boundary. Pure solvents are screened in a
// ranked table; a chosen blend can be handed to the Extraction Calculator.
// Renders inside an already-open window; never owns one.

#include "imgui.h"

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "chem/bridge.hpp"
#include "core/model.hpp"
#include "sol/kirkwood_buff.hpp"
#include "sol/solubility.hpp"
#include "sol/solvent.hpp"
#include "ui/app_state.hpp"
#include "ui/charts.hpp"
#include "ui/charts3d.hpp"
#include "ui/icons.hpp"
#include "ui/layout.hpp"
#include "ui/solubility_state.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"
#include "ui/widgets.hpp"

namespace chemcad::ui {
namespace {



// Formats a bare value across the many orders of magnitude a solubility
// prediction can span: scientific notation below 1e-3, fixed otherwise.
std::string formatSolubility(double value) {
  char buf[64];
  const double magnitude = std::fabs(value);
  if (value == 0.0) {
    return "0";
  } else if (magnitude < 1e-3 || magnitude >= 1e6) {
    std::snprintf(buf, sizeof(buf), "%.3e", value);
  } else if (magnitude < 1.0) {
    std::snprintf(buf, sizeof(buf), "%.5f", value);
  } else {
    std::snprintf(buf, sizeof(buf), "%.4g", value);
  }
  return buf;
}

// ------------------------------------------------------------------- units
// Predictions are computed in g/mL; every readout in the panel (hero card,
// axes, tooltips, screening table) goes through these two helpers so the
// units combo is honoured everywhere at once.
const std::array<const char*, 4> kUnitLabels = {"g/mL", "mg/mL", "g/100 mL", "mol/L"};

double toDisplayUnits(double gramsPerMillilitre, double molarMass, int units) {
  switch (units) {
    case 1: return gramsPerMillilitre * 1000.0;             // mg/mL
    case 2: return gramsPerMillilitre * 100.0;              // g/100 mL
    case 3: return molarMass > 1e-9 ? gramsPerMillilitre / molarMass * 1000.0 : 0.0;
    default: return gramsPerMillilitre;
  }
}

std::string formatUnits(double gramsPerMillilitre, double molarMass, int units) {
  const int clamped = std::clamp(units, 0, 3);
  return formatSolubility(toDisplayUnits(gramsPerMillilitre, molarMass, clamped)) + " " +
         kUnitLabels[static_cast<size_t>(clamped)];
}


// Slot colours report identity; interactive selection remains amber inside the controls.
ImVec4 slotColor(int index) {
  switch (index) {
    case 0: return style::col::Data;
    case 1: return style::col::DataBright;
    default: return style::col::DataDim;
  }
}

const sol::Electrolyte* activeBackground(const SolubilityState& sb);


// -------------------------------------------------------------- solute cache

void applySoluteOverrides(SolubilityState& sb) {
  if (!sb.soluteValid) return;
  sb.solute.meltingPoint = static_cast<double>(sb.meltingPointC);
  sb.solute.interactionRadius = static_cast<double>(sb.interactionRadius);
}

// Note describing how a multi-fragment sketch was interpreted, or empty when
// there is nothing surprising to say.
std::string soluteInterpretationNote(const sol::Solute& solute) {
  const sol::Ionization& ion = solute.ionization;
  if (ion.fragmentCount <= 1) return {};
  char note[224];
  if (ion.saltForm && !ion.counterIon.empty()) {
    std::snprintf(note, sizeof(note),
                  "Salt: %s counter-ion detected, modelled as the ionised form of the "
                  "%.1f g/mol skeleton.",
                  ion.counterIon.c_str(), solute.molarMass);
  } else {
    std::snprintf(note, sizeof(note),
                  "Sketch holds %d fragments; modelling the largest (%.1f g/mol). "
                  "Spectator molecules are ignored.",
                  ion.fragmentCount, solute.molarMass);
  }
  return note;
}

void recomputeSoluteFromSketch(AppState& st) {
  SolubilityState& sb = st.solubility;
  sb.soluteError.clear();
  sb.soluteNote.clear();
  sb.soluteValid = false;
  try {
    // Hand the WHOLE sketch over, not the first fragment: describeSolute picks
    // the skeleton and reads the rest as counter-ions, so a drawn
    // hydrochloride is modelled as a salt instead of as its free base, and a
    // reagent parked beside the target cannot become the solute.
    core::Molecule structure;
    for (const core::Molecule& fragment : st.doc.molecules) {
      if (fragment.empty()) continue;
      std::unordered_map<core::AtomId, core::AtomId> ids;
      ids.reserve(fragment.atomCount());
      for (const core::Atom& atom : fragment.atoms()) ids.emplace(atom.id, structure.addAtom(atom));
      for (const core::Bond& bond : fragment.bonds()) {
        const core::BondId id = structure.addBond(ids.at(bond.a), ids.at(bond.b), bond.order);
        if (core::Bond* added = structure.bond(id)) added->stereo = bond.stereo;
      }
    }
    if (structure.empty()) throw sol::SolError("Sketch is empty -- draw a structure first.");
    sb.solute = sol::describeSolute(structure);
    sb.soluteValid = true;
    sb.soluteNote = soluteInterpretationNote(sb.solute);
  } catch (const sol::SolError& err) {
    sb.soluteError = err.what();
  } catch (const chem::ChemError& err) {
    sb.soluteError = err.what();
  } catch (const std::exception& err) {
    sb.soluteError = err.what();
  }
  sb.sourceRevision = st.docRevision;
  if (sb.overrideSolute) applySoluteOverrides(sb);
  ++sb.soluteVersion;
}

void recomputeSoluteFromSmiles(AppState& st) {
  SolubilityState& sb = st.solubility;
  sb.soluteError.clear();
  sb.soluteNote.clear();
  sb.soluteValid = false;
  if (sb.manualSmiles.empty()) {
    sb.soluteError = "Enter a SMILES string.";
  } else {
    try {
      const core::Molecule mol = chem::fromSmiles(sb.manualSmiles);
      sb.solute = sol::describeSolute(mol);
      sb.soluteValid = true;
      sb.soluteNote = soluteInterpretationNote(sb.solute);
    } catch (const sol::SolError& err) {
      sb.soluteError = err.what();
    } catch (const chem::ChemError& err) {
      sb.soluteError = err.what();
    } catch (const std::exception& err) {
      sb.soluteError = err.what();
    }
  }
  if (sb.overrideSolute) applySoluteOverrides(sb);
  ++sb.soluteVersion;
}

// -------------------------------------------------------------- sweep cache

const sol::Electrolyte* activeBackground(const SolubilityState& sb) {
  if (!sb.backgroundEnabled) return nullptr;
  const std::vector<sol::Electrolyte>& all = sol::electrolytes();
  if (all.empty()) return nullptr;
  const size_t index = static_cast<size_t>(std::clamp(sb.backgroundElectrolyte, 0,
                                                      static_cast<int>(all.size()) - 1));
  return &all[index];
}

// pH the model should use: the sentinel when the solute self-buffers.
double requestedPH(const SolubilityState& sb) {
  return sb.pHAuto ? sol::kAutoPH : static_cast<double>(sb.pH);
}

std::array<int, 2> ternarySurfaceResolution(int sweepSteps) {
  const float t =
      static_cast<float>(std::clamp(sweepSteps, 2, 64) - 2) / 62.0f;
  const int ratioQuads =
      std::clamp(static_cast<int>(std::lround(24.0f + t * 24.0f)), 24, 48);
  const int cQuads =
      std::clamp(static_cast<int>(std::lround(16.0f + t * 16.0f)), 16, 32);
  return {ratioQuads, cQuads};  // at most 48 * 32 = 1536 quads
}

void recomputeSweepIfNeeded(SolubilityState& sb,
                            const std::vector<const sol::Solvent*>& chosen) {
  std::string signature;
  signature.reserve(64);
  for (const sol::Solvent* solvent : chosen) {
    signature += solvent->id;
    signature += '|';
  }
  signature += std::to_string(sb.sweepSteps);
  signature += '|';
  signature += std::to_string(static_cast<int>(sb.temperatureC * 100.0f));
  signature += '|';
  signature += std::to_string(sb.soluteVersion);
  signature += '|';
  signature += sb.pHAuto ? "auto" : std::to_string(static_cast<int>(sb.pH * 10.0f));
  if (const sol::Electrolyte* bg = activeBackground(sb)) {
    signature += '|';
    signature += bg->id;
    signature += std::to_string(static_cast<int>(sb.backgroundMolarity * 1000.0f));
  }
  if (signature == sb.sweepSignature) return;
  sb.sweepSignature = signature;
  sb.sweepPeakIndex = -1;
  sb.ternarySurface = {};
  try {
    if (chosen.size() == 3) {
      TernarySurfaceMesh mesh;
      const std::array<int, 2> resolution =
          ternarySurfaceResolution(sb.sweepSteps);
      mesh.ratioQuads = resolution[0];
      mesh.cQuads = resolution[1];
      const int ratioNodes = mesh.ratioQuads + 1;
      const int cNodes = mesh.cQuads + 1;
      mesh.nodes.reserve(static_cast<size_t>(ratioNodes) *
                         static_cast<size_t>(cNodes));

      const std::vector<const sol::Solvent*> ab{chosen[0], chosen[1]};
      const sol::Electrolyte* background = activeBackground(sb);
      const double temperature = static_cast<double>(sb.temperatureC);
      const double backgroundM = static_cast<double>(sb.backgroundMolarity);
      const double pH = requestedPH(sb);
      std::vector<sol::Component> components{
          {chosen[0], 0.0}, {chosen[1], 0.0}, {chosen[2], 0.0}};
      double minimum = DBL_MAX;
      double maximum = -DBL_MAX;
      double globalBest = -DBL_MAX;

      for (int cIndex = 0; cIndex < cNodes; ++cIndex) {
        const double fractionC =
            static_cast<double>(cIndex) / static_cast<double>(mesh.cQuads);
        const double fractionAB = 1.0 - fractionC;
        std::vector<sol::SweepPoint> row =
            sol::sweep(sb.solute, ab, mesh.ratioQuads, temperature,
                       background, backgroundM, pH);
        for (sol::SweepPoint& point : row) {
          const double ratioA = point.fractions[0];
          const double ratioB = point.fractions[1];
          const double fractionA = ratioA * fractionAB;
          const double fractionB = ratioB * fractionAB;
          point.fractions = {fractionA, fractionB, fractionC};
          components[0].volumeFraction = fractionA;
          components[1].volumeFraction = fractionB;
          components[2].volumeFraction = fractionC;
          point.prediction =
              sol::predict(sb.solute, components, temperature, background,
                           backgroundM, pH);

          TernarySurfaceNode node;
          node.point = std::move(point);
          node.cubeLinear = {static_cast<float>(ratioB), 0.0f,
                             static_cast<float>(fractionC)};
          node.cubeLog = node.cubeLinear;
          const double value = node.point.prediction.gramsPerMillilitre;
          minimum = std::min(minimum, value);
          maximum = std::max(maximum, value);
          if (value > globalBest) {
            globalBest = value;
            mesh.peakNode = static_cast<int>(mesh.nodes.size());
          }
          mesh.nodes.push_back(std::move(node));
        }
      }

      mesh.minimum = minimum;
      mesh.maximum = maximum;
      if (mesh.maximum - mesh.minimum < 1e-12) {
        mesh.maximum =
            mesh.minimum + std::max(1e-9, std::fabs(mesh.minimum) * 0.1 + 1e-9);
      }
      mesh.logFloor = std::max(maximum * 1e-6, 1e-12);
      mesh.logMinimum = std::log10(std::max(minimum, mesh.logFloor));
      mesh.logMaximum = std::log10(std::max(maximum, mesh.logFloor));
      if (mesh.logMaximum - mesh.logMinimum < 1e-9) {
        mesh.logMaximum = mesh.logMinimum + 1.0;
      }

      for (TernarySurfaceNode& node : mesh.nodes) {
        const double value = node.point.prediction.gramsPerMillilitre;
        node.cubeLinear[1] = static_cast<float>(
            (value - mesh.minimum) / (mesh.maximum - mesh.minimum));
        node.cubeLog[1] = static_cast<float>(
            (std::log10(std::max(value, mesh.logFloor)) - mesh.logMinimum) /
            (mesh.logMaximum - mesh.logMinimum));
      }


      sb.sweep.clear();
      if (mesh.peakNode >= 0) {
        sb.sweep.push_back(
            mesh.nodes[static_cast<size_t>(mesh.peakNode)].point);
        sb.sweepPeakIndex = 0;
      }
      sb.ternarySurface = std::move(mesh);
    } else {
      sb.sweep = sol::sweep(sb.solute, chosen, sb.sweepSteps,
                            static_cast<double>(sb.temperatureC),
                            activeBackground(sb),
                            static_cast<double>(sb.backgroundMolarity),
                            requestedPH(sb));
      double best = -1.0;
      for (size_t i = 0; i < sb.sweep.size(); ++i) {
        const double v = sb.sweep[i].prediction.gramsPerMillilitre;
        if (v > best) {
          best = v;
          sb.sweepPeakIndex = static_cast<int>(i);
        }
      }
    }
  } catch (const sol::SolError& err) {
    sb.sweep.clear();
    sb.ternarySurface = {};
    sb.statusMessage = std::string("Ratio sweep failed: ") + err.what();
  }
}

// ---------------------------------------------------------- screening cache

void recomputeScreeningIfNeeded(SolubilityState& sb) {
  if (!sb.soluteValid) {
    sb.screening.clear();
    sb.screeningSignature.clear();
    return;
  }
  std::string signature = std::to_string(sb.soluteVersion);
  signature += '|';
  signature += std::to_string(static_cast<int>(sb.temperatureC * 100.0f));
  signature += '|';
  signature += sb.pHAuto ? "auto" : std::to_string(static_cast<int>(sb.pH * 10.0f));
  if (const sol::Electrolyte* bg = activeBackground(sb)) {
    signature += '|';
    signature += bg->id;
    signature += std::to_string(static_cast<int>(sb.backgroundMolarity * 1000.0f));
  }
  if (signature == sb.screeningSignature) return;
  sb.screeningSignature = signature;
  try {
    sb.screening = sol::screen(sb.solute, static_cast<double>(sb.temperatureC),
                               activeBackground(sb),
                               static_cast<double>(sb.backgroundMolarity), requestedPH(sb));
  } catch (const std::exception& err) {
    sb.screening.clear();
    sb.statusMessage = std::string("Pure-solvent screen failed: ") + err.what();
  }
}

// ------------------------------------------------------------------ controls

void drawSoluteControls(AppState& st) {
  SolubilityState& sb = st.solubility;
  static bool previousUseSketch = sb.useSketch;
  static bool previousOverride = sb.overrideSolute;
  static const char* kSources[] = {"Sketch", "SMILES"};
  int source = sb.useSketch ? 0 : 1;
  if (widgets::segmented("##solute_source", kSources, 2, source)) {
    sb.useSketch = source == 0;
  }

  const bool sourceChanged = sb.useSketch != previousUseSketch;
  previousUseSketch = sb.useSketch;
  if (sb.useSketch) {
    if (sourceChanged || st.docRevision != sb.sourceRevision) recomputeSoluteFromSketch(st);
  } else {
    if (sourceChanged) {
      sb.soluteValid = false;
      sb.soluteError.clear();
    }
    const bool enter = widgets::stringInputWithHint(
        "##solute_smiles", "SMILES, e.g. CC(=O)Oc1ccccc1C(=O)O", sb.manualSmiles,
        ImGuiInputTextFlags_EnterReturnsTrue, true);
    if (enter || ImGui::IsItemDeactivatedAfterEdit()) recomputeSoluteFromSmiles(st);
  }

  if (!sb.soluteError.empty()) {
    widgets::notice(icons::Icon::Warning, sb.soluteError.c_str(), style::col::Danger);
    return;
  }
  if (!sb.soluteValid) {
    widgets::emptyState(icons::Icon::Molecule, "No solute loaded",
                        sb.useSketch ? "Draw a structure in the sketch."
                                     : "Enter a SMILES string and press Enter.");
    return;
  }

  if (!sb.soluteNote.empty()) {
    widgets::notice(icons::Icon::Info, sb.soluteNote.c_str(), style::col::Data);
  }

  char molecularWeight[32];
  char molarVolume[32];
  char logP[32];
  char meltingPoint[32];
  char hansen[64];
  std::snprintf(molecularWeight, sizeof(molecularWeight), "%.1f g/mol", sb.solute.molarMass);
  std::snprintf(molarVolume, sizeof(molarVolume), "%.1f cm3/mol", sb.solute.molarVolume);
  std::snprintf(logP, sizeof(logP), "%.2f", sb.solute.logP);
  std::snprintf(meltingPoint, sizeof(meltingPoint), "%.1f C%s", sb.solute.meltingPoint,
                sb.solute.meltingPointEstimated ? " estimated" : "");
  std::snprintf(hansen, sizeof(hansen), "%.1f / %.1f / %.1f",
                sb.solute.hansen.dispersion, sb.solute.hansen.polar,
                sb.solute.hansen.hydrogenBond);
  widgets::keyValue("Molecular weight", molecularWeight, style::col::Data);
  widgets::keyValue("Molar volume", molarVolume, style::col::Data);
  widgets::keyValue("logP", logP, style::col::Data);
  widgets::keyValue("Melting point", meltingPoint, style::col::Data);
  widgets::keyValue("Hansen dD / dP / dH", hansen, style::col::Data);

  const bool overrideChanged =
      widgets::toggle("##solute_override", "Use measured Tm / Hansen radius",
                      sb.overrideSolute);
  if (overrideChanged && sb.overrideSolute && !previousOverride) {
    sb.meltingPointC = static_cast<float>(sb.solute.meltingPoint);
    sb.interactionRadius = static_cast<float>(sb.solute.interactionRadius);
  }
  previousOverride = sb.overrideSolute;
  if (sb.overrideSolute) {
    bool changed = false;
    changed |= widgets::glyphSlider("##melting_point", icons::Icon::Thermometer,
                                    "Melting point", sb.meltingPointC, -50.0f, 300.0f,
                                    "%.1f C");
    changed |= widgets::glyphSlider("##hansen_radius", icons::Icon::Ruler,
                                    "Hansen radius", sb.interactionRadius, 1.0f, 30.0f,
                                    "%.1f MPa^0.5");
    applySoluteOverrides(sb);
    if (changed) ++sb.soluteVersion;
  }
}

void drawSolventSlot(SolubilityState& sb, int index, std::string& query) {
  static const char* kLabels[] = {"Solvent A", "Solvent B", "Solvent C"};
  std::string& id = sb.solventIds[static_cast<size_t>(index)];
  ImGui::PushID(index);
  ImGui::PushStyleColor(ImGuiCol_Text, slotColor(index));
  ImGui::TextUnformatted(kLabels[index]);
  ImGui::PopStyleColor();
  ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x);
  widgets::solventCombo("##solvent", id, query);
  widgets::glyphSlider("##ratio", icons::Icon::Droplet, "Volume parts",
                       sb.ratios[static_cast<size_t>(index)], 0.0f, 10.0f,
                       "%.2f parts",
                       "Relative volume parts of this solvent in the blend");
  ImGui::PopID();
}

void drawSolventControls(SolubilityState& sb) {
  static std::array<std::string, 3> queries;
  static const char* kCounts[] = {"Pure", "Binary", "Ternary"};
  int count = std::clamp(sb.solventCount - 1, 0, 2);
  if (widgets::segmented("##solvent_count", kCounts, 3, count)) {
    sb.solventCount = count + 1;
  }
  for (int i = 0; i < sb.solventCount; ++i) {
    drawSolventSlot(sb, i, queries[static_cast<size_t>(i)]);
  }
}

void drawConditionsControls(SolubilityState& sb) {
  widgets::glyphSlider("##temperature", icons::Icon::Thermometer, "Temperature",
                       sb.temperatureC, -20.0f, 150.0f, "%.1f C");

  const bool ionisable =
      sb.soluteValid &&
      (sb.solute.ionization.ionClass == sol::IonClass::Base ||
       sb.solute.ionization.ionClass == sol::IonClass::Acid);
  if (widgets::onlyWhen(ionisable, "pH applies only to an ionisable solute.")) {
    widgets::toggle("##ph_auto", "Self-buffered pH", sb.pHAuto,
                    "Use the pH set by the saturated solution");
    if (sb.pHAuto) {
      char value[64];
      std::snprintf(value, sizeof(value), "%.1f (pKa %.1f)",
                    sb.prediction.ionicPath ? sb.prediction.pH : 7.0,
                    sb.solute.ionization.pKa);
      widgets::keyValue("Model pH", value, style::col::Data);
    } else {
      widgets::glyphSlider("##ph", icons::Icon::Ph, "pH", sb.pH, 0.0f, 14.0f,
                           "%.1f");
    }
  }

  const sol::Salt* salt =
      sb.soluteValid && !sb.solute.canonicalSmiles.empty()
          ? sol::findSalt(sb.solute.canonicalSmiles)
          : nullptr;
  if (widgets::onlyWhen(salt != nullptr,
                        "Common-ion controls require a recognised 1:1 salt.")) {
    widgets::sectionHeader("Common-ion effect", style::col::Data);
    char ksp[64];
    std::snprintf(ksp, sizeof(ksp), "Ksp %.3g | %s / %s", salt->ksp25,
                  salt->cation.c_str(), salt->anion.c_str());
    widgets::keyValue(salt->name.c_str(), ksp, style::col::Data);
    widgets::toggle("##background_enabled", "Background electrolyte",
                    sb.backgroundEnabled);
    if (sb.backgroundEnabled) {
      const std::vector<sol::Electrolyte>& all = sol::electrolytes();
      const sol::Electrolyte* current = activeBackground(sb);
      ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x);
      if (ImGui::BeginCombo("##background_electrolyte",
                            current ? current->name.c_str() : "Select electrolyte")) {
        for (size_t i = 0; i < all.size(); ++i) {
          const bool selected = static_cast<int>(i) == sb.backgroundElectrolyte;
          if (ImGui::Selectable(all[i].name.c_str(), selected)) {
            sb.backgroundElectrolyte = static_cast<int>(i);
          }
          if (selected) ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
      }
      widgets::glyphSlider("##background_molarity", icons::Icon::Droplet,
                           "Concentration", sb.backgroundMolarity, 0.0f, 3.0f,
                           "%.2f mol/L");
      current = activeBackground(sb);
      if (current &&
          (current->cation == salt->cation || current->anion == salt->anion)) {
        widgets::notice(icons::Icon::Info,
                        "Common ion present; the model depresses solubility.",
                        style::col::Data);
      } else if (current) {
        ImGui::TextDisabled("No common ion; ionic-strength correction only.");
      }
    }
  }
}

void resolveBlend(const SolubilityState& sb,
                  std::vector<const sol::Solvent*>& chosen,
                  std::vector<sol::Component>& components,
                  std::vector<std::string>& miscibilityWarnings) {
  chosen.assign(static_cast<size_t>(sb.solventCount), nullptr);
  components.clear();
  miscibilityWarnings.clear();
  for (int i = 0; i < sb.solventCount; ++i) {
    const sol::Solvent* solvent =
        sol::findSolvent(sb.solventIds[static_cast<size_t>(i)]);
    chosen[static_cast<size_t>(i)] = solvent;
    if (solvent && sb.ratios[static_cast<size_t>(i)] > 0.0f) {
      components.push_back(
          {solvent, static_cast<double>(sb.ratios[static_cast<size_t>(i)])});
    }
  }
  for (size_t i = 0; i < chosen.size(); ++i) {
    for (size_t j = i + 1; j < chosen.size(); ++j) {
      if (chosen[i] && chosen[j] && !sol::miscibleWith(*chosen[i], *chosen[j])) {
        miscibilityWarnings.push_back(
            chosen[i]->name + " and " + chosen[j]->name +
            " form two phases; this prediction assumes one homogeneous blend.");
      }
    }
  }
}

void drawPredictionHeadline(const SolubilityState& sb, bool canPredict) {
  widgets::cardHeader(icons::Icon::Gauge, "Prediction",
                      canPredict ? "Current blend" : "Waiting for inputs",
                      style::col::Data);
  if (!canPredict) {
    widgets::emptyState(icons::Icon::Gauge, "No prediction",
                        "Load a solute and select a solvent.");
    return;
  }
  const int units = std::clamp(sb.units, 0, 3);
  const double value =
      toDisplayUnits(sb.prediction.gramsPerMillilitre, sb.solute.molarMass, units);
  const std::string formatted = formatSolubility(value);
  widgets::metric("Predicted solubility", formatted.c_str(),
                  kUnitLabels[static_cast<size_t>(units)], nullptr,
                  style::col::DataBright);
  widgets::statusDot(sb.prediction.converged ? "Converged" : "Not converged",
                     sb.prediction.converged,
                     sb.prediction.converged ? style::col::Success : style::col::Danger);
  if (sb.prediction.anchored) {
    widgets::badge("MEASURED", style::col::Success);
  } else {
    widgets::badge("MODEL", style::col::Violet);
  }
}

void drawBinarySweepPlot(SolubilityState& sb, const sol::Solvent& a,
                         const sol::Solvent& b, ImVec2 size) {
  if (sb.sweep.empty()) {
    widgets::emptyState(icons::Icon::ChartLine, "No blend response",
                        "Complete both solvent slots.");
    return;
  }
  std::vector<const sol::SweepPoint*> sorted;
  sorted.reserve(sb.sweep.size());
  for (const sol::SweepPoint& point : sb.sweep) sorted.push_back(&point);
  std::sort(sorted.begin(), sorted.end(),
            [](const sol::SweepPoint* left, const sol::SweepPoint* right) {
              return left->fractions[0] < right->fractions[0];
            });

  std::vector<double> x;
  std::vector<double> y;
  x.reserve(sorted.size());
  y.reserve(sorted.size());
  for (const sol::SweepPoint* point : sorted) {
    x.push_back(point->fractions[0] * 100.0);
    y.push_back(toDisplayUnits(point->prediction.gramsPerMillilitre,
                               sb.solute.molarMass, sb.units));
  }
  const double total = static_cast<double>(sb.ratios[0] + sb.ratios[1]);
  const double currentA = total > 0.0 ? static_cast<double>(sb.ratios[0]) / total : 0.5;
  const std::string seriesLabel = a.name + " / " + b.name;
  charts::Series series{seriesLabel.c_str(), x.data(), y.data(),
                        static_cast<int>(x.size()), style::col::Data, true, false, false};
  charts::PlotStyle plot;
  plot.xLabel = "Solvent A volume %";
  plot.yLabel = kUnitLabels[static_cast<size_t>(std::clamp(sb.units, 0, 3))];
  plot.logY = sb.logScale;
  plot.legend = true;
  plot.grid = true;
  plot.hasCursor = true;
  plot.cursorX = currentA * 100.0;
  const double hoveredX =
      charts::linePlot("##binary_response", &series, 1, size, plot);
  const bool changing = ImGui::IsItemClicked(ImGuiMouseButton_Left) ||
                        (ImGui::IsItemActive() &&
                         ImGui::IsMouseDragging(ImGuiMouseButton_Left));
  if (changing && std::isfinite(hoveredX)) {
    const float fraction = static_cast<float>(std::clamp(hoveredX / 100.0, 0.0, 1.0));
    sb.ratios[0] = fraction;
    sb.ratios[1] = 1.0f - fraction;
  }
}

void drawTernarySweepPlot(SolubilityState& sb, const sol::Solvent& a,
                          const sol::Solvent& b, const sol::Solvent& c,
                          ImVec2 size) {
  const TernarySurfaceMesh& mesh = sb.ternarySurface;
  if (mesh.nodes.empty()) {
    widgets::emptyState(icons::Icon::Grid, "No blend response",
                        "Complete all three solvent slots.");
    return;
  }

  static int view = 0;
  static const icons::Icon kViews[] = {icons::Icon::Grid, icons::Icon::Cube};
  static const char* kTooltips[] = {"Composition triangle", "3D landscape"};
  widgets::segmentedIcons("##ternary_view", kViews, kTooltips, 2, view,
                          ImGui::GetFontSize() * 6.0f);

  std::vector<double> points;
  std::vector<double> values;
  points.reserve(mesh.nodes.size() * 3);
  values.reserve(mesh.nodes.size());
  for (const TernarySurfaceNode& node : mesh.nodes) {
    points.insert(points.end(), node.point.fractions.begin(), node.point.fractions.end());
    values.push_back(toDisplayUnits(node.point.prediction.gramsPerMillilitre,
                                    sb.solute.molarMass, sb.units));
  }
  size.y = std::max(0.0f, size.y - ImGui::GetFrameHeight() -
                              ImGui::GetStyle().ItemSpacing.y);
  if (view == 0) {
    const double total = static_cast<double>(sb.ratios[0] + sb.ratios[1] + sb.ratios[2]);
    charts3d::TernaryStyle ternaryStyle;
    ternaryStyle.aLabel = a.name.c_str();
    ternaryStyle.bLabel = b.name.c_str();
    ternaryStyle.cLabel = c.name.c_str();
    ternaryStyle.low = style::col::DataDim;
    ternaryStyle.high = style::col::DataBright;
    ternaryStyle.isolines = true;
    ternaryStyle.hasMarker = total > 0.0;
    if (total > 0.0) {
      ternaryStyle.markerA = sb.ratios[0] / total;
      ternaryStyle.markerB = sb.ratios[1] / total;
      ternaryStyle.markerC = sb.ratios[2] / total;
    }
    double outA = 0.0;
    double outB = 0.0;
    double outC = 0.0;
    const int picked = charts3d::ternary(
        "##ternary_response", points.data(), values.data(),
        static_cast<int>(values.size()), size, ternaryStyle, &outA, &outB, &outC);
    if (picked == -2) {
      sb.ratios[0] = static_cast<float>(outA);
      sb.ratios[1] = static_cast<float>(outB);
      sb.ratios[2] = static_cast<float>(outC);
    }
  } else {
    static charts3d::Orbit orbit;
    charts3d::SurfaceStyle surfaceStyle;
    surfaceStyle.uLabel = "A:B ratio";
    surfaceStyle.vLabel = c.name.c_str();
    surfaceStyle.wLabel =
        kUnitLabels[static_cast<size_t>(std::clamp(sb.units, 0, 3))];
    surfaceStyle.low = style::col::DataDim;
    surfaceStyle.high = style::col::Data;
    surfaceStyle.peak = style::col::DataBright;
    surfaceStyle.logHeight = sb.logScale;
    charts3d::surface("##ternary_surface", values.data(), mesh.ratioQuads + 1,
                      mesh.cQuads + 1, size, orbit, surfaceStyle);
  }
}

const sol::SweepPoint* activeSweepPeak(const SolubilityState& sb) {
  if (sb.sweepPeakIndex < 0 ||
      static_cast<size_t>(sb.sweepPeakIndex) >= sb.sweep.size()) {
    return nullptr;
  }
  return &sb.sweep[static_cast<size_t>(sb.sweepPeakIndex)];
}

void drawBlendResponseCard(SolubilityState& sb,
                           const std::vector<const sol::Solvent*>& chosen,
                           bool canPredict, bool haveAllSlots, ImVec2 size) {
  if (!widgets::beginCard("##blend_response_card", size, style::col::BgRaised,
                          ImGuiWindowFlags_NoScrollbar |
                              ImGuiWindowFlags_NoScrollWithMouse)) {
    return;
  }
  widgets::cardHeader(icons::Icon::ChartLine, "Blend response",
                      "Drag or click the chart to set composition",
                      style::col::Data);

  if (sb.solventCount >= 2) {
    float samples = static_cast<float>(sb.sweepSteps);
    if (widgets::glyphSlider("##samples", icons::Icon::Grid, "Samples", samples,
                             2.0f, 64.0f, "%.0f")) {
      sb.sweepSteps = std::clamp(static_cast<int>(std::lround(samples)), 2, 64);
    }
  } else {
    ImGui::TextDisabled("Samples: pure solvent");
  }
  ImGui::SetNextItemWidth(layout::columnWidth(layout::measure(), 2));
  ImGui::Combo("Units", &sb.units, kUnitLabels.data(), 4);
  ImGui::SameLine(0.0f, style::metrics().gap);
  widgets::toggle("##scale", "Scale", sb.logScale,
                  "On uses a logarithmic y axis; off uses a linear y axis");

  ImVec2 plotSize = ImGui::GetContentRegionAvail();
  plotSize.x = std::max(0.0f, plotSize.x);
  plotSize.y = std::max(0.0f, plotSize.y);
  if (!sb.soluteValid) {
    widgets::emptyState(icons::Icon::Molecule, "Load a solute",
                        "Use the Solute tab to define the compound.");
  } else if (!canPredict) {
    widgets::emptyState(icons::Icon::Droplet, "Select a solvent",
                        "Complete at least one solvent slot.");
  } else if (sb.solventCount < 2) {
    widgets::emptyState(icons::Icon::ChartLine, "Pure-solvent prediction",
                        "Add a second solvent to explore blend composition.");
  } else if (!haveAllSlots) {
    widgets::emptyState(icons::Icon::Droplet, "Incomplete blend",
                        "Choose every active solvent slot.");
  } else if (sb.solventCount == 2) {
    drawBinarySweepPlot(sb, *chosen[0], *chosen[1], plotSize);
  } else {
    drawTernarySweepPlot(sb, *chosen[0], *chosen[1], *chosen[2], plotSize);
  }
  widgets::endCard();
}

void drawModelTermsCard(SolubilityState& sb,
                        const std::vector<sol::Component>& components,
                        bool canPredict, ImVec2 size) {
  if (!widgets::beginCard("##model_terms_card", size, style::col::BgSurface,
                          ImGuiWindowFlags_NoScrollbar |
                              ImGuiWindowFlags_NoScrollWithMouse)) {
    return;
  }
  widgets::cardHeader(icons::Icon::ChartLine, "Model terms",
                      "Kirkwood-Buff and Hildebrand terms", style::col::Data);
  if (!canPredict) {
    ImGui::TextDisabled("Complete the solute and blend to calculate model terms.");
    widgets::endCard();
    return;
  }

  const sol::KBResult kb =
      sol::kirkwoodBuff(sb.solute, components, sb.temperatureC);
  const sol::Mixture mixture = sol::blend(components);
  const auto hildebrand = [](const sol::Hansen& hansen) {
    return std::sqrt(hansen.dispersion * hansen.dispersion +
                     hansen.polar * hansen.polar +
                     hansen.hydrogenBond * hansen.hydrogenBond);
  };
  const double deltaGap =
      std::fabs(hildebrand(sb.solute.hansen) - hildebrand(mixture.hansen));
  const double values[] = {kb.g11, kb.g12, kb.chi, kb.lnGammaInf,
                           kb.lnIdeal, deltaGap};
  const char* labels[] = {"G11", "G12", "chi", "ln gamma", "ln ideal", "|delta|"};
  const char* units[] = {"cm3/mol", "cm3/mol", nullptr, nullptr, nullptr, "MPa^0.5"};
  std::array<charts3d::ParallelAxis, 6> axes;
  for (size_t i = 0; i < axes.size(); ++i) {
    const charts::Axis axis = charts::niceAxis(std::min(0.0, values[i]),
                                                std::max(0.0, values[i]));
    axes[i] = {labels[i], axis.min, axis.max, true, units[i]};
  }
  const charts3d::ParallelSeries series{"Current blend", values,
                                        style::col::Data, true};

  const bool details = widgets::disclosure(
      "##model_details", "Details", kb.kappaKnown ? "measured kappa" : "kappa unavailable",
      false, icons::Icon::Info, style::col::Data);
  if (details) {
    char value[64];
    std::snprintf(value, sizeof(value), "%+.2f cm3/mol", kb.g11);
    widgets::keyValue("Solvent-solvent G11", value, style::col::Data);
    std::snprintf(value, sizeof(value), "%+.2f cm3/mol", kb.g12);
    widgets::keyValue("Solute-solvent G12", value, style::col::Data);
    std::snprintf(value, sizeof(value), "%.3f", kb.chi);
    widgets::keyValue("Flory-Huggins chi", value, style::col::Data);
    std::snprintf(value, sizeof(value), "%+.3f", kb.lnGammaInf);
    widgets::keyValue("ln gamma infinity", value, style::col::Data);
    std::snprintf(value, sizeof(value), "%+.3f", kb.lnIdeal);
    widgets::keyValue("Ideal ln mole fraction", value, style::col::Data);
    std::snprintf(value, sizeof(value), "%.2f MPa^0.5", deltaGap);
    widgets::keyValue("Hildebrand difference", value, style::col::Data);
  }
  const ImVec2 chartSize(ImGui::GetContentRegionAvail().x,
                         std::max(0.0f, ImGui::GetContentRegionAvail().y));
  if (chartSize.y > 0.0f) {
    charts3d::parallelCoordinates("##model_parallel", axes.data(),
                                  static_cast<int>(axes.size()), &series, 1,
                                  chartSize);
  }
  widgets::endCard();
}

void drawPeakEvidence(AppState& st,
                      const std::vector<const sol::Solvent*>& chosen,
                      const std::vector<std::string>& miscibilityWarnings,
                      bool solventsOk, ImVec2 size) {
  SolubilityState& sb = st.solubility;
  if (!widgets::beginCard("##peak_evidence_card", size, style::col::BgSurface,
                          ImGuiWindowFlags_NoScrollbar |
                              ImGuiWindowFlags_NoScrollWithMouse)) {
    return;
  }
  widgets::cardHeader(icons::Icon::Crosshair, "Evidence",
                      "Peak, provenance and diagnostics", style::col::Data);
  const sol::SweepPoint* peak = activeSweepPeak(sb);
  if (peak && sb.solventCount >= 2) {
    const std::string peakValue =
        formatUnits(peak->prediction.gramsPerMillilitre, sb.solute.molarMass,
                    sb.units);
    widgets::keyValue("Best sampled blend", peakValue.c_str(),
                      style::col::DataBright);
    if (widgets::actionButton("##apply_peak", icons::Icon::Crosshair,
                              "Apply peak", ImVec2(ImGui::GetContentRegionAvail().x,
                                                  ImGui::GetFrameHeight()),
                              false)) {
      for (int i = 0; i < sb.solventCount; ++i) {
        sb.ratios[static_cast<size_t>(i)] =
            static_cast<float>(peak->fractions[static_cast<size_t>(i)]);
      }
    }
  } else {
    ImGui::TextDisabled(sb.solventCount < 2
                            ? "Add another solvent to search for a blend peak."
                            : "Complete the blend to calculate a peak.");
  }

  if (sb.prediction.anchored) {
    widgets::notice(icons::Icon::Book, sb.prediction.anchorNote.c_str(),
                    style::col::Success);
  } else {
    widgets::badge("MODEL", style::col::Violet);
  }
  if (sb.prediction.outsideSphere) {
    widgets::notice(icons::Icon::Warning,
                    "Outside the Hansen sphere; the prediction is extrapolated.",
                    style::col::Danger);
  } else if (!sb.prediction.converged) {
    widgets::notice(icons::Icon::Warning,
                    "Saturation composition did not converge.",
                    style::col::Danger);
  } else if (!miscibilityWarnings.empty()) {
    for (const std::string& warning : miscibilityWarnings) {
      widgets::notice(icons::Icon::Warning, warning.c_str(), style::col::Danger);
    }
  } else {
    widgets::statusDot("No model or miscibility warnings", true,
                       style::col::Success);
  }

  const bool canSend = sb.soluteValid && solventsOk && chosen.size() >= 2 &&
                       chosen[0] && chosen[1];
  if (!canSend) ImGui::BeginDisabled();
  if (widgets::actionButton("##send_extraction", icons::Icon::SepFunnel,
                            "Send to Extraction",
                            ImVec2(ImGui::GetContentRegionAvail().x,
                                   ImGui::GetFrameHeight()),
                            true)) {
    const double ratioA = std::max(0.0, static_cast<double>(sb.ratios[0]));
    const double ratioB = std::max(0.0, static_cast<double>(sb.ratios[1]));
    const double total = std::max(ratioA + ratioB, 1e-9);
    ExtractionImport& import = sb.extractionImport;
    import.pending = true;
    import.solventIdA = chosen[0]->id;
    import.solventIdB = chosen[1]->id;
    import.volumeMlA = 100.0 * ratioA / total;
    import.volumeMlB = 100.0 * ratioB / total;
    import.soluteMassMg = 100.0;
    st.tab = MainTab::Extraction;
    st.tabChangeRequested = true;
    sb.statusMessage = "Blend sent to the Extraction Calculator";
  }
  if (!canSend) ImGui::EndDisabled();
  widgets::endCard();
}

void drawPureSolvents(SolubilityState& sb, layout::Density density, ImVec2 size) {
  if (!widgets::beginCard("##pure_solvents_card", size, style::col::BgSurface,
                          ImGuiWindowFlags_NoScrollbar |
                              ImGuiWindowFlags_NoScrollWithMouse)) {
    return;
  }
  widgets::cardHeader(icons::Icon::Table, "Pure solvents",
                      "Click a bar to load Solvent A", style::col::Data);
  recomputeScreeningIfNeeded(sb);
  if (!sb.soluteValid || sb.screening.empty()) {
    ImGui::TextDisabled("Load a solute to rank pure solvents.");
    widgets::endCard();
    return;
  }

  const layout::Frame frame = layout::measure();
  const int rowsThatFit =
      std::max(1, static_cast<int>(frame.size.y / std::max(frame.row, 1.0f)) - 1);
  const int rowCount =
      std::min(rowsThatFit, static_cast<int>(sb.screening.size()));
  std::vector<std::string> labels;
  std::vector<std::string> annotations;
  std::vector<charts::BarRow> bars;
  labels.reserve(static_cast<size_t>(rowCount));
  annotations.reserve(static_cast<size_t>(rowCount));
  bars.reserve(static_cast<size_t>(rowCount));
  for (int i = 0; i < rowCount; ++i) {
    const sol::ScreenRow& row = sb.screening[static_cast<size_t>(i)];
    labels.push_back(row.solvent->name);
    annotations.push_back(
        formatUnits(row.prediction.gramsPerMillilitre, sb.solute.molarMass,
                    sb.units));
    bars.push_back({labels.back().c_str(),
                    toDisplayUnits(row.prediction.gramsPerMillilitre,
                                   sb.solute.molarMass, sb.units),
                    annotations.back().c_str(), style::col::Data,
                    row.solvent->id == sb.solventIds[0]});
  }

  int clicked = -1;
  if (density == layout::Density::Roomy &&
      layout::columnsThatFit(frame, 18.0f) >= 2) {
    const float tableWidth = layout::columnWidth(frame, 3, 2);
    const float barsWidth = layout::columnWidth(frame, 3);
    const widgets::Column columns[] = {
        {"Solvent", false, true, nullptr, 10.0f},
        {"Solubility", true, false,
         kUnitLabels[static_cast<size_t>(std::clamp(sb.units, 0, 3))], 7.0f}};
    if (widgets::beginDataTable("##pure_table", columns, 2,
                                ImVec2(tableWidth, frame.size.y))) {
      for (int i = 0; i < rowCount; ++i) {
        const sol::ScreenRow& row = sb.screening[static_cast<size_t>(i)];
        widgets::dataRow(row.solvent->id == sb.solventIds[0]
                             ? style::col::DataBright
                             : ImVec4{});
        widgets::dataCell(row.solvent->name.c_str());
        widgets::dataCell(
            formatSolubility(toDisplayUnits(row.prediction.gramsPerMillilitre,
                                             sb.solute.molarMass, sb.units))
                .c_str());
      }
      widgets::endDataTable();
    }
    ImGui::SameLine(0.0f, frame.gap);
    clicked = charts::rankedBars("##pure_bars", bars.data(), rowCount,
                                 ImVec2(barsWidth, frame.size.y));
  } else {
    clicked = charts::rankedBars("##pure_bars", bars.data(), rowCount,
                                 frame.size);
  }
  if (clicked >= 0) {
    const sol::ScreenRow& row = sb.screening[static_cast<size_t>(clicked)];
    sb.solventIds[0] = row.solvent->id;
    sb.statusMessage = "Solvent A := " + row.solvent->name;
  }
  widgets::endCard();
}

void drawEvidenceRow(AppState& st,
                     const std::vector<const sol::Solvent*>& chosen,
                     const std::vector<std::string>& miscibilityWarnings,
                     bool solventsOk, layout::Density density, ImVec2 size) {
  const layout::Frame frame = layout::measure(size);
  if (density == layout::Density::Roomy &&
      layout::columnsThatFit(frame, 20.0f) >= 2) {
    const float evidenceWidth = layout::columnWidth(frame, 2);
    drawPeakEvidence(st, chosen, miscibilityWarnings, solventsOk,
                     ImVec2(evidenceWidth, size.y));
    ImGui::SameLine(0.0f, frame.gap);
    drawPureSolvents(st.solubility, density,
                     ImVec2(evidenceWidth, size.y));
    return;
  }

  static int evidenceTab = 0;
  static const char* kLabels[] = {"Evidence", "Pure solvents"};
  static const icons::Icon kIcons[] = {icons::Icon::Crosshair, icons::Icon::Table};
  widgets::subTabs("##evidence_tabs", kLabels, kIcons, 2, evidenceTab);
  const float cardHeight =
      std::max(0.0f, size.y - ImGui::GetFrameHeight() - frame.gap);
  if (evidenceTab == 0) {
    drawPeakEvidence(st, chosen, miscibilityWarnings, solventsOk,
                     ImVec2(size.x, cardHeight));
  } else {
    drawPureSolvents(st.solubility, density, ImVec2(size.x, cardHeight));
  }
}


// One row of chrome that names the two halves of the workspace. It is drawn
// before the mode dispatch so the switch is present in both modes -- a control
// that vanishes in the view it would take you out of is a trap.
void drawSuiteModeBar(SolubilityState& sb) {
  static const char* kModeLabels[2] = {"Predict solubility", "Select solvent"};
  int mode = static_cast<int>(sb.suiteMode);

  widgets::beginToolbar("##suite_mode_bar");
  icons::draw(ImGui::GetWindowDrawList(), icons::Icon::Flask,
              ImVec2(ImGui::GetCursorScreenPos().x + ImGui::GetFontSize() * 0.5f,
                     ImGui::GetCursorScreenPos().y + ImGui::GetFrameHeight() * 0.5f),
              ImGui::GetFontSize(), style::u32(style::col::Accent));
  ImGui::Dummy(ImVec2(ImGui::GetFontSize() * 1.2f, ImGui::GetFrameHeight()));
  ImGui::SameLine();
  widgets::segmented("##suite_mode", kModeLabels, 2, mode, ImGui::GetFontSize() * 22.0f);
  sb.suiteMode = mode == 1 ? SuiteMode::Select : SuiteMode::Predict;
  ImGui::SameLine();
  widgets::helpMarker(
      "Predict answers how much of your solute dissolves in a blend you choose. "
      "Select ranks the whole solvent database for an operation -- extraction, "
      "recrystallisation, a wash -- and loads the winner back into Predict.");
  widgets::endToolbar();
}

}  // namespace

void drawSolubilitySuite(AppState& st) {
  SolubilityState& sb = st.solubility;
  drawSuiteModeBar(sb);
  if (sb.suiteMode == SuiteMode::Select) {
    drawSolventSelector(st);
    return;
  }

  const layout::Frame frame = layout::measure();
  const int columns =
      std::clamp(layout::columnsThatFit(frame, 20.0f), 2, 4);
  const float railWidth = layout::columnWidth(frame, columns);
  const float workspaceWidth =
      layout::columnWidth(frame, columns, columns - 1);

  std::vector<sol::Component> components;
  std::vector<const sol::Solvent*> chosen;
  std::vector<std::string> miscibilityWarnings;
  bool solventsOk = true;

  ImGui::PushStyleVar(ImGuiStyleVar_ChildRounding, style::metrics().radiusLg);
  ImGui::PushStyleColor(ImGuiCol_ChildBg, style::col::BgPanel);
  ImGui::BeginChild("##suite_rail", ImVec2(railWidth, frame.size.y),
                    ImGuiChildFlags_Borders,
                    ImGuiWindowFlags_NoScrollbar |
                        ImGuiWindowFlags_NoScrollWithMouse);
  ImGui::PopStyleColor();
  ImGui::PopStyleVar();

  const layout::Frame railFrame = layout::measure();
  const float railWeights[] = {0.68f, 0.32f};
  const float railMinimums[] = {railFrame.row * 8.0f, railFrame.row * 4.0f};
  float railHeights[2]{};
  layout::distribute(railFrame.size.y, railWeights, railMinimums, 2,
                     railFrame.gap, railHeights);
  const float railStartY = ImGui::GetCursorPosY();

  if (widgets::beginCard("##suite_inputs", ImVec2(railFrame.size.x, railHeights[0]),
                         style::col::BgSurface,
                         ImGuiWindowFlags_NoScrollbar |
                             ImGuiWindowFlags_NoScrollWithMouse)) {
    static int inputTab = 0;
    static const char* kInputLabels[] = {"Solute", "Solvents", "Conditions"};
    static const icons::Icon kInputIcons[] = {
        icons::Icon::Flask, icons::Icon::Droplet, icons::Icon::Thermometer};
    widgets::subTabs("##input_tabs", kInputLabels, kInputIcons, 3, inputTab);
    if (inputTab == 0) {
      drawSoluteControls(st);
    } else if (inputTab == 1) {
      drawSolventControls(sb);
    } else {
      drawConditionsControls(sb);
    }
    widgets::endCard();
  }

  try {
    resolveBlend(sb, chosen, components, miscibilityWarnings);
  } catch (const sol::SolError& error) {
    solventsOk = false;
    sb.statusMessage = std::string("Solvent database error: ") + error.what();
  }

  const bool canPredict = sb.soluteValid && solventsOk && !components.empty();
  if (canPredict) {
    try {
      sb.prediction =
          sol::predict(sb.solute, components,
                       static_cast<double>(sb.temperatureC),
                       activeBackground(sb),
                       static_cast<double>(sb.backgroundMolarity),
                       requestedPH(sb));
      sb.statusMessage =
          "Predicted " +
          formatUnits(sb.prediction.gramsPerMillilitre,
                      sb.solute.molarMass, sb.units) +
          ".";
    } catch (const sol::SolError& error) {
      sb.prediction = {};
      solventsOk = false;
      sb.statusMessage = std::string("Prediction failed: ") + error.what();
    }
  } else {
    sb.prediction = {};
  }

  const bool haveAllSlots =
      solventsOk && !chosen.empty() &&
      std::all_of(chosen.begin(), chosen.end(),
                  [](const sol::Solvent* solvent) { return solvent != nullptr; });
  if (sb.solventCount >= 2 && sb.soluteValid && haveAllSlots) {
    recomputeSweepIfNeeded(sb, chosen);
  } else {
    sb.sweep.clear();
    sb.ternarySurface = {};
    sb.sweepSignature.clear();
  }

  layout::nextRow(railStartY + railHeights[0] + railFrame.gap);
  if (widgets::beginCard("##suite_prediction",
                         ImVec2(railFrame.size.x, railHeights[1]),
                         style::col::BgRaised,
                         ImGuiWindowFlags_NoScrollbar |
                             ImGuiWindowFlags_NoScrollWithMouse)) {
    drawPredictionHeadline(sb, canPredict && solventsOk);
    widgets::endCard();
  }
  ImGui::EndChild();

  ImGui::SameLine(0.0f, frame.gap);
  ImGui::BeginChild("##suite_workspace",
                    ImVec2(workspaceWidth, frame.size.y),
                    ImGuiChildFlags_None,
                    ImGuiWindowFlags_NoScrollbar |
                        ImGuiWindowFlags_NoScrollWithMouse);
  const layout::Frame workspaceFrame = layout::measure();
  const float workspaceWeights[] = {0.52f, 0.23f, 0.25f};
  const float workspaceMinimums[] = {
      workspaceFrame.row * 10.0f,
      workspaceFrame.row * 5.0f,
      workspaceFrame.row * 6.0f};
  float workspaceHeights[3]{};
  layout::distribute(workspaceFrame.size.y, workspaceWeights,
                     workspaceMinimums, 3, workspaceFrame.gap,
                     workspaceHeights);
  const float workspaceStartY = ImGui::GetCursorPosY();

  drawBlendResponseCard(sb, chosen, canPredict && solventsOk, haveAllSlots,
                        ImVec2(workspaceFrame.size.x, workspaceHeights[0]));
  layout::nextRow(workspaceStartY + workspaceHeights[0] +
                       workspaceFrame.gap);
  drawModelTermsCard(sb, components, canPredict && solventsOk,
                     ImVec2(workspaceFrame.size.x, workspaceHeights[1]));
  layout::nextRow(workspaceStartY + workspaceHeights[0] +
                       workspaceFrame.gap + workspaceHeights[1] +
                       workspaceFrame.gap);
  drawEvidenceRow(st, chosen, miscibilityWarnings, solventsOk,
                  frame.density,
                  ImVec2(workspaceFrame.size.x, workspaceHeights[2]));
  ImGui::EndChild();
}

}  // namespace chemcad::ui

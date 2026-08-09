// Solubility Suite: solute + solvent inputs on the left rail, the prediction
// and the solubility-vs-composition graph dominant on the right. Pure
// solvents are screened in a ranked table; a chosen blend can be handed to
// the Extraction Calculator. Renders inside an already-open window; never owns one.

#include "imgui.h"

#include <algorithm>
#include <array>
#include <cfloat>
#include <cstdarg>
#include <cmath>
#include <cstdint>
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
#include "ui/solubility_state.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"
#include "ui/widgets.hpp"

namespace chemcad::ui {
namespace {

// Sized from the current font rather than fixed pixels -- the app runs at a
// 1.25 UI scale by default, and hardcoded widths clip their own labels.
inline float uiScale() { return ImGui::GetFontSize() / 13.0f; }

// Prepares a widget whose visible ImGui label is rendered to its right.
// Returns false when the row is too narrow; in that case the caption is
// wrapped above and the caller must use a hidden "##" widget label.
bool prepareTrailingLabel(const char* label) {
  const float available = ImGui::GetContentRegionAvail().x;
  const float reserved =
      ImGui::CalcTextSize(label).x + ImGui::GetStyle().ItemInnerSpacing.x;
  const float minimumWidget = ImGui::GetFontSize() * 3.0f;
  if (available - reserved >= minimumWidget) {
    ImGui::SetNextItemWidth(available - reserved);
    return true;
  }
  ImGui::TextWrapped("%s", label);
  ImGui::SetNextItemWidth(-FLT_MIN);
  return false;
}

void textDisabledWrapped(const char* format, ...) {
  va_list args;
  va_start(args, format);
  ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
  ImGui::TextWrappedV(format, args);
  ImGui::PopStyleColor();
  va_end(args);
}

bool checkboxWithResponsiveLabel(const char* label, const char* stackedLabel, bool* value) {
  const float required = ImGui::GetFrameHeight() + ImGui::GetStyle().ItemInnerSpacing.x +
                         ImGui::CalcTextSize(label).x;
  if (required <= ImGui::GetContentRegionAvail().x) return ImGui::Checkbox(label, value);
  ImGui::TextWrapped("%s", label);
  return ImGui::Checkbox(stackedLabel, value);
}

std::string ellipsizeText(const std::string& text, float maxWidth) {
  if (maxWidth <= 0.0f) return {};
  if (ImGui::CalcTextSize(text.c_str()).x <= maxWidth) return text;
  constexpr const char* kEllipsis = "...";
  const float ellipsisWidth = ImGui::CalcTextSize(kEllipsis).x;
  if (ellipsisWidth >= maxWidth) return {};
  size_t end = text.size();
  while (end > 0 &&
         ImGui::CalcTextSize(text.data(), text.data() + end).x + ellipsisWidth > maxWidth) {
    --end;
  }
  return text.substr(0, end) + kEllipsis;
}


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

// Perceptually ordered surface ramp. Keeping the stops in one palette-token
// table makes the value encoding easy to tune without touching mesh logic.
ImVec4 surfaceColour(double normalisedValue) {
  static constexpr std::array<ImVec4, 4> kStops = {
      style::col::Violet, style::col::Teal, style::col::Accent,
      style::col::AccentHover};
  const float t = std::clamp(static_cast<float>(normalisedValue), 0.0f, 1.0f);
  const float scaled = t * static_cast<float>(kStops.size() - 1);
  const size_t index = std::min(kStops.size() - 2,
                                static_cast<size_t>(scaled));
  const float local = scaled - static_cast<float>(index);
  const ImVec4& a = kStops[index];
  const ImVec4& b = kStops[index + 1];
  return ImVec4(a.x + (b.x - a.x) * local,
                a.y + (b.y - a.y) * local,
                a.z + (b.z - a.z) * local, 1.0f);
}

std::array<float, 4> cachedSurfaceColour(double normalisedValue) {
  const ImVec4 colour = surfaceColour(normalisedValue);
  return {colour.x, colour.y, colour.z, colour.w};
}

float surfaceLambertShade(const std::array<float, 3>& p00,
                          const std::array<float, 3>& p10,
                          const std::array<float, 3>& p11,
                          const std::array<float, 3>& p01) {
  const auto subtract = [](const std::array<float, 3>& left,
                           const std::array<float, 3>& right) {
    return std::array<float, 3>{left[0] - right[0],
                                left[1] - right[1],
                                left[2] - right[2]};
  };
  const auto cross = [](const std::array<float, 3>& left,
                        const std::array<float, 3>& right) {
    return std::array<float, 3>{
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0]};
  };
  const std::array<float, 3> frontNormal =
      cross(subtract(p01, p00), subtract(p10, p00));
  const std::array<float, 3> backNormal =
      cross(subtract(p11, p10), subtract(p11, p01));
  std::array<float, 3> normal{
      frontNormal[0] + backNormal[0],
      frontNormal[1] + backNormal[1],
      frontNormal[2] + backNormal[2]};
  const float length = std::sqrt(normal[0] * normal[0] +
                                 normal[1] * normal[1] +
                                 normal[2] * normal[2]);
  if (length > 1e-7f) {
    normal[0] /= length;
    normal[1] /= length;
    normal[2] /= length;
  }
  constexpr std::array<float, 3> kLight{-0.348f, 0.870f, -0.348f};
  const float lambert = std::clamp(normal[0] * kLight[0] +
                                       normal[1] * kLight[1] +
                                       normal[2] * kLight[2],
                                   0.0f, 1.0f);
  return 0.75f + lambert * 0.40f;
}

ImVec4 shadedSurfaceColour(const std::array<float, 4>& colour, float shade,
                           float alpha) {
  return ImVec4(std::clamp(colour[0] * shade, 0.0f, 1.0f),
                std::clamp(colour[1] * shade, 0.0f, 1.0f),
                std::clamp(colour[2] * shade, 0.0f, 1.0f),
                std::clamp(alpha, 0.0f, 1.0f));
}

bool barycentricPoint(ImVec2 point, ImVec2 a, ImVec2 b, ImVec2 c,
                      std::array<float, 3>& weights) {
  const float denominator =
      (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y);
  if (std::fabs(denominator) < 1e-7f) return false;
  weights[0] = ((b.y - c.y) * (point.x - c.x) +
                (c.x - b.x) * (point.y - c.y)) /
               denominator;
  weights[1] = ((c.y - a.y) * (point.x - c.x) +
                (a.x - c.x) * (point.y - c.y)) /
               denominator;
  weights[2] = 1.0f - weights[0] - weights[1];
  constexpr float kTolerance = -0.002f;
  return weights[0] >= kTolerance && weights[1] >= kTolerance &&
         weights[2] >= kTolerance;
}

// Slot identity colours are shared by mixer rows, plot markers and axes so a
// solvent keeps its colour across the panel.
ImVec4 slotColor(int index) {
  switch (index) {
    case 0: return style::col::Teal;
    case 1: return style::col::Accent;
    default: return style::col::Violet;
  }
}

const sol::Electrolyte* activeBackground(const SolubilityState& sb);

struct TimedSignature {
  std::string value;
  double startedAt = 0.0;
};

float sweepDrawProgress(const SolubilityState& sb) {
  static std::unordered_map<const SolubilityState*, TimedSignature> states;
  TimedSignature& state = states[&sb];
  if (state.value != sb.sweepSignature) {
    state.value = sb.sweepSignature;
    state.startedAt = ImGui::GetTime();
  }
  constexpr double kDrawSeconds = 0.30;
  const float linear = static_cast<float>(
      std::clamp((ImGui::GetTime() - state.startedAt) / kDrawSeconds, 0.0, 1.0));
  const float remaining = 1.0f - linear;
  return 1.0f - remaining * remaining * remaining;
}

std::string predictionSignature(const SolubilityState& sb) {
  std::string signature;
  signature.reserve(128);
  signature += std::to_string(sb.soluteVersion);
  signature += '|';
  signature += std::to_string(static_cast<int>(sb.temperatureC * 100.0f));
  signature += '|';
  for (int i = 0; i < sb.solventCount; ++i) {
    signature += sb.solventIds[static_cast<size_t>(i)];
    signature += ':';
    signature += std::to_string(
        static_cast<int>(sb.ratios[static_cast<size_t>(i)] * 10000.0f));
    signature += '|';
  }
  if (const sol::Electrolyte* background = activeBackground(sb)) {
    signature += background->id;
    signature += ':';
    signature += std::to_string(static_cast<int>(sb.backgroundMolarity * 1000.0f));
  }
  return signature;
}

float predictionFlashAlpha(const SolubilityState& sb, bool canPredict) {
  static std::unordered_map<const SolubilityState*, TimedSignature> states;
  TimedSignature& state = states[&sb];
  if (!canPredict) {
    state.value.clear();
    state.startedAt = 0.0;
    return 0.0f;
  }
  const std::string signature = predictionSignature(sb);
  if (state.value != signature) {
    state.value = signature;
    state.startedAt = ImGui::GetTime();
  }
  constexpr double kFlashSeconds = 0.40;
  const double elapsed = ImGui::GetTime() - state.startedAt;
  if (elapsed >= kFlashSeconds) return 0.0f;
  const float phase = static_cast<float>(elapsed / kFlashSeconds);
  return std::sin(phase * 3.14159265f);
}

float markerPulse(float phaseOffset = 0.0f) {
  constexpr float kCycleSeconds = 1.20f;
  const float phase =
      static_cast<float>(ImGui::GetTime()) * (2.0f * 3.14159265f / kCycleSeconds) +
      phaseOffset;
  return 0.5f + 0.5f * std::sin(phase);
}

float markerReveal(float drawProgress, float pathProgress) {
  constexpr float kRevealWindow = 0.08f;
  return std::clamp((drawProgress - pathProgress + kRevealWindow) / kRevealWindow,
                    0.0f, 1.0f);
}

void drawInteractiveBorder(const char* animationId, ImVec4 accent) {
  const bool hovered = ImGui::IsItemHovered();
  const bool held = ImGui::IsItemActive();
  ImGui::PushID(animationId);
  const float hover = widgets::hoverT(ImGui::GetID("hover"), hovered || held);
  const float press = widgets::hoverT(ImGui::GetID("press"), held);
  ImGui::PopID();

  const style::Metrics& metrics = style::metrics();
  const ImVec2 min = ImGui::GetItemRectMin();
  const ImVec2 max = ImGui::GetItemRectMax();
  const float inset = press * std::max(1.0f, metrics.hairline * 1.5f);
  ImGui::GetWindowDrawList()->AddRect(
      ImVec2(min.x + inset, min.y + inset), ImVec2(max.x - inset, max.y - inset),
      style::mix(style::col::Border, accent, hover), metrics.radiusSm, 0,
      metrics.hairline * (1.0f + 0.8f * hover + 0.5f * press));
}

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

      mesh.quads.reserve(static_cast<size_t>(mesh.ratioQuads) *
                         static_cast<size_t>(mesh.cQuads));
      for (int cIndex = 0; cIndex < mesh.cQuads; ++cIndex) {
        for (int ratioIndex = 0; ratioIndex < mesh.ratioQuads; ++ratioIndex) {
          const uint32_t n00 = static_cast<uint32_t>(
              cIndex * ratioNodes + ratioIndex);
          const uint32_t n10 = n00 + 1;
          const uint32_t n01 = n00 + static_cast<uint32_t>(ratioNodes);
          const uint32_t n11 = n01 + 1;
          TernarySurfaceQuad quad;
          quad.nodes = {n00, n10, n11, n01};
          quad.ratioIndex = ratioIndex;
          quad.cIndex = cIndex;
          const TernarySurfaceNode& p00 = mesh.nodes[n00];
          const TernarySurfaceNode& p10 = mesh.nodes[n10];
          const TernarySurfaceNode& p11 = mesh.nodes[n11];
          const TernarySurfaceNode& p01 = mesh.nodes[n01];
          const double linearT =
              (p00.cubeLinear[1] + p10.cubeLinear[1] +
               p11.cubeLinear[1] + p01.cubeLinear[1]) *
              0.25;
          const double logT =
              (p00.cubeLog[1] + p10.cubeLog[1] +
               p11.cubeLog[1] + p01.cubeLog[1]) *
              0.25;
          quad.colourLinear = cachedSurfaceColour(linearT);
          quad.colourLog = cachedSurfaceColour(logT);
          quad.shadeLinear =
              surfaceLambertShade(p00.cubeLinear, p10.cubeLinear,
                                  p11.cubeLinear, p01.cubeLinear);
          quad.shadeLog =
              surfaceLambertShade(p00.cubeLog, p10.cubeLog,
                                  p11.cubeLog, p01.cubeLog);
          quad.depth =
              (p00.cubeLinear[2] + p10.cubeLinear[2] +
               p11.cubeLinear[2] + p01.cubeLinear[2]) *
              0.25f;
          mesh.quads.push_back(std::move(quad));
        }
      }
      std::stable_sort(
          mesh.quads.begin(), mesh.quads.end(),
          [](const TernarySurfaceQuad& left,
             const TernarySurfaceQuad& right) {
            return left.depth > right.depth;
          });

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
    sb.statusMessage = std::string("Solvent screen failed: ") + err.what();
  }
}

// ------------------------------------------------------------------ widgets

// Minimum height for widgets::statCard's stacked "value above caption"
// layout (mono value line, then the dim label line beneath it), mirroring
// the padding statCard itself uses internally so a card fed this height can
// never let the label baseline land under the value's at any
// CHEMCAD_UI_SCALE -- every term below scales in lockstep with the font, so
// the inequality that keeps them apart is scale-invariant.
float statCardHeight() {
  const style::Metrics& m = style::metrics();
  const float pad = m.gap * 0.9f;
  const float valueLineH = ImGui::GetFontSize();
  const float labelLineH = ImGui::GetFontSize() * 0.82f;
  return pad * 2.0f + valueLineH + labelLineH + m.gap * 0.5f;  // + a visible gap between lines
}

// Minimum height for drawHeadlineReadout's stacked "big value / sub-value /
// caption" layout. Same reasoning as statCardHeight above.
float headlineCardHeight(int subLines) {
  const style::Metrics& m = style::metrics();
  const float pad = m.gap * 1.1f;
  const float valueFontSize = ImGui::GetFontSize() * 2.0f;
  const float subFontSize = ImGui::GetFontSize() * 0.95f;
  const float labelSize = ImGui::GetFontSize() * 0.82f;
  const float lineGap = m.gap * 0.4f;
  float h = pad * 2.0f + valueFontSize + lineGap + labelSize;
  h += static_cast<float>(subLines) * (subFontSize + lineGap);
  return h;
}

void drawHeadlineReadout(const char* label, const std::string& value,
                         const std::string& subValue, ImVec4 accent, ImVec2 size,
                         float flashAlpha) {
  const style::Metrics& m = style::metrics();
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::Dummy(size);
  const ImVec2 max(min.x + size.x, min.y + size.y);

  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddRectFilled(min, max, style::u32(style::col::BgRaised), m.radiusLg);
  dl->AddRect(min, max, style::u32(accent, 0.6f), m.radiusLg, 0, m.hairline);

  const float pad = m.gap * 1.1f;
  const float innerWidth = std::max(1.0f, size.x - pad * 2.0f);
  const auto fittedFontSize = [innerWidth](ImFont* font, float requested,
                                           const std::string& text) {
    if (text.empty()) return requested;
    // Headless runs (and any frame before the atlas is built) have no custom
    // font; measuring through a null ImFont* would segfault.
    ImFont* measureFont = font ? font : ImGui::GetFont();
    if (!measureFont) return requested;
    const float measured =
        measureFont->CalcTextSizeA(requested, FLT_MAX, 0.0f, text.c_str()).x;
    if (measured <= innerWidth || measured <= 0.0f) return requested;
    return std::max(requested * innerWidth / measured, ImGui::GetFontSize() * 0.7f);
  };
  const float valueFontSize =
      fittedFontSize(style::fonts::mono(), ImGui::GetFontSize() * 2.0f, value);
  const float subFontSize =
      fittedFontSize(style::fonts::mono(), ImGui::GetFontSize() * 0.95f, subValue);
  const float labelSize = ImGui::GetFontSize() * 0.82f;
  const float lineGap = m.gap * 0.4f;

  dl->PushClipRect(min, max, true);
  float y = min.y + pad;
  dl->AddText(style::fonts::mono(), valueFontSize, ImVec2(min.x + pad, y),
              style::u32(style::col::Text), value.c_str());
  if (flashAlpha > 0.0f) {
    dl->AddText(style::fonts::mono(), valueFontSize, ImVec2(min.x + pad, y),
                style::u32(accent, flashAlpha * 0.72f), value.c_str());
  }
  y += valueFontSize + lineGap;
  if (!subValue.empty()) {
    dl->AddText(style::fonts::mono(), subFontSize, ImVec2(min.x + pad, y), style::u32(accent),
                subValue.c_str());
  }
  dl->AddText(nullptr, labelSize, ImVec2(min.x + pad, max.y - pad - labelSize),
              style::u32(style::col::TextFaint), label);
  dl->PopClipRect();
}

// ------------------------------------------------------------ left rail
void drawSoluteCard(AppState& st) {
  SolubilityState& sb = st.solubility;
  static bool prevUseSketch = sb.useSketch;
  static bool prevOverride = sb.overrideSolute;

  if (ImGui::RadioButton("From sketch", sb.useSketch)) sb.useSketch = true;
  const float secondModeWidth =
      ImGui::GetFrameHeight() + ImGui::GetStyle().ItemInnerSpacing.x +
      ImGui::CalcTextSize("From SMILES").x;
  if (ImGui::GetContentRegionAvail().x >=
      ImGui::GetStyle().ItemSpacing.x + secondModeWidth) {
    ImGui::SameLine();
  }
  if (ImGui::RadioButton("From SMILES", !sb.useSketch)) sb.useSketch = false;

  const bool modeToggled = sb.useSketch != prevUseSketch;
  prevUseSketch = sb.useSketch;

  if (sb.useSketch) {
    if (modeToggled || st.docRevision != sb.sourceRevision) recomputeSoluteFromSketch(st);
  } else {
    if (modeToggled) {
      sb.soluteValid = false;
      sb.soluteError.clear();
    }
    ImGui::SetNextItemWidth(-FLT_MIN);
    const bool enter = widgets::stringInputWithHint(
        "##solute_smiles", "SMILES, e.g. CC(=O)Oc1ccccc1C(=O)O", sb.manualSmiles,
        ImGuiInputTextFlags_EnterReturnsTrue, true);
    if (enter || ImGui::IsItemDeactivatedAfterEdit()) recomputeSoluteFromSmiles(st);
  }

  if (!sb.soluteError.empty()) {
    ImGui::PushStyleColor(ImGuiCol_Text, style::col::Danger);
    ImGui::TextWrapped("%s", sb.soluteError.c_str());
    ImGui::PopStyleColor();
    return;
  }
  if (!sb.soluteValid) {
    textDisabledWrapped("%s", sb.useSketch ? "Draw a structure to describe the solute."
                                           : "Enter a SMILES string and press Enter.");
    return;
  }

  // How the sketch was read: which fragment is the solute, and whether a
  // counter-ion turned it into a salt. This is the difference between
  // modelling an alkaloid and modelling its hydrochloride.
  if (!sb.soluteNote.empty()) {
    ImGui::PushStyleColor(ImGuiCol_Text, style::col::Accent);
    ImGui::TextWrapped("%s", sb.soluteNote.c_str());
    ImGui::PopStyleColor();
  }

  // Two-column property grid fits the narrow rail better than five cards
  // squeezed onto one row.
  const float avail = ImGui::GetContentRegionAvail().x;
  const float spacing = ImGui::GetStyle().ItemSpacing.x;
  const float cardW = (avail - spacing) * 0.5f;
  const float cardH = statCardHeight();
  char mw[32];
  std::snprintf(mw, sizeof(mw), "%.1f", sb.solute.molarMass);
  char mv[32];
  std::snprintf(mv, sizeof(mv), "%.1f", sb.solute.molarVolume);
  char logp[32];
  std::snprintf(logp, sizeof(logp), "%.2f", sb.solute.logP);
  char hansen[64];
  std::snprintf(hansen, sizeof(hansen), "%.1f / %.1f / %.1f", sb.solute.hansen.dispersion,
               sb.solute.hansen.polar, sb.solute.hansen.hydrogenBond);
  char meltPt[32];
  if (sb.solute.meltingPointEstimated) {
    std::snprintf(meltPt, sizeof(meltPt), "%.0f (est)", sb.solute.meltingPoint);
  } else {
    std::snprintf(meltPt, sizeof(meltPt), "%.0f", sb.solute.meltingPoint);
  }
  widgets::statCard("MW (g/mol)", mw, ImVec2(cardW, cardH));
  ImGui::SameLine(0.0f, spacing);
  widgets::statCard("LOGP", logp, ImVec2(cardW, cardH));
  widgets::statCard("MOLAR VOL", mv, ImVec2(cardW, cardH));
  ImGui::SameLine(0.0f, spacing);
  widgets::statCard("TM (C)", meltPt, ImVec2(cardW, cardH));
  widgets::statCard("HANSEN  dD / dP / dH", hansen, ImVec2(avail, cardH));
  widgets::statCard("TM SOURCE",
                    sb.solute.meltingPointEstimated ? "estimated" : "literature",
                    ImVec2(avail, cardH));

  if (sb.solute.meltingPointEstimated && !sb.overrideSolute) {
    // An estimated Tm is the single biggest source of error in the
    // prediction (a Joback group-contribution guess vs. a literature value
    // can be tens of degrees off) -- make the fix impossible to miss.
    ImGui::PushStyleColor(ImGuiCol_Text, style::col::Violet);
    ImGui::TextWrapped("Melting point is a Joback estimate (%.0f C), not a measured value. "
                       "Override it below for an accurate prediction.",
                       sb.solute.meltingPoint);
    ImGui::PopStyleColor();
  }

  const bool overrideClicked = checkboxWithResponsiveLabel(
      "Override Tm / Hansen radius R0", "Enable##override_solute", &sb.overrideSolute);
  if (overrideClicked && sb.overrideSolute && !prevOverride && sb.soluteValid) {
    sb.meltingPointC = static_cast<float>(sb.solute.meltingPoint);
    sb.interactionRadius = static_cast<float>(sb.solute.interactionRadius);
  }
  prevOverride = sb.overrideSolute;
  if (sb.overrideSolute) {
    bool changed = false;
    const bool meltingLabelInline = prepareTrailingLabel("Melting point (C)");
    changed |= ImGui::SliderFloat(meltingLabelInline ? "Melting point (C)" : "##melting_point",
                                  &sb.meltingPointC, -50.0f, 300.0f, "%.1f");
    const bool radiusLabelInline = prepareTrailingLabel("Hansen radius R0");
    changed |= ImGui::SliderFloat(radiusLabelInline ? "Hansen radius R0" : "##hansen_radius",
                                  &sb.interactionRadius, 1.0f, 30.0f, "%.1f");
    if (sb.soluteValid) {
      applySoluteOverrides(sb);
      if (changed) ++sb.soluteVersion;
    }
  }
}

void drawSolventSlot(SolubilityState& sb, int index, std::string& query) {
  static const char* kSlotLabel[3] = {"Solvent A", "Solvent B", "Solvent C"};
  std::string& id = sb.solventIds[static_cast<size_t>(index)];

  // Colour chip, slot label, and combo divide exactly the row's current
  // content width. The measured label reservation keeps the combo from
  // extending beyond a narrow card.
  const float rowWidth = ImGui::GetContentRegionAvail().x;
  const float chip = ImGui::GetFontSize() * 0.9f;
  const float innerSpacing = ImGui::GetStyle().ItemInnerSpacing.x;
  const float itemSpacing = ImGui::GetStyle().ItemSpacing.x;
  const float labelWidth = ImGui::CalcTextSize(kSlotLabel[index]).x;
  const float comboWidth =
      rowWidth - chip - innerSpacing - labelWidth - itemSpacing;
  const bool comboFitsInline = comboWidth >= ImGui::GetFontSize() * 3.0f;
  const ImVec2 chipMin = ImGui::GetCursorScreenPos();
  ImGui::Dummy(ImVec2(chip, chip));
  ImGui::GetWindowDrawList()->AddRectFilled(chipMin, ImVec2(chipMin.x + chip, chipMin.y + chip),
                                            style::u32(slotColor(index)), chip * 0.25f);
  ImGui::SameLine(0.0f, innerSpacing);
  ImGui::TextUnformatted(kSlotLabel[index]);
  if (comboFitsInline) {
    ImGui::SameLine(0.0f, itemSpacing);
    ImGui::SetNextItemWidth(comboWidth);
  } else {
    ImGui::SetNextItemWidth(-FLT_MIN);
  }
  widgets::solventCombo("##solvent_pick", id, query);

  ImGui::SetNextItemWidth(-FLT_MIN);
  ImGui::SliderFloat("##ratio", &sb.ratios[static_cast<size_t>(index)], 0.0f, 10.0f,
                     "%.2f parts");
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "Relative volume parts of this solvent in the blend");
}

// Returns the resolved slots; outComponents holds the normalised blend for
// prediction. Slots that fail to resolve come back nullptr.
std::vector<const sol::Solvent*> drawSolventMixer(
    SolubilityState& sb, std::vector<sol::Component>& outComponents,
    std::vector<std::string>& miscibilityWarnings) {
  static std::array<std::string, 3> filterBuf;
  // Short item text: the trailing "Number of solvents" label already names
  // the control, and a 300 px rail cannot hold both in full.
  static const char* kCountItems[] = {"1 - pure", "2 - binary", "3 - ternary"};
  miscibilityWarnings.clear();
  int countIdx = sb.solventCount - 1;
  const bool countLabelInline = prepareTrailingLabel("Number of solvents");
  if (ImGui::Combo(countLabelInline ? "Number of solvents" : "##solvent_count", &countIdx,
                   kCountItems, 3)) {
    sb.solventCount = countIdx + 1;
  }

  std::vector<const sol::Solvent*> chosen(static_cast<size_t>(sb.solventCount), nullptr);
  for (int i = 0; i < sb.solventCount; ++i) {
    ImGui::PushID(i);
    drawSolventSlot(sb, i, filterBuf[static_cast<size_t>(i)]);
    ImGui::PopID();
    chosen[static_cast<size_t>(i)] = sol::findSolvent(sb.solventIds[static_cast<size_t>(i)]);
    const sol::Solvent* found = chosen[static_cast<size_t>(i)];
    if (found && sb.ratios[static_cast<size_t>(i)] > 0.0f) {
      outComponents.push_back({found, static_cast<double>(sb.ratios[static_cast<size_t>(i)])});
    }
  }

  if (!outComponents.empty()) {
    const sol::Mixture mixture = sol::blend(outComponents);
    ImGui::PushStyleColor(ImGuiCol_Text, style::col::TextDim);
    ImGui::TextWrapped("Blend Hansen (MPa^0.5)");
    ImGui::TextWrapped("dD %.1f | dP %.1f | dH %.1f", mixture.hansen.dispersion,
                       mixture.hansen.polar, mixture.hansen.hydrogenBond);
    ImGui::TextWrapped("Density %.3f g/mL", mixture.density);
    ImGui::PopStyleColor();
  } else {
    ImGui::TextDisabled("Select at least one solvent.");
  }

  // A blend prediction assumes one homogeneous phase; collect the pairs
  // where that assumption is physically wrong so the results workspace can
  // surface them beside the model diagnostics.
  for (size_t i = 0; i < chosen.size(); ++i) {
    for (size_t j = i + 1; j < chosen.size(); ++j) {
      if (!chosen[i] || !chosen[j]) continue;
      if (!sol::miscibleWith(*chosen[i], *chosen[j])) {
        miscibilityWarnings.push_back(
            chosen[i]->name + " and " + chosen[j]->name +
            " are immiscible. The blend forms two phases; the prediction assumes a "
            "homogeneous mixture.");
      }
    }
  }
  return chosen;
}

void drawScreeningTable(SolubilityState& sb) {
  recomputeScreeningIfNeeded(sb);
  if (!sb.soluteValid) {
    ImGui::TextDisabled("Load a solute to rank every pure solvent.");
    return;
  }
  if (sb.screening.empty()) {
    ImGui::TextDisabled("No solvents to rank.");
    return;
  }

  textDisabledWrapped("Click a row to load it as Solvent A.");

  std::vector<double> displayValues;
  std::vector<std::string> labels;
  std::vector<std::string> annotations;
  displayValues.reserve(sb.screening.size());
  labels.reserve(sb.screening.size());
  annotations.reserve(sb.screening.size());

  double bestValue = 0.0;
  double smallestPositive = 0.0;
  for (const sol::ScreenRow& row : sb.screening) {
    const double value =
        toDisplayUnits(row.prediction.gramsPerMillilitre, sb.solute.molarMass, sb.units);
    displayValues.push_back(value);
    labels.push_back(row.solvent->name);
    if (std::isfinite(value) && value > 0.0) {
      annotations.push_back(
          formatUnits(row.prediction.gramsPerMillilitre, sb.solute.molarMass, sb.units));
      bestValue = std::max(bestValue, value);
      if (smallestPositive == 0.0 || value < smallestPositive) smallestPositive = value;
    } else {
      annotations.emplace_back("--");
    }
  }

  double floorLog = 0.0;
  double logSpan = 0.0;
  if (bestValue > 0.0) {
    const double bestLog = std::log10(bestValue);
    // A two-decade log10 window. Solubility across 45 solvents spans orders of
    // magnitude, so a linear bar shows only the winner; but a wide window is no
    // better, because it compresses the handful of solvents anyone would
    // actually choose between into indistinguishable full-length bars. Two
    // decades makes the best solvent a full bar, a ten-times-worse one a half
    // bar, and anything a hundred times worse an empty one -- which is the
    // honest reading: at that point it is not a candidate.
    floorLog = std::max(std::log10(smallestPositive), bestLog - 2.0);
    logSpan = bestLog - floorLog;
  }

  std::vector<charts::BarRow> bars;
  bars.reserve(sb.screening.size());
  for (size_t i = 0; i < sb.screening.size(); ++i) {
    double scaledValue = 0.0;
    if (std::isfinite(displayValues[i]) && displayValues[i] > 0.0) {
      scaledValue = logSpan > 0.0
                        ? std::clamp((std::log10(displayValues[i]) - floorLog) / logSpan,
                                     0.0, 1.0)
                        : 1.0;
    }
    bars.push_back({labels[i].c_str(), scaledValue, annotations[i].c_str(),
                    style::col::Teal,
                    sb.screening[i].solvent->id == sb.solventIds[0]});
  }

  const float viewportHeight =
      std::max(ImGui::GetFontSize() * 7.0f, ImGui::GetContentRegionAvail().y);
  if (ImGui::BeginChild("##screen_chart_scroll", ImVec2(0.0f, viewportHeight),
                        ImGuiChildFlags_Borders)) {
    const float chartHeight =
        static_cast<float>(bars.size()) * ImGui::GetTextLineHeightWithSpacing();
    const int clicked = charts::rankedBars(
        "##screen_ranked_bars", bars.data(), static_cast<int>(bars.size()),
        ImVec2(ImGui::GetContentRegionAvail().x, chartHeight));
    if (clicked >= 0) {
      const sol::ScreenRow& row = sb.screening[static_cast<size_t>(clicked)];
      sb.solventIds[0] = row.solvent->id;
      sb.statusMessage = "Solvent A := " + row.solvent->name;
    }
  }
  ImGui::EndChild();
}

// Kirkwood-Buff integrals, the Yalkowsky ideal term and the Hildebrand
// parameters behind the current number: the model showing its working.
void drawTheoryReadout(SolubilityState& sb, const std::vector<const sol::Solvent*>& chosen) {
  std::vector<sol::Component> components;
  for (size_t i = 0; i < chosen.size() && i < sb.ratios.size(); ++i) {
    if (chosen[i]) components.push_back({chosen[i], sb.ratios[i]});
  }
  if (components.empty()) {
    ImGui::TextDisabled("Pick a solvent to see the theory readout.");
    return;
  }

  const sol::KBResult kb = sol::kirkwoodBuff(sb.solute, components, sb.temperatureC);
  const sol::Mixture mixture = sol::blend(components);
  const auto hildebrand = [](const sol::Hansen& h) {
    return std::sqrt(h.dispersion * h.dispersion + h.polar * h.polar +
                     h.hydrogenBond * h.hydrogenBond);
  };
  const double deltaSolute = hildebrand(sb.solute.hansen);
  const double deltaSolvent = hildebrand(mixture.hansen);

  const bool mono = style::pushFont(style::fonts::mono());
  ImGui::TextWrapped("Kirkwood-Buff | G11 %+.1f cm3/mol | G12 %+.1f cm3/mol | "
                     "ln gamma %+.2f",
                     kb.g11, kb.g12, kb.lnGammaInf);
  style::popFont(mono);
  if (ImGui::IsItemHovered()) {
    ImGui::SetTooltip(
        "G11 = R T kappa_T - V1 (fluctuation-theory identity; kappa_T measured)\n"
        "G12 = G11 - V1 ln gamma (Ben-Naim inversion of KB theory)\n"
        "G12 > G11 means the solute is preferentially solvated.");
  }
  textDisabledWrapped("kappa_T %.3f GPa^-1 (%s)", kb.kappaT,
                      kb.kappaKnown ? (kb.kappaSource.empty() ? "literature"
                                                              : kb.kappaSource.c_str())
                                    : "no data: incompressible limit");

  ImGui::Spacing();
  const bool mono2 = style::pushFont(style::fonts::mono());
  ImGui::TextWrapped("Yalkowsky | ln x_ideal %+.2f | x_ideal %.4g | "
                     "dS_fus %.1f J/mol K",
                     kb.lnIdeal, sb.prediction.idealMoleFraction,
                     sb.solute.entropyOfFusion);
  style::popFont(mono2);
  if (ImGui::IsItemHovered()) {
    ImGui::SetTooltip(
        "Ideal solubility: ln x = -dS_fus (Tm - T) / (R T)\n"
        "General Solubility Equation (aqueous share):\n"
        "  log S [mol/L] = 0.5 - 0.01 (Tm[C] - 25) - logP");
  }

  const bool mono3 = style::pushFont(style::fonts::mono());
  ImGui::TextWrapped("Hildebrand | solute %.1f | solvent %.1f | |d| %.1f MPa^0.5",
                     deltaSolute, deltaSolvent, std::fabs(deltaSolute - deltaSolvent));
  style::popFont(mono3);
  if (ImGui::IsItemHovered()) {
    ImGui::SetTooltip("Total Hildebrand parameter: d^2 = dD^2 + dP^2 + dH^2");
  }

  // Ionisation: pKa, the pH in force, and the Born penalty for asking an ion
  // to live in a low-dielectric solvent. Only shown when it applies.
  if (sb.prediction.ionicPath) {
    const bool mono4 = style::pushFont(style::fonts::mono());
    ImGui::TextWrapped("Ionisation | pKa ~%.1f | pH %.1f%s | ionised %.0f%% | "
                       "Born %.1f decades",
                       sb.prediction.pKa, sb.prediction.pH,
                       sb.prediction.pHSelfBuffered ? " (self)" : "",
                       sb.prediction.ionisedFraction * 100.0,
                       sb.prediction.bornPenaltyDecades);
    style::popFont(mono4);
    if (ImGui::IsItemHovered()) {
      ImGui::SetTooltip(
          "Henderson-Hasselbalch: S = S_0 (1 + 10^(pKa - pH)) for a base\n"
          "Salt crystal: the ionised branch only, so its solubility tracks\n"
          "  the medium's ability to support ions\n"
          "Born transfer: ddG = sum_i [N_A e^2 / (8 pi eps_0 r_i)] (1/eps - 1/eps_water)\n"
          "pKa is a group-contribution estimate, not a measured value.");
    }
    if (!sb.prediction.ionNote.empty()) {
      textDisabledWrapped("%s", sb.prediction.ionNote.c_str());
    }
  }

  // Honest uncertainty: measured on the anchor validation set (see
  // tests/test_kb_accuracy.cpp). Anchored pairs are measured, not estimated.
  ImGui::Spacing();
  if (sb.prediction.anchored) {
    widgets::badge("MEASURED", style::col::Success);
    textDisabledWrapped("Literature value; no model uncertainty applied.");
  } else {
    const double value = sb.prediction.gramsPerMillilitre;
    textDisabledWrapped("95%% interval %.3g - %.3g g/mL "
                        "(model sigma 0.75 log on the validation set)",
                        value / 30.0, value * 30.0);
  }
}

// ---------------------------------------------------------- hero + results

void drawResultHero(SolubilityState& sb, float flashAlpha) {
  const sol::Prediction& p = sb.prediction;
  const float avail = ImGui::GetContentRegionAvail().x;
  const float spacing = ImGui::GetStyle().ItemSpacing.x;
  const float statH = statCardHeight();
  const float rowGap = ImGui::GetStyle().ItemSpacing.y;
  const bool horizontal = avail >= ImGui::GetFontSize() * 31.0f;
  const float headlineW = horizontal ? avail * 0.34f : avail;
  const float headlineH =
      horizontal ? std::max(headlineCardHeight(1), statH * 2.0f + rowGap)
                 : headlineCardHeight(1);

  const double gPerMl = p.gramsPerMillilitre;
  const std::string headline = formatUnits(gPerMl, sb.solute.molarMass, sb.units);
  const std::string sub = sb.units == 3
                              ? formatUnits(gPerMl, sb.solute.molarMass, 0)
                              : formatUnits(gPerMl, sb.solute.molarMass, 3);
  drawHeadlineReadout("PREDICTED SOLUBILITY", headline, sub, style::col::Accent,
                      ImVec2(headlineW, headlineH), flashAlpha);

  char moleFrac[32];
  std::snprintf(moleFrac, sizeof(moleFrac), "%.4g", p.moleFraction);
  char idealFrac[32];
  std::snprintf(idealFrac, sizeof(idealFrac), "%.4g", p.idealMoleFraction);
  char ra[32];
  std::snprintf(ra, sizeof(ra), "%.2f", p.ra);
  char red[32];
  std::snprintf(red, sizeof(red), "%.2f", p.relativeEnergyDifference);
  char gamma[32];
  std::snprintf(gamma, sizeof(gamma), "%.3g", p.activityCoefficient);
  char chi[32];
  std::snprintf(chi, sizeof(chi), "%.3f", p.chi);
  const char* labels[] = {"MOLE FRAC", "IDEAL X", "Ra (MPa^0.5)",
                          "RED (Ra/R0)", "GAMMA", "CHI (F-H)"};
  const char* values[] = {moleFrac, idealFrac, ra, red, gamma, chi};

  if (horizontal) ImGui::SameLine(0.0f, spacing);
  ImGui::BeginGroup();
  const float statAvail = horizontal ? avail - headlineW - spacing : avail;
  const int columns = (!horizontal && statAvail < ImGui::GetFontSize() * 22.0f) ? 2 : 3;
  const float statW =
      std::max(1.0f, (statAvail - spacing * static_cast<float>(columns - 1)) /
                         static_cast<float>(columns));
  for (int i = 0; i < 6; ++i) {
    if (i % columns != 0) ImGui::SameLine(0.0f, spacing);
    widgets::statCard(labels[i], values[i], ImVec2(statW, statH));
  }
  ImGui::EndGroup();
}

// ------------------------------------------------------- binary ratio plot

void drawBinarySweepPlot(SolubilityState& sb, const sol::Solvent& a, const sol::Solvent& b,
                         ImVec2 canvasSize) {
  if (sb.sweep.empty()) {
    ImGui::TextDisabled("Not enough data for a ratio plot.");
    return;
  }

  std::vector<const sol::SweepPoint*> sorted;
  sorted.reserve(sb.sweep.size());
  for (const sol::SweepPoint& sp : sb.sweep) sorted.push_back(&sp);
  std::sort(sorted.begin(), sorted.end(), [](const sol::SweepPoint* x, const sol::SweepPoint* y) {
    return x->fractions[0] < y->fractions[0];
  });

  double lo = sorted.front()->prediction.gramsPerMillilitre;
  double hi = lo;
  for (const sol::SweepPoint* sp : sorted) {
    lo = std::min(lo, sp->prediction.gramsPerMillilitre);
    hi = std::max(hi, sp->prediction.gramsPerMillilitre);
  }
  if (hi - lo < 1e-12) hi = lo + std::max(1e-9, std::fabs(lo) * 0.1 + 1e-9);

  // Log mode maps through log10 with a floor six orders under the max so a
  // zero-valued endpoint cannot blow the axis up.
  const bool logY = sb.logScale && hi > 0.0;
  double yLo = lo, yHi = hi;
  if (logY) {
    const double floorV = std::max(hi * 1e-6, 1e-12);
    yLo = std::log10(std::max(lo, floorV));
    yHi = std::log10(hi);
    if (yHi - yLo < 1e-9) yHi = yLo + 1.0;
  }
  const auto mapValue = [&](double v) {
    if (!logY) return v;
    const double floorV = std::max(hi * 1e-6, 1e-12);
    return std::log10(std::max(v, floorV));
  };

  const sol::SweepPoint* peak =
      (sb.sweepPeakIndex >= 0 && static_cast<size_t>(sb.sweepPeakIndex) < sb.sweep.size())
          ? &sb.sweep[static_cast<size_t>(sb.sweepPeakIndex)]
          : nullptr;

  const double molarMass = sb.solute.molarMass;
  const std::string unitCaption =
      kUnitLabels[static_cast<size_t>(std::clamp(sb.units, 0, 3))];
  float widestAxisLabel = ImGui::CalcTextSize(unitCaption.c_str()).x;
  for (int i = 0; i <= 4; ++i) {
    const float t = static_cast<float>(i) / 4.0f;
    const double mapped = yLo + (yHi - yLo) * (1.0 - t);
    const double value = logY ? std::pow(10.0, mapped) : mapped;
    const std::string label =
        formatSolubility(toDisplayUnits(value, molarMass, sb.units));
    widestAxisLabel = std::max(widestAxisLabel, ImGui::CalcTextSize(label.c_str()).x);
  }
  const float padR =
      std::min(20.0f * uiScale(), std::max(0.0f, canvasSize.x * 0.18f));
  const float desiredPadL =
      std::max(62.0f * uiScale(), widestAxisLabel + 12.0f * uiScale());
  const float padL =
      std::min(desiredPadL, std::max(0.0f, canvasSize.x - padR - 1.0f));
  const float padB =
      std::min(46.0f * uiScale(), std::max(0.0f, canvasSize.y * 0.35f));
  const float padT =
      std::min(16.0f * uiScale(), std::max(0.0f, canvasSize.y - padB - 1.0f));

  const ImVec2 origin = ImGui::GetCursorScreenPos();
  ImGui::InvisibleButton("##binary_plot", canvasSize);
  const bool hovered = ImGui::IsItemHovered();

  // Click or drag sets the working ratio directly from the x position: the
  // graph is the ratio control, the parts sliders just display it.
  if (ImGui::IsItemClicked() ||
      (ImGui::IsItemActive() && ImGui::IsMouseDragging(ImGuiMouseButton_Left))) {
    float t = (ImGui::GetMousePos().x - (origin.x + padL)) /
              std::max(1.0f, canvasSize.x - padL - padR);
    t = std::clamp(t, 0.0f, 1.0f);
    sb.ratios[0] = t;
    sb.ratios[1] = 1.0f - t;
  }

  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddRectFilled(origin, ImVec2(origin.x + canvasSize.x, origin.y + canvasSize.y),
                    style::u32(style::col::BgSurface), style::metrics().radiusMd);
  dl->PushClipRect(origin, ImVec2(origin.x + canvasSize.x, origin.y + canvasSize.y), true);

  // Axis padding accounts for the widest formatted tick and unit caption,
  // preventing scientific-notation values from clipping at compact widths.
  const ImVec2 plotMin(origin.x + padL, origin.y + padT);
  const ImVec2 plotMax(origin.x + canvasSize.x - padR, origin.y + canvasSize.y - padB);
  dl->AddRectFilled(plotMin, plotMax, style::u32(style::col::BgPanel));
  dl->AddRect(plotMin, plotMax, style::u32(style::col::BorderStrong));

  const auto toScreen = [&](double fracA, double value) {
    const float x = plotMin.x + static_cast<float>(fracA) * (plotMax.x - plotMin.x);
    const float t = static_cast<float>((mapValue(value) - yLo) / (yHi - yLo));
    const float y = plotMax.y - t * (plotMax.y - plotMin.y);
    return ImVec2(x, y);
  };
  const auto unitLabel = [&](double v) { return formatUnits(v, molarMass, sb.units); };

  // x-axis ticks: volume % of solvent A, 0/25/50/75/100.
  for (int i = 0; i <= 4; ++i) {
    const float t = static_cast<float>(i) / 4.0f;
    const float x = plotMin.x + t * (plotMax.x - plotMin.x);
    dl->AddLine(ImVec2(x, plotMin.y), ImVec2(x, plotMax.y), style::u32(style::col::Border, 0.35f));
    char buf[8];
    std::snprintf(buf, sizeof(buf), "%.0f%%", t * 100.0f);
    const ImVec2 tickSize = ImGui::CalcTextSize(buf);
    dl->AddText(ImVec2(x - tickSize.x * 0.5f, plotMax.y + 4.0f * uiScale()),
               style::u32(style::col::TextDim), buf);
  }
  // y-axis ticks: interpolated in the mapped (linear or log) space. Ticks
  // carry bare numbers; the unit goes in one caption above the axis, so the
  // longest labels cannot clip against the rail's child window.
  for (int i = 0; i <= 4; ++i) {
    const float t = static_cast<float>(i) / 4.0f;
    const double mapped = yLo + (yHi - yLo) * (1.0 - t);
    const double value = logY ? std::pow(10.0, mapped) : mapped;
    const float y = plotMin.y + t * (plotMax.y - plotMin.y);
    dl->AddLine(ImVec2(plotMin.x, y), ImVec2(plotMax.x, y), style::u32(style::col::Border, 0.2f));
    const std::string label = formatSolubility(toDisplayUnits(value, molarMass, sb.units));
    const ImVec2 labelSize = ImGui::CalcTextSize(label.c_str());
    dl->AddText(ImVec2(plotMin.x - labelSize.x - 6.0f * uiScale(), y - labelSize.y * 0.5f),
               style::u32(style::col::TextDim), label.c_str());
  }
  // Unit caption, top-left of the axis.
  dl->AddText(ImVec2(plotMin.x - ImGui::CalcTextSize(unitCaption.c_str()).x - 6.0f * uiScale(),
                     plotMin.y - 16.0f * uiScale()),
              style::u32(style::col::TextFaint), unitCaption.c_str());

  // Axis titles: the pure endpoints of the blend (fracA=0 is 100% b).
  const float titleWidth = std::max(1.0f, (plotMax.x - plotMin.x) * 0.46f);
  const std::string leftTitle = ellipsizeText("100% " + b.name, titleWidth);
  const std::string rightTitle = ellipsizeText("100% " + a.name, titleWidth);
  const ImVec2 rightTitleSize = ImGui::CalcTextSize(rightTitle.c_str());
  dl->AddText(ImVec2(plotMin.x, plotMax.y + 18.0f * uiScale()), style::u32(style::col::TextFaint),
              leftTitle.c_str());
  dl->AddText(ImVec2(plotMax.x - rightTitleSize.x, plotMax.y + 18.0f * uiScale()),
              style::u32(style::col::TextFaint), rightTitle.c_str());

  // Shaded area under the curve -- one quad per segment rather than a single
  // convex fill, so the rise-then-fall co-solvency shape fills correctly
  // even though the region under a peaked curve is not convex.
  const float drawProgress = sweepDrawProgress(sb);
  const float revealX =
      plotMin.x + std::max(0.001f, drawProgress) * (plotMax.x - plotMin.x);
  dl->PushClipRect(plotMin, ImVec2(revealX, plotMax.y), true);
  for (size_t i = 0; i + 1 < sorted.size(); ++i) {
    const ImVec2 p0 = toScreen(sorted[i]->fractions[0], sorted[i]->prediction.gramsPerMillilitre);
    const ImVec2 p1 =
        toScreen(sorted[i + 1]->fractions[0], sorted[i + 1]->prediction.gramsPerMillilitre);
    const ImVec2 b0(p0.x, plotMax.y);
    const ImVec2 b1(p1.x, plotMax.y);
    dl->AddQuadFilled(p0, p1, b1, b0, style::u32(style::col::Teal, 0.16f));
  }

  std::vector<ImVec2> points;
  points.reserve(sorted.size());
  for (const sol::SweepPoint* sp : sorted)
    points.push_back(toScreen(sp->fractions[0], sp->prediction.gramsPerMillilitre));
  if (points.size() >= 2) {
    dl->AddPolyline(points.data(), static_cast<int>(points.size()), style::u32(style::col::Teal),
                    0, 2.4f);
  } else if (points.size() == 1) {
    dl->AddCircleFilled(points.front(), 3.5f, style::u32(style::col::Teal));
  }
  dl->PopClipRect();

  // Mark the co-solvency maximum: drop-line, filled marker, and a
  // composition + value label.
  if (peak) {
    const ImVec2 peakPt = toScreen(peak->fractions[0], peak->prediction.gramsPerMillilitre);
    const float reveal =
        markerReveal(drawProgress, static_cast<float>(peak->fractions[0]));
    if (reveal > 0.0f) {
      const float pulse = markerPulse();
      dl->AddLine(ImVec2(peakPt.x, plotMax.y), peakPt,
                  style::u32(style::col::Accent, 0.7f * reveal), 1.5f);
      dl->AddCircleFilled(peakPt, 5.5f, style::u32(style::col::Accent, reveal));
      dl->AddCircle(peakPt, 8.5f + pulse * 4.0f,
                    style::u32(style::col::Accent,
                               reveal * (0.70f - pulse * 0.48f)),
                    0, 1.5f);

      const std::string peakLabel =
          "peak " + unitLabel(peak->prediction.gramsPerMillilitre) + " at " +
          std::to_string(static_cast<int>(std::lround(peak->fractions[0] * 100.0))) +
          "% " + a.name;
      const std::string fittedPeakLabel =
          ellipsizeText(peakLabel, std::max(1.0f, plotMax.x - plotMin.x - 8.0f));
      const ImVec2 peakLabelSize = ImGui::CalcTextSize(fittedPeakLabel.c_str());
      const float labelMaxX = std::max(plotMin.x, plotMax.x - peakLabelSize.x);
      float lx = std::clamp(peakPt.x - peakLabelSize.x * 0.5f, plotMin.x, labelMaxX);
      float ly = peakPt.y - peakLabelSize.y - 10.0f * uiScale();
      if (ly < plotMin.y) ly = peakPt.y + 10.0f * uiScale();
      dl->AddRectFilled(ImVec2(lx - 4.0f, ly - 2.0f),
                        ImVec2(lx + peakLabelSize.x + 4.0f,
                               ly + peakLabelSize.y + 2.0f),
                        style::u32(style::col::BgRaised, 0.9f * reveal), 3.0f);
      dl->AddText(ImVec2(lx, ly), style::u32(style::col::Accent, reveal),
                  fittedPeakLabel.c_str());
    }
  }

  // Mark the current working ratio -- a distinct colour from the peak so
  // the two are never confused at a glance.
  const double totalRatio = static_cast<double>(sb.ratios[0]) + static_cast<double>(sb.ratios[1]);
  if (totalRatio > 1e-9) {
    const double fracA = static_cast<double>(sb.ratios[0]) / totalRatio;
    const ImVec2 marker = toScreen(fracA, sb.prediction.gramsPerMillilitre);
    const float reveal = markerReveal(drawProgress, static_cast<float>(fracA));
    if (reveal > 0.0f) {
      const float pulse = markerPulse(1.9f);
      dl->AddLine(ImVec2(marker.x, plotMin.y), ImVec2(marker.x, plotMax.y),
                  style::u32(style::col::Violet, 0.55f * reveal), 1.5f);
      dl->AddCircleFilled(marker, 4.5f, style::u32(style::col::Violet, reveal));
      dl->AddCircle(marker, 7.0f + pulse * 3.5f,
                    style::u32(style::col::Violet,
                               reveal * (0.66f - pulse * 0.44f)),
                    0, 1.35f);
      const ImVec2 currentSize = ImGui::CalcTextSize("current");
      const float currentX =
          std::clamp(marker.x + 6.0f, origin.x,
                     std::max(origin.x, origin.x + canvasSize.x - currentSize.x));
      dl->AddText(ImVec2(currentX, plotMin.y + 2.0f),
                  style::u32(style::col::Violet, reveal), "current");
    }
  }

  if (hovered) {
    const ImVec2 mouse = ImGui::GetMousePos();
    float t = (mouse.x - plotMin.x) / std::max(1.0f, plotMax.x - plotMin.x);
    t = std::clamp(t, 0.0f, 1.0f);
    size_t nearest = 0;
    double best = 1e18;
    for (size_t i = 0; i < sorted.size(); ++i) {
      const double d = std::fabs(sorted[i]->fractions[0] - t);
      if (d < best) {
        best = d;
        nearest = i;
      }
    }
    const sol::SweepPoint& sp = *sorted[nearest];
    dl->AddLine(ImVec2(mouse.x, plotMin.y), ImVec2(mouse.x, plotMax.y),
               style::u32(style::col::TextFaint, 0.5f));
    const std::string tip = a.name + " " +
                            std::to_string(static_cast<int>(std::lround(sp.fractions[0] * 100.0))) +
                            "%  /  " + b.name + " " +
                            std::to_string(static_cast<int>(std::lround(sp.fractions[1] * 100.0))) +
                            "%   ->   " + unitLabel(sp.prediction.gramsPerMillilitre);
    ImGui::SetTooltip("%s", tip.c_str());
  }

  dl->PopClipRect();
}

// ------------------------------------------------------------ ternary plot

void setTernaryComposition(SolubilityState& sb, float fractionAInPair,
                           float fractionC) {
  fractionAInPair = std::clamp(fractionAInPair, 0.0f, 1.0f);
  fractionC = std::clamp(fractionC, 0.0f, 1.0f);
  const float pairFraction = 1.0f - fractionC;
  sb.ratios[0] = pairFraction * fractionAInPair;
  sb.ratios[1] = pairFraction * (1.0f - fractionAInPair);
  sb.ratios[2] = fractionC;
}

void drawTernarySweepPlot(SolubilityState& sb, const sol::Solvent& a,
                          const sol::Solvent& b, const sol::Solvent& c,
                          ImVec2 canvasSize) {
  const TernarySurfaceMesh& mesh = sb.ternarySurface;
  if (mesh.nodes.empty() || mesh.quads.empty()) {
    ImGui::TextDisabled("Not enough data for a ternary surface.");
    return;
  }

  const bool logY = sb.logScale;
  const double yLo = logY ? mesh.logMinimum : mesh.minimum;
  const double yHi = logY ? mesh.logMaximum : mesh.maximum;
  const double molarMass = sb.solute.molarMass;
  const std::string unitCaption =
      kUnitLabels[static_cast<size_t>(std::clamp(sb.units, 0, 3))];
  const auto rawAtHeight = [&](double height) {
    const double mapped = yLo + (yHi - yLo) * height;
    return logY ? std::pow(10.0, mapped) : mapped;
  };
  const auto normaliseValue = [&](double value) {
    const double mapped =
        logY ? std::log10(std::max(value, mesh.logFloor)) : value;
    return static_cast<float>(
        std::clamp((mapped - yLo) / std::max(1e-12, yHi - yLo), 0.0, 1.0));
  };
  const auto cubePosition =
      [&](const TernarySurfaceNode& node) -> const std::array<float, 3>& {
    return logY ? node.cubeLog : node.cubeLinear;
  };

  const float fontSize = ImGui::GetFontSize();
  const float gap = style::metrics().gap;
  float widestYLabel = ImGui::CalcTextSize(unitCaption.c_str()).x;
  float widestLegendLabel = 0.0f;
  for (int i = 0; i <= 4; ++i) {
    const double raw = rawAtHeight(1.0 - static_cast<double>(i) / 4.0);
    const std::string label =
        formatSolubility(toDisplayUnits(raw, molarMass, sb.units));
    widestYLabel =
        std::max(widestYLabel, ImGui::CalcTextSize(label.c_str()).x);
    if (i % 2 == 0) {
      const std::string legendLabel =
          formatUnits(raw, molarMass, sb.units);
      widestLegendLabel = std::max(
          widestLegendLabel, ImGui::CalcTextSize(legendLabel.c_str()).x);
    }
  }

  // The depth angle is invariant; resizing changes only its length.
  constexpr float kDepthTan = 0.5317094f;  // tan(28 degrees)
  const float leftPad =
      std::min(std::max(widestYLabel + gap, fontSize * 4.2f),
               canvasSize.x * 0.34f);
  const float desiredLegend =
      widestLegendLabel + fontSize * 1.75f + gap;
  const float rightPad =
      std::min(std::max(fontSize * 4.8f, desiredLegend),
               canvasSize.x * 0.28f);
  const float topPad =
      std::min(fontSize * 1.35f, canvasSize.y * 0.12f);
  const float bottomPad =
      std::min(fontSize * 3.7f, canvasSize.y * 0.28f);
  const float innerW =
      std::max(1.0f, canvasSize.x - leftPad - rightPad);
  const float innerH =
      std::max(1.0f, canvasSize.y - topPad - bottomPad);
  float depthX =
      std::min(innerW * 0.22f, innerH * 0.24f / kDepthTan);
  depthX = std::max(0.0f, std::min(depthX, innerW * 0.42f));
  const float depthY = depthX * kDepthTan;
  const float plotW = std::max(1.0f, innerW - depthX);
  const float plotH = std::max(1.0f, innerH - depthY);

  const ImVec2 origin = ImGui::GetCursorScreenPos();
  ImGui::InvisibleButton("##ternary_cube", canvasSize);
  const bool hovered = ImGui::IsItemHovered();
  const bool clicked = ImGui::IsItemClicked(ImGuiMouseButton_Left);
  const bool dragging =
      ImGui::IsItemActive() && ImGui::IsMouseDragging(ImGuiMouseButton_Left);
  ImDrawList* dl = ImGui::GetWindowDrawList();
  const ImVec2 canvasMax(origin.x + canvasSize.x, origin.y + canvasSize.y);
  dl->AddRectFilled(origin, canvasMax, style::u32(style::col::BgSurface),
                    style::metrics().radiusMd);
  dl->PushClipRect(origin, canvasMax, true);

  const ImVec2 frontBL(origin.x + leftPad,
                       origin.y + canvasSize.y - bottomPad);
  const ImVec2 frontBR(frontBL.x + plotW, frontBL.y);
  const ImVec2 frontTL(frontBL.x, frontBL.y - plotH);
  const ImVec2 frontTR(frontBR.x, frontTL.y);
  const ImVec2 depth(depthX, -depthY);
  const auto offset = [&](ImVec2 point, float fraction) {
    return ImVec2(point.x + depth.x * fraction,
                  point.y + depth.y * fraction);
  };
  const ImVec2 backBL = offset(frontBL, 1.0f);
  const ImVec2 backBR = offset(frontBR, 1.0f);
  const ImVec2 backTL = offset(frontTL, 1.0f);
  const ImVec2 backTR = offset(frontTR, 1.0f);
  const auto screenAt = [&](float ratioB, float fractionC, float height) {
    return ImVec2(frontBL.x + plotW * ratioB + depth.x * fractionC,
                  frontBL.y - plotH * height + depth.y * fractionC);
  };
  const auto screenNode = [&](const TernarySurfaceNode& node) {
    const std::array<float, 3>& position = cubePosition(node);
    return screenAt(position[0], position[2], position[1]);
  };
  const auto surfaceValueAt = [&](float ratioB, float fractionC) {
    const float ratioCoordinate =
        std::clamp(ratioB, 0.0f, 1.0f) * mesh.ratioQuads;
    const float cCoordinate =
        std::clamp(fractionC, 0.0f, 1.0f) * mesh.cQuads;
    const int ratioIndex = std::min(
        mesh.ratioQuads - 1, static_cast<int>(std::floor(ratioCoordinate)));
    const int cIndex = std::min(
        mesh.cQuads - 1, static_cast<int>(std::floor(cCoordinate)));
    const float ratioT = ratioCoordinate - static_cast<float>(ratioIndex);
    const float cT = cCoordinate - static_cast<float>(cIndex);
    const size_t rowWidth = static_cast<size_t>(mesh.ratioQuads + 1);
    const size_t n00 =
        static_cast<size_t>(cIndex) * rowWidth +
        static_cast<size_t>(ratioIndex);
    const size_t n10 = n00 + 1;
    const size_t n01 = n00 + rowWidth;
    const size_t n11 = n01 + 1;
    const double front =
        mesh.nodes[n00].point.prediction.gramsPerMillilitre *
            (1.0f - ratioT) +
        mesh.nodes[n10].point.prediction.gramsPerMillilitre * ratioT;
    const double back =
        mesh.nodes[n01].point.prediction.gramsPerMillilitre *
            (1.0f - ratioT) +
        mesh.nodes[n11].point.prediction.gramsPerMillilitre * ratioT;
    return front * (1.0f - cT) + back * cT;
  };

  // Three faint faces and their grids establish the cube behind the surface.
  dl->AddQuadFilled(backTL, backTR, backBR, backBL,
                    style::u32(style::col::BgPanel, 0.66f));
  dl->AddQuadFilled(frontTL, backTL, backBL, frontBL,
                    style::u32(style::col::BgPanel, 0.48f));
  dl->AddQuadFilled(frontBL, frontBR, backBR, backBL,
                    style::u32(style::col::BgDeep, 0.42f));
  for (int i = 1; i < 4; ++i) {
    const float t = static_cast<float>(i) / 4.0f;
    const ImVec2 rearY0(backTL.x, backTL.y + plotH * t);
    const ImVec2 rearY1(backTR.x, backTR.y + plotH * t);
    const ImVec2 frontY(frontTL.x, frontTL.y + plotH * t);
    dl->AddLine(rearY0, rearY1,
                style::u32(style::col::Border, 0.28f));
    dl->AddLine(frontY, rearY0,
                style::u32(style::col::Border, 0.22f));
    const ImVec2 frontX(frontBL.x + plotW * t, frontBL.y);
    const ImVec2 backX(backBL.x + plotW * t, backBL.y);
    dl->AddLine(ImVec2(backTL.x + plotW * t, backTL.y), backX,
                style::u32(style::col::Border, 0.24f));
    dl->AddLine(frontX, backX,
                style::u32(style::col::Border, 0.22f));
    dl->AddLine(offset(frontBL, t), offset(frontBR, t),
                style::u32(style::col::Border, 0.18f));
  }

  struct SurfacePick {
    int quad = -1;
    float ratioA = 0.5f;
    float fractionC = 0.0f;
    double value = 0.0;
    std::array<double, 3> fractions{};
  };
  SurfacePick pick;
  const ImVec2 mouse = ImGui::GetMousePos();
  if (hovered) {
    for (size_t quadIndex = 0; quadIndex < mesh.quads.size(); ++quadIndex) {
      const TernarySurfaceQuad& quad = mesh.quads[quadIndex];
      const TernarySurfaceNode& n00 = mesh.nodes[quad.nodes[0]];
      const TernarySurfaceNode& n10 = mesh.nodes[quad.nodes[1]];
      const TernarySurfaceNode& n11 = mesh.nodes[quad.nodes[2]];
      const TernarySurfaceNode& n01 = mesh.nodes[quad.nodes[3]];
      const std::array<ImVec2, 4> points{
          screenNode(n00), screenNode(n10), screenNode(n11), screenNode(n01)};
      std::array<float, 3> weights{};
      std::array<const TernarySurfaceNode*, 3> triangle{};
      if (barycentricPoint(mouse, points[0], points[1], points[2], weights)) {
        triangle = {&n00, &n10, &n11};
      } else if (barycentricPoint(mouse, points[0], points[2], points[3],
                                  weights)) {
        triangle = {&n00, &n11, &n01};
      } else {
        continue;
      }

      pick.quad = static_cast<int>(quadIndex);
      pick.value = 0.0;
      float ratioB = 0.0f;
      pick.fractionC = 0.0f;
      pick.fractions = {};
      for (size_t i = 0; i < triangle.size(); ++i) {
        const TernarySurfaceNode& node = *triangle[i];
        const float weight = weights[i];
        pick.value +=
            node.point.prediction.gramsPerMillilitre * weight;
        ratioB += cubePosition(node)[0] * weight;
        pick.fractionC += cubePosition(node)[2] * weight;
      }
      pick.ratioA = std::clamp(1.0f - ratioB, 0.0f, 1.0f);
      pick.fractionC = std::clamp(pick.fractionC, 0.0f, 1.0f);
      const double fractionAB = 1.0 - pick.fractionC;
      pick.fractions = {
          static_cast<double>(pick.ratioA) * fractionAB,
          static_cast<double>(1.0f - pick.ratioA) * fractionAB,
          static_cast<double>(pick.fractionC)};
    }
  }

  const auto currentPairFractionA = [&]() {
    const float pairTotal = sb.ratios[0] + sb.ratios[1];
    return pairTotal > 1e-6f ? sb.ratios[0] / pairTotal : 0.5f;
  };
  const auto currentFractionC = [&]() {
    const float total = sb.ratios[0] + sb.ratios[1] + sb.ratios[2];
    return total > 1e-6f ? sb.ratios[2] / total : 0.0f;
  };
  float currentA = currentPairFractionA();
  float currentC = currentFractionC();
  if (hovered && ImGui::GetIO().MouseWheel != 0.0f) {
    const float direction = ImGui::GetIO().MouseWheel > 0.0f ? 1.0f : -1.0f;
    currentC = std::clamp(
        currentC + direction / static_cast<float>(mesh.cQuads), 0.0f, 1.0f);
    setTernaryComposition(sb, currentA, currentC);
  }
  if (clicked && pick.quad >= 0) {
    currentA = pick.ratioA;
    currentC = pick.fractionC;
    setTernaryComposition(sb, currentA, currentC);
  }
  if (dragging) {
    const float ratioB = std::clamp(
        (mouse.x - frontBL.x - depth.x * currentC) /
            std::max(1.0f, plotW),
        0.0f, 1.0f);
    currentA = 1.0f - ratioB;
    setTernaryComposition(sb, currentA, currentC);
  }

  ImGui::PushID(&sb);
  const float cubeHover =
      widgets::hoverT(ImGui::GetID("surface_hover"), hovered);
  const float cellHover =
      widgets::hoverT(ImGui::GetID("surface_cell_hover"), pick.quad >= 0);
  ImGui::PopID();

  const float drawProgress = sweepDrawProgress(sb);
  const int ratioStride = std::max(1, mesh.ratioQuads / 8);
  const int cStride = std::max(1, mesh.cQuads / 6);
  for (size_t quadIndex = 0; quadIndex < mesh.quads.size(); ++quadIndex) {
    const TernarySurfaceQuad& quad = mesh.quads[quadIndex];
    const float path =
        (static_cast<float>(quad.ratioIndex) + 0.5f) /
        static_cast<float>(mesh.ratioQuads);
    const float reveal = markerReveal(drawProgress, path);
    if (reveal <= 0.0f) continue;

    const TernarySurfaceNode& n00 = mesh.nodes[quad.nodes[0]];
    const TernarySurfaceNode& n10 = mesh.nodes[quad.nodes[1]];
    const TernarySurfaceNode& n11 = mesh.nodes[quad.nodes[2]];
    const TernarySurfaceNode& n01 = mesh.nodes[quad.nodes[3]];
    const ImVec2 p00 = screenNode(n00);
    const ImVec2 p10 = screenNode(n10);
    const ImVec2 p11 = screenNode(n11);
    const ImVec2 p01 = screenNode(n01);
    const std::array<float, 4>& cachedColour =
        logY ? quad.colourLog : quad.colourLinear;
    const float shade = logY ? quad.shadeLog : quad.shadeLinear;
    const bool activeQuad = static_cast<int>(quadIndex) == pick.quad;
    const float hoverLift = activeQuad ? cellHover * 0.08f : 0.0f;
    const ImVec4 fill =
        shadedSurfaceColour(cachedColour, shade + hoverLift,
                            reveal * (0.91f + hoverLift));
    dl->AddQuadFilled(p00, p10, p11, p01, style::u32(fill));

    const ImU32 isoline =
        style::u32(shadedSurfaceColour(cachedColour, shade,
                                       reveal * (0.20f + cubeHover * 0.07f)));
    if (quad.cIndex % cStride == 0) {
      dl->AddLine(p00, p10, isoline, style::metrics().hairline);
    }
    if (quad.cIndex == mesh.cQuads - 1) {
      dl->AddLine(p01, p11, isoline, style::metrics().hairline);
    }
    if (quad.ratioIndex % ratioStride == 0) {
      dl->AddLine(p00, p01, isoline, style::metrics().hairline);
    }
    if (quad.ratioIndex == mesh.ratioQuads - 1) {
      dl->AddLine(p10, p11, isoline, style::metrics().hairline);
    }
    if (activeQuad && cellHover > 0.002f) {
      dl->AddQuad(p00, p10, p11, p01,
                  style::u32(style::col::Text, cellHover * 0.62f),
                  1.0f + cellHover);
    }
  }

  // Current blend: exact surface marker and drop line to the cube base.
  currentA = currentPairFractionA();
  currentC = currentFractionC();
  const float currentB = 1.0f - currentA;
  const double currentSurfaceValue = surfaceValueAt(currentB, currentC);
  const float currentHeight = normaliseValue(currentSurfaceValue);
  const ImVec2 currentBase = screenAt(currentB, currentC, 0.0f);
  const ImVec2 currentMarker =
      screenAt(currentB, currentC, currentHeight);
  const float currentReveal = markerReveal(drawProgress, currentB);
  if (currentReveal > 0.0f) {
    const float pulse = markerPulse(1.9f);
    dl->AddLine(currentBase, currentMarker,
                style::u32(style::col::Violet, 0.72f * currentReveal),
                1.35f);
    dl->AddCircleFilled(currentMarker, 4.4f,
                        style::u32(style::col::Violet, currentReveal));
    dl->AddCircle(currentMarker, 6.8f + pulse * 3.5f,
                  style::u32(style::col::Violet,
                             currentReveal * (0.66f - pulse * 0.44f)),
                  0, 1.35f);
  }

  // The sampled global maximum is distinct from the current blend.
  if (mesh.peakNode >= 0 &&
      static_cast<size_t>(mesh.peakNode) < mesh.nodes.size()) {
    const TernarySurfaceNode& peak =
        mesh.nodes[static_cast<size_t>(mesh.peakNode)];
    const std::array<float, 3>& peakCube = cubePosition(peak);
    const ImVec2 marker = screenNode(peak);
    const ImVec2 peakBase =
        screenAt(peakCube[0], peakCube[2], 0.0f);
    const float reveal = markerReveal(drawProgress, peakCube[0]);
    if (reveal > 0.0f) {
      const float pulse = markerPulse();
      dl->AddLine(peakBase, marker,
                  style::u32(style::col::Accent, 0.62f * reveal), 1.25f);
      dl->AddCircleFilled(marker, 4.8f,
                          style::u32(style::col::Accent, reveal));
      dl->AddCircle(marker, 7.0f + pulse * 3.5f,
                    style::u32(style::col::Accent,
                               reveal * (0.66f - pulse * 0.44f)),
                    0, 1.35f);
      const std::string rawLabel =
          "peak " +
          formatUnits(peak.point.prediction.gramsPerMillilitre, molarMass,
                      sb.units) +
          " | A " +
          std::to_string(static_cast<int>(
              std::lround(peak.point.fractions[0] * 100.0))) +
          "% / B " +
          std::to_string(static_cast<int>(
              std::lround(peak.point.fractions[1] * 100.0))) +
          "% / C " +
          std::to_string(static_cast<int>(
              std::lround(peak.point.fractions[2] * 100.0))) +
          "%";
      const float legendBoundary = backBR.x + gap * 0.35f;
      const float labelWidth =
          std::max(1.0f, legendBoundary - frontTL.x - gap * 2.0f);
      const std::string label = ellipsizeText(rawLabel, labelWidth);
      const ImVec2 labelSize = ImGui::CalcTextSize(label.c_str());
      const float maxLabelX =
          std::max(frontTL.x + gap * 0.4f,
                   legendBoundary - labelSize.x - gap * 0.4f);
      const float labelX =
          std::clamp(marker.x + gap * 0.55f, frontTL.x + gap * 0.4f,
                     maxLabelX);
      const float labelY =
          std::clamp(marker.y - labelSize.y - gap * 0.65f,
                     frontTL.y + gap * 0.35f,
                     std::max(frontTL.y + gap * 0.35f,
                              frontBL.y - labelSize.y - gap * 0.35f));
      dl->AddRectFilled(
          ImVec2(labelX - 3.0f, labelY - 2.0f),
          ImVec2(labelX + labelSize.x + 3.0f,
                 labelY + labelSize.y + 2.0f),
          style::u32(style::col::BgRaised, 0.94f * reveal),
          style::metrics().radiusSm);
      dl->AddText(ImVec2(labelX, labelY),
                  style::u32(style::col::Accent, reveal), label.c_str());
    }
  }

  // Cube frame remains above the filled surface.
  const ImU32 frontEdge =
      style::u32(style::col::BorderStrong, 0.95f);
  const ImU32 backEdge =
      style::u32(style::col::BorderStrong, 0.58f);
  dl->AddQuad(frontTL, frontTR, frontBR, frontBL, frontEdge, 1.45f);
  dl->AddQuad(backTL, backTR, backBR, backBL, backEdge, 1.1f);
  dl->AddLine(frontTL, backTL, backEdge, 1.1f);
  dl->AddLine(frontTR, backTR, backEdge, 1.1f);
  dl->AddLine(frontBR, backBR, backEdge, 1.1f);
  dl->AddLine(frontBL, backBL, frontEdge, 1.35f);

  // Solubility axis: bare values with one unit caption.
  for (int i = 0; i <= 4; ++i) {
    const float t = static_cast<float>(i) / 4.0f;
    const std::string label = ellipsizeText(
        formatSolubility(
            toDisplayUnits(rawAtHeight(1.0 - t), molarMass, sb.units)),
        std::max(1.0f, leftPad - gap * 0.75f));
    const ImVec2 size = ImGui::CalcTextSize(label.c_str());
    const float y = frontTL.y + plotH * t;
    dl->AddLine(ImVec2(frontTL.x - 3.0f, y),
                ImVec2(frontTL.x + 3.0f, y),
                style::u32(style::col::BorderStrong));
    dl->AddText(ImVec2(frontTL.x - size.x - gap * 0.65f,
                       y - size.y * 0.5f),
                style::u32(style::col::TextDim), label.c_str());
  }
  const std::string yTitle =
      ellipsizeText("Solubility (" + unitCaption + ")",
                    std::max(1.0f, leftPad - gap * 0.5f));
  dl->AddText(ImVec2(origin.x + gap * 0.25f,
                     origin.y + gap * 0.25f),
              style::u32(style::col::TextFaint), yTitle.c_str());

  // Front-bottom A:B axis, pure A at left and pure B at right.
  const int xSubdivisions = plotW < fontSize * 16.0f ? 2 : 4;
  for (int i = 0; i <= xSubdivisions; ++i) {
    const float t =
        static_cast<float>(i) / static_cast<float>(xSubdivisions);
    const float x = frontBL.x + plotW * t;
    char label[16];
    std::snprintf(label, sizeof(label), "%.0f:%.0f",
                  (1.0f - t) * 100.0f, t * 100.0f);
    const ImVec2 size = ImGui::CalcTextSize(label);
    dl->AddLine(ImVec2(x, frontBL.y - 3.0f),
                ImVec2(x, frontBL.y + 3.0f),
                style::u32(style::col::BorderStrong));
    dl->AddText(
        ImVec2(std::clamp(x - size.x * 0.5f, origin.x,
                          std::max(origin.x, canvasMax.x - size.x)),
               frontBL.y + gap * 0.45f),
        style::u32(style::col::TextDim), label);
  }
  const std::string pairTitle =
      ellipsizeText("A " + a.name + "  ->  B " + b.name, plotW);
  const ImVec2 pairTitleSize = ImGui::CalcTextSize(pairTitle.c_str());
  dl->AddText(ImVec2(frontBL.x + (plotW - pairTitleSize.x) * 0.5f,
                     frontBL.y + fontSize + gap * 0.65f),
              style::u32(style::col::TextFaint), pairTitle.c_str());

  // Receding bottom edge: C fraction and solvent identity.
  for (int i = 1; i <= 2; ++i) {
    const float fraction = static_cast<float>(i) * 0.5f;
    const ImVec2 tick = offset(frontBL, fraction);
    char label[8];
    std::snprintf(label, sizeof(label), "%.0f%%", fraction * 100.0f);
    dl->AddLine(ImVec2(tick.x - 2.0f, tick.y - 2.0f),
                ImVec2(tick.x + 2.0f, tick.y + 2.0f),
                style::u32(slotColor(2), 0.72f));
    const ImVec2 labelSize = ImGui::CalcTextSize(label);
    dl->AddText(
        ImVec2(std::clamp(tick.x - labelSize.x - gap * 0.42f, origin.x,
                          std::max(origin.x, canvasMax.x - labelSize.x)),
               std::clamp(tick.y + 2.0f, origin.y,
                          std::max(origin.y,
                                   canvasMax.y - labelSize.y))),
        style::u32(style::col::TextDim), label);
  }
  const float cTitleWidth =
      std::max(1.0f, backBL.x - origin.x - gap);
  const std::string cTitle =
      ellipsizeText("C " + c.name, cTitleWidth);
  const ImVec2 cTitleSize = ImGui::CalcTextSize(cTitle.c_str());
  dl->AddText(
      ImVec2(std::clamp(backBL.x - cTitleSize.x - gap * 0.45f, origin.x,
                        std::max(origin.x,
                                 canvasMax.x - cTitleSize.x)),
             std::clamp(backBL.y + fontSize + gap * 0.45f, origin.y,
                        std::max(origin.y,
                                 canvasMax.y - cTitleSize.y))),
      style::u32(slotColor(2)), cTitle.c_str());

  // Compact value legend beside the cube. It always uses the unshaded ramp.
  const float legendX = backBR.x + gap * 0.55f;
  const float stripWidth =
      std::clamp(fontSize * 0.55f, 5.0f, 10.0f);
  const float legendTop = frontTL.y;
  const float legendBottom = frontBL.y;
  const float legendHeight = std::max(1.0f, legendBottom - legendTop);
  constexpr int kLegendSegments = 32;
  for (int i = 0; i < kLegendSegments; ++i) {
    const float t0 =
        static_cast<float>(i) / static_cast<float>(kLegendSegments);
    const float t1 =
        static_cast<float>(i + 1) / static_cast<float>(kLegendSegments);
    const float y0 = legendBottom - t0 * legendHeight;
    const float y1 = legendBottom - t1 * legendHeight;
    dl->AddRectFilled(ImVec2(legendX, y1),
                      ImVec2(legendX + stripWidth, y0),
                      style::u32(surfaceColour((t0 + t1) * 0.5f)));
  }
  dl->AddRect(ImVec2(legendX, legendTop),
              ImVec2(legendX + stripWidth, legendBottom),
              style::u32(style::col::BorderStrong), 0.0f, 0,
              style::metrics().hairline);
  const float legendLabelX = legendX + stripWidth + gap * 0.45f;
  const float legendLabelWidth =
      std::max(1.0f, canvasMax.x - legendLabelX - gap * 0.25f);
  const std::string legendTitle =
      ellipsizeText("VALUE (" + unitCaption + ")", legendLabelWidth);
  dl->AddText(ImVec2(legendX, std::max(origin.y, legendTop - fontSize)),
              style::u32(style::col::TextFaint), legendTitle.c_str());
  for (int i = 0; i < 3; ++i) {
    const float t = static_cast<float>(i) * 0.5f;
    const std::string label = ellipsizeText(
        formatUnits(rawAtHeight(t), molarMass, sb.units),
        legendLabelWidth);
    const ImVec2 labelSize = ImGui::CalcTextSize(label.c_str());
    const float y =
        legendBottom - t * legendHeight - labelSize.y * 0.5f;
    dl->AddLine(ImVec2(legendX + stripWidth, y + labelSize.y * 0.5f),
                ImVec2(legendX + stripWidth + 3.0f,
                       y + labelSize.y * 0.5f),
                style::u32(style::col::BorderStrong));
    dl->AddText(
        ImVec2(legendLabelX,
               std::clamp(y, origin.y,
                          std::max(origin.y,
                                   canvasMax.y - labelSize.y))),
        style::u32(style::col::TextDim), label.c_str());
  }
  const float currentLegendT = normaliseValue(currentSurfaceValue);
  const float currentLegendY =
      legendBottom - currentLegendT * legendHeight;
  dl->AddLine(ImVec2(legendX - 2.0f, currentLegendY),
              ImVec2(legendX + stripWidth + 3.0f, currentLegendY),
              style::u32(style::col::Violet), 1.6f);
  dl->AddTriangleFilled(
      ImVec2(legendX - 2.0f, currentLegendY),
      ImVec2(legendX - 6.0f, currentLegendY - 3.5f),
      ImVec2(legendX - 6.0f, currentLegendY + 3.5f),
      style::u32(style::col::Violet));

  if (pick.quad >= 0) {
    const std::string tip =
        a.name + " " +
        std::to_string(static_cast<int>(
            std::lround(pick.fractions[0] * 100.0))) +
        "%  /  " + b.name + " " +
        std::to_string(static_cast<int>(
            std::lround(pick.fractions[1] * 100.0))) +
        "%  /  " + c.name + " " +
        std::to_string(static_cast<int>(
            std::lround(pick.fractions[2] * 100.0))) +
        "%   ->   " + formatUnits(pick.value, molarMass, sb.units);
    ImGui::SetTooltip("%s", tip.c_str());
  }

  dl->AddRect(origin, canvasMax,
              style::mix(style::col::Border, style::col::Teal, cubeHover,
                         0.45f + cubeHover * 0.35f),
              style::metrics().radiusMd, 0,
              style::metrics().hairline * (1.0f + cubeHover * 0.5f));
  dl->PopClipRect();
}

const sol::SweepPoint* activeSweepPeak(const SolubilityState& sb) {
  if (sb.sweepPeakIndex < 0 ||
      static_cast<size_t>(sb.sweepPeakIndex) >= sb.sweep.size()) {
    return nullptr;
  }
  return &sb.sweep[static_cast<size_t>(sb.sweepPeakIndex)];
}

void drawPeakEvidenceCard(AppState& st, const std::vector<const sol::Solvent*>& chosen,
                          const std::vector<std::string>& miscibilityWarnings,
                          bool solventsOk, ImVec2 size) {
  SolubilityState& sb = st.solubility;
  if (!widgets::beginCard("##peak_evidence_card", size, style::col::BgSurface)) return;

  widgets::sectionHeader("Peak readout", style::col::Teal);
  const sol::SweepPoint* peak = activeSweepPeak(sb);
  if (peak && sb.solventCount >= 2) {
    widgets::badge("BEST", style::col::Accent);
    const std::string peakValue =
        formatUnits(peak->prediction.gramsPerMillilitre, sb.solute.molarMass, sb.units);
    const bool mono = style::pushFont(style::fonts::mono());
    ImGui::TextWrapped("%s", peakValue.c_str());
    style::popFont(mono);

    std::string composition;
    for (int i = 0; i < sb.solventCount && i < static_cast<int>(chosen.size()); ++i) {
      if (!chosen[static_cast<size_t>(i)]) continue;
      if (!composition.empty()) composition += "  /  ";
      composition += chosen[static_cast<size_t>(i)]->name + " ";
      composition +=
          std::to_string(static_cast<int>(std::lround(peak->fractions[static_cast<size_t>(i)] *
                                                      100.0)));
      composition += "%";
    }
    ImGui::TextWrapped("%s", composition.c_str());
    const bool applyPeak =
        widgets::ghostButton("Apply peak ratio",
                             ImVec2(ImGui::GetContentRegionAvail().x, 0.0f));
    drawInteractiveBorder("##peak_apply_feedback", style::col::Teal);
    if (applyPeak) {
      for (int i = 0; i < sb.solventCount; ++i) {
        sb.ratios[static_cast<size_t>(i)] =
            static_cast<float>(peak->fractions[static_cast<size_t>(i)]);
      }
    }
  } else if (sb.solventCount < 2) {
    textDisabledWrapped("Add a second solvent to search for a co-solvency peak.");
  } else {
    textDisabledWrapped("Complete the blend to calculate a peak.");
  }

  widgets::sectionHeader("Literature anchor", style::col::Violet);
  if (sb.prediction.anchored) {
    widgets::badge("MEASURED", style::col::Success);
    ImGui::Spacing();
    ImGui::TextWrapped("%s", sb.prediction.anchorNote.c_str());
  } else {
    widgets::badge("MODEL", style::col::Violet);
    textDisabledWrapped("No matching measured anchor.");
  }

  widgets::sectionHeader("Diagnostics", style::col::Danger);
  bool hasWarning = false;
  if (sb.prediction.outsideSphere) {
    ImGui::PushStyleColor(ImGuiCol_Text, style::col::Danger);
    ImGui::TextWrapped("Outside the Hansen sphere (RED %.2f > 1). The result is extrapolated.",
                       sb.prediction.relativeEnergyDifference);
    ImGui::PopStyleColor();
    hasWarning = true;
  }
  if (!sb.prediction.converged) {
    ImGui::PushStyleColor(ImGuiCol_Text, style::col::Danger);
    ImGui::TextWrapped("Saturation composition did not converge. Treat this result as unreliable.");
    ImGui::PopStyleColor();
    hasWarning = true;
  }
  for (const std::string& warning : miscibilityWarnings) {
    ImGui::PushStyleColor(ImGuiCol_Text, style::col::Danger);
    ImGui::TextWrapped("%s", warning.c_str());
    ImGui::PopStyleColor();
    hasWarning = true;
  }
  if (!hasWarning) {
    ImGui::PushStyleColor(ImGuiCol_Text, style::col::Success);
    ImGui::TextWrapped("No model or miscibility warnings.");
    ImGui::PopStyleColor();
  }

  widgets::sectionHeader("Workflow", style::col::Accent);
  const bool canSend = sb.soluteValid && solventsOk && chosen.size() >= 2 && chosen[0] &&
                       chosen[1];
  if (!canSend) ImGui::BeginDisabled();
  const bool sendToExtraction =
      widgets::primaryButton("Send to Extraction",
                             ImVec2(ImGui::GetContentRegionAvail().x, 0.0f));
  const bool extractionHovered = ImGui::IsItemHovered();
  drawInteractiveBorder("##extraction_handoff_feedback", style::col::AccentHover);
  if (sendToExtraction) {
    const double r0 = std::max(0.0, static_cast<double>(sb.ratios[0]));
    const double r1 = std::max(0.0, static_cast<double>(sb.ratios[1]));
    const double sum = std::max(r0 + r1, 1e-9);
    ExtractionImport& imp = sb.extractionImport;
    imp.pending = true;
    imp.solventIdA = chosen[0]->id;
    imp.solventIdB = chosen[1]->id;
    imp.volumeMlA = 100.0 * r0 / sum;
    imp.volumeMlB = 100.0 * r1 / sum;
    imp.soluteMassMg = 100.0;
    st.tab = MainTab::Extraction;
    st.tabChangeRequested = true;
    sb.statusMessage = "Blend sent to the Extraction Calculator";
  }
  if (extractionHovered)
    ImGui::SetTooltip("%s", "Load this blend into the separatory-funnel simulation");
  if (!canSend) ImGui::EndDisabled();

  widgets::endCard();
}

void drawScreeningCard(SolubilityState& sb, ImVec2 size) {
  if (!widgets::beginCard("##screening_card", size, style::col::BgSurface)) return;
  widgets::sectionHeader("Solvent screen", style::col::Accent);
  drawScreeningTable(sb);
  widgets::endCard();
}

}  // namespace

void drawSolubilitySuite(AppState& st) {
  SolubilityState& sb = st.solubility;
  const style::Metrics& metrics = style::metrics();
  const ImVec2 panelAvail = ImGui::GetContentRegionAvail();
  const float preferredRail = ImGui::GetFontSize() * 20.0f;
  const float minimumRail = ImGui::GetFontSize() * 13.0f;
  const float railW =
      std::max(1.0f, std::min({preferredRail, std::max(minimumRail, panelAvail.x * 0.34f),
                              panelAvail.x * 0.46f}));

  std::vector<sol::Component> components;
  std::vector<const sol::Solvent*> chosen;
  std::vector<std::string> miscibilityWarnings;
  bool solventsOk = true;

  // ----------------------------------------------------- left control rail
  ImGui::PushStyleVar(ImGuiStyleVar_ChildRounding, metrics.radiusLg);
  ImGui::PushStyleColor(ImGuiCol_ChildBg, style::col::BgPanel);
  ImGui::BeginChild("##suite_rail", ImVec2(railW, 0.0f), ImGuiChildFlags_Borders);
  ImGui::PopStyleColor();
  ImGui::PopStyleVar();

  widgets::sectionHeader("Solute");
  if (widgets::beginCard("##solute_summary_card", ImVec2(0.0f, 0.0f),
                         style::col::BgSurface)) {
    drawSoluteCard(st);
    widgets::endCard();
  }

  widgets::sectionHeader("Solvent blend", style::col::Teal);
  if (widgets::beginCard("##solvent_controls_card", ImVec2(0.0f, 0.0f),
                         style::col::BgSurface)) {
    try {
      chosen = drawSolventMixer(sb, components, miscibilityWarnings);
    } catch (const sol::SolError& err) {
      solventsOk = false;
      ImGui::PushStyleColor(ImGuiCol_Text, style::col::Danger);
      ImGui::TextWrapped("%s", err.what());
      ImGui::PopStyleColor();
      sb.statusMessage = std::string("Solvent database error: ") + err.what();
    }
    widgets::endCard();
  }

  widgets::sectionHeader("Conditions", style::col::Violet);
  if (widgets::beginCard("##conditions_card", ImVec2(0.0f, 0.0f),
                         style::col::BgSurface)) {
    const bool temperatureLabelInline = prepareTrailingLabel("Temperature (C)");
    ImGui::SliderFloat(temperatureLabelInline ? "Temperature (C)" : "##temperature",
                       &sb.temperatureC, -20.0f, 150.0f, "%.1f");

    // pH only exists as a control for a solute that can ionise. Showing a pH
    // slider next to naphthalene would be theatre.
    if (sb.soluteValid && (sb.solute.ionization.ionClass == sol::IonClass::Base ||
                           sb.solute.ionization.ionClass == sol::IonClass::Acid)) {
      ImGui::Spacing();
      checkboxWithResponsiveLabel("pH set by the solute itself", "Self-buffered##ph_auto",
                                  &sb.pHAuto);
      if (sb.pHAuto) {
        textDisabledWrapped("Saturated solution sets pH %.1f (pKa ~%.1f, %s).",
                            sb.prediction.ionicPath ? sb.prediction.pH : 7.0,
                            sb.solute.ionization.pKa,
                            sb.solute.ionization.siteLabel.empty()
                                ? "estimated"
                                : sb.solute.ionization.siteLabel.c_str());
      } else {
        const bool pHLabelInline = prepareTrailingLabel("pH");
        ImGui::SliderFloat(pHLabelInline ? "pH" : "##ph", &sb.pH, 0.0f, 14.0f, "%.1f");
      }
    }

    // Common-ion controls remain conditional on a recognised 1:1 salt.
    if (sb.soluteValid && !sb.solute.canonicalSmiles.empty() &&
        sol::findSalt(sb.solute.canonicalSmiles) != nullptr) {
      ImGui::Spacing();
      widgets::sectionHeader("Common ion effect");
      const sol::Salt* salt = sol::findSalt(sb.solute.canonicalSmiles);
      textDisabledWrapped("%s | Ksp %.3g | %s / %s", salt->name.c_str(), salt->ksp25,
                          salt->cation.c_str(), salt->anion.c_str());
      checkboxWithResponsiveLabel("Background electrolyte",
                                  "Enable##background_electrolyte", &sb.backgroundEnabled);
      if (sb.backgroundEnabled) {
        const std::vector<sol::Electrolyte>& all = sol::electrolytes();
        const sol::Electrolyte* current = activeBackground(sb);
        ImGui::TextUnformatted("Electrolyte");
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::BeginCombo("##bg_salt", current ? current->name.c_str() : "Select...")) {
          for (size_t i = 0; i < all.size(); ++i) {
            const bool selected = static_cast<int>(i) == sb.backgroundElectrolyte;
            if (ImGui::Selectable(all[i].name.c_str(), selected))
              sb.backgroundElectrolyte = static_cast<int>(i);
            if (selected) ImGui::SetItemDefaultFocus();
          }
          ImGui::EndCombo();
        }
        const bool concentrationLabelInline = prepareTrailingLabel("Concentration");
        ImGui::SliderFloat(concentrationLabelInline ? "Concentration" : "##background_molarity",
                           &sb.backgroundMolarity, 0.0f, 3.0f, "%.2f mol/L");
        if (current && (current->cation == salt->cation || current->anion == salt->anion)) {
          ImGui::PushStyleColor(ImGuiCol_Text, style::col::Accent);
          ImGui::TextWrapped("Common ion present -- solubility is depressed.");
          ImGui::PopStyleColor();
        } else if (current) {
          textDisabledWrapped("No common ion -- ionic-strength effect only.");
        }
      }
    }
    widgets::endCard();
  }

  if (!sb.statusMessage.empty()) {
    widgets::sectionHeader("Status", style::col::Success);
    if (widgets::beginCard("##suite_status_card", ImVec2(0.0f, 0.0f),
                           style::col::BgSurface)) {
      ImGui::TextWrapped("%s", sb.statusMessage.c_str());
      widgets::endCard();
    }
  }
  ImGui::EndChild();

  ImGui::SameLine(0.0f, metrics.gap);

  // ------------------------------------------------------ right workspace
  ImGui::BeginChild("##suite_results", ImVec2(0.0f, 0.0f), ImGuiChildFlags_None);
  const bool canPredict = sb.soluteValid && solventsOk && !components.empty();
  if (canPredict) {
    sb.prediction = sol::predict(sb.solute, components, static_cast<double>(sb.temperatureC),
                                 activeBackground(sb),
                                 static_cast<double>(sb.backgroundMolarity),
                                 requestedPH(sb));
    sb.statusMessage = "Predicted " + formatUnits(sb.prediction.gramsPerMillilitre,
                                                  sb.solute.molarMass, sb.units) +
                       ".";
  } else {
    sb.prediction = sol::Prediction{};
  }
  const float predictionFlash = predictionFlashAlpha(sb, canPredict);

  const bool haveAllSlots =
      solventsOk &&
      std::all_of(chosen.begin(), chosen.end(), [](const sol::Solvent* s) { return s != nullptr; });
  if (sb.solventCount >= 2 && sb.soluteValid && haveAllSlots) {
    recomputeSweepIfNeeded(sb, chosen);
  } else {
    sb.sweep.clear();
    sb.ternarySurface = {};
    sb.sweepSignature.clear();
  }

  const ImVec2 workspaceAvail = ImGui::GetContentRegionAvail();
  const float graphH = std::max(ImGui::GetFontSize() * 18.0f, workspaceAvail.y * 0.54f);
  if (widgets::beginCard("##composition_workspace", ImVec2(0.0f, graphH),
                         style::col::BgRaised,
                         ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    const bool headingFont = style::pushFont(style::fonts::semibold());
    ImGui::TextUnformatted("Composition response");
    style::popFont(headingFont);
    textDisabledWrapped(
        sb.solventCount == 3
            ? "Drag horizontally to set A:B. Scroll to set C, or click the surface to set both."
            : "Drag directly on the graph to set the working blend.");

    constexpr ImGuiTableFlags controlFlags =
        ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_NoSavedSettings;
    const bool stackControls =
        ImGui::GetContentRegionAvail().x < ImGui::GetFontSize() * 30.0f;
    const int controlColumns = stackControls ? 1 : 3;
    if (ImGui::BeginTable("##composition_controls", controlColumns, controlFlags)) {
      if (stackControls) {
        ImGui::TableSetupColumn("##stacked_controls", ImGuiTableColumnFlags_WidthStretch);
      } else {
        ImGui::TableSetupColumn("##grid_col", ImGuiTableColumnFlags_WidthStretch, 0.38f);
        ImGui::TableSetupColumn("##units_col", ImGuiTableColumnFlags_WidthStretch, 0.28f);
        ImGui::TableSetupColumn("##scale_col", ImGuiTableColumnFlags_WidthStretch, 0.34f);
      }
      ImGui::TableNextRow();

      ImGui::TableNextColumn();
      ImGui::TextDisabled("GRID RESOLUTION");
      if (sb.solventCount >= 2) {
        ImGui::SetNextItemWidth(-FLT_MIN);
        ImGui::SliderInt("##grid_resolution", &sb.sweepSteps, 2, 64);
      } else {
        ImGui::TextDisabled("Pure solvent");
      }

      ImGui::TableNextColumn();
      ImGui::TextDisabled("DISPLAY UNITS");
      ImGui::SetNextItemWidth(-FLT_MIN);
      ImGui::Combo("##result_units", &sb.units, kUnitLabels.data(), 4);

      ImGui::TableNextColumn();
      ImGui::TextDisabled("Y AXIS");
      ImGui::Checkbox("Log scale", &sb.logScale);
      if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("%s", "Plot solubility on a log10 scale");
      }
      ImGui::EndTable();
    }

    if (!sb.soluteValid) {
      ImGui::Spacing();
      widgets::badge("START HERE", style::col::Accent);
      ImGui::Spacing();
      const bool emptyFont = style::pushFont(style::fonts::semibold());
      ImGui::TextUnformatted("Load a solute to begin");
      style::popFont(emptyFont);
      ImGui::TextWrapped("Draw a molecular structure or switch the Solute card to SMILES. "
                         "Then choose two solvents to unlock the interactive ratio graph, or "
                         "three solvents for a colour-mapped composition surface.");
    } else if (!canPredict) {
      ImGui::Spacing();
      widgets::badge("NEEDS SOLVENT", style::col::Teal);
      ImGui::Spacing();
      ImGui::TextWrapped("Choose at least one solvent in the control rail to calculate a "
                         "prediction.");
    } else if (sb.solventCount < 2) {
      ImGui::Spacing();
      widgets::badge("PURE SOLVENT", style::col::Teal);
      ImGui::Spacing();
      ImGui::TextWrapped("The pure-solvent prediction is ready below. Add a second solvent to "
                         "explore a draggable composition curve.");
    } else if (!haveAllSlots) {
      ImGui::Spacing();
      textDisabledWrapped("Choose every active solvent slot to calculate the composition map.");
    } else {
      ImVec2 plotSize = ImGui::GetContentRegionAvail();
      plotSize.x = std::max(1.0f, plotSize.x);
      plotSize.y = std::max(1.0f, plotSize.y);
      if (sb.solventCount == 2) {
        drawBinarySweepPlot(sb, *chosen[0], *chosen[1], plotSize);
      } else {
        drawTernarySweepPlot(sb, *chosen[0], *chosen[1], *chosen[2], plotSize);
      }
    }
    widgets::endCard();
  }

  ImGui::Spacing();
  if (widgets::beginCard("##prediction_summary_card", ImVec2(0.0f, 0.0f),
                         style::col::BgSurface)) {
    widgets::sectionHeader("Current prediction", style::col::Accent);
    if (canPredict) {
      drawResultHero(sb, predictionFlash);
    } else {
      textDisabledWrapped("Prediction readouts will appear when the solute and solvent blend are "
                          "complete.");
    }
    widgets::endCard();
  }

  ImGui::Spacing();
  if (widgets::beginCard("##theory_card", ImVec2(0.0f, 0.0f), style::col::BgSurface)) {
    widgets::sectionHeader("Theory readout", style::col::Violet);
    if (canPredict) {
      drawTheoryReadout(sb, chosen);
    } else {
      textDisabledWrapped("Kirkwood-Buff integrals, the Yalkowsky ideal term and the "
                          "Hildebrand parameters appear with a complete blend.");
    }
    widgets::endCard();
  }

  ImGui::Spacing();
  const float bottomW = ImGui::GetContentRegionAvail().x;
  const float bottomGap = metrics.gap;
  const bool twoColumns = bottomW >= ImGui::GetFontSize() * 34.0f;
  const float bottomH = ImGui::GetFontSize() * 21.0f;
  const float insightsW = twoColumns ? (bottomW - bottomGap) * 0.39f : bottomW;
  const float screeningW = twoColumns ? bottomW - bottomGap - insightsW : bottomW;
  drawPeakEvidenceCard(st, chosen, miscibilityWarnings, solventsOk,
                       ImVec2(insightsW, bottomH));
  if (twoColumns) {
    ImGui::SameLine(0.0f, bottomGap);
  } else {
    ImGui::Spacing();
  }
  drawScreeningCard(sb, ImVec2(screeningW, bottomH));

  ImGui::EndChild();
}

}  // namespace chemcad::ui

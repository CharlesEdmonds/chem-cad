// Solvent Selection: operation-first solvent ranking with explicit evidence.
// This is the Solubility Suite's `Select` mode, not a panel of its own: it
// renders inside the suite's already-open window and hands its winning blend
// straight into `Predict`. Owns no GL resources.

#include "imgui.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "chem/bridge.hpp"
#include "core/model.hpp"
#include "naming/naming.hpp"
#include "sol/selection.hpp"
#include "sol/solvent.hpp"
#include "ui/app_state.hpp"
#include "ui/icons.hpp"
#include "ui/selection_state.hpp"
#include "ui/solubility_state.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"
#include "ui/widgets.hpp"

namespace chemcad::ui {
namespace {

constexpr std::array<sol::OperationKind, 6> kOperationKinds = {
    sol::OperationKind::LiquidLiquidExtraction,
    sol::OperationKind::Recrystallisation,
    sol::OperationKind::Trituration,
    sol::OperationKind::AntiSolventPrecipitation,
    sol::OperationKind::ChromatographyMobilePhase,
    sol::OperationKind::ReactionMedium,
};

constexpr std::array<const char*, 5> kSortLabels = {
    "Score", "Selectivity", "Recovery", "Greenness", "Boiling point"};

struct RankingResponse {
  std::vector<sol::SolventCandidate> candidates;
  std::string error;
};

struct ResultAnimation {
  uint64_t revision = std::numeric_limits<uint64_t>::max();
  double startedAt = 0.0;
  std::unordered_map<std::string, int> previousRanks;
  std::unordered_map<std::string, int> currentRanks;
};


std::string ellipsizeText(const std::string& text, float maximumWidth) {
  if (maximumWidth <= 0.0f) return {};
  if (ImGui::CalcTextSize(text.c_str()).x <= maximumWidth) return text;
  constexpr const char* kEllipsis = "...";
  const float ellipsisWidth = ImGui::CalcTextSize(kEllipsis).x;
  if (ellipsisWidth >= maximumWidth) return {};
  size_t end = text.size();
  while (end > 0 &&
         ImGui::CalcTextSize(text.data(), text.data() + end).x + ellipsisWidth > maximumWidth) {
    --end;
  }
  return text.substr(0, end) + kEllipsis;
}


float easeOutCubic(float value) {
  const float t = std::clamp(value, 0.0f, 1.0f);
  const float remaining = 1.0f - t;
  return 1.0f - remaining * remaining * remaining;
}

ImVec4 scoreColour(double score) {
  const float value = std::clamp(static_cast<float>(score), 0.0f, 1.0f);
  const ImVec4& a = value < 0.5f ? style::col::Danger : style::col::Accent;
  const ImVec4& b = value < 0.5f ? style::col::Accent : style::col::Success;
  const float t = value < 0.5f ? value * 2.0f : (value - 0.5f) * 2.0f;
  return ImVec4(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t,
                a.z + (b.z - a.z) * t, 1.0f);
}

ImVec4 chem21Colour(const std::string& value) {
  if (value == "recommended") return style::col::Success;
  if (value == "problematic") return style::col::Accent;
  if (value == "hazardous" || value == "highly hazardous") return style::col::Danger;
  return style::col::TextDim;
}

std::string formatSolubility(double value) {
  char buffer[64];
  const double magnitude = std::fabs(value);
  if (!std::isfinite(value)) return "infinite";
  if (value == 0.0) return "0 g/mL";
  if (magnitude < 1e-3 || magnitude >= 1e4) {
    std::snprintf(buffer, sizeof(buffer), "%.2e g/mL", value);
  } else {
    std::snprintf(buffer, sizeof(buffer), "%.4g g/mL", value);
  }
  return buffer;
}

std::string formatRatio(double value) {
  char buffer[48];
  if (!std::isfinite(value)) return "infinite";
  if (std::fabs(value) >= 1e4) {
    std::snprintf(buffer, sizeof(buffer), "%.2e x", value);
  } else {
    std::snprintf(buffer, sizeof(buffer), "%.3g x", value);
  }
  return buffer;
}

std::string formatPercent(double value) {
  char buffer[32];
  std::snprintf(buffer, sizeof(buffer), "%.1f%%", std::clamp(value, 0.0, 1.0) * 100.0);
  return buffer;
}

std::string formatBoilingPoint(const sol::SolventCandidate& candidate) {
  if (!candidate.solvent || candidate.solvent->boilingPoint <= 0.0) return "not rated";
  char buffer[64];
  if (candidate.partner && candidate.partner->boilingPoint > 0.0) {
    std::snprintf(buffer, sizeof(buffer), "%.1f / %.1f C", candidate.solvent->boilingPoint,
                  candidate.partner->boilingPoint);
  } else {
    std::snprintf(buffer, sizeof(buffer), "%.1f C", candidate.solvent->boilingPoint);
  }
  return buffer;
}

std::string candidateKey(const sol::SolventCandidate& candidate) {
  std::string key = candidate.solvent ? candidate.solvent->id : "missing";
  key += '|';
  if (candidate.partner) key += candidate.partner->id;
  key += '|';
  key += std::to_string(static_cast<int>(std::lround(candidate.partnerFraction * 10000.0)));
  return key;
}

std::string candidateTitle(const sol::SolventCandidate& candidate) {
  if (!candidate.solvent) return "Unavailable solvent";
  std::string title = candidate.solvent->name;
  if (candidate.partner) {
    title += " + ";
    title += candidate.partner->name;
    title += " ";
    title += std::to_string(static_cast<int>(std::lround(candidate.partnerFraction * 100.0)));
    title += "%";
  }
  return title;
}

bool hasTarget(const SelectionState& state) {
  return std::any_of(state.operation.species.begin(), state.operation.species.end(),
                     [](const sol::SpeciesRole& species) { return species.keep; });
}

double targetMassMg(const SelectionState& state) {
  for (const sol::SpeciesRole& species : state.operation.species) {
    if (species.keep && species.amountMg > 0.0) return species.amountMg;
  }
  return 100.0;
}

core::Molecule combinedSketch(const AppState& state) {
  core::Molecule combined;
  for (const core::Molecule& fragment : state.doc.molecules) {
    if (fragment.empty()) continue;
    std::unordered_map<core::AtomId, core::AtomId> ids;
    ids.reserve(fragment.atomCount());
    for (const core::Atom& atom : fragment.atoms()) ids.emplace(atom.id, combined.addAtom(atom));
    for (const core::Bond& bond : fragment.bonds()) {
      const core::BondId id = combined.addBond(ids.at(bond.a), ids.at(bond.b), bond.order);
      if (core::Bond* added = combined.bond(id)) added->stereo = bond.stereo;
    }
  }
  return combined;
}

void addSpecies(SelectionState& state, const core::Molecule& molecule,
                const std::string& requestedLabel) {
  if (molecule.empty()) throw sol::SolError("The structure is empty.");
  sol::SpeciesRole species;
  species.solute = sol::describeSolute(molecule);
  species.label = requestedLabel.empty() ? species.solute.name : requestedLabel;
  species.keep = true;
  species.weight = 1.0;

  SelectionSpeciesPresentation presentation;
  const chem::Properties properties = chem::computeProperties(molecule);
  presentation.formula = properties.formula;
  presentation.smiles = chem::toSmiles(molecule);
  state.operation.species.push_back(std::move(species));
  state.speciesPresentation.push_back(std::move(presentation));
  state.inputError.clear();
}

void addSketchSpecies(AppState& state) {
  try {
    const core::Molecule molecule = combinedSketch(state);
    std::string label = state.props.name;
    if (label.empty()) label = "Sketch species";
    addSpecies(state.selection, molecule, label);
    state.selection.statusMessage = "Added the current sketch";
  } catch (const std::exception& error) {
    state.selection.inputError = error.what();
  }
}

void addSmilesSpecies(AppState& state) {
  SelectionState& selection = state.selection;
  if (selection.smilesInput.empty()) {
    selection.inputError = "Enter a SMILES string first.";
    return;
  }
  try {
    const core::Molecule molecule = chem::fromSmiles(selection.smilesInput);
    addSpecies(selection, molecule, selection.smilesInput);
    selection.smilesInput.clear();
    selection.statusMessage = "Added species from SMILES";
  } catch (const std::exception& error) {
    selection.inputError = error.what();
  }
}

void addNameSpecies(AppState& state) {
  SelectionState& selection = state.selection;
  if (selection.nameInput.empty() || selection.nameLookupRunning) {
    if (selection.nameInput.empty()) selection.inputError = "Enter a chemical name first.";
    return;
  }
  const std::string name = selection.nameInput;
  const uint64_t request = ++selection.nameRequest;
  selection.nameLookupRunning = true;
  selection.inputError.clear();
  state.tasks.run<naming::Result>(
      [name] { return naming::nameToSmiles(name); },
      [&state, name, request](naming::Result result) {
        SelectionState& current = state.selection;
        if (request != current.nameRequest) return;
        current.nameLookupRunning = false;
        if (!result.ok) {
          current.inputError = result.error.empty() ? "Could not resolve that name." : result.error;
          return;
        }
        try {
          addSpecies(current, chem::fromSmiles(result.value), name);
          current.nameInput.clear();
          current.statusMessage = "Resolved and added " + name;
        } catch (const std::exception& error) {
          current.inputError = error.what();
        }
      });
}

void loadExample(AppState& state) {
  SelectionState& selection = state.selection;
  selection.operation = sol::OperationSpec{};
  selection.operation.kind = sol::OperationKind::LiquidLiquidExtraction;
  selection.operation.aqueousVolumeMl = 50.0;
  selection.operation.organicVolumeMl = 50.0;
  selection.operation.requireWaterImmiscible = true;
  selection.operation.avoidChlorinated = true;
  selection.operation.avoidPeroxideFormers = true;
  selection.speciesPresentation.clear();
  selection.inputError.clear();
  try {
    addSpecies(selection, chem::fromSmiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C"), "Caffeine");
    addSpecies(selection, chem::fromSmiles("O=C(O)c1ccccc1"), "Benzoic acid");
    selection.operation.species.back().keep = false;
    selection.operation.species.back().weight = 0.8;
    selection.statusMessage = "Loaded a caffeine / benzoic acid extraction example";
  } catch (const std::exception& error) {
    selection.inputError = error.what();
  }
}

std::string operationSignature(const sol::OperationSpec& operation) {
  std::ostringstream stream;
  stream << std::setprecision(17) << static_cast<int>(operation.kind) << '|'
         << operation.temperatureC << '|' << operation.hotTemperatureC << '|'
         << operation.coldTemperatureC << '|' << operation.aqueousVolumeMl << '|'
         << operation.organicVolumeMl << '|' << operation.pH << '|'
         << operation.requireWaterImmiscible << operation.requireWaterMiscible << '|'
         << operation.minBoilingPointC << '|' << operation.maxBoilingPointC << '|'
         << operation.avoidPeroxideFormers << operation.avoidChlorinated
         << operation.avoidAromatics << operation.excludeUnrated << '|'
         << operation.worstAcceptableClass << '|' << operation.weightSelectivity << '|'
         << operation.weightRecovery << '|' << operation.weightGreenness << '|'
         << operation.weightPracticality;
  for (const sol::SpeciesRole& species : operation.species) {
    const sol::Solute& solute = species.solute;
    stream << "||" << std::quoted(species.label) << '|' << species.keep << '|'
           << species.weight << '|' << species.amountMg << '|' << std::quoted(solute.name)
           << '|' << std::quoted(solute.canonicalSmiles) << '|'
           << solute.molarMass << '|' << solute.molarVolume << '|' << solute.meltingPoint << '|'
           << solute.entropyOfFusion << '|' << solute.meltingPointEstimated << '|'
           << solute.logP << '|' << solute.interactionRadius << '|' << solute.hansen.dispersion
           << '|' << solute.hansen.polar << '|' << solute.hansen.hydrogenBond << '|'
           << static_cast<int>(solute.ionization.ionClass) << '|' << solute.ionization.pKa << '|'
           << solute.ionization.pKaEstimated << '|' << std::quoted(solute.ionization.siteLabel)
           << '|' << solute.ionization.saltForm << '|'
           << std::quoted(solute.ionization.counterIon) << '|'
           << solute.ionization.counterIonRadiusNm << '|' << solute.ionization.netCharge << '|'
           << solute.ionization.fragmentCount;
  }
  return stream.str();
}

void requestRankingIfNeeded(AppState& state) {
  SelectionState& selection = state.selection;
  const std::string signature = operationSignature(selection.operation);
  const double now = ImGui::GetTime();

  if (selection.operation.species.empty() || !hasTarget(selection)) {
    if (selection.rankingSignature != signature || !selection.candidates.empty()) {
      selection.candidates.clear();
      selection.rankingError.clear();
      selection.rankingSignature = signature;
      ++selection.resultsRevision;
    }
    selection.observedSignature = signature;
    selection.signatureChangedAt = now;
    selection.pendingSignature.clear();
    selection.computing = false;
    return;
  }

  if (signature == selection.rankingSignature) {
    selection.observedSignature = signature;
    selection.pendingSignature.clear();
    selection.computing = false;
    return;
  }

  // Slider drags can produce dozens of distinct values. Wait for a brief
  // quiet period before queuing work so only the value the user settles on is
  // ranked; an already-running stale request is invalidated by clearing its key.
  if (signature != selection.observedSignature) {
    selection.observedSignature = signature;
    selection.signatureChangedAt = now;
    selection.pendingSignature.clear();
    selection.computing = true;
    selection.rankingError.clear();
    return;
  }
  selection.computing = true;
  constexpr double kInputSettleSeconds = 0.14;
  if (now - selection.signatureChangedAt < kInputSettleSeconds) return;
  if (signature == selection.pendingSignature) return;

  selection.pendingSignature = signature;
  selection.rankingError.clear();
  const sol::OperationSpec operation = selection.operation;
  state.tasks.run<RankingResponse>(
      [operation] {
        RankingResponse response;
        try {
          response.candidates = sol::rankSolvents(operation);
        } catch (const std::exception& error) {
          response.error = error.what();
        }
        return response;
      },
      [&state, signature](RankingResponse response) {
        SelectionState& current = state.selection;
        if (current.pendingSignature != signature) return;
        current.pendingSignature.clear();
        current.computing = false;
        current.rankingSignature = signature;
        current.rankingError = std::move(response.error);
        current.candidates = std::move(response.candidates);
        ++current.resultsRevision;
      });
}

bool drawToggleChip(const char* label, bool selected) {
  const style::Metrics& metrics = style::metrics();
  const ImGuiStyle& imguiStyle = ImGui::GetStyle();
  const ImVec2 textSize = ImGui::CalcTextSize(label);
  const ImVec2 size(textSize.x + imguiStyle.FramePadding.x * 2.0f, ImGui::GetFrameHeight());
  const ImVec2 minimum = ImGui::GetCursorScreenPos();

  ImGui::PushID(label);
  ImGui::InvisibleButton("##chip", size);
  const ImGuiID id = ImGui::GetItemID();
  const bool clicked = ImGui::IsItemClicked();
  const float hover = widgets::hoverT(id, ImGui::IsItemHovered());
  ImGui::PopID();

  const ImVec2 maximum(minimum.x + size.x, minimum.y + size.y);
  ImDrawList* draw = ImGui::GetWindowDrawList();
  const ImU32 fill = selected ? style::u32(style::col::Accent)
                              : style::mix(style::col::BgSurface, style::col::BgRaised, hover);
  const ImU32 border = selected ? style::u32(style::col::AccentActive)
                                : style::mix(style::col::Border, style::col::BorderStrong, hover);
  const ImU32 text = selected ? style::u32(style::col::OnAccent)
                              : style::u32(style::col::TextDim);
  draw->AddRectFilled(minimum, maximum, fill, metrics.radiusSm);
  draw->AddRect(minimum, maximum, border, metrics.radiusSm, 0, metrics.hairline);
  draw->AddText(ImVec2((minimum.x + maximum.x - textSize.x) * 0.5f,
                       (minimum.y + maximum.y - textSize.y) * 0.5f),
                text, label);
  return clicked;
}

void flowBeforeChip(const char* label, float& usedWidth, float availableWidth) {
  const float gap = style::metrics().gap * 0.5f;
  const float width = ImGui::CalcTextSize(label).x + ImGui::GetStyle().FramePadding.x * 2.0f;
  if (usedWidth > 0.0f && usedWidth + gap + width <= availableWidth) {
    ImGui::SameLine(0.0f, gap);
    usedWidth += gap;
  } else if (usedWidth > 0.0f) {
    usedWidth = 0.0f;
  }
  usedWidth += width;
}

bool sliderDouble(const char* visibleLabel, const char* hiddenLabel, double& value,
                  float minimum, float maximum, const char* format) {
  float temporary = static_cast<float>(value);
  const float available = ImGui::GetContentRegionAvail().x;
  const float labelWidth = ImGui::CalcTextSize(visibleLabel).x +
                           ImGui::GetStyle().ItemInnerSpacing.x;
  if (available - labelWidth >= ImGui::GetFontSize() * 4.0f) {
    ImGui::SetNextItemWidth(available - labelWidth);
    if (!ImGui::SliderFloat(visibleLabel, &temporary, minimum, maximum, format)) return false;
  } else {
    ImGui::TextWrapped("%s", visibleLabel);
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (!ImGui::SliderFloat(hiddenLabel, &temporary, minimum, maximum, format)) return false;
  }
  value = static_cast<double>(temporary);
  return true;
}

void drawOperationCard(sol::OperationKind kind, sol::OperationKind& selected, float width) {
  const style::Metrics& metrics = style::metrics();
  const float height = ImGui::GetFontSize() * 5.8f;
  const ImVec2 minimum = ImGui::GetCursorScreenPos();
  ImGui::PushID(static_cast<int>(kind));
  ImGui::InvisibleButton("##operation", ImVec2(width, height));
  const ImGuiID id = ImGui::GetItemID();
  const bool hovered = ImGui::IsItemHovered();
  const bool chosen = selected == kind;
  if (ImGui::IsItemClicked()) selected = kind;
  const float hover = widgets::hoverT(id, hovered || chosen);
  ImGui::PopID();

  const ImVec2 maximum(minimum.x + width, minimum.y + height);
  ImDrawList* draw = ImGui::GetWindowDrawList();
  const ImU32 fill = chosen ? style::mix(style::col::BgSurface, style::col::BgRaised, hover)
                            : style::mix(style::col::BgPanel, style::col::BgSurface, hover);
  const ImU32 border = chosen ? style::mix(style::col::BorderStrong, style::col::Accent, hover)
                              : style::mix(style::col::Border, style::col::BorderStrong, hover);
  draw->AddRectFilled(minimum, maximum, fill, metrics.radiusMd);
  draw->AddRect(minimum, maximum, border, metrics.radiusMd, 0,
                metrics.hairline * (chosen ? 1.8f : 1.0f));
  if (chosen) {
    draw->AddRectFilled(minimum, ImVec2(minimum.x + metrics.hairline * 3.0f, maximum.y),
                        style::u32(style::col::Accent), metrics.radiusMd);
  }

  const float padding = metrics.gap;
  const float textWidth = std::max(width - padding * 2.0f, 1.0f);
  const std::string title = ellipsizeText(sol::operationName(kind), textWidth);
  draw->AddText(ImVec2(minimum.x + padding, minimum.y + padding), style::u32(style::col::Text),
                title.c_str());
  const float purposeY = minimum.y + padding + ImGui::GetFontSize() + metrics.gap * 0.45f;
  draw->AddText(ImGui::GetFont(), ImGui::GetFontSize(),
                ImVec2(minimum.x + padding, purposeY), style::u32(style::col::TextDim),
                sol::operationDescription(kind), nullptr, textWidth);
}

void drawOperationChooser(SelectionState& state) {
  widgets::sectionHeader("Operation", style::col::Accent);
  const float available = std::max(ImGui::GetContentRegionAvail().x, 1.0f);
  const float gap = style::metrics().gap;
  int columns = 1;
  if (available >= ImGui::GetFontSize() * 42.0f) {
    columns = 3;
  } else if (available >= ImGui::GetFontSize() * 25.0f) {
    columns = 2;
  }
  const float width = std::max((available - gap * static_cast<float>(columns - 1)) /
                                   static_cast<float>(columns),
                               1.0f);
  for (size_t index = 0; index < kOperationKinds.size(); ++index) {
    if (index % static_cast<size_t>(columns) != 0) ImGui::SameLine(0.0f, gap);
    drawOperationCard(kOperationKinds[index], state.operation.kind, width);
  }
}

void drawOperationConditions(SelectionState& state) {
  widgets::sectionHeader("Working conditions", style::col::Violet);
  if (!widgets::beginCard("##operation_conditions", ImVec2(0.0f, 0.0f),
                          style::col::BgSurface)) {
    return;
  }

  float reveal = 1.0f;
  for (sol::OperationKind kind : kOperationKinds) {
    ImGui::PushID(static_cast<int>(kind));
    const float progress =
        widgets::hoverT(ImGui::GetID("condition_reveal"), state.operation.kind == kind);
    ImGui::PopID();
    if (state.operation.kind == kind) reveal = progress;
  }
  ImGui::PushID(static_cast<int>(state.operation.kind));
  ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * std::max(reveal, 0.05f));
  switch (state.operation.kind) {
    case sol::OperationKind::LiquidLiquidExtraction:
      sliderDouble("Temperature (C)", "##extraction_temperature", state.operation.temperatureC,
                   -20.0f, 120.0f, "%.1f");
      sliderDouble("Aqueous volume (mL)", "##aqueous_volume",
                   state.operation.aqueousVolumeMl, 1.0f, 1000.0f, "%.0f");
      sliderDouble("Organic volume (mL)", "##organic_volume",
                   state.operation.organicVolumeMl, 1.0f, 1000.0f, "%.0f");
      if (state.operation.pH == sol::kAutoPH) {
        if (drawToggleChip("SELF-BUFFERED PH", true)) state.operation.pH = 7.0;
        ImGui::TextWrapped("The solute sets its own saturated-solution pH.");
      } else {
        if (drawToggleChip("FIXED PH", true)) state.operation.pH = sol::kAutoPH;
        sliderDouble("pH", "##extraction_ph", state.operation.pH, 0.0f, 14.0f, "%.1f");
      }
      break;
    case sol::OperationKind::Recrystallisation:
      sliderDouble("Dissolve hot (C)", "##hot_temperature",
                   state.operation.hotTemperatureC, 20.0f, 200.0f, "%.1f");
      sliderDouble("Crystallise cold (C)", "##cold_temperature",
                   state.operation.coldTemperatureC, -20.0f, 80.0f, "%.1f");
      break;
    case sol::OperationKind::AntiSolventPrecipitation:
      sliderDouble("Addition temperature (C)", "##antisolvent_temperature",
                   state.operation.temperatureC, -20.0f, 120.0f, "%.1f");
      ImGui::TextWrapped("The engine sweeps compatible solvent / anti-solvent partners and reports "
                         "the partner fraction at the best precipitation window.");
      break;
    case sol::OperationKind::Trituration:
      sliderDouble("Wash temperature (C)", "##trituration_temperature",
                   state.operation.temperatureC, -20.0f, 120.0f, "%.1f");
      break;
    case sol::OperationKind::ChromatographyMobilePhase:
      sliderDouble("Column temperature (C)", "##chromatography_temperature",
                   state.operation.temperatureC, 0.0f, 80.0f, "%.1f");
      ImGui::TextWrapped("The polarity and selectivity window is derived from the KEEP and REMOVE "
                         "species; use the boiling-point window below for volatility limits.");
      break;
    case sol::OperationKind::ReactionMedium:
      sliderDouble("Reaction temperature (C)", "##reaction_temperature",
                   state.operation.temperatureC, -20.0f, 200.0f, "%.1f");
      break;
  }
  ImGui::PopStyleVar();
  ImGui::PopID();
  widgets::endCard();
}

bool drawSegmentButton(const char* label, bool selected, ImVec2 size) {
  ImGui::PushStyleColor(ImGuiCol_Button, selected ? style::col::Accent : style::col::BgPanel);
  ImGui::PushStyleColor(ImGuiCol_ButtonHovered,
                        selected ? style::col::AccentHover : style::col::BgRaised);
  ImGui::PushStyleColor(ImGuiCol_ButtonActive,
                        selected ? style::col::AccentActive : style::col::BgSurface);
  ImGui::PushStyleColor(ImGuiCol_Text, selected ? style::col::OnAccent : style::col::TextDim);
  const bool clicked = ImGui::Button(label, size);
  ImGui::PopStyleColor(4);
  return clicked;
}

bool drawSpeciesRow(SelectionState& state, size_t index) {
  sol::SpeciesRole& species = state.operation.species[index];
  const SelectionSpeciesPresentation& presentation = state.speciesPresentation[index];
  ImGui::PushID(static_cast<int>(index));
  bool remove = false;
  if (widgets::beginCard("##species_row", ImVec2(0.0f, 0.0f), style::col::BgRaised)) {
    const float removeSize = ImGui::GetFrameHeight();
    const float titleWidth = std::max(ImGui::GetContentRegionAvail().x - removeSize -
                                          ImGui::GetStyle().ItemSpacing.x,
                                      1.0f);
    const bool heading = style::pushFont(style::fonts::semibold());
    const std::string fitted = ellipsizeText(species.label, titleWidth);
    ImGui::TextUnformatted(fitted.c_str());
    style::popFont(heading);
    if (ImGui::IsItemHovered() && fitted != species.label) ImGui::SetTooltip("%s", species.label.c_str());
    ImGui::SameLine(0.0f, ImGui::GetStyle().ItemSpacing.x);
    remove = widgets::iconButton("##remove_species", icons::Icon::Close,
                                 ImVec2(removeSize, removeSize), false, "Remove species");

    char identity[128];
    std::snprintf(identity, sizeof(identity), "%s  |  %.2f g/mol",
                  presentation.formula.empty() ? "formula unavailable" : presentation.formula.c_str(),
                  species.solute.molarMass);
    const bool mono = style::pushFont(style::fonts::mono());
    ImGui::TextWrapped("%s", identity);
    style::popFont(mono);

    const float segmentGap = style::metrics().hairline;
    const float segmentWidth = std::max((ImGui::GetContentRegionAvail().x - segmentGap) * 0.5f, 1.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(segmentGap, ImGui::GetStyle().ItemSpacing.y));
    if (drawSegmentButton("KEEP##role", species.keep, ImVec2(segmentWidth, 0.0f))) species.keep = true;
    ImGui::SameLine(0.0f, segmentGap);
    if (drawSegmentButton("REMOVE##role", !species.keep, ImVec2(segmentWidth, 0.0f))) species.keep = false;
    ImGui::PopStyleVar();

    sliderDouble("Importance", "##species_weight", species.weight, 0.1f, 3.0f, "%.2f x");
    sliderDouble("Amount (mg, optional)", "##species_amount", species.amountMg, 0.0f, 5000.0f,
                 "%.0f");
    widgets::endCard();
  }
  ImGui::PopID();
  return remove;
}

void drawSpeciesBuilder(AppState& state) {
  SelectionState& selection = state.selection;
  widgets::sectionHeader("Species", style::col::Teal);
  if (!widgets::beginCard("##species_builder", ImVec2(0.0f, 0.0f), style::col::BgSurface)) return;

  const float available = ImGui::GetContentRegionAvail().x;
  float used = 0.0f;
  flowBeforeChip("FROM SKETCH", used, available);
  if (drawToggleChip("FROM SKETCH", selection.addMode == 0)) selection.addMode = 0;
  flowBeforeChip("SMILES", used, available);
  if (drawToggleChip("SMILES", selection.addMode == 1)) selection.addMode = 1;
  flowBeforeChip("CHEMICAL NAME", used, available);
  if (drawToggleChip("CHEMICAL NAME", selection.addMode == 2)) selection.addMode = 2;

  ImGui::Spacing();
  if (selection.addMode == 0) {
    ImGui::TextWrapped("Add the complete current sketch as one species; salt fragments are kept "
                       "with their parent structure.");
    if (widgets::primaryButton("Add current sketch", ImVec2(ImGui::GetContentRegionAvail().x, 0.0f))) {
      addSketchSpecies(state);
    }
  } else if (selection.addMode == 1) {
    ImGui::SetNextItemWidth(-FLT_MIN);
    const bool enter = widgets::stringInputWithHint("##species_smiles", "SMILES, e.g. CC(=O)Oc1ccccc1C(=O)O",
                                           selection.smilesInput,
                                           ImGuiInputTextFlags_EnterReturnsTrue, true);
    if (enter || widgets::primaryButton("Add SMILES", ImVec2(ImGui::GetContentRegionAvail().x, 0.0f))) {
      addSmilesSpecies(state);
    }
  } else {
    ImGui::SetNextItemWidth(-FLT_MIN);
    const bool enter = widgets::stringInputWithHint("##species_name", "IUPAC or common name",
                                           selection.nameInput,
                                           ImGuiInputTextFlags_EnterReturnsTrue);
    if (selection.nameLookupRunning) ImGui::BeginDisabled();
    if (enter || widgets::primaryButton(selection.nameLookupRunning ? "Resolving..." : "Resolve and add",
                                        ImVec2(ImGui::GetContentRegionAvail().x, 0.0f))) {
      addNameSpecies(state);
    }
    if (selection.nameLookupRunning) ImGui::EndDisabled();
  }

  if (!selection.inputError.empty()) {
    ImGui::PushStyleColor(ImGuiCol_Text, style::col::Danger);
    ImGui::TextWrapped("%s", selection.inputError.c_str());
    ImGui::PopStyleColor();
  }

  if (selection.speciesPresentation.size() < selection.operation.species.size()) {
    selection.speciesPresentation.resize(selection.operation.species.size());
  } else if (selection.speciesPresentation.size() > selection.operation.species.size()) {
    selection.speciesPresentation.resize(selection.operation.species.size());
  }

  ImGui::Spacing();
  const float listHeight = ImGui::GetFontSize() * 16.0f;
  ImGui::PushStyleVar(ImGuiStyleVar_ChildRounding, style::metrics().radiusMd);
  ImGui::PushStyleColor(ImGuiCol_ChildBg, style::col::BgPanel);
  if (ImGui::BeginChild("##species_list", ImVec2(0.0f, listHeight), ImGuiChildFlags_Borders)) {
    if (selection.operation.species.empty()) {
      ImGui::TextWrapped("Add every compound that matters. Mark products to isolate as KEEP and "
                         "impurities to reject as REMOVE.");
    } else {
      size_t removeIndex = selection.operation.species.size();
      for (size_t index = 0; index < selection.operation.species.size(); ++index) {
        if (drawSpeciesRow(selection, index)) removeIndex = index;
        if (index + 1 < selection.operation.species.size()) ImGui::Spacing();
      }
      if (removeIndex < selection.operation.species.size()) {
        selection.operation.species.erase(selection.operation.species.begin() +
                                           static_cast<std::ptrdiff_t>(removeIndex));
        selection.speciesPresentation.erase(selection.speciesPresentation.begin() +
                                            static_cast<std::ptrdiff_t>(removeIndex));
      }
    }
  }
  ImGui::EndChild();
  ImGui::PopStyleColor();
  ImGui::PopStyleVar();
  widgets::endCard();
}

void drawConstraintChips(SelectionState& state) {
  sol::OperationSpec& operation = state.operation;
  const float available = ImGui::GetContentRegionAvail().x;
  float used = 0.0f;

  flowBeforeChip("WATER-IMMISCIBLE", used, available);
  if (drawToggleChip("WATER-IMMISCIBLE", operation.requireWaterImmiscible)) {
    operation.requireWaterImmiscible = !operation.requireWaterImmiscible;
    if (operation.requireWaterImmiscible) operation.requireWaterMiscible = false;
  }
  flowBeforeChip("WATER-MISCIBLE", used, available);
  if (drawToggleChip("WATER-MISCIBLE", operation.requireWaterMiscible)) {
    operation.requireWaterMiscible = !operation.requireWaterMiscible;
    if (operation.requireWaterMiscible) operation.requireWaterImmiscible = false;
  }
  flowBeforeChip("NO PEROXIDE FORMERS", used, available);
  if (drawToggleChip("NO PEROXIDE FORMERS", operation.avoidPeroxideFormers)) {
    operation.avoidPeroxideFormers = !operation.avoidPeroxideFormers;
  }
  flowBeforeChip("NO CHLORINATED", used, available);
  if (drawToggleChip("NO CHLORINATED", operation.avoidChlorinated)) {
    operation.avoidChlorinated = !operation.avoidChlorinated;
  }
  flowBeforeChip("NO AROMATICS", used, available);
  if (drawToggleChip("NO AROMATICS", operation.avoidAromatics)) {
    operation.avoidAromatics = !operation.avoidAromatics;
  }
  flowBeforeChip("RATED ONLY", used, available);
  if (drawToggleChip("RATED ONLY", operation.excludeUnrated)) {
    operation.excludeUnrated = !operation.excludeUnrated;
  }
}

void drawClassSelector(sol::OperationSpec& operation) {
  static constexpr std::array<const char*, 5> kClasses = {
      "", "recommended", "problematic", "hazardous", "highly hazardous"};
  ImGui::TextColored(style::col::TextDim, "WORST ACCEPTABLE CHEM21 CLASS");
  const char* preview = operation.worstAcceptableClass.empty()
                            ? "Any class"
                            : operation.worstAcceptableClass.c_str();
  ImGui::SetNextItemWidth(-FLT_MIN);
  if (ImGui::BeginCombo("##chem21_limit", preview)) {
    for (const char* item : kClasses) {
      const bool any = item[0] == '\0';
      const std::string value = any ? std::string() : std::string(item);
      const bool selected = operation.worstAcceptableClass == value;
      ImGui::PushStyleColor(ImGuiCol_Text, any ? style::col::Text : chem21Colour(value));
      if (ImGui::Selectable(any ? "Any class" : item, selected)) {
        operation.worstAcceptableClass = value;
      }
      ImGui::PopStyleColor();
      if (selected) ImGui::SetItemDefaultFocus();
    }
    ImGui::EndCombo();
  }
  if (operation.worstAcceptableClass.empty()) {
    widgets::badge("UNRESTRICTED", style::col::TextDim);
  } else {
    widgets::badge(operation.worstAcceptableClass.c_str(),
                   chem21Colour(operation.worstAcceptableClass));
  }
}

void drawConstraints(SelectionState& state) {
  widgets::sectionHeader("Constraints", style::col::Accent);
  if (!widgets::beginCard("##constraints", ImVec2(0.0f, 0.0f), style::col::BgSurface)) return;
  drawConstraintChips(state);
  ImGui::Spacing();
  sliderDouble("Minimum boiling point (C)", "##minimum_boiling",
               state.operation.minBoilingPointC, 0.0f, 250.0f, "%.0f");
  sliderDouble("Maximum boiling point (C)", "##maximum_boiling",
               state.operation.maxBoilingPointC, 0.0f, 300.0f, "%.0f");
  ImGui::TextWrapped("A zero bound is unconstrained.");
  if (state.operation.minBoilingPointC > 0.0 &&
      state.operation.maxBoilingPointC > 0.0 &&
      state.operation.minBoilingPointC > state.operation.maxBoilingPointC) {
    ImGui::PushStyleColor(ImGuiCol_Text, style::col::Danger);
    ImGui::TextWrapped("The minimum is above the maximum; no solvent can satisfy this window.");
    ImGui::PopStyleColor();
    if (widgets::ghostButton("Swap boiling bounds")) {
      std::swap(state.operation.minBoilingPointC, state.operation.maxBoilingPointC);
    }
  }
  drawClassSelector(state.operation);
  widgets::endCard();
}

void drawWeightControl(const char* label, const char* id, double& value,
                       const char* explanation) {
  sliderDouble(label, id, value, 0.0f, 3.0f, "%.2f x");
  ImGui::PushStyleColor(ImGuiCol_Text, style::col::TextDim);
  ImGui::TextWrapped("%s", explanation);
  ImGui::PopStyleColor();
}

void drawPriorities(SelectionState& state) {
  widgets::sectionHeader("Scoring priorities", style::col::Teal);
  if (!widgets::beginCard("##priorities", ImVec2(0.0f, 0.0f), style::col::BgSurface)) return;
  drawWeightControl("Selectivity", "##weight_selectivity", state.operation.weightSelectivity,
                    "Raise this to reward separation between KEEP and REMOVE species.");
  drawWeightControl("Recovery", "##weight_recovery", state.operation.weightRecovery,
                    "Raise this to favour recovering more of every KEEP species.");
  drawWeightControl("Greenness", "##weight_greenness", state.operation.weightGreenness,
                    "Raise this to favour better CHEM21 safety, health and environment ratings.");
  drawWeightControl("Practicality", "##weight_practicality", state.operation.weightPracticality,
                    "Raise this to favour workable boiling, freezing, cost and handling properties.");
  widgets::endCard();
}

void drawBuilder(AppState& state) {
  drawOperationChooser(state.selection);
  drawOperationConditions(state.selection);
  drawSpeciesBuilder(state);
  drawConstraints(state.selection);
  drawPriorities(state.selection);
  if (!state.selection.statusMessage.empty()) {
    widgets::sectionHeader("Status", style::col::Success);
    if (widgets::beginCard("##selection_status", ImVec2(0.0f, 0.0f), style::col::BgSurface)) {
      ImGui::TextWrapped("%s", state.selection.statusMessage.c_str());
      widgets::endCard();
    }
  }
}

void drawBar(const char* id, double value, ImVec4 colour, float height, float reveal = 1.0f) {
  const float width = std::max(ImGui::GetContentRegionAvail().x, 1.0f);
  const ImVec2 minimum = ImGui::GetCursorScreenPos();
  ImGui::InvisibleButton(id, ImVec2(width, height));
  const ImVec2 maximum(minimum.x + width, minimum.y + height);
  ImDrawList* draw = ImGui::GetWindowDrawList();
  draw->AddRectFilled(minimum, maximum, style::u32(style::col::BgPanel),
                      style::metrics().radiusSm);
  const float fraction = std::clamp(static_cast<float>(value) * reveal, 0.0f, 1.0f);
  if (fraction > 0.0f) {
    draw->AddRectFilled(minimum, ImVec2(minimum.x + width * fraction, maximum.y),
                        style::u32(colour), style::metrics().radiusSm);
  }
  draw->AddRect(minimum, maximum, style::u32(style::col::Border), style::metrics().radiusSm, 0,
                style::metrics().hairline);
}

double criterionScore(const sol::SolventCandidate& candidate, const char* name) {
  for (const sol::Criterion& criterion : candidate.criteria) {
    if (criterion.name == name) return criterion.score;
  }
  return 0.0;
}

void syncResultAnimation(const SelectionState& state, ResultAnimation& animation) {
  if (animation.revision == state.resultsRevision) return;
  animation.previousRanks = std::move(animation.currentRanks);
  animation.currentRanks.clear();
  for (size_t index = 0; index < state.candidates.size(); ++index) {
    animation.currentRanks.emplace(candidateKey(state.candidates[index]),
                                   static_cast<int>(index) + 1);
  }
  animation.revision = state.resultsRevision;
  animation.startedAt = ImGui::GetTime();
}

float resultReveal(const ResultAnimation& animation) {
  constexpr double kRevealSeconds = 0.38;
  return easeOutCubic(static_cast<float>(
      std::clamp((ImGui::GetTime() - animation.startedAt) / kRevealSeconds, 0.0, 1.0)));
}

int priorRankDelta(const ResultAnimation& animation, const std::string& key, int currentRank) {
  const auto found = animation.previousRanks.find(key);
  if (found == animation.previousRanks.end()) return 0;
  return found->second - currentRank;
}

void drawCriterion(const sol::Criterion& criterion, float reveal) {
  ImGui::PushID(criterion.name.c_str());
  ImGui::TextColored(style::col::Text, "%s", criterion.name.c_str());
  ImGui::SameLine();
  ImGui::TextColored(style::col::TextDim, "%.2f  x  %.2f", criterion.score, criterion.weight);
  drawBar("##criterion_bar", criterion.score, scoreColour(criterion.score),
          std::max(style::metrics().hairline * 5.0f, ImGui::GetFontSize() * 0.36f), reveal);
  ImGui::PushStyleColor(ImGuiCol_Text, style::col::TextDim);
  ImGui::TextWrapped("%s", criterion.detail.c_str());
  ImGui::PopStyleColor();
  ImGui::PopID();
}

void drawStatGrid(const sol::SolventCandidate& candidate) {
  const std::array<std::pair<const char*, std::string>, 5> stats = {{
      {"TARGET SOLUBILITY", formatSolubility(candidate.targetSolubilityGPerMl)},
      {"CONTAMINANT SOLUBILITY", formatSolubility(candidate.contaminantSolubilityGPerMl)},
      {"SELECTIVITY", formatRatio(candidate.selectivity)},
      {"RECOVERY", formatPercent(candidate.recoveryFraction)},
      {"BOILING POINT", formatBoilingPoint(candidate)},
  }};
  const float available = std::max(ImGui::GetContentRegionAvail().x, 1.0f);
  int columns = 1;
  if (available >= ImGui::GetFontSize() * 38.0f) {
    columns = 3;
  } else if (available >= ImGui::GetFontSize() * 22.0f) {
    columns = 2;
  }
  const float gap = style::metrics().gap;
  const float cardWidth = std::max((available - gap * static_cast<float>(columns - 1)) /
                                       static_cast<float>(columns),
                                   1.0f);
  const float cardHeight = ImGui::GetFontSize() * 3.25f;
  for (size_t index = 0; index < stats.size(); ++index) {
    if (index % static_cast<size_t>(columns) != 0) ImGui::SameLine(0.0f, gap);
    widgets::statCard(stats[index].first, stats[index].second.c_str(),
                      ImVec2(cardWidth, cardHeight));
  }
  if (candidate.solvent && !candidate.solvent->chem21Class.empty()) {
    widgets::badge(candidate.solvent->chem21Class.c_str(),
                   chem21Colour(candidate.solvent->chem21Class));
  } else {
    widgets::badge("CHEM21 UNRATED", style::col::TextDim);
  }
}

void drawArithmetic(const sol::SolventCandidate& candidate) {
  const bool open = ImGui::TreeNodeEx(
      "Why this ranking", ImGuiTreeNodeFlags_SpanAvailWidth |
                              ImGuiTreeNodeFlags_NoTreePushOnOpen);
  if (!open) return;
  ImGui::Indent();
  double numerator = 0.0;
  double denominator = 0.0;
  for (const sol::Criterion& criterion : candidate.criteria) {
    const double contribution = criterion.score * criterion.weight;
    numerator += contribution;
    denominator += criterion.weight;
    ImGui::TextWrapped("%s: %.3f x %.3f = %.3f", criterion.name.c_str(),
                       criterion.score, criterion.weight, contribution);
  }
  ImGui::Separator();
  ImGui::TextWrapped("Weighted total: %.3f / %.3f = %.3f", numerator, denominator,
                     denominator > 0.0 ? numerator / denominator : 0.0);
  ImGui::TextWrapped("Engine score used for ranking: %.3f", candidate.score);
  ImGui::Unindent();
}

void sendToExtraction(AppState& state, const sol::SolventCandidate& candidate) {
  if (!candidate.solvent || !candidate.partner) return;
  const sol::Solvent* first = candidate.solvent;
  const sol::Solvent* second = candidate.partner;
  double firstFraction = 1.0 - candidate.partnerFraction;
  double secondFraction = candidate.partnerFraction;
  if (secondFraction <= 0.0 || secondFraction >= 1.0) {
    firstFraction = 0.5;
    secondFraction = 0.5;
  }
  if (second->family == "water" && first->family != "water") {
    std::swap(first, second);
    std::swap(firstFraction, secondFraction);
  }

  ExtractionImport& import = state.solubility.extractionImport;
  import.pending = true;
  import.solventIdA = first->id;
  import.solventIdB = second->id;
  import.volumeMlA = firstFraction * 100.0;
  import.volumeMlB = secondFraction * 100.0;
  import.soluteMassMg = targetMassMg(state.selection);
  state.tab = MainTab::Extraction;
  state.tabChangeRequested = true;
  state.selection.statusMessage = "Candidate sent to the Extraction Calculator";
  state.statusMessage = state.selection.statusMessage;
}

void loadIntoSuite(AppState& state, const sol::SolventCandidate& candidate) {
  if (!candidate.solvent) return;
  SolubilityState& suite = state.solubility;
  suite.solventIds[0] = candidate.solvent->id;
  suite.ratios[0] = 1.0f;
  if (candidate.partner) {
    suite.solventCount = 2;
    suite.solventIds[1] = candidate.partner->id;
    const double fraction = candidate.partnerFraction > 0.0 && candidate.partnerFraction < 1.0
                                ? candidate.partnerFraction
                                : 0.5;
    suite.ratios[0] = static_cast<float>(1.0 - fraction);
    suite.ratios[1] = static_cast<float>(fraction);
  } else {
    suite.solventCount = 1;
    suite.solventIds[1].clear();
    suite.ratios[1] = 1.0f;
  }
  suite.solventIds[2].clear();
  suite.ratios[2] = 1.0f;

  for (const sol::SpeciesRole& species : state.selection.operation.species) {
    if (!species.keep) continue;
    suite.solute = species.solute;
    suite.soluteValid = true;
    suite.soluteError.clear();
    suite.soluteNote = "Imported from the solvent ranking: " + species.label;
    suite.useSketch = false;
    suite.manualSmiles = species.solute.canonicalSmiles;
    suite.overrideSolute = false;
    ++suite.soluteVersion;
    break;
  }
  suite.sweepSignature.clear();
  suite.screeningSignature.clear();
  suite.statusMessage = "Candidate blend loaded into the prediction workspace";
  // Select and Predict are two modes of this one panel, so loading a candidate
  // switches the mode in place instead of requesting a tab change.
  suite.suiteMode = SuiteMode::Predict;
  state.selection.statusMessage = suite.statusMessage;
  state.statusMessage = suite.statusMessage;
}

bool isCompared(const SelectionState& state, const std::string& key) {
  return std::find(state.comparedCandidates.begin(), state.comparedCandidates.end(), key) !=
         state.comparedCandidates.end();
}

void toggleCompared(SelectionState& state, const std::string& key) {
  const auto found = std::find(state.comparedCandidates.begin(), state.comparedCandidates.end(), key);
  if (found != state.comparedCandidates.end()) {
    state.comparedCandidates.erase(found);
    return;
  }
  if (state.comparedCandidates.size() < 3) {
    state.comparedCandidates.push_back(key);
  } else {
    state.statusMessage = "Comparison is limited to three candidates";
  }
}

void drawCandidateCard(AppState& state, const sol::SolventCandidate& candidate,
                       int displayedRank, int engineRank, ResultAnimation& animation) {
  SelectionState& selection = state.selection;
  const std::string key = candidateKey(candidate);
  const bool pinned = isCompared(selection, key);
  ImGui::PushID(key.c_str());
  if (!widgets::beginCard("##candidate", ImVec2(0.0f, 0.0f), style::col::BgSurface)) {
    ImGui::PopID();
    return;
  }

  char rank[24];
  std::snprintf(rank, sizeof(rank), "#%d", displayedRank);
  widgets::badge(rank, displayedRank == 1 ? style::col::Accent : style::col::TextDim);
  ImGui::SameLine();
  widgets::badge(candidate.estimated ? "ESTIMATE" : "MEASURED",
                 candidate.estimated ? style::col::Violet : style::col::Success);
  const int delta = priorRankDelta(animation, key, engineRank);
  if (delta != 0 && resultReveal(animation) < 1.0f) {
    char movement[32];
    std::snprintf(movement, sizeof(movement), "%s%d", delta > 0 ? "UP " : "DOWN ",
                  std::abs(delta));
    const float movementWidth = ImGui::CalcTextSize(movement).x +
                                ImGui::GetStyle().FramePadding.x * 2.0f;
    if (ImGui::GetContentRegionAvail().x >= movementWidth) ImGui::SameLine();
    widgets::badge(movement, delta > 0 ? style::col::Success : style::col::Danger);
  }
  const bool heading = style::pushFont(style::fonts::semibold());
  const std::string fullTitle = candidateTitle(candidate);
  ImGui::TextWrapped("%s", fullTitle.c_str());
  style::popFont(heading);

  const float reveal = resultReveal(animation);
  const bool mono = style::pushFont(style::fonts::mono());
  ImGui::TextColored(scoreColour(candidate.score), "%.3f", candidate.score);
  style::popFont(mono);
  ImGui::SameLine();
  ImGui::TextColored(style::col::TextDim, "OVERALL SCORE");
  drawBar("##score_bar", candidate.score, scoreColour(candidate.score),
          std::max(style::metrics().hairline * 7.0f, ImGui::GetFontSize() * 0.55f), reveal);

  widgets::sectionHeader("Evidence", style::col::Teal);
  for (const sol::Criterion& criterion : candidate.criteria) drawCriterion(criterion, reveal);

  widgets::sectionHeader("Physical readout", style::col::Accent);
  drawStatGrid(candidate);

  if (!candidate.warnings.empty()) {
    widgets::sectionHeader("Advisories", style::col::Accent);
    for (const std::string& warning : candidate.warnings) {
      ImGui::PushStyleColor(ImGuiCol_Text, style::col::AccentHover);
      ImGui::TextWrapped("!  %s", warning.c_str());
      ImGui::PopStyleColor();
    }
  }

  drawArithmetic(candidate);
  ImGui::Spacing();

  const float available = ImGui::GetContentRegionAvail().x;
  const float gap = style::metrics().gap;
  const bool actionsInline = available >= ImGui::GetFontSize() * 42.0f;
  const int actionCount = selection.compareMode ? 3 : 2;
  const float actionWidth = actionsInline
                                ? std::max((available - gap * static_cast<float>(actionCount - 1)) /
                                               static_cast<float>(actionCount),
                                           1.0f)
                                : available;
  if (!candidate.partner) ImGui::BeginDisabled();
  if (widgets::ghostButton("Send to Extraction", ImVec2(actionWidth, 0.0f))) {
    sendToExtraction(state, candidate);
  }
  if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled) && !candidate.partner) {
    ImGui::SetTooltip("A two-solvent candidate is required for an extraction hand-off");
  }
  if (!candidate.partner) ImGui::EndDisabled();

  if (actionsInline) ImGui::SameLine(0.0f, gap);
  if (widgets::actionButton("##load_into_predict", icons::Icon::Flask,
                            "Load into Predict", ImVec2(actionWidth, 0.0f), true,
                            "Load this blend into the prediction workspace")) {
    loadIntoSuite(state, candidate);
  }
  if (selection.compareMode) {
    if (actionsInline) ImGui::SameLine(0.0f, gap);
    const bool comparisonFull = selection.comparedCandidates.size() >= 3 && !pinned;
    if (comparisonFull) ImGui::BeginDisabled();
    if (widgets::ghostButton(pinned ? "Unpin comparison" : "Pin for comparison",
                             ImVec2(actionWidth, 0.0f))) {
      toggleCompared(selection, key);
    }
    if (comparisonFull) ImGui::EndDisabled();
  }

  widgets::endCard();
  const ImVec2 cardMinimum = ImGui::GetItemRectMin();
  const ImVec2 cardMaximum = ImGui::GetItemRectMax();
  const bool hovered = ImGui::IsItemHovered();
  const float hover = widgets::hoverT(ImGui::GetID("card_hover"), hovered || pinned);
  if (hover > 0.0f) {
    ImGui::GetWindowDrawList()->AddRect(
        cardMinimum, cardMaximum,
        style::mix(style::col::BorderStrong, pinned ? style::col::Accent : style::col::Teal,
                   hover),
        style::metrics().radiusMd, 0,
        style::metrics().hairline * (1.0f + hover));
  }
  ImGui::PopID();
}

const sol::SolventCandidate* findCandidate(const SelectionState& state, const std::string& key) {
  for (const sol::SolventCandidate& candidate : state.candidates) {
    if (candidateKey(candidate) == key) return &candidate;
  }
  return nullptr;
}

void pruneCompared(SelectionState& state) {
  state.comparedCandidates.erase(
      std::remove_if(state.comparedCandidates.begin(), state.comparedCandidates.end(),
                     [&](const std::string& key) { return findCandidate(state, key) == nullptr; }),
      state.comparedCandidates.end());
}

void drawComparisonBoard(SelectionState& state, ResultAnimation& animation) {
  if (!state.compareMode) return;
  pruneCompared(state);
  widgets::sectionHeader("Pinned comparison", style::col::Violet);
  if (state.comparedCandidates.empty()) {
    if (widgets::beginCard("##comparison_empty", ImVec2(0.0f, 0.0f), style::col::BgSurface)) {
      ImGui::TextWrapped("Pin up to three candidates below to compare their trade-offs here.");
      widgets::endCard();
    }
    return;
  }

  const float available = std::max(ImGui::GetContentRegionAvail().x, 1.0f);
  const float gap = style::metrics().gap;
  const bool sideBySide = available >= ImGui::GetFontSize() * 44.0f;
  const int columns = sideBySide ? static_cast<int>(state.comparedCandidates.size()) : 1;
  const float width = std::max((available - gap * static_cast<float>(columns - 1)) /
                                   static_cast<float>(columns),
                               1.0f);
  std::string removeKey;
  for (size_t index = 0; index < state.comparedCandidates.size(); ++index) {
    const std::string& key = state.comparedCandidates[index];
    const sol::SolventCandidate* candidate = findCandidate(state, key);
    if (!candidate) continue;
    if (sideBySide && index > 0) ImGui::SameLine(0.0f, gap);
    ImGui::PushID(key.c_str());
    if (widgets::beginCard("##comparison", ImVec2(width, 0.0f), style::col::BgRaised)) {
      const float square = ImGui::GetFrameHeight();
      const float titleWidth = std::max(ImGui::GetContentRegionAvail().x - square - gap, 1.0f);
      const bool heading = style::pushFont(style::fonts::semibold());
      const std::string title = candidateTitle(*candidate);
      const std::string fitted = ellipsizeText(title, titleWidth);
      ImGui::TextUnformatted(fitted.c_str());
      style::popFont(heading);
      if (ImGui::IsItemHovered() && title != fitted) ImGui::SetTooltip("%s", title.c_str());
      ImGui::SameLine(0.0f, gap);
      if (widgets::iconButton("##unpin", icons::Icon::Close, ImVec2(square, square), false,
                              "Unpin candidate")) {
        removeKey = key;
      }
      char score[32];
      std::snprintf(score, sizeof(score), "%.3f", candidate->score);
      widgets::statCard("SCORE", score,
                        ImVec2(ImGui::GetContentRegionAvail().x, ImGui::GetFontSize() * 3.1f));
      drawBar("##comparison_score", candidate->score, scoreColour(candidate->score),
              std::max(style::metrics().hairline * 6.0f, ImGui::GetFontSize() * 0.45f),
              resultReveal(animation));
      ImGui::TextWrapped("Selectivity  %s", formatRatio(candidate->selectivity).c_str());
      ImGui::TextWrapped("Recovery  %s", formatPercent(candidate->recoveryFraction).c_str());
      ImGui::TextWrapped("Greenness  %.2f", criterionScore(*candidate, "Greenness"));
      if (candidate->solvent && !candidate->solvent->chem21Class.empty()) {
        widgets::badge(candidate->solvent->chem21Class.c_str(),
                       chem21Colour(candidate->solvent->chem21Class));
      } else {
        widgets::badge("UNRATED", style::col::TextDim);
      }
      widgets::endCard();
    }
    ImGui::PopID();
  }
  if (!removeKey.empty()) toggleCompared(state, removeKey);
}

void drawResultControls(SelectionState& state) {
  const float available = std::max(ImGui::GetContentRegionAvail().x, 1.0f);
  const bool inlineControls = available >= ImGui::GetFontSize() * 38.0f;
  if (inlineControls && ImGui::BeginTable(
                            "##result_controls", 3,
                            ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_NoSavedSettings)) {
    ImGui::TableSetupColumn("##search_column", ImGuiTableColumnFlags_WidthStretch, 0.52f);
    ImGui::TableSetupColumn("##sort_column", ImGuiTableColumnFlags_WidthStretch, 0.27f);
    ImGui::TableSetupColumn("##compare_column", ImGuiTableColumnFlags_WidthStretch, 0.21f);
    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    ImGui::SetNextItemWidth(-FLT_MIN);
    widgets::stringInputWithHint("##solvent_search", "Filter solvent names", state.resultSearch);
    ImGui::TableNextColumn();
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::Combo("##result_sort", &state.sortMode, kSortLabels.data(),
                 static_cast<int>(kSortLabels.size()));
    ImGui::TableNextColumn();
    if (drawToggleChip("COMPARE", state.compareMode)) state.compareMode = !state.compareMode;
    ImGui::EndTable();
  } else {
    ImGui::SetNextItemWidth(-FLT_MIN);
    widgets::stringInputWithHint("##solvent_search", "Filter solvent names", state.resultSearch);
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::Combo("##result_sort", &state.sortMode, kSortLabels.data(),
                 static_cast<int>(kSortLabels.size()));
    if (drawToggleChip("COMPARE UP TO 3", state.compareMode)) state.compareMode = !state.compareMode;
  }
}

std::vector<const sol::SolventCandidate*> visibleCandidates(const SelectionState& state) {
  std::vector<const sol::SolventCandidate*> result;
  result.reserve(state.candidates.size());
  for (const sol::SolventCandidate& candidate : state.candidates) {
    const std::string primary = candidate.solvent ? candidate.solvent->name : std::string();
    const std::string partner = candidate.partner ? candidate.partner->name : std::string();
    if (!widgets::containsCaseInsensitive(primary, state.resultSearch) &&
        !widgets::containsCaseInsensitive(partner, state.resultSearch)) {
      continue;
    }
    result.push_back(&candidate);
  }
  const int mode = std::clamp(state.sortMode, 0, static_cast<int>(kSortLabels.size()) - 1);
  std::stable_sort(result.begin(), result.end(), [mode](const auto* left, const auto* right) {
    switch (mode) {
      case 1: return left->selectivity > right->selectivity;
      case 2: return left->recoveryFraction > right->recoveryFraction;
      case 3: return criterionScore(*left, "Greenness") > criterionScore(*right, "Greenness");
      case 4: {
        const double leftBp =
            left->solvent && left->solvent->boilingPoint > 0.0
                ? left->solvent->boilingPoint
                : std::numeric_limits<double>::infinity();
        const double rightBp =
            right->solvent && right->solvent->boilingPoint > 0.0
                ? right->solvent->boilingPoint
                : std::numeric_limits<double>::infinity();
        return leftBp < rightBp;
      }
      default: return left->score > right->score;
    }
  });
  return result;
}

std::vector<std::string> bindingConstraints(const sol::OperationSpec& operation) {
  std::vector<std::string> bindings;
  switch (operation.kind) {
    case sol::OperationKind::LiquidLiquidExtraction:
      bindings.emplace_back("operation requires a distinct water-immiscible organic phase");
      break;
    case sol::OperationKind::Recrystallisation:
      bindings.emplace_back("solvent must stay liquid from the cold endpoint through the hot endpoint");
      break;
    case sol::OperationKind::Trituration:
      bindings.emplace_back("KEEP species must remain while REMOVE species dissolve");
      break;
    case sol::OperationKind::AntiSolventPrecipitation:
      bindings.emplace_back("primary solvent and anti-solvent must form a miscible pair");
      break;
    case sol::OperationKind::ChromatographyMobilePhase:
      bindings.emplace_back("KEEP and REMOVE species need a usable polarity separation window");
      break;
    case sol::OperationKind::ReactionMedium:
      bindings.emplace_back("every entered species must dissolve in the reaction medium");
      break;
  }
  if (operation.requireWaterImmiscible) bindings.emplace_back("water-immiscible only");
  if (operation.requireWaterMiscible) bindings.emplace_back("water-miscible only");
  if (operation.avoidPeroxideFormers) bindings.emplace_back("peroxide formers excluded");
  if (operation.avoidChlorinated) bindings.emplace_back("chlorinated solvents excluded");
  if (operation.avoidAromatics) bindings.emplace_back("aromatic solvents excluded");
  if (operation.excludeUnrated) bindings.emplace_back("CHEM21-unrated solvents excluded");
  if (!operation.worstAcceptableClass.empty()) {
    bindings.emplace_back("CHEM21 class no worse than " + operation.worstAcceptableClass);
  }
  if (operation.minBoilingPointC > 0.0) {
    bindings.emplace_back("minimum boiling point " +
                          std::to_string(static_cast<int>(operation.minBoilingPointC)) + " C");
  }
  if (operation.maxBoilingPointC > 0.0) {
    bindings.emplace_back("maximum boiling point " +
                          std::to_string(static_cast<int>(operation.maxBoilingPointC)) + " C");
  }
  return bindings;
}

std::string relaxMostRestrictive(sol::OperationSpec& operation) {
  if (operation.excludeUnrated) {
    operation.excludeUnrated = false;
    return "Include unrated solvents";
  }
  if (!operation.worstAcceptableClass.empty()) {
    operation.worstAcceptableClass.clear();
    return "Remove CHEM21 class limit";
  }
  if (operation.requireWaterImmiscible) {
    operation.requireWaterImmiscible = false;
    return "Allow water-miscible solvents";
  }
  if (operation.requireWaterMiscible) {
    operation.requireWaterMiscible = false;
    return "Allow water-immiscible solvents";
  }
  if (operation.avoidChlorinated) {
    operation.avoidChlorinated = false;
    return "Allow chlorinated solvents";
  }
  if (operation.avoidAromatics) {
    operation.avoidAromatics = false;
    return "Allow aromatic solvents";
  }
  if (operation.avoidPeroxideFormers) {
    operation.avoidPeroxideFormers = false;
    return "Allow peroxide formers";
  }
  if (operation.minBoilingPointC > 0.0 || operation.maxBoilingPointC > 0.0) {
    operation.minBoilingPointC = 0.0;
    operation.maxBoilingPointC = 0.0;
    return "Remove boiling-point window";
  }
  return {};
}

void drawComputingState() {
  if (!widgets::beginCard("##computing", ImVec2(0.0f, 0.0f), style::col::BgRaised)) return;
  widgets::badge("COMPUTING", style::col::Violet);
  ImGui::Spacing();
  const bool heading = style::pushFont(style::fonts::semibold());
  ImGui::TextUnformatted("Ranking the solvent space");
  style::popFont(heading);
  ImGui::TextWrapped("Evaluating recovery, separation, CHEM21 ratings and practical handling.");

  const float width = std::max(ImGui::GetContentRegionAvail().x, 1.0f);
  const float height = std::max(style::metrics().hairline * 7.0f, ImGui::GetFontSize() * 0.55f);
  const ImVec2 minimum = ImGui::GetCursorScreenPos();
  ImGui::Dummy(ImVec2(width, height));
  const ImVec2 maximum(minimum.x + width, minimum.y + height);
  ImDrawList* draw = ImGui::GetWindowDrawList();
  draw->AddRectFilled(minimum, maximum, style::u32(style::col::BgPanel),
                      style::metrics().radiusSm);
  const float phase = static_cast<float>(std::fmod(ImGui::GetTime() * 0.75, 1.0));
  const float segment = width * 0.28f;
  const float start = minimum.x + (width + segment) * phase - segment;
  draw->PushClipRect(minimum, maximum, true);
  draw->AddRectFilled(ImVec2(start, minimum.y), ImVec2(start + segment, maximum.y),
                      style::u32(style::col::Violet), style::metrics().radiusSm);
  draw->PopClipRect();
  widgets::endCard();
}

void drawNoResultState(SelectionState& state) {
  if (!widgets::beginCard("##no_results", ImVec2(0.0f, 0.0f), style::col::BgRaised)) return;
  widgets::badge("NO MATCH", style::col::Accent);
  ImGui::Spacing();
  const bool heading = style::pushFont(style::fonts::semibold());
  ImGui::TextUnformatted("No solvent satisfies every constraint");
  style::popFont(heading);
  ImGui::TextWrapped("The chemistry may still be workable; the current hard filters leave no "
                     "candidate for the ranking model.");
  const std::vector<std::string> bindings = bindingConstraints(state.operation);
  if (!bindings.empty()) {
    widgets::sectionHeader("Binding constraints", style::col::Accent);
    for (const std::string& binding : bindings) ImGui::TextWrapped("- %s", binding.c_str());
  }
  sol::OperationSpec relaxed = state.operation;
  const std::string relaxation = relaxMostRestrictive(relaxed);
  if (!relaxation.empty() &&
      widgets::primaryButton(("Relax: " + relaxation).c_str(),
                             ImVec2(ImGui::GetContentRegionAvail().x, 0.0f))) {
    state.operation = std::move(relaxed);
    state.statusMessage = relaxation;
  }
  widgets::endCard();
}

void drawEmptyState(AppState& state) {
  if (!widgets::beginCard("##selector_empty", ImVec2(0.0f, 0.0f), style::col::BgRaised)) return;
  widgets::badge("START HERE", style::col::Accent);
  ImGui::Spacing();
  const bool heading = style::pushFont(style::fonts::semibold());
  ImGui::TextUnformatted("Turn a real separation into a ranked solvent shortlist");
  style::popFont(heading);
  ImGui::TextWrapped("Add the compounds you need to keep and remove; ChemCAD will rank solvents "
                     "and show every weighted reason behind the order.");
  if (widgets::primaryButton("Load an extraction example",
                             ImVec2(ImGui::GetContentRegionAvail().x, 0.0f))) {
    loadExample(state);
  }
  widgets::endCard();
}

void drawResults(AppState& state) {
  SelectionState& selection = state.selection;
  widgets::sectionHeader("Ranked solvents", style::col::Accent);

  if (selection.operation.species.empty()) {
    drawEmptyState(state);
    return;
  }
  if (!hasTarget(selection)) {
    if (widgets::beginCard("##needs_target", ImVec2(0.0f, 0.0f), style::col::BgRaised)) {
      widgets::badge("NEEDS KEEP", style::col::Teal);
      ImGui::Spacing();
      ImGui::TextWrapped("At least one species must be marked KEEP before solvents can be ranked.");
      widgets::endCard();
    }
    return;
  }
  if (selection.computing) {
    drawComputingState();
    return;
  }
  if (!selection.rankingError.empty()) {
    if (widgets::beginCard("##ranking_error", ImVec2(0.0f, 0.0f), style::col::BgSurface)) {
      widgets::badge("DATABASE ERROR", style::col::Danger);
      ImGui::Spacing();
      ImGui::TextWrapped("%s", selection.rankingError.c_str());
      widgets::endCard();
    }
    return;
  }
  if (selection.candidates.empty()) {
    drawNoResultState(selection);
    return;
  }

  drawResultControls(selection);
  static std::unordered_map<const SelectionState*, ResultAnimation> animations;
  ResultAnimation& animation = animations[&selection];
  syncResultAnimation(selection, animation);
  drawComparisonBoard(selection, animation);

  const std::vector<const sol::SolventCandidate*> visible = visibleCandidates(selection);
  ImGui::TextColored(style::col::TextDim, "%zu of %zu candidates", visible.size(),
                     selection.candidates.size());
  if (visible.empty()) {
    if (widgets::beginCard("##search_empty", ImVec2(0.0f, 0.0f), style::col::BgSurface)) {
      ImGui::TextWrapped("No solvent name matches this filter.");
      if (widgets::ghostButton("Clear search")) selection.resultSearch.clear();
      widgets::endCard();
    }
    return;
  }

  for (size_t index = 0; index < visible.size(); ++index) {
    const sol::SolventCandidate& candidate = *visible[index];
    const auto original = std::find_if(
        selection.candidates.begin(), selection.candidates.end(),
        [&](const sol::SolventCandidate& value) { return &value == &candidate; });
    const int engineRank = original == selection.candidates.end()
                               ? static_cast<int>(index) + 1
                               : static_cast<int>(std::distance(selection.candidates.begin(), original)) + 1;
    drawCandidateCard(state, candidate, static_cast<int>(index) + 1, engineRank, animation);
    if (index + 1 < visible.size()) ImGui::Spacing();
  }
}

}  // namespace

void drawSolventSelector(AppState& state) {
  requestRankingIfNeeded(state);

  const style::Metrics& metrics = style::metrics();
  const ImVec2 available = ImGui::GetContentRegionAvail();
  const bool sideBySide = available.x >= ImGui::GetFontSize() * 58.0f;
  const float preferredBuilder = ImGui::GetFontSize() * 24.0f;

  if (sideBySide) {
    const float builderWidth = std::min(preferredBuilder, available.x * 0.40f);
    ImGui::PushStyleVar(ImGuiStyleVar_ChildRounding, metrics.radiusLg);
    ImGui::PushStyleColor(ImGuiCol_ChildBg, style::col::BgPanel);
    if (ImGui::BeginChild("##selector_builder", ImVec2(builderWidth, 0.0f),
                          ImGuiChildFlags_Borders)) {
      drawBuilder(state);
    }
    ImGui::EndChild();
    ImGui::PopStyleColor();
    ImGui::PopStyleVar();

    ImGui::SameLine(0.0f, metrics.gap);
    if (ImGui::BeginChild("##selector_results", ImVec2(0.0f, 0.0f))) drawResults(state);
    ImGui::EndChild();
  } else {
    const float minimumResults = ImGui::GetFontSize() * 12.0f;
    const float desiredBuilder = ImGui::GetFontSize() * 34.0f;
    const float builderHeight = std::max(
        ImGui::GetFontSize() * 14.0f,
        std::min(desiredBuilder, std::max(available.y - minimumResults - metrics.gap,
                                         ImGui::GetFontSize() * 14.0f)));
    ImGui::PushStyleVar(ImGuiStyleVar_ChildRounding, metrics.radiusLg);
    ImGui::PushStyleColor(ImGuiCol_ChildBg, style::col::BgPanel);
    if (ImGui::BeginChild("##selector_builder", ImVec2(0.0f, builderHeight),
                          ImGuiChildFlags_Borders)) {
      drawBuilder(state);
    }
    ImGui::EndChild();
    ImGui::PopStyleColor();
    ImGui::PopStyleVar();

    ImGui::Spacing();
    if (ImGui::BeginChild("##selector_results", ImVec2(0.0f, 0.0f))) drawResults(state);
    ImGui::EndChild();
  }
}

}  // namespace chemcad::ui

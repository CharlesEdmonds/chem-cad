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
#include "ui/charts3d.hpp"
#include "ui/icons.hpp"
#include "ui/layout.hpp"
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
  const double value = std::clamp(score, 0.0, 1.0);
  if (value < 0.35) return style::col::Danger;
  if (value > 0.70) return style::col::Success;
  return style::col::Data;
}

ImVec4 chem21Colour(const std::string& value) {
  if (value == "recommended") return style::col::Success;
  if (value == "problematic") return style::col::DataBright;
  if (value == "hazardous" || value == "highly hazardous") return style::col::Danger;
  return style::col::TextDim;
}

std::string formatSolubility(double value) {
  char buffer[64];
  const double magnitude = std::fabs(value);
  if (!std::isfinite(value)) return "infinite";
  if (value == 0.0) return "0";
  if (magnitude < 1e-3 || magnitude >= 1e4) {
    std::snprintf(buffer, sizeof(buffer), "%.2e", value);
  } else {
    std::snprintf(buffer, sizeof(buffer), "%.4g", value);
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
    std::snprintf(buffer, sizeof(buffer), "%.1f / %.1f", candidate.solvent->boilingPoint,
                  candidate.partner->boilingPoint);
  } else {
    std::snprintf(buffer, sizeof(buffer), "%.1f", candidate.solvent->boilingPoint);
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

bool sliderDouble(icons::Icon icon, const char* label, const char* id, double& value,
                  float minimum, float maximum, const char* format) {
  float temporary = static_cast<float>(value);
  if (!widgets::glyphSlider(id, icon, label, temporary, minimum, maximum, format)) {
    return false;
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
  const float textWidth = std::max(width - padding * 2.0f, metrics.hairline);
  const std::string title = ellipsizeText(sol::operationName(kind), textWidth);
  draw->AddText(ImVec2(minimum.x + padding, minimum.y + padding), style::u32(style::col::Text),
                title.c_str());
  const float purposeY = minimum.y + padding + ImGui::GetFontSize() + metrics.gap * 0.45f;
  draw->AddText(ImGui::GetFont(), ImGui::GetFontSize(),
                ImVec2(minimum.x + padding, purposeY), style::u32(style::col::TextDim),
                sol::operationDescription(kind), nullptr, textWidth);
}

void drawOperationChooser(SelectionState& state) {
  widgets::sectionHeader("Operation", style::col::Data);
  const layout::Frame frame = layout::measure();
  const int columns = std::min(layout::columnsThatFit(frame, 12.0f),
                               static_cast<int>(kOperationKinds.size()));
  const float width = layout::columnWidth(frame, columns);
  for (size_t index = 0; index < kOperationKinds.size(); ++index) {
    if (index % static_cast<size_t>(columns) != 0) ImGui::SameLine(0.0f, frame.gap);
    drawOperationCard(kOperationKinds[index], state.operation.kind, width);
  }
}

void drawOperationConditions(SelectionState& state) {
  widgets::sectionHeader("Working conditions", style::col::Data);
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
      sliderDouble(icons::Icon::Thermometer, "Temperature", "##extraction_temperature",
                   state.operation.temperatureC, -20.0f, 120.0f, "%.1f C");
      sliderDouble(icons::Icon::Droplet, "Aqueous volume", "##aqueous_volume",
                   state.operation.aqueousVolumeMl, 1.0f, 1000.0f, "%.0f mL");
      sliderDouble(icons::Icon::Droplet, "Organic volume", "##organic_volume",
                   state.operation.organicVolumeMl, 1.0f, 1000.0f, "%.0f mL");
      if (state.operation.pH == sol::kAutoPH) {
        if (drawToggleChip("Self-buffered pH", true)) state.operation.pH = 7.0;
        ImGui::TextWrapped("The solute sets its own saturated-solution pH.");
      } else {
        if (drawToggleChip("Fixed pH", true)) state.operation.pH = sol::kAutoPH;
        sliderDouble(icons::Icon::Ph, "pH", "##extraction_ph", state.operation.pH,
                     0.0f, 14.0f, "%.1f");
      }
      break;
    case sol::OperationKind::Recrystallisation:
      sliderDouble(icons::Icon::Flame, "Dissolve hot", "##hot_temperature",
                   state.operation.hotTemperatureC, 20.0f, 200.0f, "%.1f C");
      sliderDouble(icons::Icon::Snowflake, "Crystallise cold", "##cold_temperature",
                   state.operation.coldTemperatureC, -20.0f, 80.0f, "%.1f C");
      break;
    case sol::OperationKind::AntiSolventPrecipitation:
      sliderDouble(icons::Icon::Thermometer, "Addition temperature",
                   "##antisolvent_temperature", state.operation.temperatureC,
                   -20.0f, 120.0f, "%.1f C");
      ImGui::TextWrapped("The engine sweeps compatible solvent / anti-solvent partners and reports "
                         "the partner fraction at the best precipitation window.");
      break;
    case sol::OperationKind::Trituration:
      sliderDouble(icons::Icon::Thermometer, "Wash temperature",
                   "##trituration_temperature", state.operation.temperatureC,
                   -20.0f, 120.0f, "%.1f C");
      break;
    case sol::OperationKind::ChromatographyMobilePhase:
      sliderDouble(icons::Icon::Thermometer, "Column temperature",
                   "##chromatography_temperature", state.operation.temperatureC,
                   0.0f, 80.0f, "%.1f C");
      ImGui::TextWrapped("The polarity and selectivity window is derived from the keep and remove "
                         "species; use the boiling-point window below for volatility limits.");
      break;
    case sol::OperationKind::ReactionMedium:
      sliderDouble(icons::Icon::Thermometer, "Reaction temperature",
                   "##reaction_temperature", state.operation.temperatureC,
                   -20.0f, 200.0f, "%.1f C");
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
    const float segmentWidth = std::max(
        (ImGui::GetContentRegionAvail().x - segmentGap) * 0.5f,
        style::metrics().hairline);
    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(segmentGap, ImGui::GetStyle().ItemSpacing.y));
    if (drawSegmentButton("Keep##role", species.keep, ImVec2(segmentWidth, 0.0f))) species.keep = true;
    ImGui::SameLine(0.0f, segmentGap);
    if (drawSegmentButton("Remove##role", !species.keep, ImVec2(segmentWidth, 0.0f))) species.keep = false;
    ImGui::PopStyleVar();

    sliderDouble(icons::Icon::Balance, "Importance", "##species_weight",
                 species.weight, 0.1f, 3.0f, "%.2f x");
    sliderDouble(icons::Icon::Balance, "Amount", "##species_amount",
                 species.amountMg, 0.0f, 5000.0f, "%.0f mg");
    widgets::endCard();
  }
  ImGui::PopID();
  return remove;
}

void drawSpeciesBuilder(AppState& state) {
  SelectionState& selection = state.selection;
  widgets::sectionHeader("Species", style::col::Data);
  if (!widgets::beginCard("##species_builder", ImVec2(0.0f, 0.0f), style::col::BgSurface)) return;

  const float available = ImGui::GetContentRegionAvail().x;
  float used = 0.0f;
  flowBeforeChip("From sketch", used, available);
  if (drawToggleChip("From sketch", selection.addMode == 0)) selection.addMode = 0;
  flowBeforeChip("SMILES", used, available);
  if (drawToggleChip("SMILES", selection.addMode == 1)) selection.addMode = 1;
  flowBeforeChip("Chemical name", used, available);
  if (drawToggleChip("Chemical name", selection.addMode == 2)) selection.addMode = 2;

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
    if (selection.nameLookupRunning) {
      widgets::statusDot("Resolving chemical name", true, style::col::Data);
    }
    if (widgets::onlyWhen(!selection.nameLookupRunning,
                          "Chemical name lookup is already in progress.") &&
        (enter || widgets::primaryButton("Resolve and add",
                                         ImVec2(ImGui::GetContentRegionAvail().x, 0.0f)))) {
      addNameSpecies(state);
    }
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
  const float listHeight =
      std::max(layout::pageHeight(), ImGui::GetTextLineHeightWithSpacing());
  const layout::Frame listFrame =
      layout::measure(ImVec2(ImGui::GetContentRegionAvail().x, listHeight));
  const float listWeights[2] = {1.0f, 0.0f};
  const float listMinimums[2] = {listFrame.row * 7.0f, listFrame.control};
  float listRows[2] = {};
  layout::distribute(listHeight, listWeights, listMinimums, 2, listFrame.gap, listRows);
  static std::unordered_map<const SelectionState*, int> pages;
  int& page = pages[&selection];
  const int pageCount =
      std::max(1, static_cast<int>(selection.operation.species.size()));
  page = std::clamp(page, 0, pageCount - 1);

  ImGui::PushStyleVar(ImGuiStyleVar_ChildRounding, style::metrics().radiusMd);
  ImGui::PushStyleColor(ImGuiCol_ChildBg, style::col::BgPanel);
  if (ImGui::BeginChild("##species_list", ImVec2(0.0f, listRows[0]),
                        ImGuiChildFlags_Borders,
                        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    if (selection.operation.species.empty()) {
      ImGui::TextWrapped("Add every compound that matters. Mark products to isolate as Keep and "
                         "impurities to reject as Remove.");
    } else if (drawSpeciesRow(selection, static_cast<size_t>(page))) {
      selection.operation.species.erase(selection.operation.species.begin() + page);
      selection.speciesPresentation.erase(selection.speciesPresentation.begin() + page);
      page = std::max(page - 1, 0);
    }
  }
  ImGui::EndChild();
  if (ImGui::BeginChild("##species_pages", ImVec2(0.0f, listRows[1]),
                        ImGuiChildFlags_None,
                        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    if (page > 0 && widgets::iconButton("##previous_species", icons::Icon::ChevronLeft,
                                        ImVec2(listFrame.control, listFrame.control), false,
                                        "Previous species")) {
      --page;
    }
    if (page > 0) ImGui::SameLine(0.0f, listFrame.gap);
    ImGui::TextColored(style::col::Data, "Species %d of %d",
                       selection.operation.species.empty() ? 0 : page + 1,
                       static_cast<int>(selection.operation.species.size()));
    if (page + 1 < static_cast<int>(selection.operation.species.size())) {
      ImGui::SameLine(0.0f, listFrame.gap);
      if (widgets::iconButton("##next_species", icons::Icon::ChevronRight,
                              ImVec2(listFrame.control, listFrame.control), false,
                              "Next species")) {
        ++page;
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

  flowBeforeChip("Water-immiscible", used, available);
  if (drawToggleChip("Water-immiscible", operation.requireWaterImmiscible)) {
    operation.requireWaterImmiscible = !operation.requireWaterImmiscible;
    if (operation.requireWaterImmiscible) operation.requireWaterMiscible = false;
  }
  flowBeforeChip("Water-miscible", used, available);
  if (drawToggleChip("Water-miscible", operation.requireWaterMiscible)) {
    operation.requireWaterMiscible = !operation.requireWaterMiscible;
    if (operation.requireWaterMiscible) operation.requireWaterImmiscible = false;
  }
  flowBeforeChip("No peroxide formers", used, available);
  if (drawToggleChip("No peroxide formers", operation.avoidPeroxideFormers)) {
    operation.avoidPeroxideFormers = !operation.avoidPeroxideFormers;
  }
  flowBeforeChip("No chlorinated", used, available);
  if (drawToggleChip("No chlorinated", operation.avoidChlorinated)) {
    operation.avoidChlorinated = !operation.avoidChlorinated;
  }
  flowBeforeChip("No aromatics", used, available);
  if (drawToggleChip("No aromatics", operation.avoidAromatics)) {
    operation.avoidAromatics = !operation.avoidAromatics;
  }
  flowBeforeChip("Rated only", used, available);
  if (drawToggleChip("Rated only", operation.excludeUnrated)) {
    operation.excludeUnrated = !operation.excludeUnrated;
  }
}

void drawClassSelector(sol::OperationSpec& operation) {
  static constexpr std::array<const char*, 5> kClasses = {
      "", "recommended", "problematic", "hazardous", "highly hazardous"};
  ImGui::TextColored(style::col::TextDim, "Worst acceptable CHEM21 class");
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
    widgets::badge("Unrestricted", style::col::TextDim);
  } else {
    widgets::badge(operation.worstAcceptableClass.c_str(),
                   chem21Colour(operation.worstAcceptableClass));
  }
}

void drawConstraints(SelectionState& state) {
  widgets::sectionHeader("Constraints", style::col::Data);
  if (!widgets::beginCard("##constraints", ImVec2(0.0f, 0.0f), style::col::BgSurface)) return;
  drawConstraintChips(state);
  ImGui::Spacing();
  sliderDouble(icons::Icon::Thermometer, "Minimum boiling point", "##minimum_boiling",
               state.operation.minBoilingPointC, 0.0f, 250.0f, "%.0f C");
  sliderDouble(icons::Icon::Thermometer, "Maximum boiling point", "##maximum_boiling",
               state.operation.maxBoilingPointC, 0.0f, 300.0f, "%.0f C");
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
  sliderDouble(icons::Icon::Balance, label, id, value, 0.0f, 3.0f, "%.2f x");
  ImGui::PushStyleColor(ImGuiCol_Text, style::col::TextDim);
  ImGui::TextWrapped("%s", explanation);
  ImGui::PopStyleColor();
}

void drawPriorities(SelectionState& state) {
  widgets::sectionHeader("Scoring priorities", style::col::Data);
  if (!widgets::beginCard("##priorities", ImVec2(0.0f, 0.0f), style::col::BgSurface)) return;
  drawWeightControl("Selectivity", "##weight_selectivity", state.operation.weightSelectivity,
                    "Raise this to reward separation between keep and remove species.");
  drawWeightControl("Recovery", "##weight_recovery", state.operation.weightRecovery,
                    "Raise this to favour recovering more of every keep species.");
  drawWeightControl("Greenness", "##weight_greenness", state.operation.weightGreenness,
                    "Raise this to favour better CHEM21 safety, health and environment ratings.");
  drawWeightControl("Practicality", "##weight_practicality", state.operation.weightPracticality,
                    "Raise this to favour workable boiling, freezing, cost and handling properties.");
  widgets::endCard();
}

void drawBuilder(AppState& state) {
  SelectionState& selection = state.selection;
  static std::unordered_map<const SelectionState*, int> tabs;
  static std::unordered_map<const SelectionState*, int> conditionViews;
  int& tab = tabs[&selection];
  int& conditionView = conditionViews[&selection];

  constexpr std::array<const char*, 3> labels = {"Operation", "Conditions", "Species"};
  constexpr std::array<icons::Icon, 3> glyphs = {
      icons::Icon::Reaction, icons::Icon::Thermometer, icons::Icon::Molecule};

  const layout::Frame frame = layout::measure();
  const bool showStatus = !selection.statusMessage.empty();
  const float weights[3] = {0.0f, 1.0f, 0.0f};
  const float minimums[3] = {
      frame.control, frame.row * 5.0f, showStatus ? frame.row * 2.0f : 0.0f};
  float heights[3] = {};
  const int rowCount = showStatus ? 3 : 2;
  layout::distribute(frame.size.y, weights, minimums, rowCount, frame.gap, heights);

  if (ImGui::BeginChild("##builder_tabs", ImVec2(0.0f, heights[0]),
                        ImGuiChildFlags_None,
                        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    widgets::subTabs("##selector_tabs", labels.data(), glyphs.data(),
                     static_cast<int>(labels.size()), tab);
  }
  ImGui::EndChild();

  if (ImGui::BeginChild("##builder_content", ImVec2(0.0f, heights[1]),
                        ImGuiChildFlags_None,
                        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    if (tab == 0) {
      drawOperationChooser(selection);
    } else if (tab == 1) {
      constexpr std::array<icons::Icon, 3> conditionGlyphs = {
          icons::Icon::Thermometer, icons::Icon::Filter, icons::Icon::Balance};
      constexpr std::array<const char*, 3> conditionTips = {
          "Working conditions", "Hard constraints", "Scoring priorities"};
      widgets::segmentedIcons("##condition_view", conditionGlyphs.data(),
                              conditionTips.data(), static_cast<int>(conditionGlyphs.size()),
                              conditionView, ImGui::GetContentRegionAvail().x);
      if (conditionView == 0) {
        drawOperationConditions(selection);
      } else if (conditionView == 1) {
        drawConstraints(selection);
      } else {
        drawPriorities(selection);
      }
    } else {
      drawSpeciesBuilder(state);
    }
  }
  ImGui::EndChild();

  if (showStatus) {
    if (ImGui::BeginChild("##builder_status", ImVec2(0.0f, heights[2]),
                          ImGuiChildFlags_None,
                          ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
      widgets::notice(icons::Icon::Check, selection.statusMessage.c_str(), style::col::Success);
    }
    ImGui::EndChild();
  }
}

void drawBar(const char* id, double value, ImVec4 colour, float height, float reveal = 1.0f) {
  const float width = std::max(ImGui::GetContentRegionAvail().x, style::metrics().hairline);
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


void drawArithmetic(const sol::SolventCandidate& candidate) {
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

void drawCandidateEvidence(AppState& state, const sol::SolventCandidate& candidate,
                           int engineRank, ResultAnimation& animation) {
  SelectionState& selection = state.selection;
  const std::string key = candidateKey(candidate);
  const bool pinned = isCompared(selection, key);
  ImGui::PushID(key.c_str());
  if (!widgets::beginCard("##candidate_evidence", ImVec2(0.0f, 0.0f),
                          style::col::BgSurface)) {
    ImGui::PopID();
    return;
  }

  const bool heading = style::pushFont(style::fonts::semibold());
  ImGui::TextUnformatted(candidateTitle(candidate).c_str());
  style::popFont(heading);
  ImGui::SameLine();
  widgets::badge(candidate.estimated ? "Estimated" : "Measured",
                 candidate.estimated ? style::col::Violet : style::col::Success);

  const int delta = priorRankDelta(animation, key, engineRank);
  if (delta != 0 && resultReveal(animation) < 1.0f) {
    char movement[32];
    std::snprintf(movement, sizeof(movement), "%s %d", delta > 0 ? "Up" : "Down",
                  std::abs(delta));
    widgets::badge(movement, delta > 0 ? style::col::Success : style::col::Danger);
  }

  const float reveal = resultReveal(animation);
  for (const sol::Criterion& criterion : candidate.criteria) drawCriterion(criterion, reveal);
  if (!candidate.warnings.empty()) {
    for (const std::string& warning : candidate.warnings) {
      widgets::notice(icons::Icon::Warning, warning.c_str(), style::col::Danger);
    }
  }
  if (widgets::disclosure("##ranking_arithmetic", "Ranking arithmetic", "",
                          false, icons::Icon::Info, style::col::Accent)) {
    drawArithmetic(candidate);
  }

  const layout::Frame frame = layout::measure();
  const int actionCount = selection.compareMode ? 3 : 2;
  const int columns = std::min(layout::columnsThatFit(frame, 12.0f), actionCount);
  const float actionWidth = layout::columnWidth(frame, columns);
  int actionIndex = 0;
  auto nextAction = [&] {
    if (actionIndex > 0 && actionIndex % columns != 0) ImGui::SameLine(0.0f, frame.gap);
    ++actionIndex;
  };

  nextAction();
  if (widgets::onlyWhen(candidate.partner != nullptr,
                        "Extraction hand-off requires a two-solvent candidate.") &&
      widgets::actionButton("##send_to_extraction", icons::Icon::SepFunnel,
                            "Send to Extraction", ImVec2(actionWidth, 0.0f), false,
                            "Load this pair into the extraction calculator")) {
    sendToExtraction(state, candidate);
  }

  nextAction();
  if (widgets::actionButton("##load_into_predict", icons::Icon::Flask,
                            "Load into Predict", ImVec2(actionWidth, 0.0f), true,
                            "Load this blend into the prediction workspace")) {
    loadIntoSuite(state, candidate);
  }
  if (selection.compareMode) {
    nextAction();
    const bool comparisonFull = selection.comparedCandidates.size() >= 3 && !pinned;
    if (widgets::onlyWhen(!comparisonFull,
                          "The comparison already contains three candidates.") &&
        widgets::actionButton("##toggle_comparison", icons::Icon::ChartLine,
                              pinned ? "Unpin comparison" : "Pin for comparison",
                              ImVec2(actionWidth, 0.0f), false)) {
      toggleCompared(selection, key);
    }
  }

  widgets::endCard();
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

void drawComparisonBoard(SelectionState& state) {
  if (!state.compareMode) return;
  pruneCompared(state);
  widgets::sectionHeader("Pinned comparison", style::col::Data);
  if (!widgets::onlyWhen(!state.comparedCandidates.empty(),
                         "Pin up to three candidates from the ranked table to compare them.")) {
    return;
  }

  constexpr std::array<widgets::Column, 6> columns = {{
      {"Solvent", false, true, nullptr, 10.0f},
      {"Score", true, false, nullptr, 0.0f},
      {"Selectivity", true, false, nullptr, 0.0f},
      {"Recovery", true, false, nullptr, 0.0f},
      {"Greenness", true, false, nullptr, 0.0f},
      {"", false, false, nullptr, 0.0f},
  }};
  std::string removeKey;
  const layout::Frame frame = layout::measure();
  if (widgets::beginDataTable("##comparison_table", columns.data(),
                              static_cast<int>(columns.size()),
                              ImVec2(frame.size.x, layout::pageHeight()))) {
    for (const std::string& key : state.comparedCandidates) {
      const sol::SolventCandidate* candidate = findCandidate(state, key);
      if (!candidate) continue;
      widgets::dataRow(style::col::Accent);
      widgets::dataCell(candidateTitle(*candidate).c_str());
      widgets::dataCellf("%.3f", candidate->score);
      widgets::dataCell(formatRatio(candidate->selectivity).c_str());
      widgets::dataCell(formatPercent(candidate->recoveryFraction).c_str());
      widgets::dataCellf("%.2f", criterionScore(*candidate, "Greenness"));
      ImGui::TableNextColumn();
      ImGui::PushID(key.c_str());
      if (widgets::iconButton("##unpin", icons::Icon::Close,
                              ImVec2(frame.control, frame.control), false,
                              "Unpin candidate")) {
        removeKey = key;
      }
      ImGui::PopID();
    }
    widgets::endDataTable();
  }
  if (!removeKey.empty()) toggleCompared(state, removeKey);
}

void drawResultControls(SelectionState& state) {
  const layout::Frame frame = layout::measure();
  const int columns = std::min(layout::columnsThatFit(frame, 13.0f), 3);
  const float width = layout::columnWidth(frame, columns);

  ImGui::SetNextItemWidth(width);
  widgets::stringInputWithHint("##solvent_search", "Filter solvent names", state.resultSearch);
  if (columns > 1) ImGui::SameLine(0.0f, frame.gap);
  ImGui::SetNextItemWidth(width);
  ImGui::Combo("##result_sort", &state.sortMode, kSortLabels.data(),
               static_cast<int>(kSortLabels.size()));
  if (columns > 2) ImGui::SameLine(0.0f, frame.gap);
  if (drawToggleChip("Compare up to 3", state.compareMode)) {
    state.compareMode = !state.compareMode;
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

const sol::Solute* targetSolute(const SelectionState& state) {
  for (const sol::SpeciesRole& species : state.operation.species) {
    if (species.keep) return &species.solute;
  }
  return nullptr;
}

bool hasHansenParameters(const sol::Solute* solute) {
  return solute && std::isfinite(solute->hansen.dispersion) &&
         std::isfinite(solute->hansen.polar) &&
         std::isfinite(solute->hansen.hydrogenBond) &&
         std::isfinite(solute->interactionRadius) &&
         solute->hansen.dispersion > 0.0 && solute->interactionRadius > 0.0;
}

double candidateBoilingPoint(const sol::SolventCandidate& candidate) {
  if (!candidate.solvent || candidate.solvent->boilingPoint <= 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  return candidate.solvent->boilingPoint;
}

std::array<double, 3> candidateHansen(const sol::SolventCandidate& candidate) {
  if (!candidate.solvent) return {};
  const sol::Hansen& primary = candidate.solvent->hansen;
  if (!candidate.partner) {
    return {primary.dispersion, primary.polar, primary.hydrogenBond};
  }
  const double partnerFraction = std::clamp(candidate.partnerFraction, 0.0, 1.0);
  const double primaryFraction = 1.0 - partnerFraction;
  const sol::Hansen& partner = candidate.partner->hansen;
  return {
      primary.dispersion * primaryFraction + partner.dispersion * partnerFraction,
      primary.polar * primaryFraction + partner.polar * partnerFraction,
      primary.hydrogenBond * primaryFraction + partner.hydrogenBond * partnerFraction,
  };
}

void drawCandidateOverview(const SelectionState& state,
                           const std::vector<const sol::SolventCandidate*>& candidates,
                           std::string& selectedKey, float height) {
  static std::unordered_map<const SelectionState*, int> views;
  static std::unordered_map<const SelectionState*, charts3d::Orbit> orbits;
  int& view = views.try_emplace(&state, 1).first->second;
  charts3d::Orbit& orbit = orbits[&state];

  const layout::Frame frame =
      layout::measure(ImVec2(ImGui::GetContentRegionAvail().x, height));
  const float weights[2] = {0.0f, 1.0f};
  const float minimums[2] = {frame.control, frame.row * 7.0f};
  float rows[2] = {};
  layout::distribute(height, weights, minimums, 2, frame.gap, rows);

  const sol::Solute* solute = targetSolute(state);
  const bool hansenAvailable = hasHansenParameters(solute);
  if (ImGui::BeginChild("##overview_switch", ImVec2(0.0f, rows[0]),
                        ImGuiChildFlags_None,
                        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    if (hansenAvailable) {
      constexpr std::array<icons::Icon, 2> glyphs = {
          icons::Icon::ChartScatter, icons::Icon::ChartLine};
      constexpr std::array<const char*, 2> tips = {
          "Hansen space", "Parallel coordinates"};
      widgets::segmentedIcons("##candidate_view", glyphs.data(), tips.data(),
                              static_cast<int>(glyphs.size()), view);
    } else {
      view = 1;
      widgets::onlyWhen(false,
                        "Hansen space is unavailable because the target has no interaction sphere.");
    }
  }
  ImGui::EndChild();

  if (!ImGui::BeginChild("##overview_chart", ImVec2(0.0f, rows[1]),
                         ImGuiChildFlags_None,
                         ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    ImGui::EndChild();
    return;
  }

  std::vector<std::string> labels;
  labels.reserve(candidates.size());
  for (const sol::SolventCandidate* candidate : candidates) {
    labels.push_back(candidateTitle(*candidate));
  }

  if (view == 0 && hansenAvailable) {
    std::vector<charts3d::CloudPoint> points;
    points.reserve(candidates.size());
    for (size_t index = 0; index < candidates.size(); ++index) {
      const sol::SolventCandidate& candidate = *candidates[index];
      const std::array<double, 3> hansen = candidateHansen(candidate);
      const bool selected = candidateKey(candidate) == selectedKey;
      points.push_back({hansen[0], hansen[1], hansen[2], labels[index].c_str(),
                        selected ? style::col::Accent : style::col::Data,
                        selected ? 1.35f : 1.0f, selected});
    }
    charts3d::CloudStyle cloudStyle;
    cloudStyle.xLabel = "dD";
    cloudStyle.yLabel = "dP";
    cloudStyle.zLabel = "dH";
    cloudStyle.hasSphere = true;
    cloudStyle.sphereX = solute->hansen.dispersion;
    cloudStyle.sphereY = solute->hansen.polar;
    cloudStyle.sphereZ = solute->hansen.hydrogenBond;
    cloudStyle.sphereRadius = solute->interactionRadius;
    cloudStyle.sphereColour = style::col::DataBright;
    const int hovered = charts3d::cloud(
        "##hansen_cloud", points.data(), static_cast<int>(points.size()),
        ImVec2(frame.size.x, rows[1]), orbit, cloudStyle);
    if (hovered >= 0 && hovered < static_cast<int>(candidates.size())) {
      selectedKey = candidateKey(*candidates[static_cast<size_t>(hovered)]);
    }
  } else {
    std::vector<std::array<double, 5>> values;
    values.reserve(candidates.size());
    std::array<double, 5> minimumsByAxis;
    std::array<double, 5> maximumsByAxis;
    minimumsByAxis.fill(std::numeric_limits<double>::infinity());
    maximumsByAxis.fill(-std::numeric_limits<double>::infinity());
    for (const sol::SolventCandidate* candidate : candidates) {
      values.push_back({
          candidate->selectivity,
          candidate->recoveryFraction * 100.0,
          candidate->targetSolubilityGPerMl,
          candidate->contaminantSolubilityGPerMl,
          candidateBoilingPoint(*candidate),
      });
      for (size_t axis = 0; axis < values.back().size(); ++axis) {
        const double value = values.back()[axis];
        if (!std::isfinite(value)) continue;
        minimumsByAxis[axis] = std::min(minimumsByAxis[axis], value);
        maximumsByAxis[axis] = std::max(maximumsByAxis[axis], value);
      }
    }
    for (size_t axis = 0; axis < minimumsByAxis.size(); ++axis) {
      if (!std::isfinite(minimumsByAxis[axis])) minimumsByAxis[axis] = 0.0;
      if (!std::isfinite(maximumsByAxis[axis])) maximumsByAxis[axis] = 1.0;
    }

    const std::array<charts3d::ParallelAxis, 5> axes = {{
        {"Selectivity", minimumsByAxis[0], maximumsByAxis[0], true, "x"},
        {"Recovery", minimumsByAxis[1], maximumsByAxis[1], true, "%"},
        {"Target sol.", minimumsByAxis[2], maximumsByAxis[2], true, "g/mL"},
        {"Contam. sol.", minimumsByAxis[3], maximumsByAxis[3], false, "g/mL"},
        {"Primary bp", minimumsByAxis[4], maximumsByAxis[4], false, "C"},
    }};
    std::vector<charts3d::ParallelSeries> series;
    series.reserve(candidates.size());
    for (size_t index = 0; index < candidates.size(); ++index) {
      const bool selected = candidateKey(*candidates[index]) == selectedKey;
      series.push_back({labels[index].c_str(), values[index].data(),
                        selected ? style::col::Accent : style::col::DataDim, selected});
    }
    const int hovered = charts3d::parallelCoordinates(
        "##candidate_parallel", axes.data(), static_cast<int>(axes.size()),
        series.data(), static_cast<int>(series.size()), ImVec2(frame.size.x, rows[1]));
    if (hovered >= 0 && hovered < static_cast<int>(candidates.size())) {
      selectedKey = candidateKey(*candidates[static_cast<size_t>(hovered)]);
    }
  }
  ImGui::EndChild();
}


void drawCandidateTable(const SelectionState& state,
                        const std::vector<const sol::SolventCandidate*>& candidates,
                        std::string& selectedKey, std::string& detailKey, float height) {
  const layout::Frame frame =
      layout::measure(ImVec2(ImGui::GetContentRegionAvail().x, height));
  const bool compact = frame.density == layout::Density::Compact;
  const std::array<widgets::Column, 9> columns = {{
      {compact ? "#" : "Rank", true, false, nullptr, 0.0f},
      {"Solvent", false, true, nullptr, compact ? 7.0f : 11.0f},
      {"Score", true, false, nullptr, 0.0f},
      {compact ? "Sel." : "Selectivity", true, false, nullptr, 0.0f},
      {compact ? "Rec." : "Recovery", true, false, "%", 0.0f},
      {compact ? "T sol." : "Target solubility", true, false, "g/mL", 0.0f},
      {compact ? "C sol." : "Contaminant solubility", true, false, "g/mL", 0.0f},
      {compact ? "Bp" : "Boiling point", true, false, "C", 0.0f},
      {compact ? "Detail" : "Evidence", false, false, nullptr, 0.0f},
  }};
  static std::unordered_map<const SelectionState*, int> pages;
  int& page = pages[&state];
  const float weights[2] = {1.0f, 0.0f};
  const float minimums[2] = {frame.row * 3.0f, frame.control};
  float rows[2] = {};
  layout::distribute(height, weights, minimums, 2, frame.gap, rows);
  const int rowsPerPage = std::max(
      1, static_cast<int>(std::floor(rows[0] / std::max(frame.row, frame.em))) - 1);
  const int pageCount = std::max(
      1, (static_cast<int>(candidates.size()) + rowsPerPage - 1) / rowsPerPage);
  page = std::clamp(page, 0, pageCount - 1);
  const size_t begin = static_cast<size_t>(page * rowsPerPage);
  const size_t end = std::min(candidates.size(), begin + static_cast<size_t>(rowsPerPage));

  if (widgets::beginDataTable("##ranked_candidates", columns.data(),
                              static_cast<int>(columns.size()),
                              ImVec2(frame.size.x, rows[0]))) {
    for (size_t index = begin; index < end; ++index) {
      const sol::SolventCandidate& candidate = *candidates[index];
      const std::string key = candidateKey(candidate);
      const bool selected = key == selectedKey;
      widgets::dataRow(selected ? style::col::Accent : style::col::DataDim);
      widgets::dataCellf("%zu", index + 1);
      widgets::dataCell(candidateTitle(candidate).c_str());
      if (ImGui::IsItemClicked()) selectedKey = key;
      widgets::dataCellf("%.3f", candidate.score);
      widgets::dataCell(formatRatio(candidate.selectivity).c_str());
      widgets::dataCellf("%.1f", candidate.recoveryFraction * 100.0);
      widgets::dataCell(formatSolubility(candidate.targetSolubilityGPerMl).c_str());
      widgets::dataCell(formatSolubility(candidate.contaminantSolubilityGPerMl).c_str());
      widgets::dataCell(formatBoilingPoint(candidate).c_str());
      ImGui::TableNextColumn();
      ImGui::PushID(key.c_str());
      const char* summary = candidate.estimated ? "estimated" : "measured";
      const bool open = widgets::disclosure(
          "##evidence", "Details", summary, false, icons::Icon::Info, style::col::Accent);
      if (open) {
        detailKey = key;
      } else if (detailKey == key) {
        detailKey.clear();
      }
      ImGui::PopID();
    }
    widgets::endDataTable();
  }

  if (ImGui::BeginChild("##candidate_pages", ImVec2(0.0f, rows[1]),
                        ImGuiChildFlags_None,
                        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    if (page > 0 && widgets::iconButton(
                        "##previous_page", icons::Icon::ChevronLeft,
                        ImVec2(frame.control, frame.control), false,
                        "Previous candidates")) {
      --page;
    }
    if (page > 0) ImGui::SameLine(0.0f, frame.gap);
    ImGui::TextColored(style::col::Data, "Page %d of %d", page + 1, pageCount);
    if (page + 1 < pageCount) {
      ImGui::SameLine(0.0f, frame.gap);
      if (widgets::iconButton("##next_page", icons::Icon::ChevronRight,
                              ImVec2(frame.control, frame.control), false,
                              "Next candidates")) {
        ++page;
      }
    }
  }
  ImGui::EndChild();
}

std::vector<std::string> bindingConstraints(const sol::OperationSpec& operation) {
  std::vector<std::string> bindings;
  switch (operation.kind) {
    case sol::OperationKind::LiquidLiquidExtraction:
      bindings.emplace_back("operation requires a distinct water-immiscible organic phase");
      break;
    case sol::OperationKind::Recrystallisation:
      bindings.emplace_back(
          "solvent must stay liquid from the cold endpoint through the hot endpoint");
      break;
    case sol::OperationKind::Trituration:
      bindings.emplace_back("keep species must remain while remove species dissolve");
      break;
    case sol::OperationKind::AntiSolventPrecipitation:
      bindings.emplace_back("primary solvent and anti-solvent must form a miscible pair");
      break;
    case sol::OperationKind::ChromatographyMobilePhase:
      bindings.emplace_back(
          "keep and remove species need a usable polarity separation window");
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
    bindings.emplace_back(
        "minimum boiling point " +
        std::to_string(static_cast<int>(operation.minBoilingPointC)) + " C");
  }
  if (operation.maxBoilingPointC > 0.0) {
    bindings.emplace_back(
        "maximum boiling point " +
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
  widgets::statusDot("Computing", true, style::col::Data);
  const bool heading = style::pushFont(style::fonts::semibold());
  ImGui::TextUnformatted("Ranking the solvent space");
  style::popFont(heading);
  ImGui::TextWrapped("Evaluating recovery, separation, CHEM21 ratings and practical handling.");
  const float width =
      std::max(ImGui::GetContentRegionAvail().x, style::metrics().hairline);
  const float height =
      std::max(style::metrics().hairline * 7.0f, ImGui::GetFontSize() * 0.55f);
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
                      style::u32(style::col::Data), style::metrics().radiusSm);
  draw->PopClipRect();
  widgets::endCard();
}

void drawNoResultState(SelectionState& state) {
  if (!widgets::beginCard("##no_results", ImVec2(0.0f, 0.0f), style::col::BgRaised)) return;
  widgets::badge("No match", style::col::Danger);
  const bool heading = style::pushFont(style::fonts::semibold());
  ImGui::TextUnformatted("No solvent satisfies every constraint");
  style::popFont(heading);
  ImGui::TextWrapped("The chemistry may still be workable; the current hard filters leave no "
                     "candidate for the ranking model.");
  const std::vector<std::string> bindings = bindingConstraints(state.operation);
  if (!bindings.empty()) {
    widgets::sectionHeader("Binding constraints", style::col::Data);
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
  widgets::badge("Start here", style::col::Data);
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
  if (selection.operation.species.empty()) {
    widgets::sectionHeader("Ranked solvents", style::col::Data);
    drawEmptyState(state);
    return;
  }
  if (!hasTarget(selection)) {
    widgets::sectionHeader("Ranked solvents", style::col::Data);
    widgets::emptyState(icons::Icon::Molecule, "Target required",
                        "Mark at least one species as Keep before ranking solvents.");
    return;
  }
  if (selection.computing) {
    widgets::sectionHeader("Ranked solvents", style::col::Data);
    drawComputingState();
    return;
  }
  if (!selection.rankingError.empty()) {
    widgets::notice(icons::Icon::Warning, selection.rankingError.c_str(), style::col::Danger);
    return;
  }
  if (selection.candidates.empty()) {
    widgets::sectionHeader("Ranked solvents", style::col::Data);
    drawNoResultState(selection);
    return;
  }

  static std::unordered_map<const SelectionState*, ResultAnimation> animations;
  static std::unordered_map<const SelectionState*, std::string> selectedCandidates;
  static std::unordered_map<const SelectionState*, std::string> detailCandidates;
  ResultAnimation& animation = animations[&selection];
  std::string& selectedKey = selectedCandidates[&selection];
  std::string& detailKey = detailCandidates[&selection];
  syncResultAnimation(selection, animation);

  const std::vector<const sol::SolventCandidate*> visible = visibleCandidates(selection);
  if (visible.empty()) {
    drawResultControls(selection);
    widgets::emptyState(icons::Icon::Filter, "No matching solvent",
                        "Clear or change the name filter to restore candidates.");
    if (widgets::ghostButton("Clear search")) selection.resultSearch.clear();
    return;
  }
  const auto selectedVisible = std::find_if(
      visible.begin(), visible.end(),
      [&](const sol::SolventCandidate* candidate) {
        return candidateKey(*candidate) == selectedKey;
      });
  if (selectedVisible == visible.end()) selectedKey = candidateKey(*visible.front());
  const sol::SolventCandidate* detailCandidate =
      detailKey.empty() ? nullptr : findCandidate(selection, detailKey);
  if (detailCandidate &&
      std::find(visible.begin(), visible.end(), detailCandidate) == visible.end()) {
    detailCandidate = nullptr;
    detailKey.clear();
  }

  const layout::Frame frame = layout::measure();
  float weights[5] = {};
  float minimums[5] = {};
  weights[0] = 0.0f;
  const int controlColumns = std::min(layout::columnsThatFit(frame, 13.0f), 3);
  minimums[0] = frame.control *
                (controlColumns == 1 ? 4.8f : (controlColumns == 2 ? 3.4f : 2.4f));
  weights[1] = 1.15f;
  minimums[1] = frame.row * 8.0f;
  weights[2] = 1.0f;
  minimums[2] = frame.row * 7.0f;
  int rowCount = 3;
  int comparisonRow = -1;
  int detailRow = -1;
  if (selection.compareMode && !selection.comparedCandidates.empty()) {
    comparisonRow = rowCount++;
    weights[comparisonRow] = 0.45f;
    minimums[comparisonRow] = frame.row * 5.0f;
  }
  if (detailCandidate) {
    detailRow = rowCount++;
    weights[detailRow] = 0.8f;
    minimums[detailRow] = frame.row * 12.0f;
  }
  float heights[5] = {};
  layout::distribute(frame.size.y, weights, minimums, rowCount, frame.gap, heights);

  if (ImGui::BeginChild("##result_controls", ImVec2(0.0f, heights[0]),
                        ImGuiChildFlags_None,
                        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    widgets::sectionHeader("Ranked solvents", style::col::Data);
    drawResultControls(selection);
    ImGui::TextColored(style::col::DataDim, "%zu of %zu candidates", visible.size(),
                       selection.candidates.size());
  }
  ImGui::EndChild();
  if (ImGui::BeginChild("##result_overview", ImVec2(0.0f, heights[1]),
                        ImGuiChildFlags_None,
                        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    drawCandidateOverview(selection, visible, selectedKey, heights[1]);
  }
  ImGui::EndChild();
  if (ImGui::BeginChild("##result_table", ImVec2(0.0f, heights[2]),
                        ImGuiChildFlags_None,
                        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    drawCandidateTable(selection, visible, selectedKey, detailKey, heights[2]);
  }
  ImGui::EndChild();

  if (comparisonRow >= 0) {
    if (ImGui::BeginChild("##result_comparison", ImVec2(0.0f, heights[comparisonRow]),
                          ImGuiChildFlags_None,
                          ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
      drawComparisonBoard(selection);
    }
    ImGui::EndChild();
  }
  if (detailRow >= 0 && detailCandidate) {
    const auto original = std::find_if(
        selection.candidates.begin(), selection.candidates.end(),
        [&](const sol::SolventCandidate& value) { return &value == detailCandidate; });
    const int engineRank =
        original == selection.candidates.end()
            ? 0
            : static_cast<int>(std::distance(selection.candidates.begin(), original)) + 1;
    if (ImGui::BeginChild("##result_detail", ImVec2(0.0f, heights[detailRow]),
                          ImGuiChildFlags_None,
                          ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
      drawCandidateEvidence(state, *detailCandidate, engineRank, animation);
    }
    ImGui::EndChild();
  }
}

}  // namespace

void drawSolventSelector(AppState& state) {
  requestRankingIfNeeded(state);

  const layout::Frame frame = layout::measure();
  const int columns = std::clamp(layout::columnsThatFit(frame, 20.0f), 3, 5);
  const float builderWidth = layout::columnWidth(frame, columns);
  const float resultsWidth = layout::columnWidth(frame, columns, columns - 1);

  ImGui::PushStyleVar(ImGuiStyleVar_ChildRounding, style::metrics().radiusLg);
  ImGui::PushStyleColor(ImGuiCol_ChildBg, style::col::BgPanel);
  if (ImGui::BeginChild("##selector_builder", ImVec2(builderWidth, frame.size.y),
                        ImGuiChildFlags_Borders,
                        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    drawBuilder(state);
  }
  ImGui::EndChild();
  ImGui::PopStyleColor();
  ImGui::PopStyleVar();

  ImGui::SameLine(0.0f, frame.gap);
  if (ImGui::BeginChild("##selector_results", ImVec2(resultsWidth, frame.size.y),
                        ImGuiChildFlags_None,
                        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    drawResults(state);
  }
  ImGui::EndChild();
}

}  // namespace chemcad::ui

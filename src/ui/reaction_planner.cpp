#include "imgui.h"

#include <algorithm>
#include <cstdio>
#include <cmath>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "chem/bridge.hpp"
#include "naming/naming.hpp"
#include "rxn/engine.hpp"
#include "ui/charts.hpp"
#include "ui/icons.hpp"
#include "ui/mol_thumb.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"
#include "ui/widgets.hpp"

namespace chemcad::ui {
namespace {

// Sized from the current font rather than fixed pixels: the app runs at a 1.25
// UI scale by default, and hardcoded widths clip their own button labels.
inline float uiScale() { return ImGui::GetFontSize() / 13.0f; }
inline float materialWidth() { return ImGui::GetFontSize() * 17.7f; }
inline float materialHeight() {
  return ImGui::GetFrameHeight() * 4.0f +
         ImGui::GetTextLineHeightWithSpacing() * 7.0f +
         style::metrics().gap * 4.0f;
}
inline ImVec2 routeThumbSize() {
  return ImVec2(ImGui::GetFontSize() * 6.75f,
                ImGui::GetTextLineHeightWithSpacing() * 3.25f);
}

std::string ellipsize(const std::string& text, float maxWidth) {
  static constexpr const char* suffix = "...";
  if (maxWidth <= 0.0f) return {};
  if (ImGui::CalcTextSize(text.c_str()).x <= maxWidth) return text;
  const float suffixWidth = ImGui::CalcTextSize(suffix).x;
  if (maxWidth < suffixWidth) return {};
  const float contentWidth = maxWidth - suffixWidth;

  size_t end = text.size();
  while (end > 0 && ImGui::CalcTextSize(text.data(), text.data() + end).x > contentWidth) {
    --end;
    while (end > 0 && (static_cast<unsigned char>(text[end]) & 0xC0u) == 0x80u) --end;
  }
  return text.substr(0, end) + suffix;
}

constexpr float kRouteEnterSeconds = 0.25f;
constexpr float kSearchMinimumSeconds = 0.45f;

struct RoutePreviewCache {
  std::string signature;
  std::unordered_map<std::string, core::Molecule> molecules;
  int selectedRoute = 0;
};

struct PlannerAnimationState {
  double searchStartedAt = -1.0;
  bool awaitingResults = false;
  std::unordered_map<std::string, double> routeFirstSeen;
};

std::string formatScore(double score) {
  char buffer[32];
  const double magnitude = std::abs(score);
  if (magnitude >= 1.0e6)
    std::snprintf(buffer, sizeof(buffer), "%.3g", score);
  else
    std::snprintf(buffer, sizeof(buffer), "%.1f", score);
  return buffer;
}

std::string routeName(const rxn::Route& route, size_t routeIndex) {
  std::string name = "Route " + std::to_string(routeIndex + 1);
  if (!route.steps.empty() && !route.steps.front().reactionName.empty())
    name += " — " + route.steps.front().reactionName;
  return name;
}

void provenanceBadge(bool ai) {
  const icons::Icon icon = ai ? icons::Icon::Sparkle : icons::Icon::Book;
  const ImVec4 colour = ai ? style::col::Violet : style::col::Teal;
  const float glyphSize = ImGui::GetFontSize();
  const ImVec2 cursor = ImGui::GetCursorScreenPos();
  ImGui::Dummy(ImVec2(glyphSize, ImGui::GetFrameHeight()));
  icons::draw(ImGui::GetWindowDrawList(), icon,
              ImVec2(cursor.x + glyphSize * 0.5f,
                     cursor.y + ImGui::GetFrameHeight() * 0.5f),
              glyphSize, style::u32(colour));
  ImGui::SameLine(0.0f, style::metrics().gap * 0.5f);
  widgets::badge(ai ? "AI" : "KNOWLEDGE BASE", colour);
}

float easeOutCubic(float t) {
  t = std::clamp(t, 0.0f, 1.0f);
  const float inverse = 1.0f - t;
  return 1.0f - inverse * inverse * inverse;
}


MaterialBox* materialAt(AppState& st, bool target, int startIndex) {
  if (target) return &st.planner.target;
  if (startIndex < 0 || startIndex >= static_cast<int>(st.planner.starts.size()))
    return nullptr;
  return &st.planner.starts[static_cast<size_t>(startIndex)];
}

void rebuildPreview(MaterialBox& box) {
  box.preview.clear();
  box.previewValid = false;
  box.error.clear();
  if (box.smiles.empty()) {
    box.status = Status::Idle;
    return;
  }
  try {
    box.preview = chem::fromSmiles(box.smiles);
    box.previewValid = true;
    box.status = Status::Ok;
  } catch (const chem::ChemError& error) {
    box.status = Status::Error;
    box.error = error.what();
  } catch (const std::exception& error) {
    box.status = Status::Error;
    box.error = error.what();
  }
}

bool openInSketch(AppState& st, const std::string& smiles) {
  if (smiles.empty()) return false;
  try {
    core::Molecule molecule = chem::fromSmiles(smiles);
    st.snapshot();
    st.doc.molecules.push_back(std::move(molecule));
    st.touch();
    st.tab = MainTab::Sketch;
    st.tabChangeRequested = true;
    return true;
  } catch (const chem::ChemError& error) {
    st.planner.error = std::string("Could not open structure: ") + error.what();
    return false;
  } catch (const std::exception& error) {
    st.planner.error = std::string("Could not open structure: ") + error.what();
    return false;
  }
}

std::string joined(const std::vector<std::string>& values) {
  std::string result;
  for (const std::string& value : values) {
    if (!result.empty()) result += ", ";
    result += value;
  }
  return result;
}

void appendSignature(std::string& signature, const std::string& value) {
  signature += std::to_string(value.size());
  signature += ':';
  signature += value;
}

void appendRouteSignature(std::string& signature, const rxn::Route& route) {
  signature += "|R";
  signature += std::to_string(route.steps.size());
  for (const rxn::Step& step : route.steps) {
    signature += "|S";
    appendSignature(signature, step.reactionName);
    appendSignature(signature, step.productSmiles);
    appendSignature(signature, step.conditions);
    appendSignature(signature, step.notes);
    signature += step.source == rxn::Step::Source::LLM ? 'A' : 'K';
    for (const std::string& smiles : step.reactantSmiles)
      appendSignature(signature, smiles);
    for (const std::string& smiles : step.sideProductSmiles)
      appendSignature(signature, smiles);
    for (const std::string& reagent : step.reagents) appendSignature(signature, reagent);
  }
}

std::string routeIdentity(const rxn::Route& route) {
  std::string identity;
  identity.reserve(route.steps.size() * 128);
  appendRouteSignature(identity, route);
  return identity;
}

std::string routesSignature(const std::vector<rxn::Route>& routes) {
  std::string signature;
  signature.reserve(routes.size() * 128);
  signature += std::to_string(routes.size());
  for (const rxn::Route& route : routes) appendRouteSignature(signature, route);
  return signature;
}

void cacheSmiles(RoutePreviewCache& cache, const std::string& smiles) {
  if (smiles.empty() || cache.molecules.contains(smiles)) return;
  core::Molecule molecule;
  try {
    molecule = chem::fromSmiles(smiles);
  } catch (...) {
    // An empty graph deliberately selects the thumbnail's invalid-SMILES placeholder.
  }
  cache.molecules.emplace(smiles, std::move(molecule));
}

void rebuildRouteCache(RoutePreviewCache& cache, const std::vector<rxn::Route>& routes,
                       const std::string& signature) {
  cache.signature = signature;
  cache.molecules.clear();
  cache.selectedRoute = 0;
  for (const rxn::Route& route : routes) {
    for (const rxn::Step& step : route.steps) {
      for (const std::string& smiles : step.reactantSmiles) cacheSmiles(cache, smiles);
      cacheSmiles(cache, step.productSmiles);
      for (const std::string& smiles : step.sideProductSmiles) cacheSmiles(cache, smiles);
    }
  }
}

const core::Molecule& cachedMolecule(const RoutePreviewCache& cache,
                                     const std::string& smiles) {
  static const core::Molecule empty;
  const auto found = cache.molecules.find(smiles);
  return found == cache.molecules.end() ? empty : found->second;
}

void lookupName(AppState& st, bool target, int startIndex, MaterialBox& box) {
  if (box.nameInput.empty() || box.status == Status::Loading) return;
  const std::string requestedName = box.nameInput;
  box.status = Status::Loading;
  box.error.clear();
  st.tasks.run<naming::Result>(
      [requestedName] { return naming::nameToSmiles(requestedName); },
      [&st, target, startIndex, requestedName](naming::Result result) {
        MaterialBox* current = materialAt(st, target, startIndex);
        if (!current || current->status != Status::Loading ||
            current->nameInput != requestedName)
          return;
        if (!result.ok) {
          current->status = Status::Error;
          current->error = result.error.empty() ? "Name lookup failed." : result.error;
          return;
        }
        current->smiles = std::move(result.value);
        current->label = requestedName;
        rebuildPreview(*current);
      });
}

bool materialBoxWidget(AppState& st, MaterialBox& box, bool target, int startIndex) {
  bool remove = false;
  const float liftReserve = 3.0f * uiScale();
  const ImGuiID hoverId = ImGui::GetID("##material_hover");
  const float previousHover = ImGui::GetStateStorage()->GetFloat(hoverId, 0.0f);
  ImGui::Dummy(ImVec2(0.0f, liftReserve * (1.0f - previousHover)));

  const float availableWidth = std::max(ImGui::GetContentRegionAvail().x, 1.0f);
  const ImVec2 cardSize(std::min(materialWidth(), availableWidth), materialHeight());
  if (!widgets::beginCard("##material_box", cardSize)) return false;

  const char* heading = target ? "PRODUCT" : "STARTING MATERIAL";
  const ImVec4 headingColor = target ? style::col::Accent : style::col::TextDim;
  const float fs = ImGui::GetFontSize();
  constexpr ImGuiTableFlags headerFlags =
      ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_NoSavedSettings |
      ImGuiTableFlags_NoPadOuterX;
  const int headerColumns = target ? 1 : 2;
  if (ImGui::BeginTable("##material_header", headerColumns, headerFlags)) {
    ImGui::TableSetupColumn("##material_label", ImGuiTableColumnFlags_WidthStretch);
    if (!target)
      ImGui::TableSetupColumn("##material_remove", ImGuiTableColumnFlags_WidthFixed,
                              fs * 5.5f);
    ImGui::TableNextColumn();
    const bool pushed = style::pushFont(style::fonts::semibold());
    ImGui::TextColored(headingColor, "%s", heading);
    style::popFont(pushed);
    if (!target) {
      ImGui::TableNextColumn();
      const bool onlyStart = st.planner.starts.size() <= 1;
      ImGui::BeginDisabled(onlyStart);
      ImGui::PushStyleVar(ImGuiStyleVar_DisabledAlpha, 0.35f);
      remove = widgets::actionButton(
          "##remove", icons::Icon::Trash, "Remove", ImVec2(0.0f, 0.0f), false,
          "Remove this starting material");
      ImGui::PopStyleVar();
      ImGui::EndDisabled();
    }
    ImGui::EndTable();
  }

  const float innerWidth = std::max(ImGui::GetContentRegionAvail().x, 1.0f);
  if (moleculeThumbButton(
          "##preview", box.preview,
          ImVec2(innerWidth, ImGui::GetTextLineHeightWithSpacing() * 3.3f)) &&
      box.previewValid) {
    openInSketch(st, box.smiles);
  }
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", box.previewValid ? "Open in the Sketch canvas"
                                             : "Structure preview");

  ImGui::SetNextItemWidth(-1.0f);
  const bool smilesEnter = widgets::stringInputWithHint(
      "##smiles", "SMILES", box.smiles, ImGuiInputTextFlags_EnterReturnsTrue, true);
  const bool smilesCommitted = smilesEnter || ImGui::IsItemDeactivatedAfterEdit();
  if (smilesCommitted) rebuildPreview(box);
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "SMILES string; Enter updates the preview");

  // The reservation and the button share one width, so the pair can never
  // disagree and push the action past the card edge in a narrow flow column.
  const float lookupWidth = ImGui::CalcTextSize("Look up").x + ImGui::GetFontSize() +
                            style::metrics().gap * 2.5f;
  ImGui::SetNextItemWidth(std::max(
      1.0f, ImGui::GetContentRegionAvail().x - lookupWidth -
                ImGui::GetStyle().ItemSpacing.x));
  const bool nameEnter = widgets::stringInputWithHint(
      "##name", "Chemical name", box.nameInput, ImGuiInputTextFlags_EnterReturnsTrue);
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "IUPAC or common name; resolved via OPSIN / PubChem");
  ImGui::SameLine();
  ImGui::BeginDisabled(box.status == Status::Loading || box.nameInput.empty());
  const bool lookupClicked = widgets::actionButton(
      "##lookup", icons::Icon::Search, "Look up", ImVec2(lookupWidth, 0.0f), false,
      "Resolve an IUPAC or common name via OPSIN / PubChem");
  ImGui::EndDisabled();
  if (nameEnter || lookupClicked) lookupName(st, target, startIndex, box);

  ImGui::BeginDisabled(st.doc.empty());
  if (widgets::actionButton("##from_sketch", icons::Icon::ArrowLeft, "From sketch",
                            ImVec2(0.0f, 0.0f), false,
                            "Copy the largest sketched molecule into this box")) {
    const int largest = st.doc.largestMoleculeIndex();
    if (largest >= 0) {
      try {
        box.smiles = chem::toSmiles(st.doc.molecules[static_cast<size_t>(largest)]);
        rebuildPreview(box);
      } catch (const chem::ChemError& error) {
        box.status = Status::Error;
        box.error = error.what();
      } catch (const std::exception& error) {
        box.status = Status::Error;
        box.error = error.what();
      }
    }
  }
  ImGui::EndDisabled();
  if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
    ImGui::SetTooltip("%s", "Copy the largest sketched molecule into this box");
  const char* statusText = nullptr;
  if (box.status == Status::Loading) {
    statusText = "looking up...";
  } else if (!box.label.empty()) {
    statusText = box.label.c_str();
  }
  if (statusText && ImGui::GetContentRegionAvail().x > ImGui::GetFontSize() * 2.0f) {
    ImGui::SameLine();
    const std::string fitted =
        ellipsize(statusText, std::max(ImGui::GetContentRegionAvail().x, 1.0f));
    ImGui::TextColored(style::col::TextDim, "%s", fitted.c_str());
    if (ImGui::IsItemHovered() && fitted != statusText)
      ImGui::SetTooltip("%s", statusText);
  }
  if (!box.error.empty()) {
    widgets::notice(icons::Icon::Warning, box.error.c_str(), style::col::Danger);
  }

  widgets::endCard();
  const ImVec2 cardMin = ImGui::GetItemRectMin();
  const ImVec2 cardMax = ImGui::GetItemRectMax();
  const bool hovered = ImGui::IsItemHovered();
  const float hover = widgets::hoverT(hoverId, hovered);
  ImDrawList* draw = ImGui::GetWindowDrawList();
  const style::Metrics& metrics = style::metrics();
  draw->AddRect(cardMin, cardMax,
                style::mix(style::col::Border, headingColor, hover * 0.85f),
                metrics.radiusMd, 0, metrics.hairline * (1.0f + hover));
  if (hover > 0.0f) {
    draw->AddLine(ImVec2(cardMin.x + metrics.gap, cardMin.y),
                  ImVec2(cardMax.x - metrics.gap, cardMin.y),
                  style::u32(headingColor, hover * 0.75f),
                  metrics.hairline * (1.0f + hover));
  }
  return remove;
}

void drawArrow(ImVec2 size, ImVec4 color, float thicknessScale = 1.0f,
               bool animated = false) {
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::Dummy(size);
  ImDrawList* dl = ImGui::GetWindowDrawList();
  const ImU32 colour = style::u32(color);
  const float y = min.y + size.y * 0.5f;
  const float pad = size.x * 0.08f;
  const float left = min.x + pad;
  const float right = min.x + size.x - pad;
  const float thickness =
      style::metrics().hairline * 2.0f * thicknessScale;
  dl->AddLine(ImVec2(left, y), ImVec2(right, y), colour, thickness);
  const float head = size.y * 0.28f;
  dl->AddTriangleFilled(ImVec2(right, y), ImVec2(right - head * 1.6f, y - head),
                        ImVec2(right - head * 1.6f, y + head), colour);
  if (animated && right > left) {
    const float progress =
        std::fmod(static_cast<float>(ImGui::GetTime()) * 0.85f, 1.0f);
    const float x = left + (right - left) * progress;
    dl->AddCircleFilled(
        ImVec2(x, y),
        std::max(style::metrics().hairline * 2.0f, ImGui::GetFontSize() * 0.18f),
                        style::u32(style::col::AccentHover));
  }
}

void drawSearchingIndicator() {
  const float width = std::max(ImGui::GetContentRegionAvail().x, 1.0f);
  const float height = ImGui::GetTextLineHeightWithSpacing();
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::Dummy(ImVec2(width, height));
  ImDrawList* draw = ImGui::GetWindowDrawList();
  const char* label = "Searching routes";
  const ImVec2 labelSize = ImGui::CalcTextSize(label);
  draw->AddText(ImVec2(min.x, min.y + (height - labelSize.y) * 0.5f),
                style::u32(style::col::TextDim), label);
  const float radius =
      std::max(style::metrics().hairline * 2.0f, ImGui::GetFontSize() * 0.14f);
  const float spacing = radius * 3.0f;
  const float dotsWidth = spacing * 2.0f + radius * 2.0f;
  const float firstX =
      min.x + std::max(radius, width - dotsWidth + radius);
  for (int i = 0; i < 3; ++i) {
    const float wave = 0.5f + 0.5f * std::sin(
        static_cast<float>(ImGui::GetTime()) * 7.0f - static_cast<float>(i) * 0.9f);
    draw->AddCircleFilled(ImVec2(firstX + spacing * static_cast<float>(i),
                                min.y + height * 0.5f - wave * radius * 0.8f),
                          radius, style::u32(style::col::Accent, 0.35f + wave * 0.65f));
  }
}

void dispatchSearch(AppState& st, PlannerAnimationState& animation) {
  rxn::Request request;
  for (const MaterialBox& start : st.planner.starts) {
    if (!start.smiles.empty()) request.startSmiles.push_back(start.smiles);
  }
  if (request.startSmiles.empty()) {
    st.planner.error = "Enter at least one starting material.";
    return;
  }
  if (st.planner.target.smiles.empty()) {
    st.planner.error = "Enter a product before suggesting routes.";
    return;
  }

  request.targetSmiles = st.planner.target.smiles;
  request.maxDepth = st.planner.maxDepth;
  request.maxRoutes = st.planner.maxRoutes;
  request.allowLlm = st.planner.allowLlm && rxn::llmAvailable();
  st.planner.error.clear();
  st.planner.searching = true;
  st.planner.searched = false;
  animation.searchStartedAt = ImGui::GetTime();
  animation.awaitingResults = true;

  st.tasks.run<std::vector<rxn::Route>>(
      [request] { return rxn::suggestRoutes(request); },
      [&st](std::vector<rxn::Route> routes) {
        st.planner.routes = std::move(routes);
        st.planner.searching = false;
        st.planner.searched = true;
      });
}

void connectorControls(AppState& st, PlannerAnimationState& animation) {
  ImGui::BeginGroup();
  const float width = std::max(ImGui::GetContentRegionAvail().x, 1.0f);
  constexpr ImGuiTableFlags labelFlags =
      ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_NoSavedSettings |
      ImGuiTableFlags_NoPadOuterX;
  if (ImGui::BeginTable("##connector_labels", 2, labelFlags)) {
    ImGui::TableNextColumn();
    ImGui::TextDisabled("INPUTS");
    ImGui::TableNextColumn();
    ImGui::TextDisabled("TARGET");
    ImGui::EndTable();
  }
  drawArrow(ImVec2(width, ImGui::GetTextLineHeightWithSpacing() * 1.7f),
            style::col::Accent, 1.0f, st.planner.searching);

  // Widths must be measured INSIDE the toolbar: it insets its content by half a
  // gap on each side and tightens item spacing, so sizing against the card's
  // width pushes the trailing action past the card edge.
  widgets::beginToolbar("##route_actions");
  const float actionGap = ImGui::GetStyle().ItemSpacing.x;
  const float toolbarWidth = std::max(
      ImGui::GetContentRegionAvail().x - style::metrics().gap * 0.5f, 1.0f);
  const float clearWidth = widgets::actionButtonWidth(icons::Icon::Undo, "Clear");
  const float runWidth = std::max(ImGui::GetFrameHeight() * 3.0f,
                                  toolbarWidth - clearWidth - actionGap);
  ImGui::BeginDisabled(st.planner.searching);
  if (widgets::actionButton("##suggest_routes", icons::Icon::Play, "Suggest routes",
                            ImVec2(runWidth, 0.0f), true,
                            "Search the reaction knowledge base for routes"))
    dispatchSearch(st, animation);
  ImGui::EndDisabled();
  ImGui::SameLine(0.0f, actionGap);
  const bool nothingToClear =
      st.planner.routes.empty() && st.planner.error.empty() && !st.planner.searched;
  ImGui::BeginDisabled(st.planner.searching || nothingToClear);
  const bool clear = widgets::actionButton(
      "##clear_routes", icons::Icon::Undo, "Clear", ImVec2(clearWidth, 0.0f), false,
      "Clear the current route results");
  ImGui::EndDisabled();
  widgets::endToolbar();
  if (clear) {
    st.planner.routes.clear();
    st.planner.error.clear();
    st.planner.searched = false;
  }

  const double elapsed = animation.searchStartedAt < 0.0
                             ? kSearchMinimumSeconds
                             : ImGui::GetTime() - animation.searchStartedAt;
  if (st.planner.searching || elapsed < kSearchMinimumSeconds) {
    drawSearchingIndicator();
  } else {
    widgets::keyValue("Route source", "Knowledge base + optional AI");
  }

  const bool aiAvailable = rxn::llmAvailable();
  ImGui::BeginDisabled(!aiAvailable);
  widgets::toggle("##ai_fallback", "Use AI fallback", st.planner.allowLlm,
                  "Let the AI propose a route when the knowledge base stalls");
  ImGui::EndDisabled();
  if (!aiAvailable && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
    ImGui::SetTooltip("CHEMCAD_LLM_API_KEY is not set.");

  static constexpr const char* depthLabels[] = {"1", "2", "3", "4"};
  int depthIndex = std::clamp(st.planner.maxDepth - 1, 0, 3);
  widgets::keyValue("Maximum depth", depthLabels[depthIndex]);
  if (widgets::segmented("##route_depth", depthLabels, 4, depthIndex, width))
    st.planner.maxDepth = depthIndex + 1;
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "Maximum number of reaction steps per route");

  // The trailing ImGui label would be clipped by the narrow flow column, so the
  // value is reported above the control the same way the depth row is.
  char routesValue[16];
  std::snprintf(routesValue, sizeof(routesValue), "%d", st.planner.maxRoutes);
  widgets::keyValue("Routes", routesValue);
  ImGui::SetNextItemWidth(std::max(1.0f, ImGui::GetContentRegionAvail().x));
  ImGui::SliderInt("##route_count", &st.planner.maxRoutes, 1, 10);
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "Maximum number of routes to report");
  ImGui::EndGroup();
}

void stepArrow(float width) {
  width = std::max(width, 1.0f);
  drawArrow(ImVec2(width, ImGui::GetTextLineHeightWithSpacing()),
            style::col::TextDim);
}

void drawStep(AppState& st, const rxn::Step& step, const RoutePreviewCache& cache,
              int stepIndex) {
  ImGui::PushID(stepIndex);
  std::string heading = "Step " + std::to_string(stepIndex + 1);
  if (!step.reactionName.empty()) heading += "  ·  " + step.reactionName;
  const style::Metrics& metrics = style::metrics();
  const float tickWidth = std::max(metrics.hairline * 2.0f,
                                   ImGui::GetFontSize() * 0.16f);
  const float headingWidth = std::max(
      ImGui::GetContentRegionAvail().x - tickWidth - metrics.gap * 1.75f, 1.0f);
  const bool headingFont = style::pushFont(style::fonts::semibold());
  const std::string fittedHeading = ellipsize(heading, headingWidth);
  style::popFont(headingFont);
  widgets::sectionHeader(fittedHeading.c_str(), style::col::Teal);
  if (ImGui::IsItemHovered() && fittedHeading != heading)
    ImGui::SetTooltip("%s", heading.c_str());

  const std::string reactants = joined(step.reactantSmiles);
  const std::string reagents = joined(step.reagents);
  const std::string sideProducts = joined(step.sideProductSmiles);
  widgets::keyValue("Reaction",
                    step.reactionName.empty() ? "Unspecified" : step.reactionName.c_str());
  widgets::keyValue("Provenance",
                    step.source == rxn::Step::Source::LLM ? "AI" : "Knowledge base");
  widgets::keyValue("Reactant SMILES",
                    reactants.empty() ? "Unspecified" : reactants.c_str());
  widgets::keyValue("Reagents", reagents.empty() ? "None specified" : reagents.c_str());
  widgets::keyValue("Conditions",
                    step.conditions.empty() ? "Not specified" : step.conditions.c_str());
  widgets::keyValue("Product SMILES",
                    step.productSmiles.empty() ? "Unspecified" : step.productSmiles.c_str());
  widgets::keyValue("Side-product SMILES",
                    sideProducts.empty() ? "None predicted" : sideProducts.c_str());
  if (!step.notes.empty())
    widgets::keyValue("Mechanism / caveats", step.notes.c_str(), style::col::TextDim);

  const float rowWidth = std::max(ImGui::GetContentRegionAvail().x, 1.0f);
  const ImVec2 thumbSize(std::min(routeThumbSize().x, rowWidth), routeThumbSize().y);
  const float actionWidth =
      ImGui::CalcTextSize("Send to canvas").x + ImGui::GetFrameHeight() +
      metrics.gap;
  const float productWidth = std::max(thumbSize.x, actionWidth);
  const float rowGap = metrics.gap;
  const float minimumArrowWidth = ImGui::GetFontSize() * 6.0f;
  const bool horizontal =
      rowWidth >= thumbSize.x + productWidth + minimumArrowWidth + rowGap * 2.0f;
  const float arrowWidth =
      horizontal
          ? std::max(minimumArrowWidth,
                     rowWidth - thumbSize.x - productWidth - rowGap * 2.0f)
          : rowWidth;

  ImGui::BeginGroup();
  if (step.reactantSmiles.empty()) {
    moleculeThumbButton("##reactant_missing", cachedMolecule(cache, ""), thumbSize);
  } else {
    for (size_t i = 0; i < step.reactantSmiles.size(); ++i) {
      ImGui::PushID(static_cast<int>(i));
      moleculeThumbButton("##reactant", cachedMolecule(cache, step.reactantSmiles[i]),
                          thumbSize);
      if (ImGui::IsItemHovered() && !step.reactantSmiles[i].empty())
        ImGui::SetTooltip("%s", step.reactantSmiles[i].c_str());
      ImGui::PopID();
      if (i + 1 < step.reactantSmiles.size()) ImGui::TextUnformatted("+");
    }
  }
  ImGui::EndGroup();
  if (horizontal) ImGui::SameLine(0.0f, rowGap);
  stepArrow(arrowWidth);
  if (horizontal) ImGui::SameLine(0.0f, rowGap);

  ImGui::BeginGroup();
  moleculeThumbButton("##product", cachedMolecule(cache, step.productSmiles), thumbSize);
  if (ImGui::IsItemHovered() && !step.productSmiles.empty())
    ImGui::SetTooltip("%s", step.productSmiles.c_str());
  widgets::beginToolbar("##product_actions");
  ImGui::BeginDisabled(step.productSmiles.empty());
  if (widgets::actionButton("##copy_product", icons::Icon::Copy, "Copy SMILES",
                            ImVec2(0.0f, 0.0f), false,
                            "Copy the product SMILES to the clipboard"))
    ImGui::SetClipboardText(step.productSmiles.c_str());
  ImGui::SameLine(0.0f, metrics.gap);
  if (widgets::actionButton("##send_product", icons::Icon::ArrowRight,
                            "Send to canvas", ImVec2(0.0f, 0.0f), false,
                            "Append this product to the Sketch canvas"))
    openInSketch(st, step.productSmiles);
  ImGui::EndDisabled();
  widgets::endToolbar();
  ImGui::EndGroup();

  if (!step.sideProductSmiles.empty()) {
    ImGui::Spacing();
    for (size_t i = 0; i < step.sideProductSmiles.size(); ++i) {
      const std::string& smiles = step.sideProductSmiles[i];
      ImGui::PushID(static_cast<int>(i));
      ImGui::BeginGroup();
      const ImVec2 sideSize(ImGui::GetFontSize() * 6.3f,
                            ImGui::GetTextLineHeightWithSpacing() * 2.75f);
      moleculeThumbButton("##side_product", cachedMolecule(cache, smiles), sideSize);
      if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", smiles.c_str());
      ImGui::EndGroup();
      ImGui::PopID();
      const float contentRight =
          ImGui::GetCursorScreenPos().x + ImGui::GetContentRegionAvail().x;
      if (i + 1 < step.sideProductSmiles.size() &&
          ImGui::GetItemRectMax().x + metrics.gap + sideSize.x <= contentRight)
        ImGui::SameLine(0.0f, metrics.gap);
    }
  }
  ImGui::PopID();
}


void drawResults(AppState& st, RoutePreviewCache& cache,
                 PlannerAnimationState& animation) {
  const std::string signature = routesSignature(st.planner.routes);
  if (signature != cache.signature)
    rebuildRouteCache(cache, st.planner.routes, signature);

  if (animation.awaitingResults && !st.planner.searching && st.planner.searched) {
    animation.routeFirstSeen.clear();
    animation.awaitingResults = false;
  }

  if (!st.planner.error.empty()) {
    widgets::notice(icons::Icon::Warning, st.planner.error.c_str(),
                    style::col::Danger);
    ImGui::Spacing();
  }

  const bool hasStartingMaterial =
      std::any_of(st.planner.starts.begin(), st.planner.starts.end(),
                  [](const MaterialBox& material) { return !material.smiles.empty(); });
  if (!hasStartingMaterial) {
    widgets::emptyState(
        icons::Icon::Flask, "Add a starting material",
        "Enter a SMILES string or import the largest structure from the Sketch canvas.");
    return;
  }

  if (st.planner.target.smiles.empty()) {
    widgets::emptyState(
        icons::Icon::Molecule, "Define the target product",
        "Enter the desired product in the Target card, then choose Suggest routes.");
    return;
  }

  if (st.planner.searching && st.planner.routes.empty()) {
    widgets::emptyState(
        icons::Icon::Search, "Searching reaction space",
        "Candidate transformations are being assembled and ranked; results appear here.");
    return;
  }

  if (st.planner.searched && st.planner.routes.empty()) {
    const bool canUseAi = st.planner.allowLlm && rxn::llmAvailable();
    widgets::emptyState(
        icons::Icon::Retro, "No routes found",
        canUseAi
            ? "Increase the search depth or revise the starting material and target structures."
            : "Enable AI fallback when available, increase depth, or revise the structures.");
    return;
  }

  if (st.planner.routes.empty()) {
    widgets::emptyState(
        icons::Icon::Play, "Ready to plan",
        "Review the structures and route limits, then choose Suggest routes.");
    return;
  }

  cache.selectedRoute =
      std::clamp(cache.selectedRoute, 0,
                 static_cast<int>(st.planner.routes.size()) - 1);
  std::vector<std::string> barLabels;
  std::vector<std::string> barAnnotations;
  barLabels.reserve(st.planner.routes.size());
  barAnnotations.reserve(st.planner.routes.size());
  for (size_t i = 0; i < st.planner.routes.size(); ++i) {
    barLabels.push_back("Route " + std::to_string(i + 1));
    barAnnotations.push_back(formatScore(st.planner.routes[i].score));
  }
  std::vector<charts::BarRow> bars;
  bars.reserve(st.planner.routes.size());
  for (size_t i = 0; i < st.planner.routes.size(); ++i) {
    charts::BarRow row;
    row.label = barLabels[i].c_str();
    row.value = st.planner.routes[i].score;
    row.annotation = barAnnotations[i].c_str();
    row.accent = i == 0 ? style::col::Accent : style::col::Teal;
    row.selected = static_cast<int>(i) == cache.selectedRoute;
    bars.push_back(row);
  }
  const float rankingHeight =
      ImGui::GetTextLineHeightWithSpacing() *
      (static_cast<float>(bars.size()) + 1.0f);
  const int selected =
      charts::rankedBars("##route_ranking", bars.data(), static_cast<int>(bars.size()),
                         ImVec2(ImGui::GetContentRegionAvail().x, rankingHeight));
  if (selected >= 0) cache.selectedRoute = selected;

  double minimumScore = st.planner.routes.front().score;
  double maximumScore = minimumScore;
  for (const rxn::Route& route : st.planner.routes) {
    minimumScore = std::min(minimumScore, route.score);
    maximumScore = std::max(maximumScore, route.score);
  }
  const rxn::Route& selectedRoute =
      st.planner.routes[static_cast<size_t>(cache.selectedRoute)];
  const double relativeScore =
      maximumScore > minimumScore
          ? std::clamp((selectedRoute.score - minimumScore) /
                           (maximumScore - minimumScore),
                       0.0, 1.0)
          : 1.0;
  char confidenceBuffer[16];
  std::snprintf(confidenceBuffer, sizeof(confidenceBuffer), "%.0f%%",
                relativeScore * 100.0);
  charts::GaugeStyle gaugeStyle;
  gaugeStyle.accent = style::col::Accent;
  gaugeStyle.minLabel = "LOWER";
  gaugeStyle.maxLabel = "BEST";
  charts::gauge(
      "##route_confidence", relativeScore, confidenceBuffer,
      "Relative confidence",
      ImVec2(std::min(ImGui::GetContentRegionAvail().x,
                      ImGui::GetFontSize() * 14.0f),
             ImGui::GetTextLineHeightWithSpacing() * 6.0f),
      gaugeStyle);

  const std::string selectedScore = formatScore(selectedRoute.score);
  const std::string candidateCount = std::to_string(st.planner.routes.size());
  const std::string selectedRank =
      "#" + std::to_string(cache.selectedRoute + 1);
  constexpr ImGuiTableFlags metricFlags =
      ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_NoSavedSettings |
      ImGuiTableFlags_NoPadOuterX;
  if (ImGui::BeginTable("##result_metrics", 3, metricFlags)) {
    ImGui::TableNextColumn();
    widgets::metric("CANDIDATES", candidateCount.c_str());
    ImGui::TableNextColumn();
    widgets::metric("SELECTED", selectedRank.c_str());
    ImGui::TableNextColumn();
    widgets::metric("ROUTE SCORE", selectedScore.c_str(), nullptr, nullptr,
                    style::col::Accent);
    ImGui::EndTable();
  }

  const double now = ImGui::GetTime();
  for (size_t routeIndex = 0; routeIndex < st.planner.routes.size(); ++routeIndex) {
    const rxn::Route& route = st.planner.routes[routeIndex];
    ImGui::PushID(static_cast<int>(routeIndex));

    const std::string identity = routeIdentity(route);
    const auto [seen, inserted] = animation.routeFirstSeen.emplace(identity, now);
    (void)inserted;
    const float enter = easeOutCubic(static_cast<float>(
        (now - seen->second) / static_cast<double>(kRouteEnterSeconds)));
    ImGui::PushStyleVar(ImGuiStyleVar_Alpha,
                        ImGui::GetStyle().Alpha * (0.12f + 0.88f * enter));
    ImGui::Dummy(
        ImVec2(0.0f, (1.0f - enter) * ImGui::GetTextLineHeightWithSpacing() * 0.4f));

    const ImGuiID hoverId = ImGui::GetID("##route_hover");
    const float previousHover = ImGui::GetStateStorage()->GetFloat(hoverId, 0.0f);
    const float liftReserve = style::metrics().gap * 0.35f;
    ImGui::Dummy(ImVec2(0.0f, liftReserve * (1.0f - previousHover)));

    const bool cardOpen =
        widgets::beginCard("##route_result", ImVec2(0.0f, 0.0f), style::col::BgSurface);
    if (cardOpen) {
      const std::string title = routeName(route, routeIndex);
      const bool titleFont = style::pushFont(style::fonts::semibold());
      ImGui::TextUnformatted(title.c_str());
      style::popFont(titleFont);

      const std::string routeScore = formatScore(route.score);
      const std::string stepCount = std::to_string(route.steps.size());
      constexpr ImGuiTableFlags headlineFlags =
          ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_NoSavedSettings |
          ImGuiTableFlags_NoPadOuterX;
      if (ImGui::BeginTable("##route_headline", 3, headlineFlags)) {
        ImGui::TableNextColumn();
        widgets::metric("SCORE", routeScore.c_str(), nullptr, nullptr,
                        routeIndex == 0 ? style::col::Accent : style::col::Text);
        ImGui::TableNextColumn();
        widgets::metric("STEPS", stepCount.c_str());
        ImGui::TableNextColumn();
        provenanceBadge(route.usesLlm());
        if (routeIndex == 0) {
          ImGui::SameLine(0.0f, style::metrics().gap);
          widgets::badge("BEST", style::col::Accent);
        }
        ImGui::EndTable();
      }

      std::string summary = stepCount + (route.steps.size() == 1 ? " step" : " steps");
      summary += " · score " + routeScore;
      const bool detailsOpen =
          widgets::disclosure("##route_details", "Conditions and mechanism",
                              summary.c_str(), false, icons::Icon::Retro,
                              routeIndex == 0 ? style::col::Accent : style::col::Teal);
      if (detailsOpen) {
        ImGui::Indent(style::metrics().gap);
        const std::string bestScore = formatScore(maximumScore);
        const double routeRelative =
            maximumScore > minimumScore
                ? std::clamp((route.score - minimumScore) /
                                 (maximumScore - minimumScore),
                             0.0, 1.0)
                : 1.0;
        char routeConfidence[16];
        std::snprintf(routeConfidence, sizeof(routeConfidence), "%.0f%%",
                      routeRelative * 100.0);
        widgets::keyValue("Score target", bestScore.c_str());
        widgets::keyValue("Relative confidence", routeConfidence);
        charts::bullet(
            "##route_score_bullet", route.score, maximumScore, minimumScore,
            maximumScore,
            ImVec2(ImGui::GetContentRegionAvail().x, ImGui::GetFrameHeight() * 1.5f),
            routeIndex == 0 ? style::col::Accent : style::col::Teal);
        for (size_t stepIndex = 0; stepIndex < route.steps.size(); ++stepIndex) {
          drawStep(st, route.steps[stepIndex], cache, static_cast<int>(stepIndex));
          if (stepIndex + 1 < route.steps.size()) ImGui::Spacing();
        }
        ImGui::Unindent(style::metrics().gap);
      }
      widgets::endCard();
    }

    const ImVec2 cardMin = ImGui::GetItemRectMin();
    const ImVec2 cardMax = ImGui::GetItemRectMax();
    const bool hovered = ImGui::IsItemHovered();
    const float hover = widgets::hoverT(hoverId, hovered);
    ImGui::GetWindowDrawList()->AddRect(
        cardMin, cardMax, style::mix(style::col::Border, style::col::Teal, hover),
        style::metrics().radiusMd, 0,
        style::metrics().hairline * (1.0f + hover));

    ImGui::PopStyleVar();
    ImGui::PopID();
    ImGui::Spacing();
  }
}

void drawStartingMaterialsCard(AppState& st, float cardHeight) {
  int removeIndex = -1;
  if (widgets::beginCard("##starting_materials_card", ImVec2(0.0f, cardHeight),
                         style::col::BgSurface)) {
    widgets::cardHeader(icons::Icon::Flask, "Starting materials",
                        "Combine one or more molecular inputs", style::col::Teal);
    if (widgets::actionButton(
            "##add_starting_material", icons::Icon::Plus, "Add starting material",
            ImVec2(0.0f, 0.0f), false,
            "Add another molecular input to the reaction"))
      st.planner.starts.emplace_back();

    if (st.planner.starts.empty()) {
      widgets::emptyState(
          icons::Icon::Flask, "No starting materials",
          "Choose Add starting material, then enter SMILES or import from Sketch.");
    } else if (ImGui::BeginChild("##starting_materials_row", ImVec2(0.0f, 0.0f),
                                 ImGuiChildFlags_None,
                                 ImGuiWindowFlags_HorizontalScrollbar)) {
      for (size_t i = 0; i < st.planner.starts.size(); ++i) {
        ImGui::PushID(static_cast<int>(i));
        if (materialBoxWidget(st, st.planner.starts[i], false, static_cast<int>(i)))
          removeIndex = static_cast<int>(i);
        ImGui::PopID();
        if (i + 1 < st.planner.starts.size()) {
          ImGui::SameLine(0.0f, style::metrics().gap);
          ImGui::AlignTextToFramePadding();
          ImGui::TextColored(style::col::TextFaint, "+");
          ImGui::SameLine(0.0f, style::metrics().gap);
        }
      }
    }
    if (!st.planner.starts.empty()) ImGui::EndChild();
    widgets::endCard();
  }

  if (removeIndex >= 0 && st.planner.starts.size() > 1)
    st.planner.starts.erase(st.planner.starts.begin() + removeIndex);
}

void drawRouteCard(AppState& st, PlannerAnimationState& animation, float cardHeight) {
  if (!widgets::beginCard("##route_card", ImVec2(0.0f, cardHeight),
                          style::col::BgSurface))
    return;
  widgets::cardHeader(icons::Icon::Retro, "Route", "Search and ranking controls",
                      style::col::Accent);
  connectorControls(st, animation);
  widgets::endCard();
}

void drawTargetCard(AppState& st, float cardHeight) {
  if (!widgets::beginCard("##target_card", ImVec2(0.0f, cardHeight),
                          style::col::BgSurface))
    return;
  widgets::cardHeader(icons::Icon::Molecule, "Target",
                      "The product this route must reach", style::col::Violet);
  const float availableWidth = std::max(ImGui::GetContentRegionAvail().x, 1.0f);
  const float targetWidth = std::min(materialWidth(), availableWidth);
  const float offset = std::max((availableWidth - targetWidth) * 0.5f, 0.0f);
  ImGui::Indent(offset);
  ImGui::PushID("target");
  materialBoxWidget(st, st.planner.target, true, -1);
  ImGui::PopID();
  ImGui::Unindent(offset);
  widgets::endCard();
}

void drawResultsCard(AppState& st, RoutePreviewCache& cache,
                     PlannerAnimationState& animation) {
  if (!widgets::beginCard("##results_card", ImVec2(0.0f, 0.0f),
                          style::col::BgSurface))
    return;
  widgets::cardHeader(icons::Icon::ChartBars, "Results",
                      "Ranked retrosynthetic candidates", style::col::Violet);
  drawResults(st, cache, animation);
  widgets::endCard();
}

}  // namespace

void drawReactionPlanner(AppState& st) {
  if (st.planner.starts.empty()) st.planner.starts.emplace_back();

  static RoutePreviewCache routeCache;
  static PlannerAnimationState animation;
  const float flowCardHeight =
      materialHeight() + ImGui::GetFrameHeight() * 4.0f;
  const float availableWidth = ImGui::GetContentRegionAvail().x;
  const bool horizontalFlow = availableWidth >= 880.0f * uiScale();

  if (horizontalFlow) {
    constexpr ImGuiTableFlags flowFlags =
        ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_NoSavedSettings;
    if (ImGui::BeginTable("##reaction_flow", 3, flowFlags)) {
      // The route column carries the widest fixed content -- the run and clear
      // actions side by side -- so it gets a share close to the target card
      // rather than the narrowest one.
      ImGui::TableSetupColumn("##starts", ImGuiTableColumnFlags_WidthStretch, 1.25f);
      ImGui::TableSetupColumn("##route", ImGuiTableColumnFlags_WidthStretch, 1.05f);
      ImGui::TableSetupColumn("##target", ImGuiTableColumnFlags_WidthStretch, 1.00f);
      ImGui::TableNextColumn();
      drawStartingMaterialsCard(st, flowCardHeight);
      ImGui::TableNextColumn();
      drawRouteCard(st, animation, flowCardHeight);
      ImGui::TableNextColumn();
      drawTargetCard(st, flowCardHeight);
      ImGui::EndTable();
    }
  } else {
    drawStartingMaterialsCard(st, flowCardHeight);
    ImGui::Spacing();
    drawRouteCard(st, animation, 0.0f);
    ImGui::Spacing();
    drawTargetCard(st, 0.0f);
  }

  ImGui::Spacing();
  drawResultsCard(st, routeCache, animation);
}

}  // namespace chemcad::ui

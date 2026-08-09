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
#include "ui/layout.hpp"
#include "ui/icons.hpp"
#include "ui/mol_thumb.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"
#include "ui/widgets.hpp"

namespace chemcad::ui {
namespace {

inline float materialWidth() { return ImGui::GetFontSize() * 17.7f; }
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

constexpr float kSearchMinimumSeconds = 0.45f;

struct RoutePreviewCache {
  std::string signature;
  std::unordered_map<std::string, core::Molecule> molecules;
  int selectedRoute = 0;
};

struct PlannerAnimationState {
  double searchStartedAt = -1.0;
  bool awaitingResults = false;
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

// rxn::Step exposes no per-step score or yield contribution, so a route
// waterfall would fabricate data and is intentionally omitted.


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
  const layout::Frame outerFrame = layout::measure();
  const float liftReserve = outerFrame.gap * 0.35f;
  const ImGuiID hoverId = ImGui::GetID("##material_hover");
  const float previousHover = ImGui::GetStateStorage()->GetFloat(hoverId, 0.0f);
  ImGui::Dummy(ImVec2(0.0f, liftReserve * (1.0f - previousHover)));

  const layout::Frame frame = layout::measure();
  const ImVec2 cardSize(std::min(materialWidth(), frame.size.x),
                        std::max(layout::pageHeight(), frame.row));
  constexpr ImGuiWindowFlags cardFlags =
      ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse;
  if (!widgets::beginCard("##material_box", cardSize, style::col::BgSurface,
                          cardFlags))
    return false;

  const char* heading = target ? "Product" : "Starting material";
  const ImVec4 headingColor = target ? style::col::Violet : style::col::Teal;
  const float headerRight =
      ImGui::GetCursorPosX() + ImGui::GetContentRegionAvail().x;
  const bool pushed = style::pushFont(style::fonts::semibold());
  ImGui::TextColored(headingColor, "%s", heading);
  style::popFont(pushed);
  if (!target) {
    const float removeWidth =
        widgets::actionButtonWidth(icons::Icon::Trash, "Remove");
    ImGui::SameLine();
    ImGui::SetCursorPosX(
        std::max(ImGui::GetCursorPosX(), headerRight - removeWidth));
    const bool onlyStart = st.planner.starts.size() <= 1;
    ImGui::BeginDisabled(onlyStart);
    ImGui::PushStyleVar(ImGuiStyleVar_DisabledAlpha, 0.35f);
    remove = widgets::actionButton(
        "##remove", icons::Icon::Trash, "Remove", ImVec2(removeWidth, 0.0f),
        false, "Remove this starting material");
    ImGui::PopStyleVar();
    ImGui::EndDisabled();
  }

  const layout::Frame contentFrame = layout::measure();
  if (moleculeThumbButton(
          "##preview", box.preview,
          ImVec2(contentFrame.size.x, contentFrame.row * 3.3f)) &&
      box.previewValid) {
    openInSketch(st, box.smiles);
  }
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", box.previewValid ? "Open in the Sketch canvas"
                                             : "Structure preview");

  ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x);
  const bool smilesEnter = widgets::stringInputWithHint(
      "##smiles", "Structure (SMILES)", box.smiles,
      ImGuiInputTextFlags_EnterReturnsTrue, true);
  const bool smilesCommitted = smilesEnter || ImGui::IsItemDeactivatedAfterEdit();
  if (smilesCommitted) rebuildPreview(box);
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "SMILES string; Enter updates the preview");

  // The reservation and the button share one width, so the pair can never
  // disagree and push the action past the card edge in a narrow flow column.
  const float lookupWidth = ImGui::CalcTextSize("Look up").x + frame.em +
                            frame.gap * 2.5f;
  ImGui::SetNextItemWidth(std::max(
      frame.em, ImGui::GetContentRegionAvail().x - lookupWidth -
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
  if (statusText && layout::fits(statusText, ImGui::GetContentRegionAvail().x)) {
    ImGui::SameLine();
    ImGui::TextColored(style::col::DataDim, "%s", statusText);
  }
  if (!box.error.empty())
    widgets::notice(icons::Icon::Warning, box.error.c_str(), style::col::Danger);

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
        style::u32(style::col::DataBright));
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
  const layout::Frame frame = layout::measure();
  const char* const inputLabels[] = {"Starting materials", "Inputs"};
  const char* const targetLabels[] = {"Target product", "Target"};
  const char* inputLabel =
      layout::bestLabel(frame.size.x * 0.5f, inputLabels, 2);
  const char* targetLabel =
      layout::bestLabel(frame.size.x * 0.5f, targetLabels, 2);
  const float targetWidth = ImGui::CalcTextSize(targetLabel).x;
  ImGui::TextDisabled("%s", inputLabel);
  ImGui::SameLine();
  ImGui::SetCursorPosX(std::max(
      ImGui::GetCursorPosX(),
      ImGui::GetCursorPosX() + ImGui::GetContentRegionAvail().x - targetWidth));
  ImGui::TextDisabled("%s", targetLabel);
  drawArrow(ImVec2(frame.size.x, frame.row * 1.7f), style::col::Data,
            1.0f, st.planner.searching);

  // Widths are measured inside the toolbar because its inset is part of the
  // actual budget available to the two actions.
  widgets::beginToolbar("##route_actions");
  const layout::Frame toolbarFrame = layout::measure();
  const float actionGap = ImGui::GetStyle().ItemSpacing.x;
  const float clearWidth = widgets::actionButtonWidth(icons::Icon::Undo, "Clear");
  const float runWidth = std::max(
      toolbarFrame.control * 3.0f,
      toolbarFrame.size.x - clearWidth - actionGap - toolbarFrame.gap * 0.5f);
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
  const bool searchActive =
      st.planner.searching || elapsed < kSearchMinimumSeconds;
  widgets::statusDot(searchActive ? "Searching routes" : "Planner ready",
                     searchActive, style::col::Data);

  if (widgets::onlyWhen(rxn::llmAvailable(),
                        "AI fallback requires CHEMCAD_LLM_API_KEY.")) {
    widgets::toggle("##ai_fallback", "Use AI fallback", st.planner.allowLlm,
                    "Let the AI propose a route when the knowledge base stalls");
  }

  static constexpr const char* depthLabels[] = {"1", "2", "3", "4"};
  int depthIndex = std::clamp(st.planner.maxDepth - 1, 0, 3);
  ImGui::TextDisabled("Maximum depth");
  if (widgets::segmented("##route_depth", depthLabels, 4, depthIndex,
                         frame.size.x))
    st.planner.maxDepth = depthIndex + 1;
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "Maximum number of reaction steps per route");

  float routeCount = static_cast<float>(st.planner.maxRoutes);
  if (widgets::glyphSlider(
          "##route_count", icons::Icon::List, "Routes", routeCount, 1.0f,
          10.0f, "%.0f", "Maximum number of routes to report")) {
    st.planner.maxRoutes =
        std::clamp(static_cast<int>(std::lround(routeCount)), 1, 10);
  }
  ImGui::EndGroup();
}

void stepArrow(float width) {
  const layout::Frame frame = layout::measure();
  drawArrow(ImVec2(std::max(width, frame.em), frame.row),
            style::col::DataDim);
}

void drawStep(AppState& st, const rxn::Step& step, const RoutePreviewCache& cache,
              int stepIndex) {
  ImGui::PushID(stepIndex);
  std::string heading = "Step " + std::to_string(stepIndex + 1);
  if (!step.reactionName.empty()) heading += "  ·  " + step.reactionName;
  const style::Metrics& metrics = style::metrics();
  const layout::Frame detailFrame = layout::measure();
  const float tickWidth = std::max(metrics.hairline * 2.0f,
                                   detailFrame.em * 0.16f);
  const float headingWidth = std::max(
      detailFrame.size.x - tickWidth - detailFrame.gap * 1.75f,
      detailFrame.em);
  const bool headingFont = style::pushFont(style::fonts::semibold());
  const std::string fittedHeading = ellipsize(heading, headingWidth);
  style::popFont(headingFont);
  widgets::sectionHeader(fittedHeading.c_str(), style::col::Teal);
  if (ImGui::IsItemHovered() && fittedHeading != heading)
    ImGui::SetTooltip("%s", heading.c_str());

  const std::string reactants = joined(step.reactantSmiles);
  const std::string reagents = joined(step.reagents);
  const std::string sideProducts = joined(step.sideProductSmiles);
  widgets::keyValue(
      "Reaction",
      step.reactionName.empty() ? "Unspecified" : step.reactionName.c_str(),
      style::col::Data);
  widgets::keyValue(
      "Provenance",
      step.source == rxn::Step::Source::LLM ? "AI" : "Knowledge base",
      step.source == rxn::Step::Source::LLM ? style::col::Violet
                                            : style::col::Teal);
  widgets::keyValue("Reactant SMILES",
                    reactants.empty() ? "Unspecified" : reactants.c_str(),
                    style::col::Data);
  widgets::keyValue("Reagents",
                    reagents.empty() ? "None specified" : reagents.c_str(),
                    style::col::Data);
  widgets::keyValue(
      "Conditions",
      step.conditions.empty() ? "Not specified" : step.conditions.c_str(),
      style::col::Data);
  widgets::keyValue(
      "Product SMILES",
      step.productSmiles.empty() ? "Unspecified" : step.productSmiles.c_str(),
      style::col::Data);
  widgets::keyValue(
      "Side-product SMILES",
      sideProducts.empty() ? "None predicted" : sideProducts.c_str(),
      style::col::Data);
  if (!step.notes.empty())
    widgets::keyValue("Mechanism / caveats", step.notes.c_str(),
                      style::col::DataDim);

  const layout::Frame rowFrame = layout::measure();
  const float rowWidth = rowFrame.size.x;
  const ImVec2 thumbSize(std::min(routeThumbSize().x, rowWidth),
                         routeThumbSize().y);
  const float actionWidth =
      ImGui::CalcTextSize("Send to canvas").x + rowFrame.control + rowFrame.gap;
  const float productWidth = std::max(thumbSize.x, actionWidth);
  const float rowGap = rowFrame.gap;
  const float minimumArrowWidth = rowFrame.em * 6.0f;
  const bool horizontal = layout::columnsThatFit(rowFrame, 18.0f) >= 2;
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
    const layout::Frame sideFrame = layout::measure();
    const int sideColumns =
        layout::columnsThatFit(sideFrame, routeThumbSize().x / sideFrame.em);
    for (size_t i = 0; i < step.sideProductSmiles.size(); ++i) {
      const std::string& smiles = step.sideProductSmiles[i];
      ImGui::PushID(static_cast<int>(i));
      ImGui::BeginGroup();
      const ImVec2 sideSize(sideFrame.em * 6.3f, sideFrame.row * 2.75f);
      moleculeThumbButton("##side_product", cachedMolecule(cache, smiles), sideSize);
      if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", smiles.c_str());
      ImGui::EndGroup();
      ImGui::PopID();
      if (i + 1 < step.sideProductSmiles.size() &&
          static_cast<int>((i + 1) % static_cast<size_t>(sideColumns)) != 0)
        ImGui::SameLine(0.0f, sideFrame.gap);
    }
  }
  ImGui::PopID();
}


void drawResults(AppState& st, RoutePreviewCache& cache,
                 PlannerAnimationState& animation) {
  const std::string signature = routesSignature(st.planner.routes);
  if (signature != cache.signature)
    rebuildRouteCache(cache, st.planner.routes, signature);

  if (animation.awaitingResults && !st.planner.searching && st.planner.searched)
    animation.awaitingResults = false;

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

  double minimumScore = st.planner.routes.front().score;
  double maximumScore = minimumScore;
  for (const rxn::Route& route : st.planner.routes) {
    minimumScore = std::min(minimumScore, route.score);
    maximumScore = std::max(maximumScore, route.score);
  }

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
    row.accent = static_cast<int>(i) == cache.selectedRoute
                     ? style::col::Accent
                     : style::col::Data;
    row.selected = static_cast<int>(i) == cache.selectedRoute;
    bars.push_back(row);
  }

  const layout::Frame resultFrame = layout::measure(
      ImVec2(ImGui::GetContentRegionAvail().x, layout::pageHeight()));
  const float resultStartY = ImGui::GetCursorPosY();
  const float sectionWeights[] = {1.0f, 1.35f};
  const float sectionMinimums[] = {resultFrame.row * 5.0f,
                                   resultFrame.row * 5.0f};
  float sectionHeights[2]{};
  layout::distribute(resultFrame.size.y, sectionWeights, sectionMinimums, 2,
                     resultFrame.gap, sectionHeights);

  charts::GaugeStyle gaugeStyle;
  gaugeStyle.accent = style::col::Data;
  gaugeStyle.minLabel = "Lower";
  gaugeStyle.maxLabel = "Best";

  const int chartColumns =
      std::min(2, layout::columnsThatFit(resultFrame, 20.0f));
  if (chartColumns == 2) {
    const float chartWidth = layout::columnWidth(resultFrame, 2);
    const int selected = charts::rankedBars(
        "##route_ranking", bars.data(), static_cast<int>(bars.size()),
        ImVec2(chartWidth, sectionHeights[0]));
    if (selected >= 0) cache.selectedRoute = selected;

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
    ImGui::SameLine(0.0f, resultFrame.gap);
    charts::gauge("##route_confidence", relativeScore, confidenceBuffer,
                  "Relative confidence",
                  ImVec2(chartWidth, sectionHeights[0]), gaugeStyle);
  } else {
    const float chartWeights[] = {1.0f, 1.0f};
    const float chartMinimums[] = {resultFrame.row * 2.0f,
                                   resultFrame.row * 2.0f};
    float chartHeights[2]{};
    layout::distribute(sectionHeights[0], chartWeights, chartMinimums, 2,
                       resultFrame.gap, chartHeights);
    const int selected = charts::rankedBars(
        "##route_ranking", bars.data(), static_cast<int>(bars.size()),
        ImVec2(resultFrame.size.x, chartHeights[0]));
    if (selected >= 0) cache.selectedRoute = selected;

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
    layout::nextRow(resultStartY + chartHeights[0] + resultFrame.gap);
    charts::gauge("##route_confidence", relativeScore, confidenceBuffer,
                  "Relative confidence",
                  ImVec2(resultFrame.size.x, chartHeights[1]), gaugeStyle);
  }

  layout::nextRow(resultStartY + sectionHeights[0] + resultFrame.gap);
  const bool listOpen = ImGui::BeginChild(
      "##route_results_list", ImVec2(0.0f, sectionHeights[1]),
      ImGuiChildFlags_None, ImGuiWindowFlags_AlwaysVerticalScrollbar);
  if (listOpen) {
    const widgets::Column columns[] = {
        {"Route name", false, true, nullptr, 12.0f},
        {"Steps", true, false, nullptr, 4.0f},
        {"Score", true, false, nullptr, 5.0f},
        {"Source", false, false, nullptr, 8.0f},
    };
    if (widgets::beginDataTable("##routes_table", columns, 4,
                                ImVec2(0.0f, 0.0f))) {
      for (size_t routeIndex = 0; routeIndex < st.planner.routes.size();
           ++routeIndex) {
        const rxn::Route& route = st.planner.routes[routeIndex];
        const std::string name = routeName(route, routeIndex);
        const std::string source =
            route.usesLlm() ? "AI-assisted" : "Knowledge base";
        widgets::dataRow(static_cast<int>(routeIndex) == cache.selectedRoute
                             ? style::col::Accent
                             : ImVec4(0.0f, 0.0f, 0.0f, 0.0f));
        widgets::dataCell(name.c_str());
        widgets::dataCellf("%zu", route.steps.size());
        widgets::dataCellf("%.1f", route.score);
        widgets::dataCell(source.c_str());
      }
      widgets::endDataTable();
    }

    for (size_t routeIndex = 0; routeIndex < st.planner.routes.size();
         ++routeIndex) {
      const rxn::Route& route = st.planner.routes[routeIndex];
      ImGui::PushID(static_cast<int>(routeIndex));
      const std::string title = routeName(route, routeIndex);
      if (widgets::disclosure(
              "##route_details", title.c_str(), "Steps, conditions and provenance",
              false, icons::Icon::Retro, style::col::Accent)) {
        ImGui::Indent(resultFrame.gap);
        for (size_t stepIndex = 0; stepIndex < route.steps.size(); ++stepIndex) {
          drawStep(st, route.steps[stepIndex], cache,
                   static_cast<int>(stepIndex));
          if (stepIndex + 1 < route.steps.size()) ImGui::Spacing();
        }
        ImGui::Unindent(resultFrame.gap);
      }
      ImGui::PopID();
    }
  }
  ImGui::EndChild();
}

void drawStartingMaterialsCard(AppState& st, ImVec2 cardSize) {
  static int selectedMaterial = 0;
  int removeIndex = -1;
  constexpr ImGuiWindowFlags cardFlags =
      ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse;
  if (widgets::beginCard("##starting_materials_card", cardSize,
                         style::col::BgSurface, cardFlags)) {
    widgets::cardHeader(icons::Icon::Flask, "Starting materials",
                        "Combine one or more molecular inputs", style::col::Teal);
    if (widgets::actionButton(
            "##add_starting_material", icons::Icon::Plus, "Add starting material",
            ImVec2(0.0f, 0.0f), false,
            "Add another molecular input to the reaction")) {
      st.planner.starts.emplace_back();
      selectedMaterial = static_cast<int>(st.planner.starts.size()) - 1;
    }

    if (st.planner.starts.empty()) {
      widgets::emptyState(
          icons::Icon::Flask, "No starting materials",
          "Choose Add starting material, then enter SMILES or import from Sketch.");
    } else {
      selectedMaterial = std::clamp(
          selectedMaterial, 0, static_cast<int>(st.planner.starts.size()) - 1);
      if (st.planner.starts.size() > 1) {
        std::vector<std::string> tabLabels;
        std::vector<const char*> tabPointers;
        tabLabels.reserve(st.planner.starts.size());
        tabPointers.reserve(st.planner.starts.size());
        for (size_t i = 0; i < st.planner.starts.size(); ++i)
          tabLabels.push_back("Material " + std::to_string(i + 1));
        for (const std::string& label : tabLabels)
          tabPointers.push_back(label.c_str());
        widgets::subTabs("##material_tabs", tabPointers.data(), nullptr,
                         static_cast<int>(tabPointers.size()), selectedMaterial);
      }
      ImGui::PushID(selectedMaterial);
      if (materialBoxWidget(
              st, st.planner.starts[static_cast<size_t>(selectedMaterial)],
              false, selectedMaterial))
        removeIndex = selectedMaterial;
      ImGui::PopID();
    }
    widgets::endCard();
  }

  if (removeIndex >= 0 && st.planner.starts.size() > 1) {
    st.planner.starts.erase(st.planner.starts.begin() + removeIndex);
    selectedMaterial = std::clamp(
        selectedMaterial, 0, static_cast<int>(st.planner.starts.size()) - 1);
  }
}

void drawRouteCard(AppState& st, PlannerAnimationState& animation,
                   ImVec2 cardSize) {
  constexpr ImGuiWindowFlags cardFlags =
      ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse;
  if (!widgets::beginCard("##route_card", cardSize, style::col::BgSurface,
                          cardFlags))
    return;
  widgets::cardHeader(icons::Icon::Retro, "Route", "Search and ranking controls",
                      style::col::Accent);
  connectorControls(st, animation);
  widgets::endCard();
}

void drawTargetCard(AppState& st, ImVec2 cardSize) {
  constexpr ImGuiWindowFlags cardFlags =
      ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse;
  if (!widgets::beginCard("##target_card", cardSize, style::col::BgSurface,
                          cardFlags))
    return;
  widgets::cardHeader(icons::Icon::Molecule, "Target product",
                      "The product this route must reach", style::col::Violet);
  const layout::Frame frame = layout::measure();
  const float targetWidth = std::min(materialWidth(), frame.size.x);
  const float offset = std::max((frame.size.x - targetWidth) * 0.5f, 0.0f);
  ImGui::Indent(offset);
  ImGui::PushID("target");
  materialBoxWidget(st, st.planner.target, true, -1);
  ImGui::PopID();
  ImGui::Unindent(offset);
  widgets::endCard();
}

void drawResultsCard(AppState& st, RoutePreviewCache& cache,
                     PlannerAnimationState& animation, float cardHeight) {
  constexpr ImGuiWindowFlags cardFlags =
      ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse;
  if (!widgets::beginCard("##results_card", ImVec2(0.0f, cardHeight),
                          style::col::BgSurface, cardFlags))
    return;
  widgets::cardHeader(icons::Icon::ChartBars, "Results",
                      "Ranked retrosynthetic candidates", style::col::Data);
  drawResults(st, cache, animation);
  widgets::endCard();
}

}  // namespace

void drawReactionPlanner(AppState& st) {
  if (st.planner.starts.empty()) st.planner.starts.emplace_back();

  static RoutePreviewCache routeCache;
  static PlannerAnimationState animation;
  static int flowTab = 0;

  const layout::Frame page = layout::measure(
      ImVec2(ImGui::GetContentRegionAvail().x, layout::pageHeight()));
  const float pageStartY = ImGui::GetCursorPosY();
  const float pageWeights[] = {1.0f, 0.75f};
  const float pageMinimums[] = {page.row * 18.0f, page.row * 10.0f};
  float pageHeights[2]{};
  layout::distribute(page.size.y, pageWeights, pageMinimums, 2, page.gap,
                     pageHeights);

  const layout::Frame flowFrame =
      layout::measure(ImVec2(page.size.x, pageHeights[0]));
  const int flowColumns =
      std::min(3, layout::columnsThatFit(flowFrame, 18.0f));
  if (flowColumns == 3) {
    const float cardWidth = layout::columnWidth(flowFrame, 3);
    const ImVec2 cardSize(cardWidth, pageHeights[0]);
    drawStartingMaterialsCard(st, cardSize);
    ImGui::SameLine(0.0f, flowFrame.gap);
    drawRouteCard(st, animation, cardSize);
    ImGui::SameLine(0.0f, flowFrame.gap);
    drawTargetCard(st, cardSize);
  } else {
    constexpr ImGuiWindowFlags flowFlags =
        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse;
    const bool flowOpen = ImGui::BeginChild(
        "##reaction_flow", ImVec2(0.0f, pageHeights[0]),
        ImGuiChildFlags_None, flowFlags);
    if (flowOpen) {
      static constexpr const char* labels[] = {
          "Starting materials", "Route", "Target product"};
      static constexpr icons::Icon glyphs[] = {
          icons::Icon::Flask, icons::Icon::Retro, icons::Icon::Molecule};
      widgets::subTabs("##flow_tabs", labels, glyphs, 3, flowTab);
      const ImVec2 cardSize(0.0f, layout::pageHeight());
      if (flowTab == 0)
        drawStartingMaterialsCard(st, cardSize);
      else if (flowTab == 1)
        drawRouteCard(st, animation, cardSize);
      else
        drawTargetCard(st, cardSize);
    }
    ImGui::EndChild();
  }

  layout::nextRow(pageStartY + pageHeights[0] + page.gap);
  drawResultsCard(st, routeCache, animation, pageHeights[1]);
}

}  // namespace chemcad::ui

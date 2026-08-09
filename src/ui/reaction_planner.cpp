#include "imgui.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "chem/bridge.hpp"
#include "naming/naming.hpp"
#include "rxn/engine.hpp"
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
inline float materialWidth() { return 230.0f * uiScale(); }
inline float materialHeight() { return 228.0f * uiScale(); }
inline ImVec2 routeThumbSize() { return ImVec2(88.0f * uiScale(), 66.0f * uiScale()); }

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
};

struct PlannerAnimationState {
  double searchStartedAt = -1.0;
  bool awaitingResults = false;
  std::unordered_map<std::string, double> routeFirstSeen;
};

float easeOutCubic(float t) {
  t = std::clamp(t, 0.0f, 1.0f);
  const float inverse = 1.0f - t;
  return 1.0f - inverse * inverse * inverse;
}

int resizeStringInput(ImGuiInputTextCallbackData* data) {
  if (data->EventFlag != ImGuiInputTextFlags_CallbackResize) return 0;
  auto* value = static_cast<std::string*>(data->UserData);
  value->resize(static_cast<size_t>(data->BufTextLen));
  data->Buf = value->data();
  return 0;
}

bool stringInputWithHint(const char* label, const char* hint, std::string& value,
                         ImGuiInputTextFlags flags, bool mono = false) {
  flags |= ImGuiInputTextFlags_CallbackResize;
  const bool pushed = mono ? style::pushFont(style::fonts::mono()) : false;
  const bool result = ImGui::InputTextWithHint(label, hint, value.data(),
                                               value.capacity() + 1, flags,
                                               resizeStringInput, &value);
  style::popFont(pushed);
  return result;
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
                              fs * 1.35f);
    ImGui::TableNextColumn();
    const bool pushed = style::pushFont(style::fonts::semibold());
    ImGui::TextColored(headingColor, "%s", heading);
    style::popFont(pushed);
    if (!target) {
      ImGui::TableNextColumn();
      const bool onlyStart = st.planner.starts.size() <= 1;
      ImGui::BeginDisabled(onlyStart);
      ImGui::PushStyleVar(ImGuiStyleVar_DisabledAlpha, 0.35f);
      remove = widgets::iconButton("##remove", icons::Icon::Close,
                                   ImVec2(fs * 1.25f, fs * 1.25f), false,
                                   "Remove this starting material");
      ImGui::PopStyleVar();
      ImGui::EndDisabled();
    }
    ImGui::EndTable();
  }

  const float innerWidth = std::max(ImGui::GetContentRegionAvail().x, 1.0f);
  if (moleculeThumbButton("##preview", box.preview,
                          ImVec2(innerWidth, 67.0f * uiScale())) &&
      box.previewValid) {
    openInSketch(st, box.smiles);
  }
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", box.previewValid ? "Open in the Sketch canvas"
                                             : "Structure preview");

  ImGui::SetNextItemWidth(-1.0f);
  const bool smilesEnter = stringInputWithHint("##smiles", "SMILES", box.smiles,
                                               ImGuiInputTextFlags_EnterReturnsTrue, true);
  const bool smilesCommitted = smilesEnter || ImGui::IsItemDeactivatedAfterEdit();
  if (smilesCommitted) rebuildPreview(box);
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "SMILES string; Enter updates the preview");

  const float lookupWidth = ImGui::CalcTextSize("Look up").x +
                            ImGui::GetStyle().FramePadding.x * 2.0f +
                            ImGui::GetStyle().ItemSpacing.x;
  ImGui::SetNextItemWidth(
      std::max(1.0f, ImGui::GetContentRegionAvail().x - lookupWidth));
  const bool nameEnter = stringInputWithHint("##name", "Chemical name", box.nameInput,
                                             ImGuiInputTextFlags_EnterReturnsTrue);
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "IUPAC or common name; resolved via OPSIN / PubChem");
  ImGui::SameLine();
  ImGui::BeginDisabled(box.status == Status::Loading || box.nameInput.empty());
  const bool lookupClicked = widgets::ghostButton("Look up");
  ImGui::EndDisabled();
  if (nameEnter || lookupClicked) lookupName(st, target, startIndex, box);

  ImGui::BeginDisabled(st.doc.empty());
  if (widgets::ghostButton("From sketch")) {
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
  ImVec4 statusColor = style::col::TextDim;
  if (box.status == Status::Loading) {
    statusText = "looking up...";
  } else if (!box.error.empty()) {
    statusText = "Invalid input";
    statusColor = style::col::Danger;
  } else if (!box.label.empty()) {
    statusText = box.label.c_str();
  }
  if (statusText && ImGui::GetContentRegionAvail().x > ImGui::GetFontSize() * 2.0f) {
    ImGui::SameLine();
    const std::string fitted =
        ellipsize(statusText, std::max(ImGui::GetContentRegionAvail().x, 1.0f));
    ImGui::TextColored(statusColor, "%s", fitted.c_str());
    if (ImGui::IsItemHovered() && fitted != statusText)
      ImGui::SetTooltip("%s", statusText);
    else if (ImGui::IsItemHovered() && !box.error.empty())
      ImGui::SetTooltip("%s", box.error.c_str());
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
      std::max(1.6f, style::metrics().hairline * 2.0f) * thicknessScale;
  dl->AddLine(ImVec2(left, y), ImVec2(right, y), colour, thickness);
  const float head = size.y * 0.28f;
  dl->AddTriangleFilled(ImVec2(right, y), ImVec2(right - head * 1.6f, y - head),
                        ImVec2(right - head * 1.6f, y + head), colour);
  if (animated && right > left) {
    const float progress =
        std::fmod(static_cast<float>(ImGui::GetTime()) * 0.85f, 1.0f);
    const float x = left + (right - left) * progress;
    dl->AddCircleFilled(ImVec2(x, y), std::max(2.5f, ImGui::GetFontSize() * 0.18f),
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
  const float radius = std::max(2.0f, ImGui::GetFontSize() * 0.14f);
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
  drawArrow(ImVec2(width, 34.0f * uiScale()), style::col::Accent, 1.0f,
            st.planner.searching);

  ImGui::BeginDisabled(st.planner.searching);
  if (widgets::primaryButton("Suggest Routes", ImVec2(width, 0.0f)))
    dispatchSearch(st, animation);
  ImGui::EndDisabled();
  if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
    ImGui::SetTooltip("%s", "Search the reaction knowledge base for routes");

  const double elapsed = animation.searchStartedAt < 0.0
                             ? kSearchMinimumSeconds
                             : ImGui::GetTime() - animation.searchStartedAt;
  if (st.planner.searching || elapsed < kSearchMinimumSeconds) {
    drawSearchingIndicator();
  } else {
    ImGui::PushTextWrapPos();
    ImGui::TextColored(style::col::TextDim, "%s",
                       "Knowledge base with optional AI fallback");
    ImGui::PopTextWrapPos();
  }

  const bool aiAvailable = rxn::llmAvailable();
  ImGui::BeginDisabled(!aiAvailable);
  ImGui::Checkbox("Use AI fallback", &st.planner.allowLlm);
  ImGui::EndDisabled();
  if (!aiAvailable && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
    ImGui::SetTooltip("CHEMCAD_LLM_API_KEY is not set.");
  else if (aiAvailable && ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "Let the AI propose a route when the knowledge base stalls");

  const float depthLabel = ImGui::CalcTextSize("Depth").x +
                           ImGui::GetStyle().ItemInnerSpacing.x;
  ImGui::SetNextItemWidth(
      std::max(1.0f, ImGui::GetContentRegionAvail().x - depthLabel));
  ImGui::SliderInt("Depth", &st.planner.maxDepth, 1, 4);
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "Maximum number of reaction steps per route");

  const float routesLabel = ImGui::CalcTextSize("Routes").x +
                            ImGui::GetStyle().ItemInnerSpacing.x;
  ImGui::SetNextItemWidth(
      std::max(1.0f, ImGui::GetContentRegionAvail().x - routesLabel));
  ImGui::SliderInt("Routes", &st.planner.maxRoutes, 1, 10);
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "Maximum number of routes to report");
  ImGui::EndGroup();
}

void stepArrow(const rxn::Step& step, float width) {
  width = std::max(width, 1.0f);
  ImGui::BeginGroup();
  const std::string reagents = joined(step.reagents);
  ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + width);
  if (reagents.empty())
    ImGui::TextDisabled("no added reagent");
  else
    ImGui::TextWrapped("%s", reagents.c_str());
  ImGui::PopTextWrapPos();
  drawArrow(ImVec2(width, 22.0f * uiScale()), style::col::TextDim);
  ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + width);
  if (step.conditions.empty())
    ImGui::TextDisabled("conditions not specified");
  else
    ImGui::TextColored(style::col::TextDim, "%s", step.conditions.c_str());
  ImGui::PopTextWrapPos();
  ImGui::EndGroup();
}

void drawStep(AppState& st, const rxn::Step& step, const RoutePreviewCache& cache,
              int stepIndex) {
  ImGui::PushID(stepIndex);
  std::string heading = "Step " + std::to_string(stepIndex + 1);
  if (!step.reactionName.empty()) heading += "  ·  " + step.reactionName;
  const style::Metrics& metrics = style::metrics();
  const float tickWidth = std::max(2.0f, ImGui::GetFontSize() * 0.16f);
  const float headingWidth = std::max(
      ImGui::GetContentRegionAvail().x - tickWidth - metrics.gap * 1.75f, 1.0f);
  const bool headingFont = style::pushFont(style::fonts::semibold());
  const std::string fittedHeading = ellipsize(heading, headingWidth);
  style::popFont(headingFont);
  widgets::sectionHeader(fittedHeading.c_str(), style::col::Teal);
  if (ImGui::IsItemHovered() && fittedHeading != heading)
    ImGui::SetTooltip("%s", heading.c_str());

  const float rowWidth = std::max(ImGui::GetContentRegionAvail().x, 1.0f);
  const ImVec2 thumbSize(std::min(routeThumbSize().x, rowWidth), routeThumbSize().y);
  const float openButtonWidth = ImGui::CalcTextSize("Open in Sketch").x +
                                ImGui::GetStyle().FramePadding.x * 2.0f;
  const float productWidth = std::max(thumbSize.x, openButtonWidth);
  const float rowGap = style::metrics().gap;
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
    ImGui::TextDisabled("unspecified reactant");
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
  stepArrow(step, arrowWidth);
  if (horizontal) ImGui::SameLine(0.0f, rowGap);

  ImGui::BeginGroup();
  moleculeThumbButton("##product", cachedMolecule(cache, step.productSmiles), thumbSize);
  if (ImGui::IsItemHovered() && !step.productSmiles.empty())
    ImGui::SetTooltip("%s", step.productSmiles.c_str());
  if (widgets::ghostButton("Open in Sketch")) openInSketch(st, step.productSmiles);
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "Append this product to the Sketch canvas");
  ImGui::EndGroup();

  ImGui::Spacing();
  ImGui::TextColored(style::col::Accent, "Side products:");
  ImGui::SameLine();
  if (step.sideProductSmiles.empty()) {
    ImGui::TextDisabled("none predicted");
  } else {
    ImGui::NewLine();
    for (size_t i = 0; i < step.sideProductSmiles.size(); ++i) {
      const std::string& smiles = step.sideProductSmiles[i];
      ImGui::PushID(static_cast<int>(i));
      ImGui::BeginGroup();
      moleculeThumbButton("##side_product", cachedMolecule(cache, smiles),
                          ImVec2(82.0f * uiScale(), 56.0f * uiScale()));
      ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + 92.0f * uiScale());
      const bool mono = style::pushFont(style::fonts::mono());
      ImGui::TextUnformatted(smiles.c_str());
      style::popFont(mono);
      ImGui::PopTextWrapPos();
      ImGui::EndGroup();
      ImGui::PopID();
      const float contentRight =
          ImGui::GetCursorScreenPos().x + ImGui::GetContentRegionAvail().x;
      if (i + 1 < step.sideProductSmiles.size() &&
          ImGui::GetItemRectMax().x + style::metrics().gap +
                  92.0f * uiScale() <=
              contentRight)
        ImGui::SameLine(0.0f, style::metrics().gap);
    }
  }

  if (!step.notes.empty()) {
    ImGui::Spacing();
    ImGui::PushTextWrapPos();
    ImGui::TextColored(style::col::TextDim, "%s", step.notes.c_str());
    ImGui::PopTextWrapPos();
  }
  ImGui::PopID();
}

void drawEmptyState(const char* title, const char* detail, ImVec4 accent) {
  if (!widgets::beginCard("##empty_state", ImVec2(0.0f, 0.0f), style::col::BgDeep))
    return;
  const bool pushed = style::pushFont(style::fonts::semibold());
  ImGui::TextColored(accent, "%s", title);
  style::popFont(pushed);
  ImGui::PushTextWrapPos();
  ImGui::TextColored(style::col::TextDim, "%s", detail);
  ImGui::PopTextWrapPos();
  widgets::endCard();
}

void drawResults(AppState& st, RoutePreviewCache& cache,
                 PlannerAnimationState& animation) {
  const std::string signature = routesSignature(st.planner.routes);
  if (signature != cache.signature) rebuildRouteCache(cache, st.planner.routes, signature);

  if (animation.awaitingResults && !st.planner.searching && st.planner.searched) {
    animation.routeFirstSeen.clear();
    animation.awaitingResults = false;
  }

  if (!st.planner.error.empty()) {
    drawEmptyState("CHECK THE ROUTE INPUTS", st.planner.error.c_str(), style::col::Danger);
    ImGui::Spacing();
  }

  if (st.planner.target.smiles.empty()) {
    drawEmptyState(
        "DEFINE A TARGET",
        "Add one or more starting materials, enter the desired product in the Target "
        "card, then use Suggest Routes to search the reaction knowledge base.",
        style::col::Accent);
    return;
  }

  if (st.planner.searching && st.planner.routes.empty()) {
    drawEmptyState("SEARCH IN PROGRESS",
                   "Candidate transformations are being assembled and ranked. Results "
                   "will appear here as soon as the search completes.",
                   style::col::Violet);
    return;
  }

  if (st.planner.searched && st.planner.routes.empty()) {
    const bool canUseAi = st.planner.allowLlm && rxn::llmAvailable();
    drawEmptyState(
        "NO ROUTES FOUND",
        canUseAi
            ? "No route matched the current starting materials, target, and depth. "
              "Try a deeper search or revise the structures."
            : "No knowledge-base route matched. Enable AI fallback when available, "
              "increase the search depth, or revise the structures.",
        style::col::TextDim);
    return;
  }

  if (st.planner.routes.empty()) {
    drawEmptyState("READY TO SEARCH",
                   "Review the structures and route limits above, then choose Suggest "
                   "Routes to generate ranked candidates.",
                   style::col::Teal);
    return;
  }

  ImGui::TextDisabled("%d ranked candidate%s", static_cast<int>(st.planner.routes.size()),
                      st.planner.routes.size() == 1 ? "" : "s");
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
    ImGui::Dummy(ImVec2(0.0f, (1.0f - enter) * 8.0f * uiScale()));

    const ImGuiID hoverId = ImGui::GetID("##route_hover");
    const float previousHover = ImGui::GetStateStorage()->GetFloat(hoverId, 0.0f);
    const float liftReserve = 3.0f * uiScale();
    ImGui::Dummy(ImVec2(0.0f, liftReserve * (1.0f - previousHover)));

    const bool cardOpen =
        widgets::beginCard("##route_result", ImVec2(0.0f, 0.0f), style::col::BgSurface);
    if (cardOpen) {
      std::string title = "Route " + std::to_string(routeIndex + 1);
      const bool open = ImGui::CollapsingHeader(
          title.c_str(), routeIndex == 0 ? ImGuiTreeNodeFlags_DefaultOpen : 0);
      const ImVec2 headerMin = ImGui::GetItemRectMin();
      const ImVec2 headerMax = ImGui::GetItemRectMax();
      const bool headerHovered = ImGui::IsItemHovered();
      const float rowHover = widgets::hoverT(ImGui::GetItemID(), headerHovered);
      ImGui::GetWindowDrawList()->AddRect(
          headerMin, headerMax,
          style::mix(style::col::Border, style::col::Teal, rowHover),
          style::metrics().radiusSm, 0,
          style::metrics().hairline * (1.0f + rowHover));
      if (headerHovered)
        ImGui::SetTooltip("%s",
                          routeIndex == 0 ? "Highest-ranked route" : "Alternative route");

      if (open) {
        std::string stepsText = std::to_string(route.steps.size()) + " step";
        if (route.steps.size() != 1) stepsText += 's';
        widgets::badge(stepsText.c_str(), style::col::TextDim);
        ImGui::SameLine();
        widgets::badge(route.usesLlm() ? "AI" : "KB",
                       route.usesLlm() ? style::col::Violet : style::col::Teal);
        if (routeIndex == 0) {
          ImGui::SameLine();
          widgets::badge("BEST", style::col::Accent);
        }

        for (size_t stepIndex = 0; stepIndex < route.steps.size(); ++stepIndex) {
          drawStep(st, route.steps[stepIndex], cache, static_cast<int>(stepIndex));
          if (stepIndex + 1 < route.steps.size()) ImGui::Spacing();
        }
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
    widgets::sectionHeader("STARTING MATERIALS", style::col::Teal);
    if (widgets::ghostButton("+ Add starting material")) st.planner.starts.emplace_back();
    if (ImGui::IsItemHovered())
      ImGui::SetTooltip("%s", "Reactions can combine several starting materials");
    const float noteWidth =
        ImGui::CalcTextSize("Combine one or more inputs").x + style::metrics().gap;
    if (ImGui::GetContentRegionAvail().x > noteWidth) {
      ImGui::SameLine(0.0f, style::metrics().gap);
      ImGui::TextDisabled("Combine one or more inputs");
    }

    if (ImGui::BeginChild("##starting_materials_row", ImVec2(0.0f, 0.0f),
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
    ImGui::EndChild();
    widgets::endCard();
  }

  if (removeIndex >= 0 && st.planner.starts.size() > 1)
    st.planner.starts.erase(st.planner.starts.begin() + removeIndex);
}

void drawRouteCard(AppState& st, PlannerAnimationState& animation, float cardHeight) {
  if (!widgets::beginCard("##route_card", ImVec2(0.0f, cardHeight),
                          style::col::BgSurface))
    return;
  widgets::sectionHeader("ROUTE", style::col::Accent);
  connectorControls(st, animation);
  widgets::endCard();
}

void drawTargetCard(AppState& st, float cardHeight) {
  if (!widgets::beginCard("##target_card", ImVec2(0.0f, cardHeight),
                          style::col::BgSurface))
    return;
  widgets::sectionHeader("TARGET", style::col::Violet);
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
  widgets::sectionHeader("RESULTS", style::col::Violet);
  drawResults(st, cache, animation);
  widgets::endCard();
}

}  // namespace

void drawReactionPlanner(AppState& st) {
  if (st.planner.starts.empty()) st.planner.starts.emplace_back();

  static RoutePreviewCache routeCache;
  static PlannerAnimationState animation;
  const float flowCardHeight = materialHeight() + 92.0f * uiScale();
  const float availableWidth = ImGui::GetContentRegionAvail().x;
  const bool horizontalFlow = availableWidth >= 880.0f * uiScale();

  if (horizontalFlow) {
    constexpr ImGuiTableFlags flowFlags =
        ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_NoSavedSettings;
    if (ImGui::BeginTable("##reaction_flow", 3, flowFlags)) {
      ImGui::TableSetupColumn("##starts", ImGuiTableColumnFlags_WidthStretch, 1.40f);
      ImGui::TableSetupColumn("##route", ImGuiTableColumnFlags_WidthStretch, 0.85f);
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

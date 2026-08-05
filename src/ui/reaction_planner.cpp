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
// Width that fits `label` inside a regular framed button.
inline float buttonWidthFor(const char* label) {
  return ImGui::CalcTextSize(label).x + ImGui::GetStyle().FramePadding.x * 2.0f +
         ImGui::GetStyle().ItemSpacing.x;
}

struct RoutePreviewCache {
  std::string signature;
  std::unordered_map<std::string, core::Molecule> molecules;
};

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

std::string routesSignature(const std::vector<rxn::Route>& routes) {
  std::string signature;
  signature.reserve(routes.size() * 128);
  signature += std::to_string(routes.size());
  for (const rxn::Route& route : routes) {
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
  const ImVec2 cardSize(materialWidth(), materialHeight());
  if (!widgets::beginCard("##material_box", cardSize)) return false;

  const float innerWidth = ImGui::GetContentRegionAvail().x;
  const char* heading = target ? "PRODUCT" : "STARTING MATERIAL";
  const ImVec4 headingColor = target ? style::col::Accent : style::col::TextDim;

  // Header row: accent tick + label, close chip on the right for starts.
  const float fs = ImGui::GetFontSize();
  {
    const ImVec2 pos = ImGui::GetCursorScreenPos();
    ImGui::Dummy(ImVec2(innerWidth, fs));
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const float tickW = std::max(2.0f, fs * 0.16f);
    const float tickH = fs * 0.85f;
    const float tickY = pos.y + (fs - tickH) * 0.5f + 1.0f;
    dl->AddRectFilled(ImVec2(pos.x, tickY), ImVec2(pos.x + tickW, tickY + tickH),
                      style::u32(headingColor), tickW * 0.5f);
    const bool pushed = style::pushFont(style::fonts::semibold());
    dl->AddText(ImVec2(pos.x + tickW + style::metrics().gap * 0.6f, pos.y),
                style::u32(headingColor), heading);
    style::popFont(pushed);
  }
  if (!target) {
    const float chip = fs * 1.25f;
    ImGui::SameLine();
    ImGui::SetCursorPosX(ImGui::GetCursorPosX() + innerWidth - chip -
                         ImGui::GetStyle().WindowPadding.x * 0.5f);
    ImGui::SetCursorPosY(ImGui::GetCursorPosY() - fs * 0.9f);
    const bool onlyStart = st.planner.starts.size() <= 1;
    ImGui::BeginDisabled(onlyStart);
    ImGui::PushStyleVar(ImGuiStyleVar_DisabledAlpha, 0.35f);
    remove = widgets::iconButton("##remove", icons::Icon::Close, ImVec2(chip, chip),
                                 false, "Remove this starting material");
    ImGui::PopStyleVar();
    ImGui::EndDisabled();
  }

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

  // Reserve exactly what the adjacent "Look up" button needs so neither clips.
  const float lookupWidth = ImGui::CalcTextSize("Look up").x +
                            ImGui::GetStyle().FramePadding.x * 2.0f +
                            ImGui::GetStyle().ItemSpacing.x * 2.0f;
  ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x - lookupWidth);
  const bool nameEnter = stringInputWithHint("##name", "Chemical name", box.nameInput,
                                             ImGuiInputTextFlags_EnterReturnsTrue);
  const bool nameCommitted = nameEnter;
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "IUPAC or common name; resolved via OPSIN / PubChem");
  ImGui::SameLine();
  ImGui::BeginDisabled(box.status == Status::Loading || box.nameInput.empty());
  const bool lookupClicked = widgets::ghostButton("Look up");
  ImGui::EndDisabled();
  if (nameCommitted || lookupClicked) lookupName(st, target, startIndex, box);

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
  ImGui::SameLine();
  if (box.status == Status::Loading) {
    ImGui::TextDisabled("looking up...");
  } else if (!box.error.empty()) {
    ImGui::TextColored(style::col::Danger, "%s", box.error.c_str());
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", box.error.c_str());
  } else if (!box.label.empty()) {
    ImGui::TextDisabled("%s", box.label.c_str());
  }

  widgets::endCard();
  return remove;
}

void drawArrow(ImVec2 size, ImVec4 color, float thicknessScale = 1.0f) {
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::Dummy(size);
  ImDrawList* dl = ImGui::GetWindowDrawList();
  const ImU32 colour = style::u32(color);
  const float y = min.y + size.y * 0.5f;
  const float pad = size.x * 0.08f;
  const float left = min.x + pad;
  const float right = min.x + size.x - pad;
  const float t = std::max(1.6f, style::metrics().hairline * 2.0f) * thicknessScale;
  dl->AddLine(ImVec2(left, y), ImVec2(right, y), colour, t);
  const float head = size.y * 0.28f;
  dl->AddTriangleFilled(ImVec2(right, y), ImVec2(right - head * 1.6f, y - head),
                        ImVec2(right - head * 1.6f, y + head), colour);
}

void dispatchSearch(AppState& st) {
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

  st.tasks.run<std::vector<rxn::Route>>(
      [request] { return rxn::suggestRoutes(request); },
      [&st](std::vector<rxn::Route> routes) {
        st.planner.routes = std::move(routes);
        st.planner.searching = false;
        st.planner.searched = true;
      });
}

void connectorControls(AppState& st) {
  // The column must fit its widest control ("Use AI fallback" checkbox and the
  // Suggest button), otherwise their labels clip at non-default UI scales.
  const float buttonWidth = buttonWidthFor("Suggest Routes");
  const float checkboxWidth = ImGui::CalcTextSize("Use AI fallback").x +
                              ImGui::GetFrameHeight() + ImGui::GetStyle().ItemSpacing.x * 2.0f;
  const float kWidth = std::max({190.0f * uiScale(), buttonWidth, checkboxWidth});
  ImGui::BeginGroup();
  drawArrow(ImVec2(kWidth, 42.0f * uiScale()), style::col::Accent);
  ImGui::SetCursorPosX(ImGui::GetCursorPosX() + (kWidth - buttonWidth) * 0.5f);
  ImGui::BeginDisabled(st.planner.searching);
  if (widgets::primaryButton("Suggest Routes", ImVec2(buttonWidth, 0.0f)))
    dispatchSearch(st);
  ImGui::EndDisabled();
  if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
    ImGui::SetTooltip("%s", "Search the reaction knowledge base for routes");
  if (st.planner.searching) {
    ImGui::SetCursorPosX(ImGui::GetCursorPosX() + (kWidth - ImGui::CalcTextSize("searching...").x) * 0.5f);
    ImGui::TextDisabled("searching...");
  }

  const bool aiAvailable = rxn::llmAvailable();
  ImGui::BeginDisabled(!aiAvailable);
  ImGui::Checkbox("Use AI fallback", &st.planner.allowLlm);
  ImGui::EndDisabled();
  if (!aiAvailable && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
    ImGui::SetTooltip("CHEMCAD_LLM_API_KEY is not set.");
  else if (aiAvailable && ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "Let the AI propose a route when the knowledge base stalls");

  ImGui::SetNextItemWidth(kWidth - ImGui::CalcTextSize("Routes").x -
                          ImGui::GetStyle().ItemSpacing.x * 2.0f);
  ImGui::SliderInt("Depth", &st.planner.maxDepth, 1, 4);
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "Maximum number of reaction steps per route");
  ImGui::SetNextItemWidth(kWidth - ImGui::CalcTextSize("Routes").x -
                          ImGui::GetStyle().ItemSpacing.x * 2.0f);
  ImGui::SliderInt("Routes", &st.planner.maxRoutes, 1, 10);
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "Maximum number of routes to report");
  ImGui::EndGroup();
}

void stepArrow(const rxn::Step& step) {
  const float kWidth = 185.0f * uiScale();
  ImGui::BeginGroup();
  const std::string reagents = joined(step.reagents);
  ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + kWidth);
  if (reagents.empty())
    ImGui::TextDisabled("no added reagent");
  else
    ImGui::TextWrapped("%s", reagents.c_str());
  ImGui::PopTextWrapPos();
  drawArrow(ImVec2(kWidth, 22.0f * uiScale()), style::col::TextDim);
  ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + kWidth);
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
  widgets::sectionHeader(heading.c_str(), style::col::Teal);

  ImGui::BeginGroup();
  if (step.reactantSmiles.empty()) {
    ImGui::TextDisabled("unspecified reactant");
  } else {
    for (size_t i = 0; i < step.reactantSmiles.size(); ++i) {
      ImGui::PushID(static_cast<int>(i));
      moleculeThumbButton("##reactant", cachedMolecule(cache, step.reactantSmiles[i]),
                          routeThumbSize());
      if (ImGui::IsItemHovered() && !step.reactantSmiles[i].empty())
        ImGui::SetTooltip("%s", step.reactantSmiles[i].c_str());
      ImGui::PopID();
      if (i + 1 < step.reactantSmiles.size()) {
        ImGui::SameLine();
        ImGui::AlignTextToFramePadding();
        ImGui::TextUnformatted("+");
        ImGui::SameLine();
      }
    }
  }
  ImGui::EndGroup();
  ImGui::SameLine(0.0f, 12.0f);
  stepArrow(step);
  ImGui::SameLine(0.0f, 12.0f);

  ImGui::BeginGroup();
  moleculeThumbButton("##product", cachedMolecule(cache, step.productSmiles), routeThumbSize());
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
      if (i + 1 < step.sideProductSmiles.size()) ImGui::SameLine(0.0f, 12.0f);
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

void drawResults(AppState& st, RoutePreviewCache& cache) {
  const std::string signature = routesSignature(st.planner.routes);
  if (signature != cache.signature) rebuildRouteCache(cache, st.planner.routes, signature);

  if (!st.planner.error.empty())
    ImGui::TextColored(style::col::Danger, "%s", st.planner.error.c_str());

  if (st.planner.searched && st.planner.routes.empty()) {
    ImGui::Spacing();
    ImGui::TextUnformatted("No routes found.");
    if (!st.planner.allowLlm || !rxn::llmAvailable()) {
      ImGui::TextDisabled(
          "Enable AI fallback or set CHEMCAD_LLM_API_KEY to let the model suggest a route.");
    }
    return;
  }
  if (st.planner.routes.empty()) return;

  widgets::sectionHeader("Suggested Routes", style::col::Teal);
  ImGui::BeginChild("##route_results", ImVec2(0.0f, 0.0f), ImGuiChildFlags_None,
                    ImGuiWindowFlags_HorizontalScrollbar);
  for (size_t routeIndex = 0; routeIndex < st.planner.routes.size(); ++routeIndex) {
    const rxn::Route& route = st.planner.routes[routeIndex];
    ImGui::PushID(static_cast<int>(routeIndex));
    std::string title = "Route " + std::to_string(routeIndex + 1);
    const bool open = ImGui::CollapsingHeader(
        title.c_str(), routeIndex == 0 ? ImGuiTreeNodeFlags_DefaultOpen : 0);
    if (ImGui::IsItemHovered())
      ImGui::SetTooltip("%s", routeIndex == 0 ? "Highest-ranked route" : "Alternative route");

    ImGui::SameLine();
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

    if (open) {
      if (widgets::beginCard("##route_card", ImVec2(0.0f, 0.0f))) {
        for (size_t stepIndex = 0; stepIndex < route.steps.size(); ++stepIndex) {
          drawStep(st, route.steps[stepIndex], cache, static_cast<int>(stepIndex));
          if (stepIndex + 1 < route.steps.size()) ImGui::Spacing();
        }
      }
      widgets::endCard();
    }
    ImGui::PopID();
    ImGui::Spacing();
  }
  ImGui::EndChild();
}

}  // namespace

void drawReactionPlanner(AppState& st) {
  if (st.planner.starts.empty()) st.planner.starts.emplace_back();

  if (widgets::ghostButton("+ Add starting material")) st.planner.starts.emplace_back();
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "Reactions can combine several starting materials");
  ImGui::SameLine();
  ImGui::TextDisabled("Build a route from one or more inputs to one product");

  int removeIndex = -1;
  ImGui::BeginChild("##reaction_materials", ImVec2(0.0f, materialHeight() + 20.0f),
                    ImGuiChildFlags_None, ImGuiWindowFlags_HorizontalScrollbar);
  for (size_t i = 0; i < st.planner.starts.size(); ++i) {
    ImGui::PushID(static_cast<int>(i));
    if (materialBoxWidget(st, st.planner.starts[i], false, static_cast<int>(i)))
      removeIndex = static_cast<int>(i);
    ImGui::PopID();
    ImGui::SameLine(0.0f, 9.0f);
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted("+");
    ImGui::SameLine(0.0f, 9.0f);
  }

  ImGui::PushID("route_controls");
  connectorControls(st);
  ImGui::PopID();
  ImGui::SameLine(0.0f, 10.0f);
  ImGui::PushID("target");
  materialBoxWidget(st, st.planner.target, true, -1);
  ImGui::PopID();
  ImGui::EndChild();

  if (removeIndex >= 0 && st.planner.starts.size() > 1)
    st.planner.starts.erase(st.planner.starts.begin() + removeIndex);

  static RoutePreviewCache routeCache;
  drawResults(st, routeCache);
}

}  // namespace chemcad::ui

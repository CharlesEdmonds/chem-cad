#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include "imgui.h"

#include "chem/bridge.hpp"
#include "naming/naming.hpp"
#include "ui/ui.hpp"

namespace chemcad::ui {
namespace {

using Clock = std::chrono::steady_clock;

struct BuildNameState {
  std::array<char, 256> input{};
  Status status = Status::Idle;
  std::string error;
};

BuildNameState& buildNameState() {
  static BuildNameState state;
  return state;
}

core::Molecule flatten(const core::Document& document) {
  core::Molecule result;
  for (const auto& fragment : document.molecules) {
    std::unordered_map<core::AtomId, core::AtomId> ids;
    ids.reserve(fragment.atomCount());
    for (const auto& atom : fragment.atoms()) ids.emplace(atom.id, result.addAtom(atom));
    for (const auto& bond : fragment.bonds()) {
      const core::BondId id = result.addBond(ids.at(bond.a), ids.at(bond.b), bond.order);
      if (core::Bond* added = result.bond(id)) added->stereo = bond.stereo;
    }
  }
  return result;
}

void offsetForAppend(const core::Document& document, core::Molecule& molecule) {
  if (molecule.empty()) return;

  bool hasExisting = false;
  float existingMaxX = -std::numeric_limits<float>::infinity();
  for (const auto& fragment : document.molecules) {
    for (const auto& atom : fragment.atoms()) {
      hasExisting = true;
      existingMaxX = std::max(existingMaxX, atom.pos.x);
    }
  }
  if (!hasExisting) return;

  float newMinX = std::numeric_limits<float>::infinity();
  for (const auto& atom : molecule.atoms()) newMinX = std::min(newMinX, atom.pos.x);
  const float dx = existingMaxX - newMinX + 2.0f;
  for (auto& atom : molecule.mutableAtoms()) atom.pos.x += dx;
}

void recomputeProperties(AppState& st) {
  const std::string previousSmiles = st.props.smiles;
  st.props.chemError.clear();
  st.props.formula.clear();
  st.props.smiles.clear();
  st.props.mw = 0;
  st.props.logP = 0;
  st.props.rings = 0;

  try {
    core::Molecule molecule = flatten(st.doc);
    if (!molecule.empty()) {
      const chem::Properties properties = chem::computeProperties(molecule);
      st.props.formula = properties.formula;
      st.props.mw = properties.mw;
      st.props.logP = properties.logP;
      st.props.rings = properties.rings;
      st.props.smiles = chem::toSmiles(molecule);
    }
  } catch (const chem::ChemError& e) {
    st.props.chemError = e.what();
  } catch (const std::exception& e) {
    st.props.chemError = e.what();
  }

  st.props.computedForRevision = st.docRevision;
  if (st.props.smiles != previousSmiles) {
    st.props.name.clear();
    st.props.nameError.clear();
    st.props.nameStatus = Status::Idle;
  } else if (!st.props.smiles.empty() &&
             (st.props.nameStatus == Status::Ok || st.props.nameStatus == Status::Error)) {
    // Geometry-only edits preserve a completed lookup without issuing another request.
    st.props.nameRequestedForRevision = st.docRevision;
  } else if (st.props.smiles.empty()) {
    st.props.name.clear();
    st.props.nameError.clear();
    st.props.nameStatus = Status::Idle;
  }
}

void requestAutomaticName(AppState& st) {
  const uint64_t revision = st.docRevision;
  const std::string smiles = st.props.smiles;
  st.props.nameRequestedForRevision = revision;
  st.props.nameStatus = Status::Loading;
  st.props.name.clear();
  st.props.nameError.clear();
  st.tasks.run<naming::Result>(
      [smiles] { return naming::smilesToName(smiles); },
      [&st, smiles, revision](naming::Result result) {
        if (st.props.smiles != smiles) return;
        if (st.props.computedForRevision != st.docRevision) {
          if (st.props.nameRequestedForRevision == revision) {
            st.props.nameStatus = Status::Idle;
            st.props.nameRequestedForRevision = ~0ull;
          }
          return;
        }
        st.props.nameRequestedForRevision = st.docRevision;
        if (result.ok) {
          st.props.name = std::move(result.value);
          st.props.nameError.clear();
          st.props.nameStatus = Status::Ok;
        } else {
          st.props.name.clear();
          st.props.nameError =
              result.error.empty() ? "Name lookup was unavailable." : std::move(result.error);
          st.props.nameStatus = Status::Error;
        }
      });
}

void submitBuild(AppState& st) {
  BuildNameState& build = buildNameState();
  const std::string name = build.input.data();
  if (name.empty() || build.status == Status::Loading) return;

  build.status = Status::Loading;
  build.error.clear();
  st.tasks.run<naming::Result>(
      [name] { return naming::nameToSmiles(name); },
      [&st, name](naming::Result result) {
        BuildNameState& state = buildNameState();
        if (!result.ok) {
          state.status = Status::Error;
          state.error =
              result.error.empty() ? "Could not resolve that name." : std::move(result.error);
          return;
        }
        try {
          core::Molecule molecule = chem::fromSmiles(result.value);
          if (molecule.empty()) throw chem::ChemError("name resolved to an empty structure");
          offsetForAppend(st.doc, molecule);
          st.snapshot();
          st.doc.molecules.push_back(std::move(molecule));
          st.touch();
          st.sel.clear();
          state.status = Status::Ok;
          state.error.clear();
          state.input.fill('\0');
          st.statusMessage = "Built structure for " + name;
        } catch (const std::exception& e) {
          state.status = Status::Error;
          state.error = e.what();
        }
      });
}

}  // namespace

void drawPropertiesPanel(AppState& st) {
  const Clock::time_point now = Clock::now();
  const bool propertyDelayPassed =
      now - st.props.lastEdit >= std::chrono::milliseconds(250);
  if (st.props.computedForRevision != st.docRevision && propertyDelayPassed) {
    recomputeProperties(st);
  }

  ImGui::SeparatorText("Properties");
  if (!st.props.chemError.empty()) {
    ImGui::TextColored(ImVec4(0.95f, 0.35f, 0.32f, 1.0f), "%s", st.props.chemError.c_str());
  } else if (st.props.computedForRevision != st.docRevision) {
    ImGui::TextDisabled("Updating...");
  } else if (st.props.smiles.empty()) {
    ImGui::TextDisabled("Draw a structure to see its properties.");
  } else {
    ImGui::Text("Formula: %s", st.props.formula.c_str());
    ImGui::Text("MW: %.2f", st.props.mw);
    ImGui::Text("cLogP: %.2f", st.props.logP);
    ImGui::Text("Rings: %d", st.props.rings);

    std::vector<char> smiles(st.props.smiles.begin(), st.props.smiles.end());
    smiles.push_back('\0');
    const float copyWidth = ImGui::CalcTextSize("Copy").x + ImGui::GetStyle().FramePadding.x * 2.0f;
    ImGui::SetNextItemWidth(std::max(50.0f, ImGui::GetContentRegionAvail().x - copyWidth -
                                               ImGui::GetStyle().ItemSpacing.x));
    ImGui::InputText("##canonical_smiles", smiles.data(), smiles.size(),
                     ImGuiInputTextFlags_ReadOnly | ImGuiInputTextFlags_AutoSelectAll);
    ImGui::SameLine();
    if (ImGui::Button("Copy")) {
      st.clipboardSmiles = st.props.smiles;
      ImGui::SetClipboardText(st.props.smiles.c_str());
      st.statusMessage = "Canonical SMILES copied";
    }
  }

  ImGui::SeparatorText("Name");
  ImGui::Checkbox("Auto-name", &st.props.autoName);
  const bool nameDelayPassed = now - st.props.lastEdit >= std::chrono::milliseconds(1500);
  if (st.props.autoName && !st.props.smiles.empty() &&
      st.props.computedForRevision == st.docRevision &&
      st.props.nameRequestedForRevision != st.docRevision && nameDelayPassed) {
    requestAutomaticName(st);
  }

  if (st.props.nameStatus == Status::Loading) {
    ImGui::TextDisabled("Looking up...");
  } else if (st.props.nameStatus == Status::Ok) {
    ImGui::TextWrapped("%s", st.props.name.c_str());
  } else if (st.props.nameStatus == Status::Error) {
    ImGui::TextColored(ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled), "%s",
                       st.props.nameError.c_str());
  }

  ImGui::Spacing();
  ImGui::TextUnformatted("Name to structure");
  BuildNameState& build = buildNameState();
  ImGui::SetNextItemWidth(-1.0f);
  const bool enter = ImGui::InputTextWithHint(
      "##name_to_structure", "e.g. acetylsalicylic acid", build.input.data(),
      build.input.size(), ImGuiInputTextFlags_EnterReturnsTrue);
  if (build.status == Status::Loading) ImGui::BeginDisabled();
  if (enter || ImGui::Button("Build")) submitBuild(st);
  if (build.status == Status::Loading) {
    ImGui::EndDisabled();
    ImGui::SameLine();
    ImGui::TextDisabled("Resolving...");
  }
  if (build.status == Status::Error && !build.error.empty()) {
    ImGui::TextColored(ImVec4(0.95f, 0.35f, 0.32f, 1.0f), "%s", build.error.c_str());
  }
}

}  // namespace chemcad::ui

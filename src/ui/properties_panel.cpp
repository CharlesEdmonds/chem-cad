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
#include "ui/icons.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"
#include "ui/widgets.hpp"

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

  if (!st.props.chemError.empty()) {
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(style::col::Danger, "%s", st.props.chemError.c_str());
    ImGui::PopTextWrapPos();
  } else if (st.props.computedForRevision != st.docRevision) {
    ImGui::TextDisabled("Updating...");
  } else if (st.props.smiles.empty()) {
    ImGui::TextWrapped("Draw a structure to see its properties.");
  } else {
    // Identity dashboard: 2x2 stat grid, mono values.
    const float avail = std::max(ImGui::GetContentRegionAvail().x, 1.0f);
    const float spacing = ImGui::GetStyle().ItemSpacing.x;
    const bool twoColumns = avail >= ImGui::GetFontSize() * 14.0f;
    const float cardW = twoColumns ? std::max((avail - spacing) * 0.5f, 1.0f) : avail;
    const float cardH = ImGui::GetFontSize() * 3.1f;

    char mw[32];
    std::snprintf(mw, sizeof(mw), "%.2f", st.props.mw);
    char logp[32];
    std::snprintf(logp, sizeof(logp), "%.2f", st.props.logP);
    char rings[16];
    std::snprintf(rings, sizeof(rings), "%d", st.props.rings);

    const float valueWidth =
        std::max(cardW - style::metrics().gap * 1.8f, 1.0f);
    const bool formulaFont = style::pushFont(style::fonts::mono());
    const std::string fittedFormula = ellipsize(st.props.formula, valueWidth);
    style::popFont(formulaFont);
    widgets::statCard("FORMULA", fittedFormula.c_str(), ImVec2(cardW, cardH));
    if (ImGui::IsItemHovered() && fittedFormula != st.props.formula)
      ImGui::SetTooltip("%s", st.props.formula.c_str());
    if (twoColumns) ImGui::SameLine(0.0f, spacing);
    widgets::statCard("MW g/mol", mw, ImVec2(cardW, cardH));
    widgets::statCard("CLOGP", logp, ImVec2(cardW, cardH));
    if (twoColumns) ImGui::SameLine(0.0f, spacing);
    widgets::statCard("RINGS", rings, ImVec2(cardW, cardH));

    ImGui::Spacing();
    std::vector<char> smiles(st.props.smiles.begin(), st.props.smiles.end());
    smiles.push_back('\0');
    const float chip = ImGui::GetFrameHeight();
    ImGui::SetNextItemWidth(std::max(1.0f, ImGui::GetContentRegionAvail().x - chip -
                                               ImGui::GetStyle().ItemSpacing.x));
    const bool mono = style::pushFont(style::fonts::mono());
    ImGui::InputText("##canonical_smiles", smiles.data(), smiles.size(),
                     ImGuiInputTextFlags_ReadOnly | ImGuiInputTextFlags_AutoSelectAll);
    style::popFont(mono);
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", "Canonical SMILES");
    ImGui::SameLine();
    if (widgets::iconButton("##copy_smiles", icons::Icon::Copy, ImVec2(chip, chip), false,
                            "Copy canonical SMILES")) {
      st.clipboardSmiles = st.props.smiles;
      ImGui::SetClipboardText(st.props.smiles.c_str());
      st.statusMessage = "Canonical SMILES copied";
    }
  }

  widgets::sectionHeader("Name");
  ImGui::Checkbox("Auto-name", &st.props.autoName);
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "Resolve the sketched structure's IUPAC name automatically");
  const bool nameDelayPassed = now - st.props.lastEdit >= std::chrono::milliseconds(1500);
  if (st.props.autoName && !st.props.smiles.empty() &&
      st.props.computedForRevision == st.docRevision &&
      st.props.nameRequestedForRevision != st.docRevision && nameDelayPassed) {
    requestAutomaticName(st);
  }

  if (st.props.nameStatus == Status::Loading) {
    ImGui::TextDisabled("Looking up...");
  } else if (st.props.nameStatus == Status::Ok) {
    const bool pushed = style::pushFont(style::fonts::semibold());
    ImGui::PushTextWrapPos();
    ImGui::TextUnformatted(st.props.name.c_str());
    ImGui::PopTextWrapPos();
    style::popFont(pushed);
  } else if (st.props.nameStatus == Status::Error) {
    ImGui::TextWrapped("%s", st.props.nameError.c_str());
  }

  widgets::sectionHeader("Build from name");
  BuildNameState& build = buildNameState();
  ImGui::SetNextItemWidth(-1.0f);
  const bool enter = ImGui::InputTextWithHint(
      "##name_to_structure", "e.g. acetylsalicylic acid", build.input.data(),
      build.input.size(), ImGuiInputTextFlags_EnterReturnsTrue);
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "IUPAC or common name; Enter builds the structure");
  if (build.status == Status::Loading) ImGui::BeginDisabled();
  if (enter || widgets::primaryButton("Build")) submitBuild(st);
  if (build.status == Status::Loading) {
    ImGui::EndDisabled();
    const float statusWidth = ImGui::CalcTextSize("Resolving...").x;
    if (ImGui::GetContentRegionAvail().x >=
        statusWidth + ImGui::GetStyle().ItemSpacing.x)
      ImGui::SameLine();
    ImGui::TextDisabled("Resolving...");
  }
  if (build.status == Status::Error && !build.error.empty()) {
    ImGui::TextWrapped("%s", build.error.c_str());
  }
}

}  // namespace chemcad::ui

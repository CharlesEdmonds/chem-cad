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
#include "ui/charts.hpp"
#include "ui/element_data.hpp"
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

struct DocumentSummary {
  std::size_t explicitAtoms = 0;
  std::size_t bonds = 0;
  int implicitHydrogens = 0;
  int formalCharge = 0;
  std::unordered_map<uint8_t, int> elementCounts;

  std::size_t totalAtoms() const {
    return explicitAtoms + static_cast<std::size_t>(std::max(implicitHydrogens, 0));
  }
};

struct CompositionEntry {
  uint8_t atomicNumber = 0;
  int count = 0;
  std::string symbol;
  std::string name;
  ImVec4 colour = style::col::Teal;
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

  std::size_t end = text.size();
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

DocumentSummary summariseDocument(const core::Document& document) {
  DocumentSummary summary;
  for (const auto& fragment : document.molecules) {
    summary.bonds += fragment.bondCount();
    for (const auto& atom : fragment.atoms()) {
      ++summary.explicitAtoms;
      summary.formalCharge += atom.charge;
      ++summary.elementCounts[atom.atomicNumber];
      const int implicit = chem::implicitHCount(fragment, atom.id);
      if (implicit > 0) {
        summary.implicitHydrogens += implicit;
        summary.elementCounts[1] += implicit;
      }
    }
  }
  return summary;
}
const DocumentSummary& cachedDocumentSummary(const AppState& st) {
  struct Cache {
    const AppState* owner = nullptr;
    uint64_t revision = ~0ull;
    DocumentSummary summary;
  };
  static Cache cache;
  if (cache.owner != &st || cache.revision != st.docRevision) {
    cache.owner = &st;
    cache.revision = st.docRevision;
    cache.summary = summariseDocument(st.doc);
  }
  return cache.summary;
}

ImVec4 compositionColour(std::size_t index) {
  switch (index % 6) {
    case 0:
      return style::col::Accent;
    case 1:
      return style::col::Teal;
    case 2:
      return style::col::Violet;
    case 3:
      return style::col::Success;
    case 4:
      return style::col::Danger;
    default:
      return style::col::AccentHover;
  }
}

std::vector<CompositionEntry> makeComposition(const DocumentSummary& summary) {
  std::vector<CompositionEntry> entries;
  entries.reserve(summary.elementCounts.size());
  for (const auto& [atomicNumber, count] : summary.elementCounts) {
    if (count <= 0) continue;
    CompositionEntry entry;
    entry.atomicNumber = atomicNumber;
    entry.count = count;
    if (const ElementData* element = findElement(atomicNumber)) {
      entry.symbol = element->symbol;
      entry.name = element->name;
    } else {
      const char* symbol = chem::symbolFor(atomicNumber);
      entry.symbol = symbol && *symbol ? symbol : ("Z" + std::to_string(atomicNumber));
      entry.name = "Atomic number " + std::to_string(atomicNumber);
    }
    entries.push_back(std::move(entry));
  }
  std::sort(entries.begin(), entries.end(), [](const CompositionEntry& a,
                                                const CompositionEntry& b) {
    if (a.count != b.count) return a.count > b.count;
    return a.atomicNumber < b.atomicNumber;
  });
  for (std::size_t i = 0; i < entries.size(); ++i) entries[i].colour = compositionColour(i);
  return entries;
}

std::string countText(std::size_t count) { return std::to_string(count); }

std::string signedCharge(int charge) {
  if (charge > 0) return "+" + std::to_string(charge);
  return std::to_string(charge);
}

void drawHeadline(const AppState& st, const DocumentSummary& summary, bool propertiesReady) {
  if (!widgets::beginCard("##property_headline", ImVec2(0.0f, 0.0f),
                          style::col::BgRaised)) {
    return;
  }
  widgets::cardHeader(icons::Icon::Molecule, "Molecular overview",
                      "The current structure at a glance", style::col::Accent);

  const float fs = ImGui::GetFontSize();
  const float available = std::max(ImGui::GetContentRegionAvail().x, fs);
  int columns = 2;
  if (available >= fs * 36.0f) {
    columns = 5;
  } else if (available >= fs * 23.0f) {
    columns = 3;
  }

  char mw[32];
  std::snprintf(mw, sizeof(mw), "%.2f", st.props.mw);
  const std::string atoms = countText(summary.totalAtoms());
  const std::string bonds = countText(summary.bonds);
  const std::string charge = signedCharge(summary.formalCharge);

  if (ImGui::BeginTable("##headline_metrics", columns,
                        ImGuiTableFlags_SizingStretchSame |
                            ImGuiTableFlags_NoSavedSettings)) {
    for (int metricIndex = 0; metricIndex < 5; ++metricIndex) {
      if (metricIndex % columns == 0) ImGui::TableNextRow();
      ImGui::TableSetColumnIndex(metricIndex % columns);
      switch (metricIndex) {
        case 0: {
          const bool mono = style::pushFont(style::fonts::mono());
          const std::string formula = propertiesReady
                                          ? ellipsize(st.props.formula,
                                                      ImGui::GetContentRegionAvail().x)
                                          : "--";
          style::popFont(mono);
          widgets::metric("FORMULA", formula.c_str(), nullptr, nullptr,
                          style::col::Accent);
          if (propertiesReady && formula != st.props.formula && ImGui::IsItemHovered())
            ImGui::SetTooltip("%s", st.props.formula.c_str());
          break;
        }
        case 1:
          widgets::metric("MOLAR MASS", propertiesReady ? mw : "--", "g/mol", nullptr,
                          style::col::Teal);
          break;
        case 2:
          widgets::metric("ATOMS", atoms.c_str(), nullptr, nullptr, style::col::Violet);
          break;
        case 3:
          widgets::metric("BONDS", bonds.c_str());
          break;
        default:
          widgets::metric("CHARGE", charge.c_str(), "e", nullptr,
                          summary.formalCharge == 0 ? style::col::Text : style::col::Accent);
          break;
      }
    }
    ImGui::EndTable();
  }
  widgets::endCard();
}

void drawCompositionGroup(const DocumentSummary& summary,
                          const std::vector<CompositionEntry>& composition) {
  char summaryText[64];
  std::snprintf(summaryText, sizeof(summaryText), "%zu elements / %zu atoms",
                composition.size(), summary.totalAtoms());
  if (!widgets::disclosure("##composition_disclosure", "Composition", summaryText, true,
                           icons::Icon::Atom, style::col::Teal)) {
    return;
  }

  ImGui::Indent(style::metrics().gap);
  widgets::cardHeader(icons::Icon::Atom, "Element composition",
                      "Includes implicit hydrogens", style::col::Teal);

  std::vector<charts::StackSegment> segments;
  segments.reserve(composition.size());
  for (const auto& entry : composition) {
    segments.push_back(
        charts::StackSegment{entry.symbol.c_str(), static_cast<double>(entry.count), entry.colour});
  }
  const std::string total = countText(summary.totalAtoms());
  const float chartHeight = std::max(ImGui::GetFontSize() * 11.0f,
                                     ImGui::GetTextLineHeightWithSpacing() * 8.0f);
  charts::donut("##element_composition", segments.data(), static_cast<int>(segments.size()),
                total.c_str(), "ATOMS", ImVec2(ImGui::GetContentRegionAvail().x, chartHeight));

  for (const auto& entry : composition) {
    const std::string key = entry.symbol + "  " + entry.name;
    const std::string value = std::to_string(entry.count) +
                              (entry.count == 1 ? " atom" : " atoms");
    widgets::keyValue(key.c_str(), value.c_str(), entry.colour);
  }
  ImGui::Unindent(style::metrics().gap);
}

void drawPhysicochemicalGroup(const AppState& st) {
  char summary[96];
  std::snprintf(summary, sizeof(summary), "%.2f g/mol / cLogP %.2f", st.props.mw,
                st.props.logP);
  if (!widgets::disclosure("##physchem_disclosure", "Physicochemical", summary, true,
                           icons::Icon::Balance, style::col::Accent)) {
    return;
  }

  ImGui::Indent(style::metrics().gap);
  widgets::cardHeader(icons::Icon::Balance, "Mass and partitioning",
                      "Calculated from the complete structure", style::col::Accent);
  char mw[32];
  char logp[32];
  std::snprintf(mw, sizeof(mw), "%.2f g/mol", st.props.mw);
  std::snprintf(logp, sizeof(logp), "%.2f", st.props.logP);
  widgets::keyValue("Molecular weight", mw, style::col::Teal);
  widgets::keyValue("Calculated logP", logp, style::col::Violet);
  ImGui::Unindent(style::metrics().gap);
}

void drawStructureGroup(const AppState& st, const DocumentSummary& summary,
                        bool propertiesReady) {
  char summaryText[64];
  if (propertiesReady) {
    std::snprintf(summaryText, sizeof(summaryText), "%zu bonds / %d rings", summary.bonds,
                  st.props.rings);
  } else {
    std::snprintf(summaryText, sizeof(summaryText), "%zu bonds / rings updating",
                  summary.bonds);
  }
  if (!widgets::disclosure("##structure_disclosure", "Structure", summaryText, false,
                           icons::Icon::Ruler, style::col::Violet)) {
    return;
  }

  ImGui::Indent(style::metrics().gap);
  widgets::cardHeader(icons::Icon::Ruler, "Topology and geometry",
                      "Counts from the current sketch", style::col::Violet);
  const std::string explicitAtoms = countText(summary.explicitAtoms);
  const std::string implicitHydrogens = std::to_string(summary.implicitHydrogens);
  const std::string bonds = countText(summary.bonds);
  const std::string rings = propertiesReady ? std::to_string(st.props.rings) : "--";
  const std::string charge = signedCharge(summary.formalCharge) + " e";
  widgets::keyValue("Explicit atoms", explicitAtoms.c_str());
  widgets::keyValue("Implicit hydrogens", implicitHydrogens.c_str());
  widgets::keyValue("Bonds", bonds.c_str());
  widgets::keyValue("Rings", rings.c_str());
  widgets::keyValue("Formal charge", charge.c_str(),
                    summary.formalCharge == 0 ? style::col::Text : style::col::Accent);
  ImGui::Unindent(style::metrics().gap);
}

void drawIdentityGroup(AppState& st, Clock::time_point now) {
  const bool nameDelayPassed = now - st.props.lastEdit >= std::chrono::milliseconds(1500);
  if (st.props.autoName && !st.props.smiles.empty() &&
      st.props.computedForRevision == st.docRevision &&
      st.props.nameRequestedForRevision != st.docRevision && nameDelayPassed) {
    requestAutomaticName(st);
  }
  const char* identitySummary = st.props.nameStatus == Status::Ok && !st.props.name.empty()
                                    ? st.props.name.c_str()
                                    : "Canonical SMILES and naming";
  if (!widgets::disclosure("##identity_disclosure", "Identity", identitySummary, true,
                           icons::Icon::Molecule, style::col::Teal)) {
    return;
  }

  ImGui::Indent(style::metrics().gap);
  widgets::cardHeader(icons::Icon::Molecule, "Chemical identity",
                      "Canonical representation and resolved name", style::col::Teal);

  std::vector<char> smiles(st.props.smiles.begin(), st.props.smiles.end());
  smiles.push_back('\0');
  const float actionWidth = ImGui::GetFontSize() * 5.5f;
  ImGui::SetNextItemWidth(std::max(ImGui::GetContentRegionAvail().x - actionWidth -
                                      ImGui::GetStyle().ItemSpacing.x,
                                  ImGui::GetFontSize()));
  const bool mono = style::pushFont(style::fonts::mono());
  ImGui::InputText("##canonical_smiles", smiles.data(), smiles.size(),
                   ImGuiInputTextFlags_ReadOnly | ImGuiInputTextFlags_AutoSelectAll);
  style::popFont(mono);
  if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", st.props.smiles.c_str());
  ImGui::SameLine();
  if (widgets::actionButton("##copy_smiles", icons::Icon::Copy, "Copy",
                            ImVec2(actionWidth, ImGui::GetFrameHeight()), false,
                            "Copy canonical SMILES")) {
    st.clipboardSmiles = st.props.smiles;
    ImGui::SetClipboardText(st.props.smiles.c_str());
    st.statusMessage = "Canonical SMILES copied";
  }

  widgets::toggle("##auto_name", "Resolve IUPAC name automatically", st.props.autoName,
                  "Looks up the current canonical structure after editing settles");

  if (st.props.nameStatus == Status::Loading) {
    widgets::notice(icons::Icon::Search, "Resolving the systematic name...",
                    style::col::Violet);
  } else if (st.props.nameStatus == Status::Ok) {
    widgets::keyValue("Resolved name", st.props.name.c_str(), style::col::Accent);
  } else if (st.props.nameStatus == Status::Error) {
    widgets::notice(icons::Icon::Warning, st.props.nameError.c_str(), style::col::Danger);
  }
  ImGui::Unindent(style::metrics().gap);
}

void drawBuildGroup(AppState& st, bool defaultOpen) {
  if (!widgets::disclosure("##build_name_disclosure", "Build from name",
                           "Chemical name to structure", defaultOpen, icons::Icon::Molecule,
                           style::col::Accent)) {
    return;
  }

  ImGui::Indent(style::metrics().gap);
  widgets::cardHeader(icons::Icon::Molecule, "Build from a chemical name",
                      "Accepts systematic or common names", style::col::Accent);
  BuildNameState& build = buildNameState();
  ImGui::SetNextItemWidth(-1.0f);
  const bool enter = ImGui::InputTextWithHint(
      "##name_to_structure", "e.g. acetylsalicylic acid", build.input.data(),
      build.input.size(), ImGuiInputTextFlags_EnterReturnsTrue);
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "Press Enter to append the resolved structure to the sketch");

  if (build.status == Status::Loading) ImGui::BeginDisabled();
  const bool buildClicked = widgets::actionButton(
      "##build_structure", icons::Icon::Molecule, "Build structure",
      ImVec2(ImGui::GetFontSize() * 10.0f, ImGui::GetFrameHeight()), true,
      "Resolve this name and append its structure");
  if (build.status == Status::Loading) ImGui::EndDisabled();
  if (enter || buildClicked) submitBuild(st);

  if (build.status == Status::Loading) {
    widgets::notice(icons::Icon::Search, "Resolving the chemical name...",
                    style::col::Violet);
  } else if (build.status == Status::Error && !build.error.empty()) {
    widgets::notice(icons::Icon::Warning, build.error.c_str(), style::col::Danger);
  }
  ImGui::Unindent(style::metrics().gap);
}

}  // namespace

void drawPropertiesPanel(AppState& st) {
  const Clock::time_point now = Clock::now();
  const bool propertyDelayPassed =
      now - st.props.lastEdit >= std::chrono::milliseconds(250);
  if (st.props.computedForRevision != st.docRevision && propertyDelayPassed) {
    recomputeProperties(st);
  }

  const DocumentSummary& summary = cachedDocumentSummary(st);
  if (summary.explicitAtoms == 0) {
    widgets::emptyState(
        icons::Icon::Molecule, "No structure selected",
        "Sketch a molecule, or build one from a chemical name, to inspect its properties.");
    ImGui::Spacing();
    drawBuildGroup(st, true);
    return;
  }

  const bool propertiesReady = st.props.computedForRevision == st.docRevision &&
                               st.props.chemError.empty() && !st.props.smiles.empty();
  drawHeadline(st, summary, propertiesReady);

  if (!st.props.chemError.empty()) {
    ImGui::Spacing();
    widgets::notice(icons::Icon::Warning, st.props.chemError.c_str(), style::col::Danger);
  } else if (st.props.computedForRevision != st.docRevision) {
    ImGui::Spacing();
    widgets::notice(icons::Icon::Timer, "Updating calculated properties...",
                    style::col::Violet);
  }

  const std::vector<CompositionEntry> composition = makeComposition(summary);
  ImGui::Spacing();
  drawCompositionGroup(summary, composition);
  if (propertiesReady) {
    ImGui::Spacing();
    drawPhysicochemicalGroup(st);
  }
  ImGui::Spacing();
  drawStructureGroup(st, summary, propertiesReady);
  if (propertiesReady) {
    ImGui::Spacing();
    drawIdentityGroup(st, now);
  }
  ImGui::Spacing();
  drawBuildGroup(st, false);
}

}  // namespace chemcad::ui

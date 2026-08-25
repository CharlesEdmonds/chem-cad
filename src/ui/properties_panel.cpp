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
#include "ui/charts3d.hpp"
#include "ui/element_data.hpp"
#include "ui/icons.hpp"
#include "ui/layout.hpp"
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
  static constexpr ImVec4 colours[] = {
      style::col::DataBright, style::col::Data, style::col::DataDim,
      style::col::Teal};
  return colours[index % (sizeof(colours) / sizeof(colours[0]))];
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
void advanceVerticalGap(float gap) {
  layout::nextRow(ImGui::GetCursorPosY() + gap -
                       ImGui::GetStyle().ItemSpacing.y);
}


void drawHeadline(const AppState& st, const DocumentSummary& summary,
                  bool propertiesReady, float height,
                  const layout::Frame& page) {
  if (!widgets::beginCard("##property_headline", ImVec2(0.0f, height),
                          style::col::BgRaised)) {
    return;
  }
  widgets::cardHeader(icons::Icon::Molecule, "Molecular overview",
                      "Current structure at a glance", style::col::Data);

  char mw[32];
  std::snprintf(mw, sizeof(mw), "%.2f", st.props.mw);
  const std::string atoms = countText(summary.totalAtoms());
  const std::string bonds = countText(summary.bonds);
  const std::string charge = signedCharge(summary.formalCharge);
  const float gap = page.gap * 0.5f;
  const float width =
      std::max((ImGui::GetContentRegionAvail().x - gap * 4.0f) / 5.0f, 0.0f);
  const float metricHeight = std::max(ImGui::GetContentRegionAvail().y, 0.0f);
  const std::string formula =
      propertiesReady ? ellipsize(st.props.formula, width) : "--";

  static constexpr const char* captions[] = {
      "Formula", "Molar mass", "Atoms", "Bonds", "Charge"};
  const char* values[] = {
      formula.c_str(), propertiesReady ? mw : "--", atoms.c_str(), bonds.c_str(),
      charge.c_str()};
  const char* units[] = {nullptr, "g/mol", nullptr, nullptr, "e"};
  const ImVec4 accents[] = {
      style::col::DataBright, style::col::Data, style::col::DataDim,
      style::col::DataDim,
      summary.formalCharge == 0 ? style::col::DataDim : style::col::Danger};

  for (int index = 0; index < 5; ++index) {
    ImGui::PushID(index);
    if (ImGui::BeginChild("##headline_metric", ImVec2(width, metricHeight),
                          ImGuiChildFlags_None,
                          ImGuiWindowFlags_NoScrollbar |
                              ImGuiWindowFlags_NoScrollWithMouse)) {
      widgets::metric(captions[index], values[index], units[index], nullptr,
                      accents[index]);
      if (index == 0 && propertiesReady && formula != st.props.formula &&
          ImGui::IsItemHovered()) {
        ImGui::SetTooltip("%s", st.props.formula.c_str());
      }
    }
    ImGui::EndChild();
    ImGui::PopID();
    if (index != 4) {
      ImGui::SameLine(0.0f, gap);
    }
  }
  widgets::endCard();
}

void drawCompositionTable(const std::vector<CompositionEntry>& composition,
                          ImVec2 size, const layout::Frame& frame) {
  static constexpr widgets::Column columns[] = {
      {"Element", false, true, nullptr, 7.0f},
      {"Count", true, false, "atoms", 4.0f}};
  if (!widgets::beginDataTable("##composition_table", columns, 2, size)) return;

  const int availableRows =
      std::max(static_cast<int>(size.y / std::max(frame.row, frame.em)) - 1, 1);
  const int compositionRows = static_cast<int>(composition.size());
  const int visibleRows =
      compositionRows <= availableRows
          ? compositionRows
          : std::max(availableRows - 1, 0);
  for (int index = 0; index < visibleRows; ++index) {
    const CompositionEntry& entry = composition[static_cast<std::size_t>(index)];
    const std::string label = entry.symbol + "  " + entry.name;
    widgets::dataRow(entry.colour);
    widgets::dataCell(label.c_str());
    widgets::dataCellf("%d", entry.count);
  }
  if (visibleRows < static_cast<int>(composition.size())) {
    int omittedAtoms = 0;
    for (std::size_t index = static_cast<std::size_t>(visibleRows);
         index < composition.size(); ++index) {
      omittedAtoms += composition[index].count;
    }
    const std::string label =
        "+" + std::to_string(composition.size() -
                             static_cast<std::size_t>(visibleRows)) +
        " more";
    widgets::dataRow(style::col::DataDim);
    widgets::dataCell(label.c_str());
    widgets::dataCellf("%d", omittedAtoms);
  }
  widgets::endDataTable();
}

void drawCompositionGroup(const DocumentSummary& summary,
                          const std::vector<CompositionEntry>& composition) {
  char summaryText[64];
  std::snprintf(summaryText, sizeof(summaryText), "%zu elements / %zu atoms",
                composition.size(), summary.totalAtoms());
  if (!widgets::disclosure("##composition_disclosure", "Element composition",
                           summaryText, true, icons::Icon::Atom,
                           style::col::Accent)) {
    return;
  }

  ImGui::Indent(style::metrics().gap);
  const layout::Frame body = layout::measure();
  std::vector<charts::StackSegment> segments;
  segments.reserve(composition.size());
  for (const auto& entry : composition) {
    segments.push_back(charts::StackSegment{
        entry.symbol.c_str(), static_cast<double>(entry.count), entry.colour});
  }
  const std::string total = countText(summary.totalAtoms());

  if (body.wide) {
    const float chartWidth = layout::columnWidth(body, 2);
    if (ImGui::BeginChild("##composition_chart",
                          ImVec2(chartWidth, body.size.y), ImGuiChildFlags_None,
                          ImGuiWindowFlags_NoScrollbar |
                              ImGuiWindowFlags_NoScrollWithMouse)) {
      charts::donut("##element_composition", segments.data(),
                    static_cast<int>(segments.size()), total.c_str(), "ATOMS",
                    ImGui::GetContentRegionAvail());
    }
    ImGui::EndChild();
    ImGui::SameLine(0.0f, body.gap);
    if (ImGui::BeginChild("##composition_counts",
                          ImVec2(layout::columnWidth(body, 2), body.size.y),
                          ImGuiChildFlags_None,
                          ImGuiWindowFlags_NoScrollbar |
                              ImGuiWindowFlags_NoScrollWithMouse)) {
      drawCompositionTable(composition, ImGui::GetContentRegionAvail(), body);
    }
    ImGui::EndChild();
  } else {
    const float weights[] = {1.1f, 0.9f};
    const float minimums[] = {0.0f, 0.0f};
    float rows[2]{};
    layout::distribute(body.size.y, weights, minimums, 2, body.gap, rows);
    if (ImGui::BeginChild("##composition_chart",
                          ImVec2(body.size.x, rows[0]), ImGuiChildFlags_None,
                          ImGuiWindowFlags_NoScrollbar |
                              ImGuiWindowFlags_NoScrollWithMouse)) {
      charts::donut("##element_composition", segments.data(),
                    static_cast<int>(segments.size()), total.c_str(), "ATOMS",
                    ImGui::GetContentRegionAvail());
    }
    ImGui::EndChild();
    advanceVerticalGap(body.gap);
    if (ImGui::BeginChild("##composition_counts",
                          ImVec2(body.size.x, rows[1]), ImGuiChildFlags_None,
                          ImGuiWindowFlags_NoScrollbar |
                              ImGuiWindowFlags_NoScrollWithMouse)) {
      drawCompositionTable(composition, ImGui::GetContentRegionAvail(), body);
    }
    ImGui::EndChild();
  }
  ImGui::Unindent(style::metrics().gap);
}

void drawPhysicochemicalGroup(const AppState& st, bool propertiesReady) {
  char summary[96];
  std::snprintf(summary, sizeof(summary), "%.2f g/mol / cLogP %.2f",
                st.props.mw, st.props.logP);
  if (!widgets::disclosure("##physchem_disclosure", "Physicochemical",
                           propertiesReady ? summary : "Updating", true,
                           icons::Icon::Balance, style::col::Accent)) {
    return;
  }
  ImGui::Indent(style::metrics().gap);
  if (widgets::onlyWhen(propertiesReady,
                        "Calculated properties are still updating")) {
    char mw[32];
    char logp[32];
    std::snprintf(mw, sizeof(mw), "%.2f g/mol", st.props.mw);
    std::snprintf(logp, sizeof(logp), "%.2f", st.props.logP);
    widgets::keyValue("Molecular weight", mw, style::col::DataBright);
    widgets::keyValue("Calculated logP", logp, style::col::Data);
  }
  ImGui::Unindent(style::metrics().gap);
}

void drawConformerCloud(const AppState& st, ImVec2 size,
                        const layout::Frame& frame) {
  const auto& atoms = st.viewer3d.model.atoms;
  std::vector<std::string> labels;
  labels.reserve(atoms.size());
  for (const auto& atom : atoms) {
    if (const ElementData* element = findElement(atom.atomicNumber)) {
      labels.push_back(element->symbol);
    } else {
      const char* symbol = chem::symbolFor(atom.atomicNumber);
      labels.emplace_back(symbol && *symbol ? symbol : "?");
    }
  }

  std::vector<charts3d::CloudPoint> points;
  points.reserve(atoms.size());
  for (std::size_t index = 0; index < atoms.size(); ++index) {
    const auto& atom = atoms[index];
    points.push_back(charts3d::CloudPoint{
        atom.x, atom.y, atom.z, labels[index].c_str(), style::col::Data,
        1.0f, false});
  }
  static charts3d::Orbit orbit;
  charts3d::CloudStyle cloudStyle;
  cloudStyle.xLabel = "x (A)";
  cloudStyle.yLabel = "y (A)";
  cloudStyle.zLabel = "z (A)";
  cloudStyle.showAxes = true;
  cloudStyle.showLabels = frame.density != layout::Density::Compact;
  cloudStyle.hasSphere = false;
  charts3d::cloud("##conformer_cloud", points.data(),
                  static_cast<int>(points.size()), size, orbit, cloudStyle);
}

void drawStructureTable(const AppState& st, const DocumentSummary& summary,
                        bool propertiesReady, ImVec2 size) {
  static constexpr widgets::Column columns[] = {
      {"Measure", false, true, nullptr, 10.0f},
      {"Count", true, false, nullptr, 4.0f}};
  if (!widgets::beginDataTable("##structure_counts", columns, 2, size)) return;
  const std::string explicitAtoms = countText(summary.explicitAtoms);
  const std::string implicitHydrogens = std::to_string(summary.implicitHydrogens);
  const std::string bonds = countText(summary.bonds);
  const std::string rings =
      propertiesReady ? std::to_string(st.props.rings) : "--";
  const std::string charge = signedCharge(summary.formalCharge);
  const char* labels[] = {
      "Explicit atoms", "Implicit hydrogens", "Bonds", "Rings", "Formal charge"};
  const char* values[] = {
      explicitAtoms.c_str(), implicitHydrogens.c_str(), bonds.c_str(),
      rings.c_str(), charge.c_str()};
  for (int index = 0; index < 5; ++index) {
    widgets::dataRow(index == 4 && summary.formalCharge != 0
                         ? style::col::Danger
                         : style::col::DataDim);
    widgets::dataCell(labels[index]);
    widgets::dataCell(values[index]);
  }
  widgets::endDataTable();
}

void drawStructureGroup(const AppState& st, const DocumentSummary& summary,
                        bool propertiesReady) {
  char summaryText[64];
  if (propertiesReady) {
    std::snprintf(summaryText, sizeof(summaryText), "%zu bonds / %d rings",
                  summary.bonds, st.props.rings);
  } else {
    std::snprintf(summaryText, sizeof(summaryText),
                  "%zu bonds / rings updating", summary.bonds);
  }
  if (!widgets::disclosure("##structure_disclosure", "Structure", summaryText,
                           true, icons::Icon::Ruler, style::col::Accent)) {
    return;
  }

  ImGui::Indent(style::metrics().gap);
  const layout::Frame body = layout::measure();
  const bool hasConformer =
      st.viewer3d.hasModel &&
      st.viewer3d.sourceRevision == st.docRevision &&
      !st.viewer3d.model.atoms.empty();
  if (hasConformer) {
    const float weights[] = {0.42f, 0.58f};
    const float minimums[] = {body.row * 6.0f, 0.0f};
    float rows[2]{};
    layout::distribute(body.size.y, weights, minimums, 2, body.gap, rows);
    if (ImGui::BeginChild("##structure_table_region",
                          ImVec2(body.size.x, rows[0]), ImGuiChildFlags_None,
                          ImGuiWindowFlags_NoScrollbar |
                              ImGuiWindowFlags_NoScrollWithMouse)) {
      drawStructureTable(st, summary, propertiesReady,
                         ImGui::GetContentRegionAvail());
    }
    ImGui::EndChild();
    advanceVerticalGap(body.gap);
    if (ImGui::BeginChild("##conformer_region", ImVec2(body.size.x, rows[1]),
                          ImGuiChildFlags_None,
                          ImGuiWindowFlags_NoScrollbar |
                              ImGuiWindowFlags_NoScrollWithMouse)) {
      const ImVec2 cloudSize = ImGui::GetContentRegionAvail();
      widgets::hudFrame(ImGui::GetCursorScreenPos(),
                        ImVec2(ImGui::GetCursorScreenPos().x + cloudSize.x,
                               ImGui::GetCursorScreenPos().y + cloudSize.y),
                        style::col::Data);
      drawConformerCloud(st, cloudSize, body);
    }
    ImGui::EndChild();
  } else {
    drawStructureTable(st, summary, propertiesReady,
                       ImGui::GetContentRegionAvail());
  }
  ImGui::Unindent(style::metrics().gap);
}

void updateAutomaticName(AppState& st, Clock::time_point now) {
  const bool nameDelayPassed =
      now - st.props.lastEdit >= std::chrono::milliseconds(1500);
  if (st.props.autoName && !st.props.smiles.empty() &&
      st.props.computedForRevision == st.docRevision &&
      st.props.nameRequestedForRevision != st.docRevision &&
      nameDelayPassed) {
    requestAutomaticName(st);
  }
}

void drawIdentityGroup(AppState& st, bool propertiesReady) {
  const char* identitySummary =
      st.props.nameStatus == Status::Ok && !st.props.name.empty()
          ? st.props.name.c_str()
          : "Canonical SMILES and naming";
  if (!widgets::disclosure("##identity_disclosure", "Identity",
                           propertiesReady ? identitySummary : "Updating",
                           true, icons::Icon::Molecule,
                           style::col::Accent)) {
    return;
  }

  ImGui::Indent(style::metrics().gap);
  if (!widgets::onlyWhen(propertiesReady,
                         "Canonical identity is still updating")) {
    ImGui::Unindent(style::metrics().gap);
    return;
  }

  std::vector<char> smiles(st.props.smiles.begin(), st.props.smiles.end());
  smiles.push_back('\0');
  const float actionWidth =
      widgets::actionButtonWidth(icons::Icon::Copy, "Copy");
  ImGui::SetNextItemWidth(std::max(
      ImGui::GetContentRegionAvail().x - actionWidth -
          ImGui::GetStyle().ItemSpacing.x,
      ImGui::GetFontSize()));
  const bool mono = style::pushFont(style::fonts::mono());
  ImGui::InputText("##canonical_smiles", smiles.data(), smiles.size(),
                   ImGuiInputTextFlags_ReadOnly |
                       ImGuiInputTextFlags_AutoSelectAll);
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

  widgets::toggle("##auto_name", "Resolve IUPAC name automatically",
                  st.props.autoName,
                  "Looks up the current canonical structure after editing settles");
  if (st.props.nameStatus == Status::Loading) {
    widgets::statusDot("Resolving systematic name", true,
                       style::col::Violet);
  } else if (st.props.nameStatus == Status::Ok) {
    widgets::keyValue("Resolved name", st.props.name.c_str(),
                      style::col::Violet);
  } else if (st.props.nameStatus == Status::Error) {
    widgets::notice(icons::Icon::Warning, st.props.nameError.c_str(),
                    style::col::Danger);
  }
  ImGui::Unindent(style::metrics().gap);
}

void drawBuildGroup(AppState& st, bool defaultOpen) {
  if (!widgets::disclosure("##build_name_disclosure", "Build from name",
                           "Chemical name to structure", defaultOpen,
                           icons::Icon::Molecule, style::col::Accent)) {
    return;
  }

  ImGui::Indent(style::metrics().gap);
  BuildNameState& build = buildNameState();
  ImGui::SetNextItemWidth(-1.0f);
  const bool enter = ImGui::InputTextWithHint(
      "##name_to_structure", "e.g. acetylsalicylic acid", build.input.data(),
      build.input.size(), ImGuiInputTextFlags_EnterReturnsTrue);
  if (ImGui::IsItemHovered()) {
    ImGui::SetTooltip(
        "%s", "Press Enter to append the resolved structure to the sketch");
  }

  if (build.status == Status::Loading) ImGui::BeginDisabled();
  const float actionWidth =
      widgets::actionButtonWidth(icons::Icon::Molecule, "Build structure", true);
  const bool buildClicked = widgets::actionButton(
      "##build_structure", icons::Icon::Molecule, "Build structure",
      ImVec2(actionWidth, ImGui::GetFrameHeight()), true,
      "Resolve this name and append its structure");
  if (build.status == Status::Loading) ImGui::EndDisabled();
  if (enter || buildClicked) submitBuild(st);

  if (build.status == Status::Loading) {
    widgets::statusDot("Resolving chemical name", true, style::col::Violet);
  } else if (build.status == Status::Error && !build.error.empty()) {
    widgets::notice(icons::Icon::Warning, build.error.c_str(),
                    style::col::Danger);
  }
  ImGui::Unindent(style::metrics().gap);
}

}  // namespace

void drawPropertiesPanel(AppState& st) {
  static int secondaryTab = 0;
  const Clock::time_point now = Clock::now();
  const bool propertyDelayPassed =
      now - st.props.lastEdit >= std::chrono::milliseconds(250);
  if (st.props.computedForRevision != st.docRevision &&
      propertyDelayPassed) {
    recomputeProperties(st);
  }

  const layout::Frame page = layout::measure();
  const float budget = std::min(page.size.y, layout::pageHeight());
  const DocumentSummary& summary = cachedDocumentSummary(st);
  if (summary.explicitAtoms == 0) {
    // The build group is a fixed stack -- disclosure, input, button -- so it
    // takes exactly what it needs and the empty state absorbs the rest. Sharing
    // the band by weight instead clipped the input in a short dock node.
    const float buildBand = page.control * 3.0f + page.gap * 2.0f;
    const float weights[] = {1.0f, 0.0f};
    const float minimums[] = {page.row * 2.0f, buildBand};
    float rows[2]{};
    layout::distribute(budget, weights, minimums, 2, page.gap, rows);
    if (ImGui::BeginChild("##properties_empty",
                          ImVec2(page.size.x, rows[0]), ImGuiChildFlags_None,
                          ImGuiWindowFlags_NoScrollbar |
                              ImGuiWindowFlags_NoScrollWithMouse)) {
      widgets::emptyState(
          icons::Icon::Molecule, "No structure selected",
          "Sketch a molecule, or build one from a chemical name, to inspect its properties.");
    }
    ImGui::EndChild();
    advanceVerticalGap(page.gap);
    if (ImGui::BeginChild("##properties_empty_build",
                          ImVec2(page.size.x, rows[1]), ImGuiChildFlags_None,
                          ImGuiWindowFlags_NoScrollbar |
                              ImGuiWindowFlags_NoScrollWithMouse)) {
      drawBuildGroup(st, true);
    }
    ImGui::EndChild();
    return;
  }

  const bool propertiesReady =
      st.props.computedForRevision == st.docRevision &&
      st.props.chemError.empty() && !st.props.smiles.empty();
  if (propertiesReady) updateAutomaticName(st, now);

  const bool hasNotice =
      !st.props.chemError.empty() ||
      st.props.computedForRevision != st.docRevision;
  float rows[4]{};
  if (hasNotice) {
    const float weights[] = {0.0f, 0.0f, 0.0f, 1.0f};
    const float minimums[] = {
        page.row * 4.4f, page.control, page.control, 0.0f};
    layout::distribute(budget, weights, minimums, 4, page.gap, rows);
  } else {
    const float weights[] = {0.0f, 0.0f, 1.0f};
    const float minimums[] = {page.row * 4.4f, page.control, 0.0f};
    float compactRows[3]{};
    layout::distribute(budget, weights, minimums, 3, page.gap,
                       compactRows);
    rows[0] = compactRows[0];
    rows[2] = compactRows[1];
    rows[3] = compactRows[2];
  }

  drawHeadline(st, summary, propertiesReady, rows[0], page);
  advanceVerticalGap(page.gap);

  if (rows[1] > 0.0f) {
    if (ImGui::BeginChild("##properties_notice",
                          ImVec2(page.size.x, rows[1]), ImGuiChildFlags_None,
                          ImGuiWindowFlags_NoScrollbar |
                              ImGuiWindowFlags_NoScrollWithMouse)) {
      if (!st.props.chemError.empty()) {
        widgets::notice(icons::Icon::Warning, st.props.chemError.c_str(),
                        style::col::Danger);
      } else {
        widgets::statusDot("Updating calculated properties", true,
                           style::col::Data);
      }
    }
    ImGui::EndChild();
    advanceVerticalGap(page.gap);
  }

  if (ImGui::BeginChild("##properties_tabs",
                        ImVec2(page.size.x, rows[2]), ImGuiChildFlags_None,
                        ImGuiWindowFlags_NoScrollbar |
                            ImGuiWindowFlags_NoScrollWithMouse)) {
    static constexpr const char* labels[] = {
        "Elements", "Physical", "Structure", "Identity", "Build"};
    static constexpr icons::Icon glyphs[] = {
        icons::Icon::Atom, icons::Icon::Balance, icons::Icon::Ruler,
        icons::Icon::Molecule, icons::Icon::Sparkle};
    widgets::subTabs("##property_secondary_tabs", labels, glyphs, 5,
                     secondaryTab);
  }
  ImGui::EndChild();
  advanceVerticalGap(page.gap);

  if (ImGui::BeginChild("##properties_secondary",
                        ImVec2(page.size.x, rows[3]), ImGuiChildFlags_None,
                        ImGuiWindowFlags_NoScrollbar |
                            ImGuiWindowFlags_NoScrollWithMouse)) {
    const std::vector<CompositionEntry> composition = makeComposition(summary);
    switch (secondaryTab) {
      case 0:
        drawCompositionGroup(summary, composition);
        break;
      case 1:
        drawPhysicochemicalGroup(st, propertiesReady);
        break;
      case 2:
        drawStructureGroup(st, summary, propertiesReady);
        break;
      case 3:
        drawIdentityGroup(st, propertiesReady);
        break;
      default:
        drawBuildGroup(st, true);
        break;
    }
  }
  ImGui::EndChild();
}

}  // namespace chemcad::ui

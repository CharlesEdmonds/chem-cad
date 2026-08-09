#include <cmath>
#include <cstdio>
#include <cfloat>
#include <filesystem>
#include <string>

#include "imgui.h"
#include "imgui_internal.h"

#include "ui/charts.hpp"
#include "ui/element_data.hpp"
#include "ui/icons.hpp"
#include "ui/layout.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"
#include "ui/widgets.hpp"

namespace chemcad::ui {
namespace {

const char* bondOrderName(core::BondOrder order) {
  switch (order) {
    case core::BondOrder::Single:
      return "single";
    case core::BondOrder::Double:
      return "double";
    case core::BondOrder::Triple:
      return "triple";
    case core::BondOrder::Aromatic:
      return "aromatic";
  }
  return "unknown";
}

std::string hoverDescription(AppState& st) {
  if (st.hoverAtom.valid()) {
    if (const core::Atom* atom = st.atomAt(st.hoverAtom)) {
      const ElementData* element = findElement(atom->atomicNumber);
      std::string text = "Atom: ";
      text += element ? element->symbol : ("Z=" + std::to_string(atom->atomicNumber));
      if (atom->charge != 0) {
        text += atom->charge > 0 ? " +" : " ";
        text += std::to_string(atom->charge);
      }
      return text;
    }
  }
  if (st.hoverBond.valid()) {
    if (const core::Bond* bond = st.bondAt(st.hoverBond)) {
      return std::string("Bond: ") + bondOrderName(bond->order);
    }
  }
  return "No item hovered";
}

std::string toolHint(const AppState& st) {
  switch (st.tool) {
    case Tool::Bond:
      return "Bond: click an atom to attach CH3; drag to place or connect carbon (M preset)";
    case Tool::Chain:
      return "Chain: drag from an atom or empty space to draw a carbon chain";
    case Tool::RingTemplate:
      return "Ring: click the canvas or an atom to stamp the selected ring";
    case Tool::Atom:
      return "Atom: click to place or retype the selected element";
    case Tool::Eraser:
      return "Eraser: click an atom or bond to remove it";
    case Tool::ChargePlus:
      return "Positive charge: click an atom to increase formal charge";
    case Tool::ChargeMinus:
      return "Negative charge: click an atom to decrease formal charge";
    case Tool::Select:
      return "Select: click or drag around atoms and bonds";
  }
  return "Choose a tool, then draw in the Sketch canvas";
}

const char* toolName(Tool tool) {
  switch (tool) {
    case Tool::Select:
      return "Select";
    case Tool::Eraser:
      return "Eraser";
    case Tool::Bond:
      return "Bond";
    case Tool::Chain:
      return "Chain";
    case Tool::RingTemplate:
      return "Ring";
    case Tool::Atom:
      return "Atom";
    case Tool::ChargePlus:
      return "Charge +";
    case Tool::ChargeMinus:
      return "Charge -";
  }
  return "Tool";
}

icons::Icon ringIcon(RingKind ring) {
  switch (ring) {
    case RingKind::Cyclopropane:
      return icons::Icon::RingCyclopropane;
    case RingKind::Cyclobutane:
      return icons::Icon::RingCyclobutane;
    case RingKind::Cyclopentane:
      return icons::Icon::RingCyclopentane;
    case RingKind::Cyclohexane:
      return icons::Icon::RingCyclohexane;
    case RingKind::Cycloheptane:
      return icons::Icon::RingCycloheptane;
    case RingKind::Cyclooctane:
      return icons::Icon::RingCyclooctane;
    case RingKind::Benzene:
      return icons::Icon::RingBenzene;
    case RingKind::Cyclopentadiene:
      return icons::Icon::RingCyclopentadiene;
    case RingKind::Naphthalene:
      return icons::Icon::RingNaphthalene;
  }
  return icons::Icon::RingBenzene;
}

icons::Icon toolIcon(const AppState& st) {
  switch (st.tool) {
    case Tool::Select:
      return icons::Icon::Select;
    case Tool::Eraser:
      return icons::Icon::Eraser;
    case Tool::Chain:
      return icons::Icon::Chain;
    case Tool::RingTemplate:
      return ringIcon(st.currentRing);
    case Tool::Atom:
      return icons::Icon::Atom;
    case Tool::ChargePlus:
      return icons::Icon::ChargePlus;
    case Tool::ChargeMinus:
      return icons::Icon::ChargeMinus;
    case Tool::Bond:
      switch (st.currentStereo) {
        case core::BondStereo::Wedge:
          return icons::Icon::StereoWedge;
        case core::BondStereo::Hash:
          return icons::Icon::StereoHash;
        case core::BondStereo::Wavy:
          return icons::Icon::StereoWavy;
        case core::BondStereo::None:
          break;
      }
      switch (st.currentOrder) {
        case core::BondOrder::Single:
          return icons::Icon::BondSingle;
        case core::BondOrder::Double:
          return icons::Icon::BondDouble;
        case core::BondOrder::Triple:
          return icons::Icon::BondTriple;
        case core::BondOrder::Aromatic:
          return icons::Icon::BondAromatic;
      }
  }
  return icons::Icon::Bond;
}

std::string documentName(const AppState& st) {
  if (st.projectPath.empty()) return "Untitled";
  const std::filesystem::path path(st.projectPath);
  const std::string name = path.filename().string();
  return name.empty() ? "Untitled" : name;
}

void glyphText(icons::Icon icon, const char* text, ImVec4 colour = style::col::Data) {
  const style::Metrics& m = style::metrics();
  const float iconSize = m.iconSize * 0.66f;
  const ImVec2 measured = ImGui::CalcTextSize(text);
  const float width = iconSize + m.gap * 0.55f + measured.x;
  const float height = ImGui::GetTextLineHeight();
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::Dummy(ImVec2(width, height));
  ImDrawList* dl = ImGui::GetWindowDrawList();
  icons::draw(dl, icon, ImVec2(min.x + iconSize * 0.5f, min.y + height * 0.5f),
              iconSize, style::u32(colour));
  dl->AddText(ImVec2(min.x + iconSize + m.gap * 0.55f, min.y), style::u32(colour), text);
}

void compactMetric(icons::Icon icon, const char* caption, const char* value) {
  const style::Metrics& m = style::metrics();
  const float iconSize = m.iconSize * 0.58f;
  const ImVec2 valueSize = ImGui::CalcTextSize(value);
  const float smallSize = ImGui::GetFontSize() * 0.68f;
  const ImVec2 captionSize =
      style::fonts::semibold()->CalcTextSizeA(smallSize, FLT_MAX, 0.0f, caption);
  const float width =
      iconSize + m.gap * 0.45f + valueSize.x + m.gap * 0.35f + captionSize.x;
  const float height = ImGui::GetTextLineHeight();
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::Dummy(ImVec2(width, height));

  ImDrawList* dl = ImGui::GetWindowDrawList();
  icons::draw(dl, icon, ImVec2(min.x + iconSize * 0.5f, min.y + height * 0.5f),
              iconSize, style::u32(style::col::DataDim));
  const float valueX = min.x + iconSize + m.gap * 0.45f;
  const bool mono = style::pushFont(style::fonts::mono());
  dl->AddText(ImVec2(valueX, min.y), style::u32(style::col::DataBright), value);
  style::popFont(mono);
  dl->AddText(style::fonts::semibold(), smallSize,
              ImVec2(valueX + valueSize.x + m.gap * 0.35f,
                     min.y + (height - smallSize) * 0.5f),
              style::u32(style::col::DataDim), caption);
}

void drawBadgeCluster(const AppState& st) {
  const style::Metrics& m = style::metrics();
  const std::string formula = st.props.formula.empty() ? "--" : st.props.formula;
  char zoom[16];
  std::snprintf(zoom, sizeof(zoom), "%.0f%%", std::round(st.cam.zoom * 100.0f));

  const bool mono = style::pushFont(style::fonts::mono());
  widgets::badge(formula.c_str(), style::col::DataDim);
  style::popFont(mono);
  ImGui::SameLine(0.0f, m.gap * 0.5f);
  widgets::badge(zoom, style::col::DataDim);
  if (st.tasks.busy()) {
    ImGui::SameLine(0.0f, m.gap * 0.5f);
    widgets::badge("Working", style::col::Data);
  }
}

}  // namespace

void drawStatusBar(AppState& st) {
  static charts::Trace frameTimes(120);
  frameTimes.push(static_cast<double>(ImGui::GetIO().DeltaTime) * 1000.0);

  size_t atomCount = 0;
  size_t bondCount = 0;
  for (const core::Molecule& molecule : st.doc.molecules) {
    atomCount += molecule.atomCount();
    bondCount += molecule.bondCount();
  }
  const size_t selectionCount = st.sel.atoms.size() + st.sel.bonds.size();

  char atoms[16];
  char bonds[16];
  char selection[16];
  std::snprintf(atoms, sizeof(atoms), "%zu", atomCount);
  std::snprintf(bonds, sizeof(bonds), "%zu", bondCount);
  std::snprintf(selection, sizeof(selection), "%zu", selectionCount);
  char frameTime[24];
  std::snprintf(frameTime, sizeof(frameTime), "%.1f ms", frameTimes.latest());

  ImGuiViewport* viewport = ImGui::GetMainViewport();
  const style::Metrics& m = style::metrics();
  const float height = ImGui::GetTextLineHeightWithSpacing() * 2.0f + m.gap * 1.15f;
  constexpr ImGuiWindowFlags flags =
      ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoDocking |
      ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoNavFocus |
      ImGuiWindowFlags_NoBringToFrontOnFocus;

  if (ImGui::BeginViewportSideBar("##status", viewport, ImGuiDir_Down, height, flags)) {
    const ImVec2 wmin = ImGui::GetWindowPos();
    ImGui::GetWindowDrawList()->AddLine(
        wmin, ImVec2(wmin.x + ImGui::GetWindowWidth(), wmin.y),
        style::u32(style::col::GridLine), m.hairline);

    const layout::Frame frame = layout::measure();
    const float columnWeights[]{0.28f, 0.34f, 0.38f};
    float columnWidths[3]{};
    layout::distribute(frame.size.x, columnWeights, nullptr, 3, frame.gap,
                       columnWidths);
    constexpr ImGuiWindowFlags childFlags =
        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse;

    ImGui::BeginChild("##status_identity", ImVec2(columnWidths[0], frame.size.y),
                      false, childFlags);
    const std::string identity =
        std::string(toolName(st.tool)) + "  ·  " + documentName(st) + (st.dirty ? " *" : "");
    glyphText(toolIcon(st), identity.c_str());
    const std::string hover = hoverDescription(st);
    ImGui::TextColored(style::col::DataDim, "%s", hover.c_str());
    ImGui::EndChild();

    ImGui::SameLine(0.0f, frame.gap);
    ImGui::BeginChild("##status_message", ImVec2(columnWidths[1], frame.size.y),
                      false, childFlags);
    const std::string idleMessage = toolHint(st);
    const char* message =
        st.statusMessage.empty() ? idleMessage.c_str() : st.statusMessage.c_str();
    ImGui::TextColored(style::col::Data, "%s", message);
    if (ImGui::IsItemHovered() &&
        ImGui::CalcTextSize(message).x > ImGui::GetContentRegionAvail().x) {
      ImGui::SetTooltip("%s", message);
    }
    widgets::statusDot(st.solubility.funnelRunning ? "Fluid solver running"
                                                   : "Fluid solver paused",
                       st.solubility.funnelRunning, style::col::Data);
    ImGui::EndChild();

    ImGui::SameLine(0.0f, frame.gap);
    ImGui::BeginChild("##status_metrics", ImVec2(columnWidths[2], frame.size.y),
                      false, childFlags);
    compactMetric(icons::Icon::Atom, "Atoms", atoms);
    ImGui::SameLine(0.0f, frame.gap);
    compactMetric(icons::Icon::Bond, "Bonds", bonds);
    ImGui::SameLine(0.0f, frame.gap);
    compactMetric(icons::Icon::Select, "Selected", selection);

    compactMetric(icons::Icon::ChartLine, "Frame time", frameTime);
    ImGui::SameLine(0.0f, frame.gap);
    charts::SparklineStyle sparkStyle;
    sparkStyle.accent = style::col::Data;
    sparkStyle.fill = true;
    sparkStyle.showLatest = true;
    sparkStyle.autoFloor = false;
    sparkStyle.floorValue = 0.0;
    sparkStyle.ceilingValue = 0.0;
    charts::sparkline("##frame_time", frameTimes,
                      ImVec2(frame.em * 5.2f, ImGui::GetTextLineHeight()),
                      sparkStyle);
    if (ImGui::IsItemHovered()) {
      ImGui::SetTooltip("Frame time %.2f ms", frameTimes.latest());
    }
    if (frame.density != layout::Density::Compact) {
      ImGui::SameLine(0.0f, frame.gap);
      drawBadgeCluster(st);
    }
    ImGui::EndChild();
  }
  ImGui::End();
}

}  // namespace chemcad::ui

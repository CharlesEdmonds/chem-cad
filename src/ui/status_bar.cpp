#include <cmath>
#include <cstdio>
#include <string>

#include "imgui.h"
#include "imgui_internal.h"

#include "ui/element_data.hpp"
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

// Thin vertical hairline between segments.
void divider() {
  ImGui::SameLine(0.0f, 10.0f);
  const ImVec2 pos = ImGui::GetCursorScreenPos();
  const float h = ImGui::GetFontSize();
  ImGui::Dummy(ImVec2(1.0f, h));
  ImGui::GetWindowDrawList()->AddLine(ImVec2(pos.x + 0.5f, pos.y + 2.0f),
                                      ImVec2(pos.x + 0.5f, pos.y + h - 2.0f),
                                      style::u32(style::col::BorderStrong));
  ImGui::SameLine(0.0f, 10.0f);
}

}  // namespace

void drawStatusBar(AppState& st) {
  ImGuiViewport* viewport = ImGui::GetMainViewport();
  const float height = ImGui::GetFrameHeight() + 2.0f;
  constexpr ImGuiWindowFlags flags =
      ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDecoration |
      ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoSavedSettings |
      ImGuiWindowFlags_NoNavFocus | ImGuiWindowFlags_NoBringToFrontOnFocus;

  if (ImGui::BeginViewportSideBar("##status", viewport, ImGuiDir_Down, height, flags)) {
    // Hairline separating the bar from the workspace.
    const ImVec2 wmin = ImGui::GetWindowPos();
    ImGui::GetWindowDrawList()->AddLine(
        wmin, ImVec2(wmin.x + ImGui::GetWindowWidth(), wmin.y),
        style::u32(style::col::Border), style::metrics().hairline);

    if (ImGui::BeginMenuBar()) {
      const std::string idleMessage = toolHint(st);
      const char* message =
          st.statusMessage.empty() ? idleMessage.c_str() : st.statusMessage.c_str();
      ImGui::TextUnformatted(message);
      divider();

      const std::string hover = hoverDescription(st);
      ImGui::TextDisabled("%s", hover.c_str());

      // Right-aligned cluster: formula pill, zoom pill, busy indicator.
      const std::string formula =
          st.props.formula.empty() ? "--" : st.props.formula;
      char zoom[16];
      std::snprintf(zoom, sizeof(zoom), "%.0f%%", std::round(st.cam.zoom * 100.0f));

      const float padX = style::metrics().gap * 0.8f;
      const float spacing = ImGui::GetStyle().ItemSpacing.x;
      auto pillWidth = [&](const char* text) {
        return ImGui::CalcTextSize(text).x + padX * 2.0f;
      };
      float total = pillWidth(formula.c_str()) + pillWidth(zoom) + spacing;
      const bool busy = st.tasks.busy();
      if (busy) total += pillWidth("working") + spacing;

      ImGui::SetCursorPosX(ImGui::GetWindowWidth() - total -
                           ImGui::GetStyle().WindowPadding.x - 2.0f);
      const bool mono = style::pushFont(style::fonts::mono());
      widgets::badge(formula.c_str(), style::col::TextDim);
      style::popFont(mono);
      ImGui::SameLine(0.0f, spacing);
      widgets::badge(zoom, style::col::TextDim);
      if (busy) {
        ImGui::SameLine(0.0f, spacing);
        widgets::badge("working", style::col::Accent);
      }
      ImGui::EndMenuBar();
    }
  }
  ImGui::End();
}

}  // namespace chemcad::ui

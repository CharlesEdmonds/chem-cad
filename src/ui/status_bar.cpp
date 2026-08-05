#include <cmath>
#include <string>

#include "imgui.h"
#include "imgui_internal.h"

#include "ui/element_data.hpp"
#include "ui/ui.hpp"

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

void divider() {
  ImGui::SameLine();
  ImGui::TextDisabled("|");
  ImGui::SameLine();
}

}  // namespace

void drawStatusBar(AppState& st) {
  ImGuiViewport* viewport = ImGui::GetMainViewport();
  const float height = ImGui::GetFrameHeight();
  constexpr ImGuiWindowFlags flags =
      ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDecoration |
      ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoSavedSettings |
      ImGuiWindowFlags_NoNavFocus | ImGuiWindowFlags_NoBringToFrontOnFocus;

  if (ImGui::BeginViewportSideBar("##status", viewport, ImGuiDir_Down, height, flags)) {
    if (ImGui::BeginMenuBar()) {
      const char* message = st.statusMessage.empty()
                                ? "Hint: choose a tool, then draw in the Sketch canvas"
                                : st.statusMessage.c_str();
      ImGui::TextUnformatted(message);
      divider();

      const std::string hover = hoverDescription(st);
      ImGui::TextUnformatted(hover.c_str());
      divider();

      if (st.props.formula.empty()) {
        ImGui::TextDisabled("Formula: --");
      } else {
        ImGui::Text("Formula: %s", st.props.formula.c_str());
      }
      divider();

      ImGui::Text("Zoom: %.0f%%", std::round(st.cam.zoom * 100.0f));
      if (st.tasks.busy()) {
        divider();
        ImGui::TextDisabled("working...");
      }
      ImGui::EndMenuBar();
    }
  }
  ImGui::End();
}

}  // namespace chemcad::ui

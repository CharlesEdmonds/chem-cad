#include "ui/mol_thumb.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <string>

#include "ui/theme.hpp"

namespace chemcad::ui {
namespace {

constexpr std::array<const char*, 119> kSymbols = {
    "?",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg",
    "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
    "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
    "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir",
    "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv",
    "Ts", "Og"};

ImU32 elementColour(uint8_t z) {
  switch (z) {
    case 1:
      return IM_COL32(225, 225, 225, 255);
    case 6:
      return IM_COL32(235, 235, 235, 255);
    case 7:
      return IM_COL32(105, 150, 255, 255);
    case 8:
      return IM_COL32(255, 100, 100, 255);
    case 9:
    case 17:
      return IM_COL32(100, 225, 120, 255);
    case 15:
      return IM_COL32(255, 165, 70, 255);
    case 16:
      return IM_COL32(245, 210, 70, 255);
    case 35:
      return IM_COL32(185, 90, 70, 255);
    case 53:
      return IM_COL32(175, 105, 215, 255);
    default:
      return IM_COL32(180, 205, 220, 255);
  }
}

void dashedLine(ImDrawList* dl, ImVec2 a, ImVec2 b, ImU32 colour, float thickness,
                float dash = 4.0f, float gap = 3.0f) {
  const float dx = b.x - a.x;
  const float dy = b.y - a.y;
  const float length = std::sqrt(dx * dx + dy * dy);
  if (length < 0.01f) return;
  const float ux = dx / length;
  const float uy = dy / length;
  for (float at = 0.0f; at < length; at += dash + gap) {
    const float end = std::min(at + dash, length);
    dl->AddLine(ImVec2(a.x + ux * at, a.y + uy * at),
                ImVec2(a.x + ux * end, a.y + uy * end), colour, thickness);
  }
}

void dashedRect(ImDrawList* dl, ImVec2 min, ImVec2 max, ImU32 colour) {
  dashedLine(dl, min, ImVec2(max.x, min.y), colour, 1.0f);
  dashedLine(dl, ImVec2(max.x, min.y), max, colour, 1.0f);
  dashedLine(dl, max, ImVec2(min.x, max.y), colour, 1.0f);
  dashedLine(dl, ImVec2(min.x, max.y), min, colour, 1.0f);
}

std::string atomLabel(const core::Atom& atom) {
  const char* symbol =
      atom.atomicNumber < kSymbols.size() ? kSymbols[atom.atomicNumber] : kSymbols[0];
  std::string label;
  if (atom.isotope != 0) label += std::to_string(atom.isotope);
  label += symbol;
  if (atom.charge != 0) {
    const int magnitude = std::abs(static_cast<int>(atom.charge));
    if (magnitude > 1) label += std::to_string(magnitude);
    label += atom.charge > 0 ? '+' : '-';
  }
  return label;
}

}  // namespace

void drawMoleculeThumb(ImDrawList* dl, const core::Molecule& mol, ImVec2 min, ImVec2 max) {
  if (!dl || max.x <= min.x || max.y <= min.y) return;

  if (mol.empty()) {
    const ImU32 muted = ImGui::GetColorU32(ImGuiCol_TextDisabled);
    dashedRect(dl, ImVec2(min.x + 4.0f, min.y + 4.0f),
               ImVec2(max.x - 4.0f, max.y - 4.0f), muted);
    const ImVec2 textSize = ImGui::CalcTextSize("?");
    dl->AddText(ImVec2((min.x + max.x - textSize.x) * 0.5f,
                       (min.y + max.y - textSize.y) * 0.5f),
                muted, "?");
    return;
  }

  float loX = mol.atoms().front().pos.x;
  float hiX = loX;
  float loY = mol.atoms().front().pos.y;
  float hiY = loY;
  for (const core::Atom& atom : mol.atoms()) {
    loX = std::min(loX, atom.pos.x);
    hiX = std::max(hiX, atom.pos.x);
    loY = std::min(loY, atom.pos.y);
    hiY = std::max(hiY, atom.pos.y);
  }

  constexpr float kPadding = 8.0f;
  const float availableW = std::max(1.0f, max.x - min.x - 2.0f * kPadding);
  const float availableH = std::max(1.0f, max.y - min.y - 2.0f * kPadding);
  const float spanX = hiX - loX;
  const float spanY = hiY - loY;
  float scale = 1.0f;
  if (spanX > 0.001f || spanY > 0.001f) {
    const float sx = spanX > 0.001f ? availableW / spanX : availableH / spanY;
    const float sy = spanY > 0.001f ? availableH / spanY : availableW / spanX;
    scale = std::min(sx, sy);
  }

  const float sourceCx = (loX + hiX) * 0.5f;
  const float sourceCy = (loY + hiY) * 0.5f;
  const ImVec2 targetCentre((min.x + max.x) * 0.5f, (min.y + max.y) * 0.5f);
  const auto screenPos = [&](const core::Atom& atom) {
    return ImVec2(targetCentre.x + (atom.pos.x - sourceCx) * scale,
                  targetCentre.y + (atom.pos.y - sourceCy) * scale);
  };

  const ImU32 bondColour = ImGui::GetColorU32(ImGuiCol_Text);
  const ImU32 aromaticColour = ImGui::GetColorU32(ImGuiCol_TextDisabled);
  constexpr float kOffset = 2.5f;
  for (const core::Bond& bond : mol.bonds()) {
    const core::Atom* a = mol.atom(bond.a);
    const core::Atom* b = mol.atom(bond.b);
    if (!a || !b) continue;
    const ImVec2 pa = screenPos(*a);
    const ImVec2 pb = screenPos(*b);
    const float dx = pb.x - pa.x;
    const float dy = pb.y - pa.y;
    const float length = std::sqrt(dx * dx + dy * dy);
    if (length < 0.01f) continue;
    const ImVec2 normal(-dy / length * kOffset, dx / length * kOffset);
    const auto shifted = [&](ImVec2 p, float direction) {
      return ImVec2(p.x + normal.x * direction, p.y + normal.y * direction);
    };

    switch (bond.order) {
      case core::BondOrder::Single:
        dl->AddLine(pa, pb, bondColour, 1.6f);
        break;
      case core::BondOrder::Double:
        dl->AddLine(shifted(pa, -1.0f), shifted(pb, -1.0f), bondColour, 1.4f);
        dl->AddLine(shifted(pa, 1.0f), shifted(pb, 1.0f), bondColour, 1.4f);
        break;
      case core::BondOrder::Triple:
        dl->AddLine(pa, pb, bondColour, 1.3f);
        dl->AddLine(shifted(pa, -1.5f), shifted(pb, -1.5f), bondColour, 1.3f);
        dl->AddLine(shifted(pa, 1.5f), shifted(pb, 1.5f), bondColour, 1.3f);
        break;
      case core::BondOrder::Aromatic:
        dl->AddLine(shifted(pa, -0.7f), shifted(pb, -0.7f), bondColour, 1.4f);
        dashedLine(dl, shifted(pa, 1.0f), shifted(pb, 1.0f), aromaticColour, 1.2f,
                   3.0f, 2.5f);
        break;
    }
  }

  const ImU32 labelBackground = ImGui::GetColorU32(ImGuiCol_ChildBg, 0.92f);
  for (const core::Atom& atom : mol.atoms()) {
    const bool showCarbon =
        atom.atomicNumber != 6 || atom.charge != 0 || mol.degree(atom.id) == 0;
    if (!showCarbon) continue;
    const std::string label = atomLabel(atom);
    const ImVec2 textSize = ImGui::CalcTextSize(label.c_str());
    const ImVec2 centre = screenPos(atom);
    const ImVec2 textPos(centre.x - textSize.x * 0.5f, centre.y - textSize.y * 0.5f);
    dl->AddRectFilled(ImVec2(textPos.x - 1.5f, textPos.y - 0.5f),
                      ImVec2(textPos.x + textSize.x + 1.5f,
                             textPos.y + textSize.y + 0.5f),
                      labelBackground, 2.0f);
    dl->AddText(textPos, elementColour(atom.atomicNumber), label.c_str());
  }
}

bool moleculeThumbButton(const char* id, const core::Molecule& mol, ImVec2 size) {
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::InvisibleButton(id, size);
  const bool clicked = ImGui::IsItemClicked();
  const bool hovered = ImGui::IsItemHovered();
  const ImVec2 max(min.x + size.x, min.y + size.y);
  ImDrawList* dl = ImGui::GetWindowDrawList();
  const float radius = style::metrics().radiusMd;
  const ImU32 background =
      ImGui::GetColorU32(hovered ? ImGuiCol_FrameBgHovered : ImGuiCol_FrameBg);
  const ImU32 border =
      hovered ? style::u32(style::col::Accent, 0.75f) : style::u32(style::col::Border);
  dl->AddRectFilled(min, max, background, radius);
  dl->AddRect(min, max, border, radius, 0, hovered ? 1.6f : 1.0f);
  drawMoleculeThumb(dl, mol, ImVec2(min.x + 2.0f, min.y + 2.0f),
                    ImVec2(max.x - 2.0f, max.y - 2.0f));
  return clicked;
}

}  // namespace chemcad::ui

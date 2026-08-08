#include "ui/canvas_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <optional>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

#include "chem/bridge.hpp"
#include "core/sprout.hpp"
#include "ui/theme.hpp"

namespace chemcad::ui::canvas {
namespace {

struct LabelData {
  bool visible = false;
  std::string isotope;
  std::string main;
  std::string subscript;
  std::string charge;
  float fontSize = 0.0f;
  float smallSize = 0.0f;
  float isotopeWidth = 0.0f;
  float mainWidth = 0.0f;
  float subscriptWidth = 0.0f;
  float chargeWidth = 0.0f;
  ImVec2 halfExtent{};
};

using LabelMap = std::unordered_map<uint64_t, LabelData>;

uint64_t atomKey(int mol, core::AtomId id) {
  return (static_cast<uint64_t>(static_cast<uint32_t>(mol)) << 32) | id;
}

ImVec2 toIm(core::Vec2 p) { return {p.x, p.y}; }

ImU32 atomColor(uint8_t z) {
  switch (z) {
    case 6: return IM_COL32(224, 228, 234, 255);
    case 7: return IM_COL32(90, 145, 255, 255);
    case 8: return IM_COL32(246, 82, 82, 255);
    case 9:
    case 17: return IM_COL32(92, 226, 112, 255);
    case 16: return IM_COL32(245, 211, 74, 255);
    case 15: return IM_COL32(255, 156, 70, 255);
    case 35: return IM_COL32(150, 66, 49, 255);
    case 53: return IM_COL32(164, 104, 224, 255);
    default: return IM_COL32(193, 158, 229, 255);
  }
}

ImVec2 textSize(ImFont* font, float size, const std::string& text) {
  if (text.empty()) return {};
  return font->CalcTextSizeA(size, 10000.0f, 0.0f, text.c_str());
}

int implicitHydrogens(AppState& st, Runtime& rt, int molIndex, const core::Molecule& mol,
                      core::AtomId atomId) {
  if (rt.hydrogenRevision != st.docRevision) {
    rt.hydrogenRevision = st.docRevision;
    rt.hydrogenCounts.clear();
  }
  const uint64_t key = atomKey(molIndex, atomId);
  if (const auto found = rt.hydrogenCounts.find(key); found != rt.hydrogenCounts.end()) {
    return found->second;
  }
  int result = 0;
  try {
    result = std::max(0, chem::implicitHCount(mol, atomId));
  } catch (const std::exception& error) {
    st.statusMessage = std::string("Hydrogen display: ") + error.what();
  }
  rt.hydrogenCounts.emplace(key, result);
  return result;
}

std::string chargeText(int charge) {
  if (charge == 0) return {};
  std::string text;
  const int magnitude = std::abs(charge);
  if (magnitude > 1) text = std::to_string(magnitude);
  text += charge > 0 ? "+" : "-";
  return text;
}

LabelData makeLabel(AppState& st, Runtime& rt, int molIndex, const core::Molecule& mol,
                    const core::Atom& atom) {
  LabelData label;
  const int degree = mol.degree(atom.id);
  // Skeletal notation: a plain carbon vertex is an unlabelled line end or
  // corner. Only carbons carrying a charge, an isotope or no bonds at all get
  // a written symbol; terminal CH3 groups stay implicit like every other
  // carbon.
  if (atom.atomicNumber == 6 && atom.charge == 0 && atom.isotope == 0 && degree > 0) {
    return label;
  }

  int hydrogens = implicitHydrogens(st, rt, molIndex, mol, atom.id);
  if (atom.atomicNumber == 6 && degree == 0 && atom.explicitH < 0 && hydrogens == 0)
    hydrogens = 4;
  label.visible = true;

  std::string symbol = "?";
  try {
    if (const char* value = chem::symbolFor(atom.atomicNumber); value && *value) symbol = value;
  } catch (const std::exception& error) {
    st.statusMessage = std::string("Element display: ") + error.what();
  }

  label.fontSize = std::clamp(ImGui::GetFontSize() * std::sqrt(st.cam.zoom), 12.0f, 32.0f);
  label.smallSize = label.fontSize * 0.68f;
  if (atom.isotope != 0) label.isotope = std::to_string(atom.isotope);
  label.main = symbol;
  if (hydrogens > 0) label.main += "H";
  if (hydrogens > 1) label.subscript = std::to_string(hydrogens);
  label.charge = chargeText(atom.charge);

  ImFont* font = ImGui::GetFont();
  label.isotopeWidth = textSize(font, label.smallSize, label.isotope).x;
  label.mainWidth = textSize(font, label.fontSize, label.main).x;
  label.subscriptWidth = textSize(font, label.smallSize, label.subscript).x;
  label.chargeWidth = textSize(font, label.smallSize, label.charge).x;
  const float totalWidth = label.isotopeWidth + label.mainWidth + label.subscriptWidth +
                           label.chargeWidth;
  label.halfExtent = {totalWidth * 0.5f + 1.0f, label.fontSize * 0.62f + label.smallSize * 0.2f};
  return label;
}

LabelMap buildLabels(AppState& st, Runtime& rt) {
  LabelMap labels;
  for (int mi = 0; mi < static_cast<int>(st.doc.molecules.size()); ++mi) {
    const core::Molecule& mol = st.doc.molecules[static_cast<size_t>(mi)];
    for (const core::Atom& atom : mol.atoms()) {
      LabelData label = makeLabel(st, rt, mi, mol, atom);
      if (label.visible) labels.emplace(atomKey(mi, atom.id), std::move(label));
    }
  }
  return labels;
}

std::optional<core::Vec2> smallestCycleCentroid(const core::Molecule& mol,
                                                 const core::Bond& excluded) {
  std::queue<core::AtomId> pending;
  std::unordered_map<core::AtomId, core::AtomId> parent;
  pending.push(excluded.a);
  parent.emplace(excluded.a, core::kInvalidAtom);
  while (!pending.empty() && !parent.contains(excluded.b)) {
    const core::AtomId current = pending.front();
    pending.pop();
    for (core::BondId bondId : mol.incidentBonds(current)) {
      if (bondId == excluded.id) continue;
      const core::Bond* bond = mol.bond(bondId);
      if (!bond) continue;
      const core::AtomId next = bond->a == current ? bond->b : bond->a;
      if (parent.contains(next)) continue;
      parent.emplace(next, current);
      pending.push(next);
    }
  }
  if (!parent.contains(excluded.b)) return std::nullopt;

  core::Vec2 centroid{};
  int count = 0;
  for (core::AtomId id = excluded.b; id != core::kInvalidAtom; id = parent[id]) {
    if (const core::Atom* atom = mol.atom(id)) {
      centroid.x += atom->pos.x;
      centroid.y += atom->pos.y;
      ++count;
    }
  }
  if (count == 0) return std::nullopt;
  centroid.x /= static_cast<float>(count);
  centroid.y /= static_cast<float>(count);
  return centroid;
}

void shortenForLabels(ImVec2& a, ImVec2& b, const LabelData* aLabel, const LabelData* bLabel) {
  const float dx = b.x - a.x;
  const float dy = b.y - a.y;
  const float length = std::sqrt(dx * dx + dy * dy);
  if (length < 1e-3f) return;
  const ImVec2 direction{dx / length, dy / length};
  auto trim = [&](const LabelData* label) {
    if (!label) return 0.0f;
    return std::abs(direction.x) * label->halfExtent.x +
           std::abs(direction.y) * label->halfExtent.y + 2.0f;
  };
  float trimA = trim(aLabel);
  float trimB = trim(bLabel);
  const float total = trimA + trimB;
  if (total > length * 0.8f && total > 0.0f) {
    const float factor = length * 0.8f / total;
    trimA *= factor;
    trimB *= factor;
  }
  a.x += direction.x * trimA;
  a.y += direction.y * trimA;
  b.x -= direction.x * trimB;
  b.y -= direction.y * trimB;
}

void dashedLine(ImDrawList* draw, ImVec2 a, ImVec2 b, ImU32 color, float thickness) {
  constexpr int segments = 5;
  for (int i = 0; i < segments; ++i) {
    const float t0 = static_cast<float>(i) / static_cast<float>(segments);
    const float t1 = std::min(1.0f, t0 + 0.11f);
    draw->AddLine({a.x + (b.x - a.x) * t0, a.y + (b.y - a.y) * t0},
                  {a.x + (b.x - a.x) * t1, a.y + (b.y - a.y) * t1}, color, thickness);
  }
}

void drawBond(ImDrawList* draw, const AppState& st, const CanvasRect& rect, int molIndex,
              const core::Molecule& mol, const core::Bond& bond, const LabelMap& labels) {
  const core::Atom* atomA = mol.atom(bond.a);
  const core::Atom* atomB = mol.atom(bond.b);
  if (!atomA || !atomB) return;
  ImVec2 a = toIm(st.cam.worldToScreen(atomA->pos, rect.origin));
  ImVec2 b = toIm(st.cam.worldToScreen(atomB->pos, rect.origin));
  const auto aIt = labels.find(atomKey(molIndex, atomA->id));
  const auto bIt = labels.find(atomKey(molIndex, atomB->id));
  shortenForLabels(a, b, aIt == labels.end() ? nullptr : &aIt->second,
                   bIt == labels.end() ? nullptr : &bIt->second);

  const float dx = b.x - a.x;
  const float dy = b.y - a.y;
  const float length = std::sqrt(dx * dx + dy * dy);
  if (length < 1e-3f) return;
  const ImVec2 perpendicular{-dy / length, dx / length};
  const float thickness = std::max(1.2f, 2.0f * st.cam.zoom);
  const BondRef ref{molIndex, bond.id};
  const bool selected = st.sel.contains(ref);
  const ImU32 normal = IM_COL32(224, 228, 234, 255);
  const ImU32 accent = style::u32(style::col::Accent);
  const ImU32 color = selected ? accent : normal;

  float inwardSign = 0.0f;
  if (const auto centroid = smallestCycleCentroid(mol, bond)) {
    const ImVec2 c = toIm(st.cam.worldToScreen(*centroid, rect.origin));
    const ImVec2 midpoint{(a.x + b.x) * 0.5f, (a.y + b.y) * 0.5f};
    inwardSign = ((c.x - midpoint.x) * perpendicular.x +
                   (c.y - midpoint.y) * perpendicular.y) >= 0.0f
                      ? 1.0f
                      : -1.0f;
  }
  auto shifted = [&](ImVec2 p, float worldOffset) {
    const float offset = worldOffset * st.cam.scale();
    return ImVec2{p.x + perpendicular.x * offset, p.y + perpendicular.y * offset};
  };

  if (bond.stereo == core::BondStereo::Wedge) {
    const float halfWidth = 0.14f * st.cam.scale();
    draw->AddTriangleFilled(a, {b.x + perpendicular.x * halfWidth, b.y + perpendicular.y * halfWidth},
                            {b.x - perpendicular.x * halfWidth, b.y - perpendicular.y * halfWidth},
                            color);
  } else if (bond.stereo == core::BondStereo::Hash) {
    constexpr int strokes = 6;
    for (int i = 0; i < strokes; ++i) {
      const float t = (static_cast<float>(i) + 1.0f) / (static_cast<float>(strokes) + 1.0f);
      const ImVec2 center{a.x + dx * t, a.y + dy * t};
      const float halfWidth = 0.14f * st.cam.scale() * t;
      draw->AddLine({center.x - perpendicular.x * halfWidth,
                     center.y - perpendicular.y * halfWidth},
                    {center.x + perpendicular.x * halfWidth,
                     center.y + perpendicular.y * halfWidth}, color, thickness);
    }
  } else if (bond.stereo == core::BondStereo::Wavy) {
    // Squiggly bond: stereochemistry unknown or a mixture (MDL "either").
    const float amplitude = 0.05f * st.cam.scale();
    constexpr int segments = 18;
    ImVec2 points[segments + 1];
    for (int i = 0; i <= segments; ++i) {
      const float t = static_cast<float>(i) / segments;
      const float wave = std::sin(t * 3.0f * 6.2831853f) * amplitude;
      points[i] = ImVec2{a.x + dx * t + perpendicular.x * wave,
                         a.y + dy * t + perpendicular.y * wave};
    }
    draw->AddPolyline(points, segments + 1, color, 0, thickness);
  } else {
    switch (bond.order) {
      case core::BondOrder::Single:
        draw->AddLine(a, b, color, thickness);
        break;
      case core::BondOrder::Double:
        if (inwardSign != 0.0f) {
          draw->AddLine(shifted(a, 0.015f * inwardSign), shifted(b, 0.015f * inwardSign),
                        color, thickness);
          draw->AddLine(shifted(a, 0.165f * inwardSign), shifted(b, 0.165f * inwardSign),
                        color, thickness);
        } else {
          draw->AddLine(shifted(a, -0.09f), shifted(b, -0.09f), color, thickness);
          draw->AddLine(shifted(a, 0.09f), shifted(b, 0.09f), color, thickness);
        }
        break;
      case core::BondOrder::Triple:
        draw->AddLine(a, b, color, thickness);
        draw->AddLine(shifted(a, -0.12f), shifted(b, -0.12f), color, thickness);
        draw->AddLine(shifted(a, 0.12f), shifted(b, 0.12f), color, thickness);
        break;
      case core::BondOrder::Aromatic: {
        draw->AddLine(a, b, color, thickness);
        const float sign = inwardSign == 0.0f ? 1.0f : inwardSign;
        ImVec2 innerA = shifted(a, 0.15f * sign);
        ImVec2 innerB = shifted(b, 0.15f * sign);
        innerA = {innerA.x + (innerB.x - innerA.x) * 0.18f,
                  innerA.y + (innerB.y - innerA.y) * 0.18f};
        innerB = {innerB.x + (innerA.x - innerB.x) * 0.18f,
                  innerB.y + (innerA.y - innerB.y) * 0.18f};
        dashedLine(draw, innerA, innerB, color, thickness * 0.85f);
        break;
      }
    }
  }

  if (st.hoverBond == ref) {
    draw->AddLine(a, b, style::u32(style::col::Teal, 0.73f), thickness + 3.0f);
  }
}

void drawLabel(ImDrawList* draw, const AppState& st, const CanvasRect& rect, int molIndex,
               const core::Atom& atom, const LabelData& label, ImU32 background) {
  const ImVec2 center = toIm(st.cam.worldToScreen(atom.pos, rect.origin));
  const AtomRef ref{molIndex, atom.id};
  const bool selected = st.sel.contains(ref);
  const bool hovered = st.hoverAtom == ref;
  const float radius = std::max(label.halfExtent.x, label.halfExtent.y) + 3.0f;
  draw->AddCircleFilled(center, radius, background, 24);

  ImFont* font = ImGui::GetFont();
  const float totalWidth = label.isotopeWidth + label.mainWidth + label.subscriptWidth +
                           label.chargeWidth;
  float x = center.x - totalWidth * 0.5f;
  const float baselineY = center.y - label.fontSize * 0.5f;
  const ImU32 color = selected ? style::u32(style::col::Accent) : atomColor(atom.atomicNumber);
  if (!label.isotope.empty()) {
    draw->AddText(font, label.smallSize, {x, baselineY - label.smallSize * 0.35f}, color,
                  label.isotope.c_str());
    x += label.isotopeWidth;
  }
  draw->AddText(font, label.fontSize, {x, baselineY}, color, label.main.c_str());
  x += label.mainWidth;
  if (!label.subscript.empty()) {
    draw->AddText(font, label.smallSize, {x, baselineY + label.fontSize * 0.48f}, color,
                  label.subscript.c_str());
    x += label.subscriptWidth;
  }
  if (!label.charge.empty()) {
    draw->AddText(font, label.smallSize, {x, baselineY - label.smallSize * 0.35f}, color,
                  label.charge.c_str());
  }
  if (selected || hovered) {
    draw->AddCircle(center, radius + 2.0f, selected ? style::u32(style::col::Accent, 0.90f)
                                                  : style::u32(style::col::Teal, 0.86f),
                    24, selected ? 2.2f : 1.4f);
  }
}

void drawUnlabelledAtomFeedback(ImDrawList* draw, const AppState& st, const CanvasRect& rect,
                                int molIndex, const core::Atom& atom) {
  const AtomRef ref{molIndex, atom.id};
  const bool selected = st.sel.contains(ref);
  const bool hovered = st.hoverAtom == ref;
  if (!selected && !hovered) return;
  const ImVec2 center = toIm(st.cam.worldToScreen(atom.pos, rect.origin));
  if (selected) draw->AddCircleFilled(center, 4.0f, style::u32(style::col::Accent, 0.86f), 16);
  draw->AddCircle(center, 8.0f, hovered ? style::u32(style::col::Teal, 0.90f)
                                      : style::u32(style::col::Accent, 0.90f),
                  20, selected ? 2.2f : 1.4f);
}

void drawGesturePreview(ImDrawList* draw, const AppState& st, const Runtime& rt,
                        const CanvasRect& rect) {
  const ImU32 preview = style::u32(style::col::Teal, 0.71f);
  if (rt.gesture == Gesture::Bond && rt.dragged && rt.downAtom.valid()) {
    const core::Atom* atom = nullptr;
    if (rt.downAtom.mol >= 0 && rt.downAtom.mol < static_cast<int>(st.doc.molecules.size())) {
      atom = st.doc.molecules[static_cast<size_t>(rt.downAtom.mol)].atom(rt.downAtom.id);
    }
    if (atom) {
      const core::Vec2 direction = core::snapAngle(
          {rt.currentWorld.x - atom->pos.x, rt.currentWorld.y - atom->pos.y});
      const core::Vec2 end{atom->pos.x + direction.x * core::kBondLength,
                           atom->pos.y + direction.y * core::kBondLength};
      draw->AddLine(toIm(st.cam.worldToScreen(atom->pos, rect.origin)),
                    toIm(st.cam.worldToScreen(end, rect.origin)), preview,
                    std::max(1.5f, 2.0f * st.cam.zoom));
      draw->AddCircle(toIm(st.cam.worldToScreen(end, rect.origin)), 5.0f, preview, 16, 1.5f);
    }
  } else if (rt.gesture == Gesture::Chain) {
    const std::vector<core::Vec2> points = makeChainPreview(st, rt);
    for (size_t i = 1; i < points.size(); ++i) {
      draw->AddLine(toIm(st.cam.worldToScreen(points[i - 1], rect.origin)),
                    toIm(st.cam.worldToScreen(points[i], rect.origin)), preview,
                    std::max(1.5f, 2.0f * st.cam.zoom));
    }
  } else if (rt.gesture == Gesture::Marquee && rt.dragged) {
    draw->AddRect(rt.downScreen, rt.currentScreen, style::u32(style::col::Teal, 0.90f), 0.0f, 0, 1.5f);
    draw->AddRectFilled(rt.downScreen, rt.currentScreen, style::u32(style::col::Teal, 0.10f));
  }
}

void drawRingGhost(ImDrawList* draw, const AppState& st, const CanvasRect& rect) {
  if (st.tool != Tool::RingTemplate || !ImGui::IsItemHovered()) return;
  const ImVec2 mouse = ImGui::GetIO().MousePos;
  const core::Vec2 cursor = st.cam.screenToWorld({mouse.x, mouse.y}, rect.origin);
  const RingGeometry geometry = makeRingGeometry(st, cursor);
  for (size_t i = 0; i < geometry.edges.size(); ++i) {
    const auto [a, b] = geometry.edges[i];
    const ImVec2 from = toIm(st.cam.worldToScreen(geometry.positions[static_cast<size_t>(a)],
                                                  rect.origin));
    const ImVec2 to = toIm(st.cam.worldToScreen(geometry.positions[static_cast<size_t>(b)],
                                                rect.origin));
    draw->AddLine(from, to, IM_COL32(126, 205, 255, 130), 1.6f);
    if (geometry.orders[i] == core::BondOrder::Double ||
        geometry.orders[i] == core::BondOrder::Aromatic) {
      const ImVec2 delta{to.x - from.x, to.y - from.y};
      const float len = std::sqrt(delta.x * delta.x + delta.y * delta.y);
      if (len > 0.1f) {
        const ImVec2 perp{-delta.y / len * 5.0f, delta.x / len * 5.0f};
        draw->AddLine({from.x + perp.x, from.y + perp.y}, {to.x + perp.x, to.y + perp.y},
                      IM_COL32(126, 205, 255, 100), 1.2f);
      }
    }
  }
}

}  // namespace

void render(AppState& st, Runtime& rt, const CanvasRect& rect) {
  ImDrawList* draw = ImGui::GetWindowDrawList();
  const ImU32 background = IM_COL32(35, 39, 47, 255);
  draw->AddRectFilled(rect.min, rect.max, background);

  const LabelMap labels = buildLabels(st, rt);
  for (int mi = 0; mi < static_cast<int>(st.doc.molecules.size()); ++mi) {
    const core::Molecule& mol = st.doc.molecules[static_cast<size_t>(mi)];
    for (const core::Bond& bond : mol.bonds()) drawBond(draw, st, rect, mi, mol, bond, labels);
  }
  for (int mi = 0; mi < static_cast<int>(st.doc.molecules.size()); ++mi) {
    const core::Molecule& mol = st.doc.molecules[static_cast<size_t>(mi)];
    for (const core::Atom& atom : mol.atoms()) {
      const auto found = labels.find(atomKey(mi, atom.id));
      if (found != labels.end()) {
        drawLabel(draw, st, rect, mi, atom, found->second, background);
      } else {
        drawUnlabelledAtomFeedback(draw, st, rect, mi, atom);
      }
    }
  }
  drawRingGhost(draw, st, rect);
  drawGesturePreview(draw, st, rt, rect);
}

}  // namespace chemcad::ui::canvas

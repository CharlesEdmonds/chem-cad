#include "ui/canvas_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>

#include "chem/bridge.hpp"
#include "core/sprout.hpp"

namespace chemcad::ui::canvas {
namespace {

constexpr float kPi = 3.14159265358979323846f;

float distanceSquared(core::Vec2 a, core::Vec2 b) {
  const float dx = a.x - b.x;
  const float dy = a.y - b.y;
  return dx * dx + dy * dy;
}

float segmentDistance(core::Vec2 p, core::Vec2 a, core::Vec2 b) {
  const float dx = b.x - a.x;
  const float dy = b.y - a.y;
  const float denom = dx * dx + dy * dy;
  if (denom <= 1e-8f) return std::sqrt(distanceSquared(p, a));
  const float t = std::clamp(((p.x - a.x) * dx + (p.y - a.y) * dy) / denom, 0.0f, 1.0f);
  return std::sqrt(distanceSquared(p, {a.x + t * dx, a.y + t * dy}));
}

void clearTransientRefs(AppState& st) {
  st.sel.clear();
  st.hoverAtom = {};
  st.hoverBond = {};
}

void markChanged(AppState& st) {
  st.touch();
  st.hoverAtom = {};
  st.hoverBond = {};
}

core::BondOrder nextOrder(core::BondOrder order) {
  switch (order) {
    case core::BondOrder::Single: return core::BondOrder::Double;
    case core::BondOrder::Double: return core::BondOrder::Triple;
    default: return core::BondOrder::Single;
  }
}

core::BondId addStyledBond(core::Molecule& mol, core::AtomId a, core::AtomId b,
                           core::BondOrder order, core::BondStereo stereo) {
  const core::BondId id = mol.addBond(a, b, order);
  if (core::Bond* bond = mol.bond(id)) bond->stereo = stereo;
  return id;
}

const core::Atom* atomAt(const AppState& st, AtomRef ref) {
  if (!ref.valid() || ref.mol >= static_cast<int>(st.doc.molecules.size())) return nullptr;
  return st.doc.molecules[static_cast<size_t>(ref.mol)].atom(ref.id);
}

const core::Bond* bondAt(const AppState& st, BondRef ref) {
  if (!ref.valid() || ref.mol >= static_cast<int>(st.doc.molecules.size())) return nullptr;
  return st.doc.molecules[static_cast<size_t>(ref.mol)].bond(ref.id);
}

AtomRef nearestAtomInMolecule(const AppState& st, int molIndex, core::Vec2 p, float radius,
                              core::AtomId exclude = core::kInvalidAtom) {
  AtomRef result;
  if (molIndex < 0 || molIndex >= static_cast<int>(st.doc.molecules.size())) return result;
  float best = radius * radius;
  for (const core::Atom& atom : st.doc.molecules[static_cast<size_t>(molIndex)].atoms()) {
    if (atom.id == exclude) continue;
    const float d = distanceSquared(atom.pos, p);
    if (d <= best) {
      best = d;
      result = {molIndex, atom.id};
    }
  }
  return result;
}

void removeEmptyFragments(AppState& st) {
  std::erase_if(st.doc.molecules, [](const core::Molecule& mol) { return mol.empty(); });
}

void eraseHovered(AppState& st) {
  if (st.hoverAtom.valid()) {
    const AtomRef ref = st.hoverAtom;
    core::Molecule* mol = st.molecule(ref.mol);
    if (!mol || !mol->atom(ref.id)) return;
    st.snapshot();
    mol->removeAtom(ref.id);
    removeEmptyFragments(st);
    clearTransientRefs(st);
    st.touch();
    return;
  }
  if (st.hoverBond.valid()) {
    const BondRef ref = st.hoverBond;
    core::Molecule* mol = st.molecule(ref.mol);
    if (!mol || !mol->bond(ref.id)) return;
    st.snapshot();
    mol->removeBond(ref.id);
    clearTransientRefs(st);
    st.touch();
  }
}

void deleteSelection(AppState& st) {
  if (st.sel.empty()) return;
  st.snapshot();
  for (int molIndex = static_cast<int>(st.doc.molecules.size()) - 1; molIndex >= 0; --molIndex) {
    core::Molecule& mol = st.doc.molecules[static_cast<size_t>(molIndex)];
    for (const BondRef& ref : st.sel.bonds) {
      if (ref.mol == molIndex) mol.removeBond(ref.id);
    }
    for (const AtomRef& ref : st.sel.atoms) {
      if (ref.mol == molIndex) mol.removeAtom(ref.id);
    }
  }
  removeEmptyFragments(st);
  clearTransientRefs(st);
  st.touch();
}

void undoOrRedo(AppState& st, bool redo) {
  const bool changed = redo ? st.undo.redo(st.doc) : st.undo.undo(st.doc);
  if (!changed) return;
  clearTransientRefs(st);
  st.touch();
}

uint64_t refKey(int mol, core::AtomId id) {
  return (static_cast<uint64_t>(static_cast<uint32_t>(mol)) << 32) | id;
}

core::Molecule moleculeForClipboard(const AppState& st) {
  core::Molecule result;
  std::unordered_map<uint64_t, core::AtomId> ids;
  const bool wholeDocument = st.sel.empty();

  auto selectedAtom = [&](int mol, core::AtomId id) {
    if (wholeDocument) return true;
    if (st.sel.contains(AtomRef{mol, id})) return true;
    for (const BondRef& ref : st.sel.bonds) {
      if (ref.mol != mol) continue;
      const core::Bond* bond = st.doc.molecules[static_cast<size_t>(mol)].bond(ref.id);
      if (bond && (bond->a == id || bond->b == id)) return true;
    }
    return false;
  };

  for (int mi = 0; mi < static_cast<int>(st.doc.molecules.size()); ++mi) {
    const core::Molecule& mol = st.doc.molecules[static_cast<size_t>(mi)];
    for (const core::Atom& atom : mol.atoms()) {
      if (!selectedAtom(mi, atom.id)) continue;
      core::Atom copy = atom;
      copy.id = core::kInvalidAtom;
      ids.emplace(refKey(mi, atom.id), result.addAtom(copy));
    }
  }
  for (int mi = 0; mi < static_cast<int>(st.doc.molecules.size()); ++mi) {
    const core::Molecule& mol = st.doc.molecules[static_cast<size_t>(mi)];
    for (const core::Bond& bond : mol.bonds()) {
      const auto a = ids.find(refKey(mi, bond.a));
      const auto b = ids.find(refKey(mi, bond.b));
      if (a == ids.end() || b == ids.end()) continue;
      const core::BondId newId = result.addBond(a->second, b->second, bond.order);
      if (core::Bond* copy = result.bond(newId)) copy->stereo = bond.stereo;
    }
  }
  return result;
}

void copySelection(AppState& st) {
  try {
    core::Molecule copy = moleculeForClipboard(st);
    if (copy.empty()) {
      st.statusMessage = "Nothing to copy";
      return;
    }
    st.clipboardSmiles = chem::toSmiles(copy);
    ImGui::SetClipboardText(st.clipboardSmiles.c_str());
    st.statusMessage = "Copied SMILES";
  } catch (const std::exception& error) {
    st.statusMessage = std::string("Copy failed: ") + error.what();
  }
}

void pasteClipboard(AppState& st) {
  const char* text = ImGui::GetClipboardText();
  const std::string smiles = text && *text ? text : st.clipboardSmiles;
  if (smiles.empty()) {
    st.statusMessage = "Clipboard does not contain SMILES";
    return;
  }
  try {
    core::Molecule pasted = chem::fromSmiles(smiles);
    float dx = 0.6f;
    float dy = -0.6f;
    if (!st.doc.empty() && !pasted.empty()) {
      float existingMaxX = -std::numeric_limits<float>::infinity();
      float pastedMinX = std::numeric_limits<float>::infinity();
      for (const core::Molecule& mol : st.doc.molecules) {
        for (const core::Atom& atom : mol.atoms()) existingMaxX = std::max(existingMaxX, atom.pos.x);
      }
      for (const core::Atom& atom : pasted.atoms()) pastedMinX = std::min(pastedMinX, atom.pos.x);
      if (std::isfinite(existingMaxX) && std::isfinite(pastedMinX)) {
        dx = existingMaxX - pastedMinX + 1.5f;
      }
    }
    for (core::Atom& atom : pasted.mutableAtoms()) {
      atom.pos.x += dx;
      atom.pos.y += dy;
    }
    st.snapshot();
    st.doc.molecules.push_back(std::move(pasted));
    st.touch();
    st.statusMessage = "Pasted structure";
  } catch (const std::exception& error) {
    st.statusMessage = std::string("Paste failed: ") + error.what();
  }
}

void retypeHoveredAtom(AppState& st, uint8_t atomicNumber) {
  core::Atom* atom = st.atomAt(st.hoverAtom);
  if (!atom || atom->atomicNumber == atomicNumber) return;
  st.snapshot();
  atom = st.atomAt(st.hoverAtom);
  if (atom) atom->atomicNumber = atomicNumber;
  markChanged(st);
}

void setHoveredBondOrder(AppState& st, core::BondOrder order) {
  core::Bond* bond = st.bondAt(st.hoverBond);
  if (!bond || bond->order == order) return;
  st.snapshot();
  bond = st.bondAt(st.hoverBond);
  if (bond) bond->order = order;
  markChanged(st);
}

void toggleHoveredStereo(AppState& st, core::BondStereo desired) {
  core::Bond* bond = st.bondAt(st.hoverBond);
  if (!bond) return;
  st.snapshot();
  bond = st.bondAt(st.hoverBond);
  if (bond) bond->stereo = bond->stereo == desired ? core::BondStereo::None : desired;
  markChanged(st);
}

std::vector<core::Vec2> regularRing(int count, core::Vec2 center, float phase) {
  const float radius = core::kBondLength / (2.0f * std::sin(kPi / static_cast<float>(count)));
  std::vector<core::Vec2> points;
  points.reserve(static_cast<size_t>(count));
  for (int i = 0; i < count; ++i) {
    const float angle = phase + 2.0f * kPi * static_cast<float>(i) / static_cast<float>(count);
    points.push_back({center.x + radius * std::cos(angle), center.y + radius * std::sin(angle)});
  }
  return points;
}

core::Vec2 meanPosition(const core::Molecule& mol, core::AtomId skipA = core::kInvalidAtom,
                        core::AtomId skipB = core::kInvalidAtom) {
  core::Vec2 mean{};
  int count = 0;
  for (const core::Atom& atom : mol.atoms()) {
    if (atom.id == skipA || atom.id == skipB) continue;
    mean.x += atom.pos.x;
    mean.y += atom.pos.y;
    ++count;
  }
  if (count > 0) {
    mean.x /= static_cast<float>(count);
    mean.y /= static_cast<float>(count);
  }
  return mean;
}

std::vector<core::Vec2> fusedPolygon(core::Vec2 a, core::Vec2 b, int count, int sign) {
  std::vector<core::Vec2> points;
  points.reserve(static_cast<size_t>(count));
  points.push_back(a);
  points.push_back(b);
  float angle = std::atan2(b.y - a.y, b.x - a.x);
  for (int i = 2; i < count; ++i) {
    angle += static_cast<float>(sign) * 2.0f * kPi / static_cast<float>(count);
    points.push_back({points.back().x + std::cos(angle), points.back().y + std::sin(angle)});
  }
  return points;
}

core::Vec2 pointsMean(const std::vector<core::Vec2>& points) {
  core::Vec2 mean{};
  if (points.empty()) return mean;
  for (core::Vec2 p : points) {
    mean.x += p.x;
    mean.y += p.y;
  }
  mean.x /= static_cast<float>(points.size());
  mean.y /= static_cast<float>(points.size());
  return mean;
}

int ringSize(RingKind kind) {
  switch (kind) {
    case RingKind::Cyclopropane: return 3;
    case RingKind::Cyclobutane: return 4;
    case RingKind::Cyclopentane:
    case RingKind::Cyclopentadiene: return 5;
    case RingKind::Cyclohexane:
    case RingKind::Benzene: return 6;
    case RingKind::Cycloheptane: return 7;
    case RingKind::Cyclooctane: return 8;
    case RingKind::Naphthalene: return 10;
  }
  return 6;
}

void fillSimpleEdges(RingGeometry& geometry, RingKind kind) {
  const int count = static_cast<int>(geometry.positions.size());
  for (int i = 0; i < count; ++i) {
    geometry.edges.emplace_back(i, (i + 1) % count);
    core::BondOrder order = core::BondOrder::Single;
    if (kind == RingKind::Benzene) order = core::BondOrder::Aromatic;
    if (kind == RingKind::Cyclopentadiene && (i == 0 || i == 2)) {
      order = core::BondOrder::Double;
    }
    geometry.orders.push_back(order);
  }
}

struct Transform2D {
  core::Vec2 from;
  core::Vec2 to;
  core::Vec2 targetFrom;
  core::Vec2 targetTo;
  bool mirror = false;

  core::Vec2 apply(core::Vec2 p) const {
    const core::Vec2 u{to.x - from.x, to.y - from.y};
    const core::Vec2 v{targetTo.x - targetFrom.x, targetTo.y - targetFrom.y};
    const float sourceLen = std::max(core::length(u), 1e-6f);
    const float targetLen = core::length(v);
    const float x = ((p.x - from.x) * u.x + (p.y - from.y) * u.y) / sourceLen;
    float y = (-(p.x - from.x) * u.y + (p.y - from.y) * u.x) / sourceLen;
    if (mirror) y = -y;
    const core::Vec2 targetUnit = core::normalize(v);
    const core::Vec2 targetPerp{-targetUnit.y, targetUnit.x};
    const float scale = targetLen / sourceLen;
    return {targetFrom.x + targetUnit.x * x * scale + targetPerp.x * y * scale,
            targetFrom.y + targetUnit.y * x * scale + targetPerp.y * y * scale};
  }
};

RingGeometry naphthaleneGeometry() {
  RingGeometry geometry;
  constexpr float h = 0.8660254037844386f;
  geometry.positions = {{-h, 1.0f}, {-2.0f * h, 0.5f}, {-2.0f * h, -0.5f},
                        {-h, -1.0f}, {0.0f, -0.5f}, {h, -1.0f},
                        {2.0f * h, -0.5f}, {2.0f * h, 0.5f}, {h, 1.0f},
                        {0.0f, 0.5f}};
  for (int i = 0; i < 10; ++i) geometry.edges.emplace_back(i, (i + 1) % 10);
  geometry.edges.emplace_back(9, 4);
  geometry.orders.assign(geometry.edges.size(), core::BondOrder::Aromatic);
  geometry.reuse.resize(geometry.positions.size());
  return geometry;
}

void transformNaphthalene(RingGeometry& geometry, const Transform2D& transform) {
  for (core::Vec2& p : geometry.positions) p = transform.apply(p);
}

void stampRing(AppState& st, core::Vec2 cursor) {
  RingGeometry geometry = makeRingGeometry(st, cursor);
  if (geometry.positions.empty()) return;

  int molIndex = -1;
  for (const AtomRef& ref : geometry.reuse) {
    if (ref.valid()) {
      molIndex = ref.mol;
      break;
    }
  }
  st.snapshot();
  if (molIndex < 0) {
    st.doc.molecules.emplace_back();
    molIndex = static_cast<int>(st.doc.molecules.size()) - 1;
  }
  core::Molecule& mol = st.doc.molecules[static_cast<size_t>(molIndex)];
  std::vector<core::AtomId> ids(geometry.positions.size(), core::kInvalidAtom);
  for (size_t i = 0; i < geometry.positions.size(); ++i) {
    if (i < geometry.reuse.size() && geometry.reuse[i].valid()) {
      ids[i] = geometry.reuse[i].id;
    } else {
      core::Atom atom;
      atom.atomicNumber = 6;
      atom.pos = geometry.positions[i];
      ids[i] = mol.addAtom(atom);
    }
  }
  for (size_t i = 0; i < geometry.edges.size(); ++i) {
    const auto [from, to] = geometry.edges[i];
    const core::BondOrder order = geometry.orders[i];
    const core::BondId existing = mol.bondBetween(ids[static_cast<size_t>(from)],
                                                   ids[static_cast<size_t>(to)]);
    if (existing != core::kInvalidBond) {
      if (core::Bond* bond = mol.bond(existing)) bond->order = order;
    } else {
      mol.addBond(ids[static_cast<size_t>(from)], ids[static_cast<size_t>(to)], order);
    }
  }
  markChanged(st);
}

void finishBondGesture(AppState& st, Runtime& rt) {
  const bool drag = rt.dragged;
  if (rt.downBond.valid() && !drag) {
    core::Bond* bond = st.bondAt(rt.downBond);
    if (!bond) return;
    st.snapshot();
    bond = st.bondAt(rt.downBond);
    if (bond) bond->order = nextOrder(bond->order);
    markChanged(st);
    return;
  }

  if (rt.downAtom.valid()) {
    core::Molecule* mol = st.molecule(rt.downAtom.mol);
    const core::Atom* from = mol ? mol->atom(rt.downAtom.id) : nullptr;
    if (!mol || !from) return;
    const core::Vec2 fromPos = from->pos;
    core::Vec2 end;
    if (drag) {
      const core::Vec2 snapped = core::snapAngle(
          {rt.currentWorld.x - fromPos.x, rt.currentWorld.y - fromPos.y});
      end = {fromPos.x + snapped.x * core::kBondLength,
             fromPos.y + snapped.y * core::kBondLength};
    } else {
      end = core::sproutPosition(*mol, rt.downAtom.id);
    }
    AtomRef target = drag ? nearestAtomInMolecule(st, rt.downAtom.mol, end, 0.35f,
                                                  rt.downAtom.id)
                          : AtomRef{};
    if (target.valid() && mol->bondBetween(rt.downAtom.id, target.id) != core::kInvalidBond) {
      st.statusMessage = "Those atoms are already bonded";
      return;
    }
    st.snapshot();
    core::AtomId targetId = target.id;
    if (!target.valid()) {
      core::Atom atom;
      atom.atomicNumber = 6;
      atom.pos = end;
      targetId = mol->addAtom(atom);
    }
    addStyledBond(*mol, rt.downAtom.id, targetId, st.currentOrder, st.currentStereo);
    markChanged(st);
    return;
  }

  if (!drag && !rt.downBond.valid()) {
    st.snapshot();
    st.doc.molecules.emplace_back();
    core::Molecule& mol = st.doc.molecules.back();
    core::Atom first;
    first.atomicNumber = 6;
    first.pos = rt.downWorld;
    const core::AtomId a = mol.addAtom(first);
    core::Atom second;
    second.atomicNumber = 6;
    second.pos = core::sproutPosition(mol, a);
    const core::AtomId b = mol.addAtom(second);
    addStyledBond(mol, a, b, st.currentOrder, st.currentStereo);
    markChanged(st);
  }
}

void finishChainGesture(AppState& st, Runtime& rt) {
  const std::vector<core::Vec2> points = makeChainPreview(st, rt);
  if (points.size() < 2) return;
  int molIndex = rt.downAtom.valid() ? rt.downAtom.mol : -1;
  st.snapshot();
  core::AtomId previous = core::kInvalidAtom;
  size_t firstNew = 0;
  if (molIndex >= 0) {
    previous = rt.downAtom.id;
    firstNew = 1;
  } else {
    st.doc.molecules.emplace_back();
    molIndex = static_cast<int>(st.doc.molecules.size()) - 1;
  }
  core::Molecule& mol = st.doc.molecules[static_cast<size_t>(molIndex)];
  for (size_t i = firstNew; i < points.size(); ++i) {
    core::Atom atom;
    atom.atomicNumber = 6;
    atom.pos = points[i];
    const core::AtomId id = mol.addAtom(atom);
    if (previous != core::kInvalidAtom) mol.addBond(previous, id, core::BondOrder::Single);
    previous = id;
  }
  markChanged(st);
}

void startSelectionGesture(AppState& st, Runtime& rt, const ImGuiIO& io) {
  rt.shiftAtStart = io.KeyShift;
  if (st.hoverAtom.valid()) {
    if (!rt.shiftAtStart && !st.sel.contains(st.hoverAtom)) st.sel.clear();
    if (!st.sel.contains(st.hoverAtom)) st.sel.atoms.push_back(st.hoverAtom);
    rt.gesture = Gesture::Move;
    rt.moveOrigins.clear();
    for (const AtomRef& ref : st.sel.atoms) {
      if (const core::Atom* atom = atomAt(st, ref)) rt.moveOrigins.push_back({ref, atom->pos});
    }
    return;
  }
  if (st.hoverBond.valid()) {
    if (!rt.shiftAtStart) st.sel.clear();
    if (!st.sel.contains(st.hoverBond)) st.sel.bonds.push_back(st.hoverBond);
    rt.gesture = Gesture::None;
    return;
  }
  if (!rt.shiftAtStart) st.sel.clear();
  rt.gesture = Gesture::Marquee;
}

void updateMoveGesture(AppState& st, Runtime& rt) {
  const core::Vec2 delta{rt.currentWorld.x - rt.downWorld.x,
                         rt.currentWorld.y - rt.downWorld.y};
  for (const MoveOrigin& origin : rt.moveOrigins) {
    if (core::Atom* atom = st.atomAt(origin.ref)) {
      atom->pos = {origin.pos.x + delta.x, origin.pos.y + delta.y};
    }
  }
}

void finishMarquee(AppState& st, Runtime& rt) {
  const float minX = std::min(rt.downWorld.x, rt.currentWorld.x);
  const float maxX = std::max(rt.downWorld.x, rt.currentWorld.x);
  const float minY = std::min(rt.downWorld.y, rt.currentWorld.y);
  const float maxY = std::max(rt.downWorld.y, rt.currentWorld.y);
  auto inside = [&](core::Vec2 p) {
    return p.x >= minX && p.x <= maxX && p.y >= minY && p.y <= maxY;
  };
  for (int mi = 0; mi < static_cast<int>(st.doc.molecules.size()); ++mi) {
    const core::Molecule& mol = st.doc.molecules[static_cast<size_t>(mi)];
    for (const core::Atom& atom : mol.atoms()) {
      const AtomRef ref{mi, atom.id};
      if (inside(atom.pos) && !st.sel.contains(ref)) st.sel.atoms.push_back(ref);
    }
    for (const core::Bond& bond : mol.bonds()) {
      const core::Atom* a = mol.atom(bond.a);
      const core::Atom* b = mol.atom(bond.b);
      const BondRef ref{mi, bond.id};
      if (a && b && inside(a->pos) && inside(b->pos) && !st.sel.contains(ref)) {
        st.sel.bonds.push_back(ref);
      }
    }
  }
}

void handleImmediateTool(AppState& st, core::Vec2 cursor) {
  switch (st.tool) {
    case Tool::Eraser:
      eraseHovered(st);
      break;
    case Tool::Atom:
      if (st.hoverAtom.valid()) {
        retypeHoveredAtom(st, st.currentElement);
      } else {
        st.snapshot();
        st.doc.molecules.emplace_back();
        core::Atom atom;
        atom.atomicNumber = st.currentElement;
        atom.pos = cursor;
        st.doc.molecules.back().addAtom(atom);
        markChanged(st);
      }
      break;
    case Tool::ChargePlus:
    case Tool::ChargeMinus: {
      core::Atom* atom = st.atomAt(st.hoverAtom);
      if (!atom) break;
      const int delta = st.tool == Tool::ChargePlus ? 1 : -1;
      const int charge = std::clamp(static_cast<int>(atom->charge) + delta, -4, 4);
      if (charge == atom->charge) break;
      st.snapshot();
      atom = st.atomAt(st.hoverAtom);
      if (atom) atom->charge = static_cast<int8_t>(charge);
      markChanged(st);
      break;
    }
    case Tool::RingTemplate:
      stampRing(st, cursor);
      break;
    default:
      break;
  }
}

void handleKeyboard(AppState& st) {
  ImGuiIO& io = ImGui::GetIO();
  if (io.WantTextInput) return;
  const bool eligible = ImGui::IsWindowHovered(ImGuiHoveredFlags_ChildWindows) ||
                        ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows);
  if (!eligible) return;

  if (io.KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_Z, false)) {
    undoOrRedo(st, io.KeyShift);
    return;
  }
  if (io.KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_Y, false)) {
    undoOrRedo(st, true);
    return;
  }
  if (io.KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_C, false)) {
    copySelection(st);
    return;
  }
  if (io.KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_V, false)) {
    pasteClipboard(st);
    return;
  }
  if (ImGui::IsKeyPressed(ImGuiKey_Escape, false)) {
    st.tool = Tool::Select;
    st.sel.clear();
    return;
  }
  if (ImGui::IsKeyPressed(ImGuiKey_Delete, false) ||
      ImGui::IsKeyPressed(ImGuiKey_Backspace, false)) {
    deleteSelection(st);
    return;
  }

  if (!io.KeyCtrl && st.hoverAtom.valid()) {
    constexpr std::array<std::pair<ImGuiKey, uint8_t>, 10> elements{{
        {ImGuiKey_C, 6}, {ImGuiKey_N, 7}, {ImGuiKey_O, 8}, {ImGuiKey_S, 16},
        {ImGuiKey_P, 15}, {ImGuiKey_F, 9}, {ImGuiKey_B, 5}, {ImGuiKey_I, 53},
        {ImGuiKey_L, 17}, {ImGuiKey_R, 35},
    }};
    for (const auto& [key, atomicNumber] : elements) {
      if (ImGui::IsKeyPressed(key, false)) {
        retypeHoveredAtom(st, atomicNumber);
        return;
      }
    }
  }
  if (st.hoverBond.valid()) {
    if (ImGui::IsKeyPressed(ImGuiKey_1, false)) {
      setHoveredBondOrder(st, core::BondOrder::Single);
    } else if (ImGui::IsKeyPressed(ImGuiKey_2, false)) {
      setHoveredBondOrder(st, core::BondOrder::Double);
    } else if (ImGui::IsKeyPressed(ImGuiKey_3, false)) {
      setHoveredBondOrder(st, core::BondOrder::Triple);
    } else if (ImGui::IsKeyPressed(ImGuiKey_W, false)) {
      toggleHoveredStereo(st, io.KeyShift ? core::BondStereo::Hash : core::BondStereo::Wedge);
    }
  }
}

}  // namespace

void hitTest(AppState& st, const CanvasRect& rect, bool canvasHovered) {
  st.hoverAtom = {};
  st.hoverBond = {};
  if (!canvasHovered) return;
  const ImVec2 mouse = ImGui::GetIO().MousePos;
  const core::Vec2 cursor = st.cam.screenToWorld({mouse.x, mouse.y}, rect.origin);

  float bestAtom = 0.30f * 0.30f;
  for (int mi = 0; mi < static_cast<int>(st.doc.molecules.size()); ++mi) {
    const core::Molecule& mol = st.doc.molecules[static_cast<size_t>(mi)];
    for (const core::Atom& atom : mol.atoms()) {
      const float distance = distanceSquared(cursor, atom.pos);
      if (distance <= bestAtom) {
        bestAtom = distance;
        st.hoverAtom = {mi, atom.id};
      }
    }
  }
  if (st.hoverAtom.valid()) return;

  float bestBond = 0.15f;
  for (int mi = 0; mi < static_cast<int>(st.doc.molecules.size()); ++mi) {
    const core::Molecule& mol = st.doc.molecules[static_cast<size_t>(mi)];
    for (const core::Bond& bond : mol.bonds()) {
      const core::Atom* a = mol.atom(bond.a);
      const core::Atom* b = mol.atom(bond.b);
      if (!a || !b) continue;
      const float distance = segmentDistance(cursor, a->pos, b->pos);
      if (distance <= bestBond) {
        bestBond = distance;
        st.hoverBond = {mi, bond.id};
      }
    }
  }
}

std::vector<core::Vec2> makeChainPreview(const AppState& st, const Runtime& rt) {
  core::Vec2 start = rt.downWorld;
  if (const core::Atom* atom = atomAt(st, rt.downAtom)) start = atom->pos;
  const core::Vec2 raw{rt.currentWorld.x - start.x, rt.currentWorld.y - start.y};
  const float distance = core::length(raw);
  if (distance < 0.3f) return {start};
  const core::Vec2 axis = core::snapAngle(raw);
  const float alongPerBond = std::cos(kPi / 6.0f) * core::kBondLength;
  const int count = std::clamp(static_cast<int>(std::floor(distance / alongPerBond + 0.5f)), 1, 100);
  const float baseAngle = core::angleOf(axis);
  std::vector<core::Vec2> points;
  points.reserve(static_cast<size_t>(count + 1));
  points.push_back(start);
  for (int i = 0; i < count; ++i) {
    const float angle = baseAngle + (i % 2 == 0 ? kPi / 6.0f : -kPi / 6.0f);
    const core::Vec2 direction = core::fromAngle(angle, core::kBondLength);
    points.push_back({points.back().x + direction.x, points.back().y + direction.y});
  }
  return points;
}

RingGeometry makeRingGeometry(const AppState& st, core::Vec2 cursor) {
  if (st.currentRing == RingKind::Naphthalene) {
    RingGeometry geometry = naphthaleneGeometry();
    if (const core::Bond* bond = bondAt(st, st.hoverBond)) {
      const core::Molecule& mol = st.doc.molecules[static_cast<size_t>(st.hoverBond.mol)];
      const core::Atom* a = mol.atom(bond->a);
      const core::Atom* b = mol.atom(bond->b);
      if (a && b) {
        Transform2D first{{-0.8660254f, 1.0f}, {-1.7320508f, 0.5f}, a->pos, b->pos, false};
        Transform2D second = first;
        second.mirror = true;
        RingGeometry candidate = geometry;
        transformNaphthalene(candidate, first);
        RingGeometry mirrored = geometry;
        transformNaphthalene(mirrored, second);
        const core::Vec2 rest = meanPosition(mol, bond->a, bond->b);
        const core::Vec2 midpoint{(a->pos.x + b->pos.x) * 0.5f, (a->pos.y + b->pos.y) * 0.5f};
        const core::Vec2 away{midpoint.x - rest.x, midpoint.y - rest.y};
        const core::Vec2 c1 = pointsMean(candidate.positions);
        const core::Vec2 c2 = pointsMean(mirrored.positions);
        geometry = ((c1.x - midpoint.x) * away.x + (c1.y - midpoint.y) * away.y >=
                    (c2.x - midpoint.x) * away.x + (c2.y - midpoint.y) * away.y)
                       ? std::move(candidate)
                       : std::move(mirrored);
        geometry.reuse[0] = {st.hoverBond.mol, bond->a};
        geometry.reuse[1] = {st.hoverBond.mol, bond->b};
        return geometry;
      }
    }
    if (const core::Atom* atom = atomAt(st, st.hoverAtom)) {
      const core::Molecule& mol = st.doc.molecules[static_cast<size_t>(st.hoverAtom.mol)];
      const core::Vec2 rest = meanPosition(mol, atom->id);
      core::Vec2 away = core::normalize({atom->pos.x - rest.x, atom->pos.y - rest.y});
      if (core::length(away) < 0.1f) away = {0.0f, 1.0f};
      const core::Vec2 baseFrom{-0.8660254f, 1.0f};
      const core::Vec2 baseCenter{};
      Transform2D transform{baseFrom, baseCenter, atom->pos,
                            {atom->pos.x + away.x * core::length(baseCenter),
                             atom->pos.y + away.y * core::length(baseCenter)},
                            false};
      // The anchor-to-centre vector defines orientation; its length must remain non-zero.
      transform.to = baseCenter;
      transform.targetTo = {atom->pos.x + away.x * core::length(
                                {baseCenter.x - baseFrom.x, baseCenter.y - baseFrom.y}),
                            atom->pos.y + away.y * core::length(
                                {baseCenter.x - baseFrom.x, baseCenter.y - baseFrom.y})};
      transformNaphthalene(geometry, transform);
      geometry.reuse[0] = st.hoverAtom;
      return geometry;
    }
    for (core::Vec2& p : geometry.positions) {
      p.x += cursor.x;
      p.y += cursor.y;
    }
    return geometry;
  }

  const int count = ringSize(st.currentRing);
  RingGeometry geometry;
  geometry.reuse.resize(static_cast<size_t>(count));
  if (const core::Bond* bond = bondAt(st, st.hoverBond)) {
    const core::Molecule& mol = st.doc.molecules[static_cast<size_t>(st.hoverBond.mol)];
    const core::Atom* a = mol.atom(bond->a);
    const core::Atom* b = mol.atom(bond->b);
    if (a && b) {
      std::vector<core::Vec2> left = fusedPolygon(a->pos, b->pos, count, 1);
      std::vector<core::Vec2> right = fusedPolygon(a->pos, b->pos, count, -1);
      const core::Vec2 rest = meanPosition(mol, bond->a, bond->b);
      const core::Vec2 midpoint{(a->pos.x + b->pos.x) * 0.5f, (a->pos.y + b->pos.y) * 0.5f};
      core::Vec2 away{midpoint.x - rest.x, midpoint.y - rest.y};
      if (mol.atomCount() <= 2 || core::length(away) < 0.05f) away = {0.0f, 1.0f};
      const core::Vec2 leftMean = pointsMean(left);
      const core::Vec2 rightMean = pointsMean(right);
      geometry.positions = ((leftMean.x - midpoint.x) * away.x +
                                (leftMean.y - midpoint.y) * away.y >=
                            (rightMean.x - midpoint.x) * away.x +
                                (rightMean.y - midpoint.y) * away.y)
                               ? std::move(left)
                               : std::move(right);
      geometry.reuse[0] = {st.hoverBond.mol, bond->a};
      geometry.reuse[1] = {st.hoverBond.mol, bond->b};
      fillSimpleEdges(geometry, st.currentRing);
      return geometry;
    }
  }
  if (const core::Atom* atom = atomAt(st, st.hoverAtom)) {
    const core::Molecule& mol = st.doc.molecules[static_cast<size_t>(st.hoverAtom.mol)];
    const core::Vec2 rest = meanPosition(mol, atom->id);
    core::Vec2 away = core::normalize({atom->pos.x - rest.x, atom->pos.y - rest.y});
    if (core::length(away) < 0.1f) away = {0.0f, 1.0f};
    const float radius = core::kBondLength / (2.0f * std::sin(kPi / static_cast<float>(count)));
    const core::Vec2 center{atom->pos.x + away.x * radius, atom->pos.y + away.y * radius};
    geometry.positions = regularRing(count, center, std::atan2(-away.y, -away.x));
    geometry.positions[0] = atom->pos;
    geometry.reuse[0] = st.hoverAtom;
  } else {
    geometry.positions = regularRing(count, cursor, kPi * 0.5f);
  }
  fillSimpleEdges(geometry, st.currentRing);
  return geometry;
}

void handleInput(AppState& st, Runtime& rt, const CanvasRect& rect, bool canvasHovered,
                 bool canvasActive) {
  ImGuiIO& io = ImGui::GetIO();
  const core::Vec2 mouseWorld = st.cam.screenToWorld({io.MousePos.x, io.MousePos.y}, rect.origin);
  rt.currentWorld = mouseWorld;
  rt.currentScreen = io.MousePos;

  if (canvasHovered && io.MouseWheel != 0.0f) {
    st.cam.zoomAt(std::pow(1.1f, io.MouseWheel), {io.MousePos.x, io.MousePos.y}, rect.origin);
  }

  const bool spacePan = ImGui::IsKeyDown(ImGuiKey_Space) && ImGui::IsMouseDown(ImGuiMouseButton_Left);
  const bool middlePan = ImGui::IsMouseDown(ImGuiMouseButton_Middle);
  if ((canvasHovered || canvasActive) && (spacePan || middlePan)) {
    const float scale = st.cam.scale();
    if (scale > 0.0f) {
      st.cam.pan.x -= io.MouseDelta.x / scale;
      st.cam.pan.y += io.MouseDelta.y / scale;
    }
  }

  handleKeyboard(st);

  if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) && canvasHovered &&
      !ImGui::IsKeyDown(ImGuiKey_Space)) {
    rt.downWorld = mouseWorld;
    rt.currentWorld = mouseWorld;
    rt.downScreen = io.MousePos;
    rt.currentScreen = io.MousePos;
    rt.downAtom = st.hoverAtom;
    rt.downBond = st.hoverBond;
    rt.dragged = false;
    rt.moveOrigins.clear();
    switch (st.tool) {
      case Tool::Bond: rt.gesture = Gesture::Bond; break;
      case Tool::Chain: rt.gesture = Gesture::Chain; break;
      case Tool::Select: startSelectionGesture(st, rt, io); break;
      case Tool::Eraser:
      case Tool::Atom:
      case Tool::ChargePlus:
      case Tool::ChargeMinus:
      case Tool::RingTemplate:
        handleImmediateTool(st, mouseWorld);
        rt.gesture = Gesture::None;
        break;
    }
  }

  if (rt.gesture != Gesture::None && ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
    const float dx = io.MousePos.x - rt.downScreen.x;
    const float dy = io.MousePos.y - rt.downScreen.y;
    if (!rt.dragged && dx * dx + dy * dy >= io.MouseDragThreshold * io.MouseDragThreshold) {
      rt.dragged = true;
      if (rt.gesture == Gesture::Move && !rt.moveOrigins.empty()) st.snapshot();
    }
    if (rt.gesture == Gesture::Move && rt.dragged) updateMoveGesture(st, rt);
  }

  if (ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
    switch (rt.gesture) {
      case Gesture::Bond: finishBondGesture(st, rt); break;
      case Gesture::Chain: finishChainGesture(st, rt); break;
      case Gesture::Marquee:
        if (rt.dragged) finishMarquee(st, rt);
        break;
      case Gesture::Move:
        if (rt.dragged && !rt.moveOrigins.empty()) {
          st.touch();
          st.hoverAtom = {};
          st.hoverBond = {};
        }
        break;
      case Gesture::None: break;
    }
    rt.gesture = Gesture::None;
    rt.moveOrigins.clear();
    rt.downAtom = {};
    rt.downBond = {};
    rt.dragged = false;
  }
}

}  // namespace chemcad::ui::canvas

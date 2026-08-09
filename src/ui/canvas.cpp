#include <algorithm>
#include <cmath>
#include <cstdio>

#include "imgui.h"

#include "ui/canvas_internal.hpp"
#include "ui/icons.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"
#include "ui/widgets.hpp"

namespace chemcad::ui {
namespace {

// One ghost button of the zoom overlay.
bool overlayButton(const char* id, ImVec2 pos, ImVec2 size, const char* tooltip) {
  ImGui::SetCursorScreenPos(pos);
  ImGui::InvisibleButton(id, size);
  const bool clicked = ImGui::IsItemClicked();
  if (ImGui::IsItemHovered() && tooltip) ImGui::SetTooltip("%s", tooltip);
  return clicked;
}

// Bottom-right canvas chrome: zoom out / zoom % / zoom in / fit. Submitted
// after the canvas's own button, so they win hover where they overlap.
void drawZoomOverlay(AppState& st, ImVec2 rectMin, ImVec2 rectMax) {
  const style::Metrics& m = style::metrics();
  const float h = ImGui::GetFontSize() * 1.7f;
  const float gap = m.gap * 0.5f;
  const float margin = m.gap * 1.2f;
  const core::Vec2 centre{(rectMin.x + rectMax.x) * 0.5f, (rectMin.y + rectMax.y) * 0.5f};
  const core::Vec2 origin{rectMin.x, rectMin.y};

  char pct[8];
  std::snprintf(pct, sizeof(pct), "%.0f%%", std::round(st.cam.zoom * 100.0f));
  const float pctW = ImGui::CalcTextSize(pct).x + m.gap * 1.2f;

  const float totalW = h * 3.0f + pctW + gap * 3.0f;
  const float y = rectMax.y - margin - h;
  float x = rectMax.x - margin - totalW;

  ImDrawList* dl = ImGui::GetWindowDrawList();
  const auto chrome = [&](ImVec2 p, ImVec2 s, bool hovered) {
    dl->AddRectFilled(p, ImVec2(p.x + s.x, p.y + s.y),
                      style::u32(style::col::BgRaised, hovered ? 0.95f : 0.72f), m.radiusSm);
    dl->AddRect(p, ImVec2(p.x + s.x, p.y + s.y), style::u32(style::col::BorderStrong, 0.8f),
                m.radiusSm, 0, m.hairline);
  };

  if (overlayButton("##zoom_out", ImVec2(x, y), ImVec2(h, h), "Zoom out"))
    st.cam.zoomAt(1.0f / 1.2f, centre, origin);
  chrome(ImVec2(x, y), ImVec2(h, h), ImGui::IsItemHovered());
  icons::draw(dl, icons::Icon::Minus, ImVec2(x + h * 0.5f, y + h * 0.5f), h * 0.34f,
              style::u32(style::col::Text));
  x += h + gap;

  if (overlayButton("##zoom_reset", ImVec2(x, y), ImVec2(pctW, h), "Reset zoom (Ctrl+0)"))
    st.cam.zoom = 1.0f;
  chrome(ImVec2(x, y), ImVec2(pctW, h), ImGui::IsItemHovered());
  {
    const ImVec2 w = ImGui::CalcTextSize(pct);
    dl->AddText(ImVec2(x + (pctW - w.x) * 0.5f, y + (h - w.y) * 0.5f), style::u32(style::col::Text),
                pct);
  }
  x += pctW + gap;

  if (overlayButton("##zoom_in", ImVec2(x, y), ImVec2(h, h), "Zoom in (mouse wheel)"))
    st.cam.zoomAt(1.2f, centre, origin);
  chrome(ImVec2(x, y), ImVec2(h, h), ImGui::IsItemHovered());
  icons::draw(dl, icons::Icon::Plus, ImVec2(x + h * 0.5f, y + h * 0.5f), h * 0.34f,
              style::u32(style::col::Text));
  x += h + gap;

  if (overlayButton("##zoom_fit", ImVec2(x, y), ImVec2(h, h), "Fit to window (Ctrl+F)"))
    st.cam.fit(st.doc, st.canvasSize);
  chrome(ImVec2(x, y), ImVec2(h, h), ImGui::IsItemHovered());
  icons::draw(dl, icons::Icon::ZoomFit, ImVec2(x + h * 0.5f, y + h * 0.5f), h * 0.34f,
              style::u32(style::col::Text));
}

constexpr float kAtomEraseRadius = 0.30f;
constexpr float kBondEraseRadius = 0.15f;

struct EraseDrag {
  bool active = false;
  bool dragged = false;
  core::Vec2 downWorld;
  core::Vec2 previousWorld;
  ImVec2 downScreen{};
  Selection pending;
};

float distanceSquared(core::Vec2 a, core::Vec2 b) {
  const float dx = a.x - b.x;
  const float dy = a.y - b.y;
  return dx * dx + dy * dy;
}

float pointSegmentDistanceSquared(core::Vec2 point, core::Vec2 a, core::Vec2 b) {
  const float dx = b.x - a.x;
  const float dy = b.y - a.y;
  const float lengthSquared = dx * dx + dy * dy;
  if (lengthSquared <= 1e-8f) return distanceSquared(point, a);
  const float t =
      std::clamp(((point.x - a.x) * dx + (point.y - a.y) * dy) / lengthSquared, 0.0f, 1.0f);
  return distanceSquared(point, {a.x + dx * t, a.y + dy * t});
}

float cross(core::Vec2 a, core::Vec2 b, core::Vec2 c) {
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

bool onSegment(core::Vec2 point, core::Vec2 a, core::Vec2 b) {
  constexpr float epsilon = 1e-6f;
  return std::abs(cross(a, b, point)) <= epsilon &&
         point.x >= std::min(a.x, b.x) - epsilon &&
         point.x <= std::max(a.x, b.x) + epsilon &&
         point.y >= std::min(a.y, b.y) - epsilon &&
         point.y <= std::max(a.y, b.y) + epsilon;
}

bool segmentsIntersect(core::Vec2 a, core::Vec2 b, core::Vec2 c, core::Vec2 d) {
  constexpr float epsilon = 1e-6f;
  const float abC = cross(a, b, c);
  const float abD = cross(a, b, d);
  const float cdA = cross(c, d, a);
  const float cdB = cross(c, d, b);
  if (((abC > epsilon && abD < -epsilon) || (abC < -epsilon && abD > epsilon)) &&
      ((cdA > epsilon && cdB < -epsilon) || (cdA < -epsilon && cdB > epsilon))) {
    return true;
  }
  return (std::abs(abC) <= epsilon && onSegment(c, a, b)) ||
         (std::abs(abD) <= epsilon && onSegment(d, a, b)) ||
         (std::abs(cdA) <= epsilon && onSegment(a, c, d)) ||
         (std::abs(cdB) <= epsilon && onSegment(b, c, d));
}

float segmentDistanceSquared(core::Vec2 a, core::Vec2 b, core::Vec2 c, core::Vec2 d) {
  if (segmentsIntersect(a, b, c, d)) return 0.0f;
  return std::min({pointSegmentDistanceSquared(a, c, d),
                   pointSegmentDistanceSquared(b, c, d),
                   pointSegmentDistanceSquared(c, a, b),
                   pointSegmentDistanceSquared(d, a, b)});
}

void addPendingBond(EraseDrag& drag, BondRef ref) {
  if (ref.valid() && !drag.pending.contains(ref)) drag.pending.bonds.push_back(ref);
}

void addPendingAtom(const AppState& st, EraseDrag& drag, AtomRef ref) {
  if (!ref.valid() || drag.pending.contains(ref) ||
      ref.mol >= static_cast<int>(st.doc.molecules.size())) {
    return;
  }
  const core::Molecule& mol = st.doc.molecules[static_cast<size_t>(ref.mol)];
  if (!mol.atom(ref.id)) return;
  drag.pending.atoms.push_back(ref);
  // Atom erasure removes incident bonds in Molecule::removeAtom. Include them
  // in the preview so the highlighted set is exactly what release removes.
  for (core::BondId id : mol.incidentBonds(ref.id)) addPendingBond(drag, {ref.mol, id});
}

void collectClickTarget(const AppState& st, EraseDrag& drag) {
  // hitTest gives atoms precedence over bonds; preserving that ordering keeps a
  // zero-distance gesture identical to the original click eraser.
  if (st.hoverAtom.valid()) {
    addPendingAtom(st, drag, st.hoverAtom);
  } else {
    addPendingBond(drag, st.hoverBond);
  }
}

void collectEraseSweep(const AppState& st, EraseDrag& drag, core::Vec2 from, core::Vec2 to) {
  const float atomRadiusSquared = kAtomEraseRadius * kAtomEraseRadius;
  const float bondRadiusSquared = kBondEraseRadius * kBondEraseRadius;
  for (int mi = 0; mi < static_cast<int>(st.doc.molecules.size()); ++mi) {
    const core::Molecule& mol = st.doc.molecules[static_cast<size_t>(mi)];
    for (const core::Atom& atom : mol.atoms()) {
      if (pointSegmentDistanceSquared(atom.pos, from, to) <= atomRadiusSquared) {
        addPendingAtom(st, drag, {mi, atom.id});
      }
    }
    for (const core::Bond& bond : mol.bonds()) {
      const core::Atom* a = mol.atom(bond.a);
      const core::Atom* b = mol.atom(bond.b);
      if (a && b && segmentDistanceSquared(from, to, a->pos, b->pos) <= bondRadiusSquared) {
        addPendingBond(drag, {mi, bond.id});
      }
    }
  }
}

void clearEraseDrag(EraseDrag& drag) {
  drag.active = false;
  drag.dragged = false;
  drag.pending.clear();
}

void commitErase(AppState& st, EraseDrag& drag) {
  if (drag.pending.empty()) {
    clearEraseDrag(drag);
    return;
  }

  // One snapshot for the whole sweep. Bonds incident to a pending atom stay
  // in the preview but are deleted by Molecule::removeAtom, the same primitive
  // the click eraser has always used.
  st.snapshot();
  for (int mi = static_cast<int>(st.doc.molecules.size()) - 1; mi >= 0; --mi) {
    core::Molecule& mol = st.doc.molecules[static_cast<size_t>(mi)];
    for (const BondRef& ref : drag.pending.bonds) {
      if (ref.mol != mi) continue;
      const core::Bond* bond = mol.bond(ref.id);
      if (!bond) continue;
      if (drag.pending.contains(AtomRef{mi, bond->a}) ||
          drag.pending.contains(AtomRef{mi, bond->b})) {
        continue;
      }
      mol.removeBond(ref.id);
    }
    for (const AtomRef& ref : drag.pending.atoms) {
      if (ref.mol == mi) mol.removeAtom(ref.id);
    }
  }
  st.doc.molecules.erase(
      std::remove_if(st.doc.molecules.begin(), st.doc.molecules.end(),
                     [](const core::Molecule& mol) { return mol.empty(); }),
      st.doc.molecules.end());
  st.sel.clear();
  st.hoverAtom = {};
  st.hoverBond = {};
  st.touch();
  clearEraseDrag(drag);
}

void handleEraseInput(AppState& st, EraseDrag& drag, const canvas::CanvasRect& rect,
                      bool canvasHovered) {
  ImGuiIO& io = ImGui::GetIO();
  if (drag.active && st.tool != Tool::Eraser) {
    clearEraseDrag(drag);
    return;
  }

  if (drag.active &&
      (ImGui::IsKeyPressed(ImGuiKey_Escape, false) ||
       ImGui::IsMouseClicked(ImGuiMouseButton_Right))) {
    clearEraseDrag(drag);
    if (ImGui::IsKeyPressed(ImGuiKey_Escape, false)) {
      st.tool = Tool::Select;
      st.sel.clear();
    }
    return;
  }

  if (!drag.active && st.tool == Tool::Eraser && canvasHovered &&
      ImGui::IsMouseClicked(ImGuiMouseButton_Left) &&
      !ImGui::IsKeyDown(ImGuiKey_Space)) {
    drag.active = true;
    drag.dragged = false;
    drag.downScreen = io.MousePos;
    drag.downWorld = st.cam.screenToWorld({io.MousePos.x, io.MousePos.y}, rect.origin);
    drag.previousWorld = drag.downWorld;
    drag.pending.clear();
    collectClickTarget(st, drag);
  }
  if (!drag.active) return;

  const core::Vec2 current =
      st.cam.screenToWorld({io.MousePos.x, io.MousePos.y}, rect.origin);
  const float dx = io.MousePos.x - drag.downScreen.x;
  const float dy = io.MousePos.y - drag.downScreen.y;
  if (!drag.dragged &&
      dx * dx + dy * dy >= io.MouseDragThreshold * io.MouseDragThreshold) {
    drag.dragged = true;
  }
  if (drag.dragged &&
      (ImGui::IsMouseDown(ImGuiMouseButton_Left) ||
       ImGui::IsMouseReleased(ImGuiMouseButton_Left))) {
    // Sweep the whole frame-to-frame segment so low frame rates cannot tunnel
    // through a small atom or bond hit region.
    collectEraseSweep(st, drag, drag.previousWorld, current);
    drag.previousWorld = current;
  }

  if (ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
    commitErase(st, drag);
  } else if (!ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
    // Covers focus loss as well as an off-canvas release missed by the item.
    clearEraseDrag(drag);
  }
}

void drawEraseOverlay(const AppState& st, const EraseDrag& drag,
                      const canvas::CanvasRect& rect, bool canvasHovered, float animation) {
  ImDrawList* draw = ImGui::GetWindowDrawList();
  const float highlightAlpha = 0.62f + animation * 0.28f;
  const ImU32 danger = style::u32(style::col::Danger, highlightAlpha);
  const float baseThickness = std::max(1.2f, 2.0f * st.cam.zoom);

  for (const BondRef& ref : drag.pending.bonds) {
    if (ref.mol < 0 || ref.mol >= static_cast<int>(st.doc.molecules.size())) continue;
    const core::Molecule& mol = st.doc.molecules[static_cast<size_t>(ref.mol)];
    const core::Bond* bond = mol.bond(ref.id);
    if (!bond) continue;
    const core::Atom* a = mol.atom(bond->a);
    const core::Atom* b = mol.atom(bond->b);
    if (!a || !b) continue;
    const core::Vec2 screenA = st.cam.worldToScreen(a->pos, rect.origin);
    const core::Vec2 screenB = st.cam.worldToScreen(b->pos, rect.origin);
    draw->AddLine({screenA.x, screenA.y}, {screenB.x, screenB.y}, danger,
                  baseThickness + 1.5f + animation * 1.5f);
  }

  for (const AtomRef& ref : drag.pending.atoms) {
    if (ref.mol < 0 || ref.mol >= static_cast<int>(st.doc.molecules.size())) continue;
    const core::Molecule& mol = st.doc.molecules[static_cast<size_t>(ref.mol)];
    const core::Atom* atom = mol.atom(ref.id);
    if (!atom) continue;
    const core::Vec2 screen = st.cam.worldToScreen(atom->pos, rect.origin);
    const bool labelled = atom->atomicNumber != 6 || atom->charge != 0 || atom->isotope != 0 ||
                          mol.degree(atom->id) == 0;
    const float fontSize =
        std::clamp(ImGui::GetFontSize() * std::sqrt(st.cam.zoom), 12.0f, 32.0f);
    const float radius = labelled ? fontSize * 0.95f : 8.0f;
    draw->AddCircle({screen.x, screen.y}, radius, danger, 24, 1.4f + animation * 0.8f);
  }

  if (st.tool == Tool::Eraser && (canvasHovered || drag.active)) {
    const ImVec2 mouse = ImGui::GetIO().MousePos;
    const float cursorRadius =
        std::max(style::metrics().iconSize * 0.45f, kAtomEraseRadius * st.cam.scale());
    draw->AddCircle(mouse, cursorRadius, style::u32(style::col::Danger, 0.48f), 28,
                    style::metrics().hairline);
    icons::draw(draw, icons::Icon::Eraser, mouse, style::metrics().iconSize * 0.58f,
                style::u32(style::col::Danger, 0.88f));
  }
}

}  // namespace

void drawCanvas(AppState& st) {
  static canvas::Runtime runtime;
  static EraseDrag eraseDrag;

  ImVec2 size = ImGui::GetContentRegionAvail();
  size.x = std::max(size.x, 1.0f);
  size.y = std::max(size.y, 1.0f);
  const ImVec2 origin = ImGui::GetCursorScreenPos();
  ImGui::InvisibleButton("##canvas", size,
                         ImGuiButtonFlags_MouseButtonLeft |
                             ImGuiButtonFlags_MouseButtonRight |
                             ImGuiButtonFlags_MouseButtonMiddle);
  const ImVec2 rectMin = ImGui::GetItemRectMin();
  const ImVec2 rectMax = ImGui::GetItemRectMax();
  const bool hovered = ImGui::IsItemHovered();
  const bool active = ImGui::IsItemActive();

  st.canvasOrigin = {rectMin.x, rectMin.y};
  st.canvasSize = {size.x, size.y};

  // Content appearing in a previously empty sketch (build-from-name, MOL
  // import, project open) lands wherever the source put it; frame it. An
  // intentional zoom/pan on a non-empty doc is never second-guessed.
  static bool wasEmpty = true;
  const bool isEmpty = st.doc.empty();
  if (wasEmpty && !isEmpty) st.cam.fit(st.doc, st.canvasSize);
  wasEmpty = isEmpty;

  const canvas::CanvasRect rect{{rectMin.x, rectMin.y},
                                {size.x, size.y},
                                rectMin,
                                rectMax};

  canvas::hitTest(st, rect, hovered);
  const bool eraserPress =
      st.tool == Tool::Eraser && hovered &&
      ImGui::IsMouseClicked(ImGuiMouseButton_Left) &&
      !ImGui::IsKeyDown(ImGuiKey_Space);
  // The legacy input handler commits Eraser on mouse-down. Own eraser frames
  // here so both a click and a sweep use the same release-time transaction.
  if (!eraseDrag.active && !eraserPress) {
    canvas::handleInput(st, runtime, rect, hovered, active);
  }
  handleEraseInput(st, eraseDrag, rect, hovered);
  canvas::hitTest(st, rect, hovered);
  if (eraseDrag.active) {
    st.hoverAtom = {};
    st.hoverBond = {};
  }
  const float eraseHighlight =
      widgets::hoverT(ImGui::GetID("##erase_highlight"), !eraseDrag.pending.empty());

  ImDrawList* draw = ImGui::GetWindowDrawList();
  draw->PushClipRect(rectMin, rectMax, true);
  canvas::render(st, runtime, rect);
  drawEraseOverlay(st, eraseDrag, rect, hovered, eraseHighlight);
  draw->PopClipRect();

  drawZoomOverlay(st, rectMin, rectMax);
}

}  // namespace chemcad::ui

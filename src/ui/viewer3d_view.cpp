// 3D molecule viewer: a turntable ball-and-stick / licorice / space-filling
// renderer drawn entirely through ImDrawList (weak perspective + painter's
// algorithm, no GL resources), re-embedded from the sketch on every
// docRevision change. Renders inside an already-open window; never owns one.

#include "imgui.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "chem/bridge.hpp"
#include "chem/embed3d.hpp"
#include "core/model.hpp"
#include "ui/app_state.hpp"
#include "ui/theme.hpp"
#include "ui/viewer3d_state.hpp"
#include "ui/widgets.hpp"

namespace chemcad::ui {
namespace {

constexpr float kPi = 3.14159265f;

// ------------------------------------------------------------ element table
ImVec4 elementColor(uint8_t z) {
  switch (z) {
    case 1:  return {0.90f, 0.90f, 0.92f, 1.0f};  // H
    case 6:  return {0.52f, 0.52f, 0.56f, 1.0f};  // C
    case 7:  return {0.24f, 0.46f, 0.95f, 1.0f};  // N
    case 8:  return {0.88f, 0.22f, 0.22f, 1.0f};  // O
    case 9:  return {0.55f, 0.85f, 0.38f, 1.0f};  // F
    case 15: return {0.95f, 0.58f, 0.16f, 1.0f};  // P
    case 16: return {0.92f, 0.82f, 0.26f, 1.0f};  // S
    case 17: return {0.36f, 0.80f, 0.38f, 1.0f};  // Cl
    case 35: return {0.68f, 0.26f, 0.18f, 1.0f};  // Br
    case 53: return {0.66f, 0.32f, 0.78f, 1.0f};  // I
    default: return {0.78f, 0.50f, 0.88f, 1.0f};  // anything exotic
  }
}

// Covalent radii in Angstrom; ball-and-stick draws a fraction of these.
float covalentRadius(uint8_t z) {
  switch (z) {
    case 1:  return 0.31f;
    case 6:  return 0.76f;
    case 7:  return 0.71f;
    case 8:  return 0.66f;
    case 9:  return 0.57f;
    case 15: return 1.07f;
    case 16: return 1.05f;
    case 17: return 1.02f;
    case 35: return 1.20f;
    case 53: return 1.39f;
    default: return 1.00f;
  }
}

float vdwRadius(uint8_t z) {
  switch (z) {
    case 1:  return 1.20f;
    case 6:  return 1.70f;
    case 7:  return 1.55f;
    case 8:  return 1.52f;
    case 9:  return 1.47f;
    case 15: return 1.80f;
    case 16: return 1.80f;
    case 17: return 1.75f;
    case 35: return 1.85f;
    case 53: return 1.98f;
    default: return 1.80f;
  }
}

// ------------------------------------------------------------- model sync
void syncModel(AppState& st) {
  Viewer3DState& vs = st.viewer3d;
  if (vs.sourceRevision == st.docRevision && (vs.hasModel || !vs.errorMessage.empty())) return;
  vs.sourceRevision = st.docRevision;
  vs.hasModel = false;
  vs.errorMessage.clear();
  try {
    const core::Molecule* found = nullptr;
    for (const core::Molecule& mol : st.doc.molecules) {
      if (!mol.empty()) {
        found = &mol;
        break;
      }
    }
    if (!found) {
      vs.errorMessage = "Sketch is empty.";
      return;
    }
    vs.model = chem::embed3D(*found);
    const chem::Properties props = chem::computeProperties(*found);
    vs.formula = props.formula;
    vs.molWeight = props.mw;
    vs.hasModel = true;
  } catch (const std::exception& err) {
    vs.errorMessage = err.what();
  }
}

// ------------------------------------------------------------- projection
struct Projected {
  float x = 0.0f;     // screen-space
  float y = 0.0f;
  float depth = 0.0f; // rotated z: larger = further from camera
  float persp = 1.0f; // weak-perspective scale factor applied
};

// Turntable rotation (yaw about Y, then pitch about X) followed by a weak
// perspective divide. Depth cueing does the rest of the 3D read.
Projected project(const chem::Atom3D& a, float yawRad, float pitchRad, ImVec2 centre,
                  float scale, float perspK) {
  const float cy = std::cos(yawRad), sy = std::sin(yawRad);
  const float cp = std::cos(pitchRad), sp = std::sin(pitchRad);
  const float x1 = a.x * cy + a.z * sy;
  const float z1 = -a.x * sy + a.z * cy;
  const float y2 = a.y * cp - z1 * sp;
  const float z2 = a.y * sp + z1 * cp;
  const float persp = 1.0f / (1.0f + z2 * perspK);
  return Projected{centre.x + x1 * scale * persp, centre.y - y2 * scale * persp, z2, persp};
}

ImVec4 shadeAtom(ImVec4 base, float fog) {
  // Distant atoms sink toward the background so depth reads without GL.
  const ImVec4 bg = style::col::BgDeep;
  return ImVec4(base.x + (bg.x - base.x) * fog, base.y + (bg.y - base.y) * fog,
                base.z + (bg.z - base.z) * fog, 1.0f);
}

void drawSphere(ImDrawList* dl, ImVec2 p, float r, ImVec4 color) {
  dl->AddCircleFilled(p, r, style::u32(color), 24);
  // Fake Lambert shading: darker rim at the bottom-right, specular dot at
  // the top-left. Three draw calls read as a lit sphere at molecule sizes.
  const ImVec4 dark(color.x * 0.55f, color.y * 0.55f, color.z * 0.55f, 1.0f);
  dl->AddCircle(p, r, style::u32(dark, 0.75f), 24, std::max(1.0f, r * 0.10f));
  const ImVec4 hi(color.x + (1.0f - color.x) * 0.55f, color.y + (1.0f - color.y) * 0.55f,
                  color.z + (1.0f - color.z) * 0.55f, 1.0f);
  dl->AddCircleFilled(ImVec2(p.x - r * 0.34f, p.y - r * 0.38f), r * 0.34f, style::u32(hi, 0.85f),
                      16);
}

void drawBondSegment(ImDrawList* dl, ImVec2 from, ImVec2 to, float thick, ImVec4 color) {
  dl->AddLine(from, to, style::u32(color), std::max(1.0f, thick));
}

// Parallel offset lines for multiple bonds. count includes the central one.
void drawBond(ImDrawList* dl, ImVec2 pa, ImVec2 pb, float thick, ImVec4 colorA, ImVec4 colorB,
              int order) {
  const ImVec2 mid((pa.x + pb.x) * 0.5f, (pa.y + pb.y) * 0.5f);
  const auto half = [&](ImVec2 a, ImVec2 b, ImVec4 c, float ox, float oy) {
    drawBondSegment(dl, ImVec2(a.x + ox, a.y + oy), ImVec2(b.x + ox, b.y + oy), thick, c);
  };
  float px = 0.0f, py = 0.0f;
  const float dx = pb.x - pa.x, dy = pb.y - pa.y;
  const float len = std::sqrt(dx * dx + dy * dy);
  if (len > 1e-3f) {
    px = -dy / len;
    py = dx / len;
  }
  const float off = thick * 1.15f;
  if (order == 2) {
    half(pa, mid, colorA, px * off * 0.5f, py * off * 0.5f);
    half(mid, pb, colorB, px * off * 0.5f, py * off * 0.5f);
    half(pa, mid, colorA, -px * off * 0.5f, -py * off * 0.5f);
    half(mid, pb, colorB, -px * off * 0.5f, -py * off * 0.5f);
  } else if (order == 3) {
    half(pa, mid, colorA, 0.0f, 0.0f);
    half(mid, pb, colorB, 0.0f, 0.0f);
    half(pa, mid, colorA, px * off, py * off);
    half(mid, pb, colorB, px * off, py * off);
    half(pa, mid, colorA, -px * off, -py * off);
    half(mid, pb, colorB, -px * off, -py * off);
  } else if (order == 4) {
    // Aromatic: solid line plus a thinner offset inner line.
    half(pa, mid, colorA, 0.0f, 0.0f);
    half(mid, pb, colorB, 0.0f, 0.0f);
    drawBondSegment(dl, ImVec2(pa.x + px * off * 0.6f, pa.y + py * off * 0.6f),
                    ImVec2(mid.x + px * off * 0.6f, mid.y + py * off * 0.6f), thick * 0.45f,
                    colorA);
    drawBondSegment(dl, ImVec2(mid.x + px * off * 0.6f, mid.y + py * off * 0.6f),
                    ImVec2(pb.x + px * off * 0.6f, pb.y + py * off * 0.6f), thick * 0.45f,
                    colorB);
  } else {
    half(pa, mid, colorA, 0.0f, 0.0f);
    half(mid, pb, colorB, 0.0f, 0.0f);
  }
}

// ---------------------------------------------------------------- drawing
void drawViewerCanvas(AppState& st, ImVec2 min, ImVec2 max, bool compact = false) {
  Viewer3DState& vs = st.viewer3d;
  ImDrawList* dl = ImGui::GetWindowDrawList();
  const style::Metrics& m = style::metrics();

  // Stage: inset panel with a soft elliptical ground shadow under the model.
  dl->AddRectFilled(min, max, style::u32(style::col::BgDeep), m.radiusMd);
  const ImVec2 centre((min.x + max.x) * 0.5f, (min.y + max.y) * 0.5f);

  if (!vs.hasModel) {
    const char* message = vs.errorMessage.empty() ? "Draw a molecule in the Sketch tab."
                                                  : vs.errorMessage.c_str();
    const ImVec2 textSize = ImGui::CalcTextSize(message);
    dl->AddText(ImVec2(centre.x - textSize.x * 0.5f, centre.y - textSize.y * 0.5f),
                style::u32(style::col::TextFaint), message);
    return;
  }

  const chem::Embedded3D& model = vs.model;
  const float radius = model.radius;
  const float fit = std::min(max.x - min.x, max.y - min.y) * 0.42f / radius;
  const float scale = fit * vs.zoom;
  const float perspK = 0.18f / radius;
  const float yaw = vs.yawDeg * kPi / 180.0f;
  const float pitch = vs.pitchDeg * kPi / 180.0f;

  std::vector<Projected> pts(model.atoms.size());
  for (size_t i = 0; i < model.atoms.size(); ++i)
    pts[i] = project(model.atoms[i], yaw, pitch, centre, scale, perspK);

  // Ground shadow: an ellipse under the model that widens as it spins.
  const float shadowW = radius * scale * 0.9f;
  const float shadowH = shadowW * 0.22f;
  const float shadowY = centre.y + radius * scale * 0.75f;
  for (int i = 0; i < 6; ++i) {
    const float t = 1.0f - static_cast<float>(i) / 6.0f;
    dl->AddEllipseFilled(ImVec2(centre.x, shadowY), ImVec2(shadowW * t, shadowH * t),
                         style::u32(style::col::Text, 0.016f + 0.012f * t));
  }

  // Painter's algorithm: bonds first (they sit behind atom balls), then
  // atoms, both sorted far-to-near.
  std::vector<const chem::Bond3D*> bonds;
  bonds.reserve(model.bonds.size());
  for (const chem::Bond3D& b : model.bonds) bonds.push_back(&b);
  std::sort(bonds.begin(), bonds.end(), [&](const chem::Bond3D* x, const chem::Bond3D* y) {
    return (pts[static_cast<size_t>(x->a)].depth + pts[static_cast<size_t>(x->b)].depth) >
           (pts[static_cast<size_t>(y->a)].depth + pts[static_cast<size_t>(y->b)].depth);
  });

  const bool spacefill = vs.style == 2;
  const bool licorice = vs.style == 1;
  const float bondThick =
      spacefill ? 0.0f
                : (licorice ? 0.20f : 0.13f) * scale * 2.0f *
                      (pts.empty() ? 1.0f : 1.0f);

  if (!spacefill) {
    for (const chem::Bond3D* b : bonds) {
      const Projected& pa = pts[static_cast<size_t>(b->a)];
      const Projected& pb = pts[static_cast<size_t>(b->b)];
      const float fogA = std::clamp((pa.depth / radius + 1.0f) * 0.5f * 0.4f, 0.0f, 0.45f);
      const float fogB = std::clamp((pb.depth / radius + 1.0f) * 0.5f * 0.4f, 0.0f, 0.45f);
      const float thick = std::max(1.0f, bondThick * (pa.persp + pb.persp) * 0.5f);
      drawBond(dl, ImVec2(pa.x, pa.y), ImVec2(pb.x, pb.y), thick,
               shadeAtom(elementColor(model.atoms[static_cast<size_t>(b->a)].atomicNumber), fogA),
               shadeAtom(elementColor(model.atoms[static_cast<size_t>(b->b)].atomicNumber), fogB),
               b->order);
    }
  }

  std::vector<size_t> order(model.atoms.size());
  for (size_t i = 0; i < order.size(); ++i) order[i] = i;
  std::sort(order.begin(), order.end(),
            [&](size_t x, size_t y) { return pts[x].depth > pts[y].depth; });

  const chem::Atom3D* hoveredAtom = nullptr;
  float hoveredDist = 1e9f;
  for (size_t idx : order) {
    const chem::Atom3D& atom = model.atoms[idx];
    const Projected& p = pts[idx];
    const float baseRadius = spacefill ? vdwRadius(atom.atomicNumber)
                                       : covalentRadius(atom.atomicNumber) *
                                             (licorice ? 0.22f : 0.42f);
    const float r = std::max(2.0f, baseRadius * scale * p.persp);
    const float fog = std::clamp((p.depth / radius + 1.0f) * 0.5f * 0.4f, 0.0f, 0.45f);
    drawSphere(dl, ImVec2(p.x, p.y), r, shadeAtom(elementColor(atom.atomicNumber), fog));

    const ImVec2 mouse = ImGui::GetMousePos();
    const float dist = std::sqrt((mouse.x - p.x) * (mouse.x - p.x) +
                                 (mouse.y - p.y) * (mouse.y - p.y));
    if (dist < r + 4.0f && dist < hoveredDist) {
      hoveredDist = dist;
      hoveredAtom = &atom;
    }
  }

  if (!compact && hoveredAtom && ImGui::IsItemHovered()) {
    ImGui::SetTooltip("%s", chem::symbolFor(hoveredAtom->atomicNumber));
  }

  // Caption: identity strip, bottom-left.
  char caption[96];
  std::snprintf(caption, sizeof(caption), "%s   %.2f g/mol   %zu atoms", vs.formula.c_str(),
                vs.molWeight, model.atoms.size());
  dl->AddText(style::fonts::mono(), ImGui::GetFontSize() * 0.9f,
              ImVec2(min.x + m.gap, max.y - m.gap - ImGui::GetFontSize()),
              style::u32(style::col::TextDim), caption);

  // Interaction hint, bottom-right (the corner overlay has no room for it).
  if (!compact) {
    const char* hint = "drag to rotate  |  wheel to zoom  |  double-click resets";
    const ImVec2 hintSize = ImGui::CalcTextSize(hint);
    dl->AddText(ImVec2(max.x - m.gap - hintSize.x, max.y - m.gap - hintSize.y),
                style::u32(style::col::TextFaint), hint);
  }
}

}  // namespace

void drawViewer3D(AppState& st) {
  Viewer3DState& vs = st.viewer3d;
  syncModel(st);

  if (vs.autoRotate) {
    vs.yawDeg = std::fmod(vs.yawDeg + ImGui::GetIO().DeltaTime * 25.0f, 360.0f);
  }

  // Control row: packs on one line when the panel is wide, stacks when the
  // docked preview column narrows.
  static const char* kStyles[] = {"Ball and stick", "Licorice", "Space-filling"};
  const float rowW = 170.0f + ImGui::CalcTextSize("Auto-rotate").x + 30.0f + 90.0f;
  const bool oneRow = ImGui::GetContentRegionAvail().x >= rowW;
  ImGui::SetNextItemWidth(std::min(170.0f, ImGui::GetContentRegionAvail().x));
  ImGui::Combo("##v3d_style", &vs.style, kStyles, 3);
  if (oneRow) ImGui::SameLine();
  ImGui::Checkbox("Auto-rotate", &vs.autoRotate);
  if (oneRow) ImGui::SameLine();
  if (widgets::ghostButton("Reset view")) {
    vs.yawDeg = 35.0f;
    vs.pitchDeg = -18.0f;
    vs.zoom = 1.0f;
    vs.autoRotate = false;
  }

  ImVec2 size = ImGui::GetContentRegionAvail();
  size.x = std::max(size.x, 60.0f);
  size.y = std::max(size.y, 200.0f);
  ImGui::InvisibleButton("##v3d_canvas", size);
  const ImVec2 min = ImGui::GetItemRectMin();
  const ImVec2 max = ImGui::GetItemRectMax();

  if (ImGui::IsItemActive() && ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
    const ImVec2 delta = ImGui::GetIO().MouseDelta;
    vs.yawDeg = std::fmod(vs.yawDeg + delta.x * 0.45f + 360.0f, 360.0f);
    vs.pitchDeg = std::clamp(vs.pitchDeg + delta.y * 0.45f, -89.0f, 89.0f);
    vs.autoRotate = false;  // grabbing the model always wins over the spinner
  }
  if (ImGui::IsItemHovered()) {
    const float wheel = ImGui::GetIO().MouseWheel;
    if (wheel != 0.0f) vs.zoom = std::clamp(vs.zoom * (1.0f + wheel * 0.12f), 0.25f, 6.0f);
    if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
      vs.yawDeg = 35.0f;
      vs.pitchDeg = -18.0f;
      vs.zoom = 1.0f;
    }
  }

  drawViewerCanvas(st, min, max);
}

}  // namespace chemcad::ui

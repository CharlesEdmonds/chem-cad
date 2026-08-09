#include "ui/icons.hpp"

#include <algorithm>
#include <cmath>

namespace chemcad::ui::icons {
namespace {

// All glyphs are defined in a normalized grid: h is half the requested size,
// coordinates are fractions of h, so (±0.4, ±0.4) is a comfortable inset.

ImVec2 at(ImVec2 c, float h, float x, float y) {
  return ImVec2(c.x + x * h, c.y + y * h);
}

void stroke(ImDrawList* dl, ImVec2 a, ImVec2 b, ImU32 col, float t) {
  dl->AddLine(a, b, col, t);
}

void dashed(ImDrawList* dl, ImVec2 a, ImVec2 b, ImU32 col, float t, float dash,
            float gap) {
  const float dx = b.x - a.x;
  const float dy = b.y - a.y;
  const float length = std::sqrt(dx * dx + dy * dy);
  if (length < 0.01f) return;
  const float ux = dx / length;
  const float uy = dy / length;
  for (float s = 0.0f; s < length; s += dash + gap) {
    const float e = std::min(s + dash, length);
    dl->AddLine(ImVec2(a.x + ux * s, a.y + uy * s), ImVec2(a.x + ux * e, a.y + uy * e),
                col, t);
  }
}

void polyline(ImDrawList* dl, ImVec2 c, float h, const ImVec2* pts, int count,
              ImU32 col, float t, bool closed) {
  ImVec2 scaled[8];
  for (int i = 0; i < count; ++i) scaled[i] = at(c, h, pts[i].x, pts[i].y);
  dl->AddPolyline(scaled, count, col, closed ? ImDrawFlags_Closed : 0, t);
}

// Regular polygon with a vertex straight up, matching the orientation
// ChemDraw uses for ring templates.
void polygonVerts(ImVec2 c, float radius, int sides, ImVec2* out) {
  for (int i = 0; i < sides; ++i) {
    const float angle =
        (-90.0f + 360.0f / static_cast<float>(sides) * static_cast<float>(i)) *
        3.14159265f / 180.0f;
    out[i] = ImVec2(c.x + radius * std::cos(angle), c.y + radius * std::sin(angle));
  }
}

void regularPolygon(ImDrawList* dl, ImVec2 c, float radius, int sides, ImU32 col, float t) {
  ImVec2 pts[8];
  polygonVerts(c, radius, sides, pts);
  dl->AddPolyline(pts, sides, col, ImDrawFlags_Closed, t);
}

void hexagon(ImDrawList* dl, ImVec2 c, float radius, ImU32 col, float t) {
  regularPolygon(dl, c, radius, 6, col, t);
}

// Short stroke parallel to edge (a,b), pulled towards the ring centre: the
// inner line of a double bond in a ring glyph.
void innerBondLine(ImDrawList* dl, ImVec2 c, ImVec2 a, ImVec2 b, float offset, ImU32 col,
                   float t) {
  const ImVec2 mid((a.x + b.x) * 0.5f, (a.y + b.y) * 0.5f);
  float nx = c.x - mid.x;
  float ny = c.y - mid.y;
  const float len = std::sqrt(nx * nx + ny * ny);
  if (len < 1e-4f) return;
  nx = nx / len * offset;
  ny = ny / len * offset;
  const float sx = (b.x - a.x) * 0.5f * 0.64f;
  const float sy = (b.y - a.y) * 0.5f * 0.64f;
  dl->AddLine(ImVec2(mid.x + nx - sx, mid.y + ny - sy), ImVec2(mid.x + nx + sx, mid.y + ny + sy),
              col, t);
}

void arc(ImDrawList* dl, ImVec2 c, float radius, float start, float end, ImU32 col,
         float t, int segments = 16) {
  dl->PathArcTo(c, radius, start, end, segments);
  dl->PathStroke(col, 0, t);
}

void arrow(ImDrawList* dl, ImVec2 c, float h, float x0, float y0, float x1, float y1,
           ImU32 col, float t, float head = 0.18f) {
  const ImVec2 a = at(c, h, x0, y0);
  const ImVec2 b = at(c, h, x1, y1);
  const float dx = b.x - a.x;
  const float dy = b.y - a.y;
  const float length = std::sqrt(dx * dx + dy * dy);
  if (length < 0.01f) return;
  const float ux = dx / length;
  const float uy = dy / length;
  const float wing = head * h;
  stroke(dl, a, b, col, t);
  stroke(dl, b, ImVec2(b.x - ux * wing - uy * wing * 0.72f,
                       b.y - uy * wing + ux * wing * 0.72f),
         col, t);
  stroke(dl, b, ImVec2(b.x - ux * wing + uy * wing * 0.72f,
                       b.y - uy * wing - ux * wing * 0.72f),
         col, t);
}

void fourPointStar(ImDrawList* dl, ImVec2 c, float h, float scale, ImU32 col,
                   float t) {
  const ImVec2 pts[] = {{0.0f, -0.42f}, {0.09f, -0.09f}, {0.42f, 0.0f},
                        {0.09f, 0.09f}, {0.0f, 0.42f},   {-0.09f, 0.09f},
                        {-0.42f, 0.0f}, {-0.09f, -0.09f}};
  polyline(dl, c, h * scale, pts, 8, col, t, true);
}

}  // namespace

void draw(ImDrawList* dl, Icon icon, ImVec2 centre, float size, ImU32 color,
          float thickness) {
  if (!dl || size <= 0.0f || icon == Icon::None) return;
  // Glyph coordinates are normalized around ±0.4. Using `size` as the unit
  // lets the drawing occupy ~80% of the requested box instead of only ~40%.
  const float h = size;
  const float t = thickness > 0.0f ? thickness : std::max(1.4f, size * 0.078f);

  switch (icon) {
    case Icon::None:
      return;
    case Icon::Select: {
      const ImVec2 pts[] = {{-0.26f, -0.42f}, {-0.26f, 0.24f}, {-0.07f, 0.11f},
                            {0.07f, 0.40f},   {0.20f, 0.32f},  {0.07f, 0.05f},
                            {0.28f, 0.05f}};
      polyline(dl, centre, h, pts, 7, color, t * 0.9f, true);
      return;
    }
    case Icon::Eraser: {
      const ImVec2 body[] = {{-0.36f, -0.04f}, {0.00f, -0.40f},
                             {0.38f, -0.02f},  {0.02f, 0.36f}};
      polyline(dl, centre, h, body, 4, color, t, true);
      stroke(dl, at(centre, h, -0.06f, 0.28f), at(centre, h, 0.20f, 0.02f), color,
             t * 0.8f);
      stroke(dl, at(centre, h, 0.08f, 0.44f), at(centre, h, 0.44f, 0.44f), color, t);
      return;
    }
    case Icon::Bond:
      stroke(dl, at(centre, h, -0.36f, 0.28f), at(centre, h, 0.36f, -0.28f), color,
             t * 1.15f);
      return;
    case Icon::Chain: {
      const ImVec2 pts[] = {{-0.42f, 0.18f}, {-0.14f, -0.20f},
                            {0.14f, 0.18f},  {0.42f, -0.20f}};
      polyline(dl, centre, h, pts, 4, color, t * 1.1f, false);
      return;
    }
    case Icon::RingCyclopropane:
      regularPolygon(dl, centre, h * 0.42f, 3, color, t);
      return;
    case Icon::RingCyclobutane:
      regularPolygon(dl, centre, h * 0.42f, 4, color, t);
      return;
    case Icon::RingCyclopentane:
      regularPolygon(dl, centre, h * 0.42f, 5, color, t);
      return;
    case Icon::RingCyclohexane:
      regularPolygon(dl, centre, h * 0.42f, 6, color, t);
      return;
    case Icon::RingCycloheptane:
      regularPolygon(dl, centre, h * 0.42f, 7, color, t);
      return;
    case Icon::RingCyclooctane:
      regularPolygon(dl, centre, h * 0.42f, 8, color, t);
      return;
    case Icon::RingBenzene:
      regularPolygon(dl, centre, h * 0.42f, 6, color, t);
      dl->AddCircle(centre, h * 0.20f, color, 24, t * 0.9f);
      return;
    case Icon::RingCyclopentadiene: {
      // Apex is the sp3 CH2; the two double bonds flank the bottom edge.
      ImVec2 v[5];
      polygonVerts(centre, h * 0.42f, 5, v);
      dl->AddPolyline(v, 5, color, ImDrawFlags_Closed, t);
      innerBondLine(dl, centre, v[1], v[2], h * 0.10f, color, t * 0.9f);
      innerBondLine(dl, centre, v[3], v[4], h * 0.10f, color, t * 0.9f);
      return;
    }
    case Icon::RingNaphthalene: {
      // Two hexagons fused along the vertical centre edge, one aromatic
      // circle per ring.
      const float r = h * 0.22f;
      const float dx = r * 0.8660254f;  // apothem: shared edge lands on x = cx
      regularPolygon(dl, ImVec2(centre.x - dx, centre.y), r, 6, color, t * 0.95f);
      regularPolygon(dl, ImVec2(centre.x + dx, centre.y), r, 6, color, t * 0.95f);
      dl->AddCircle(ImVec2(centre.x - dx, centre.y), r * 0.44f, color, 20, t * 0.85f);
      dl->AddCircle(ImVec2(centre.x + dx, centre.y), r * 0.44f, color, 20, t * 0.85f);
      return;
    }
    case Icon::Atom:
      dl->AddCircle(centre, h * 0.40f, color, 24, t);
      dl->AddCircleFilled(centre, h * 0.11f, color, 12);
      return;
    case Icon::ChargePlus:
      dl->AddCircle(centre, h * 0.42f, color, 24, t);
      stroke(dl, at(centre, h, -0.18f, 0.0f), at(centre, h, 0.18f, 0.0f), color, t);
      stroke(dl, at(centre, h, 0.0f, -0.18f), at(centre, h, 0.0f, 0.18f), color, t);
      return;
    case Icon::ChargeMinus:
      dl->AddCircle(centre, h * 0.42f, color, 24, t);
      stroke(dl, at(centre, h, -0.18f, 0.0f), at(centre, h, 0.18f, 0.0f), color, t);
      return;
    case Icon::BondSingle:
      stroke(dl, at(centre, h, -0.40f, 0.0f), at(centre, h, 0.40f, 0.0f), color,
             t * 1.1f);
      return;
    case Icon::BondDouble:
      stroke(dl, at(centre, h, -0.40f, -0.13f), at(centre, h, 0.40f, -0.13f), color,
             t);
      stroke(dl, at(centre, h, -0.40f, 0.13f), at(centre, h, 0.40f, 0.13f), color, t);
      return;
    case Icon::BondTriple:
      stroke(dl, at(centre, h, -0.40f, -0.22f), at(centre, h, 0.40f, -0.22f), color,
             t * 0.9f);
      stroke(dl, at(centre, h, -0.40f, 0.0f), at(centre, h, 0.40f, 0.0f), color,
             t * 0.9f);
      stroke(dl, at(centre, h, -0.40f, 0.22f), at(centre, h, 0.40f, 0.22f), color,
             t * 0.9f);
      return;
    case Icon::BondAromatic:
      stroke(dl, at(centre, h, -0.40f, -0.12f), at(centre, h, 0.40f, -0.12f), color,
             t);
      dashed(dl, at(centre, h, -0.40f, 0.14f), at(centre, h, 0.40f, 0.14f), color,
             t * 0.9f, h * 0.16f, h * 0.12f);
      return;
    case Icon::StereoNone:
      stroke(dl, at(centre, h, -0.40f, 0.0f), at(centre, h, 0.40f, 0.0f), color, t);
      return;
    case Icon::StereoWedge:
      dl->AddTriangleFilled(at(centre, h, -0.34f, -0.16f), at(centre, h, -0.34f, 0.16f),
                            at(centre, h, 0.40f, 0.0f), color);
      return;
    case Icon::StereoHash: {
      const float xs[] = {-0.30f, -0.11f, 0.08f, 0.27f};
      const float hs[] = {0.22f, 0.15f, 0.09f, 0.04f};
      for (int i = 0; i < 4; ++i) {
        stroke(dl, at(centre, h, xs[i], -hs[i]), at(centre, h, xs[i], hs[i]), color,
               t * 0.95f);
      }
      return;
    }
    case Icon::StereoWavy: {
      // Squiggly bond for unknown or mixed stereochemistry.
      constexpr int segments = 16;
      ImVec2 points[segments + 1];
      for (int i = 0; i <= segments; ++i) {
        const float u = static_cast<float>(i) / segments;
        points[i] = ImVec2(centre.x - 0.32f * h + u * 0.64f * h,
                           centre.y + std::sin(u * 3.0f * 6.2831853f) * 0.10f * h);
      }
      dl->AddPolyline(points, segments + 1, color, 0, t);
      return;
    }
    case Icon::Flask: {
      stroke(dl, at(centre, h, -0.09f, -0.42f), at(centre, h, -0.09f, -0.16f),
             color, t);
      stroke(dl, at(centre, h, 0.09f, -0.42f), at(centre, h, 0.09f, -0.16f),
             color, t);
      stroke(dl, at(centre, h, -0.13f, -0.42f), at(centre, h, 0.13f, -0.42f),
             color, t);
      stroke(dl, at(centre, h, -0.09f, -0.16f), at(centre, h, -0.23f, 0.00f),
             color, t);
      stroke(dl, at(centre, h, 0.09f, -0.16f), at(centre, h, 0.23f, 0.00f),
             color, t);
      arc(dl, at(centre, h, 0.0f, 0.16f), h * 0.28f, -0.61f, 3.75f, color, t,
          20);
      stroke(dl, at(centre, h, -0.25f, 0.20f), at(centre, h, 0.25f, 0.20f),
             color, t * 0.8f);
      return;
    }
    case Icon::Beaker: {
      const ImVec2 vessel[] = {{-0.32f, -0.34f}, {-0.26f, 0.38f},
                               {0.26f, 0.38f},   {0.32f, -0.34f}};
      polyline(dl, centre, h, vessel, 4, color, t, false);
      const ImVec2 lip[] = {{-0.36f, -0.34f}, {0.36f, -0.34f},
                            {0.42f, -0.27f}};
      polyline(dl, centre, h, lip, 3, color, t, false);
      stroke(dl, at(centre, h, 0.12f, -0.08f), at(centre, h, 0.29f, -0.08f),
             color, t * 0.8f);
      stroke(dl, at(centre, h, 0.10f, 0.12f), at(centre, h, 0.27f, 0.12f),
             color, t * 0.8f);
      return;
    }
    case Icon::SepFunnel: {
      dl->AddRect(at(centre, h, -0.12f, -0.44f), at(centre, h, 0.12f, -0.35f),
                  color, h * 0.03f, 0, t);
      stroke(dl, at(centre, h, -0.06f, -0.35f), at(centre, h, -0.06f, -0.27f),
             color, t);
      stroke(dl, at(centre, h, 0.06f, -0.35f), at(centre, h, 0.06f, -0.27f),
             color, t);
      dl->PathLineTo(at(centre, h, -0.06f, -0.27f));
      dl->PathBezierCubicCurveTo(at(centre, h, -0.38f, -0.16f),
                                 at(centre, h, -0.27f, 0.14f),
                                 at(centre, h, 0.0f, 0.27f));
      dl->PathBezierCubicCurveTo(at(centre, h, 0.27f, 0.14f),
                                 at(centre, h, 0.38f, -0.16f),
                                 at(centre, h, 0.06f, -0.27f));
      dl->PathStroke(color, 0, t);
      stroke(dl, at(centre, h, 0.0f, 0.27f), at(centre, h, 0.0f, 0.44f), color,
             t);
      stroke(dl, at(centre, h, -0.17f, 0.31f), at(centre, h, 0.17f, 0.31f),
             color, t);
      stroke(dl, at(centre, h, -0.09f, 0.25f), at(centre, h, 0.09f, 0.37f),
             color, t * 0.85f);
      return;
    }
    case Icon::Droplet: {
      dl->PathLineTo(at(centre, h, 0.0f, -0.44f));
      dl->PathBezierCubicCurveTo(at(centre, h, -0.08f, -0.20f),
                                 at(centre, h, -0.30f, 0.04f),
                                 at(centre, h, -0.28f, 0.20f));
      dl->PathBezierCubicCurveTo(at(centre, h, -0.25f, 0.48f),
                                 at(centre, h, 0.25f, 0.48f),
                                 at(centre, h, 0.28f, 0.20f));
      dl->PathBezierCubicCurveTo(at(centre, h, 0.30f, 0.04f),
                                 at(centre, h, 0.08f, -0.20f),
                                 at(centre, h, 0.0f, -0.44f));
      dl->PathStroke(color, ImDrawFlags_Closed, t);
      return;
    }
    case Icon::Vial: {
      dl->AddRect(at(centre, h, -0.25f, -0.26f), at(centre, h, 0.25f, 0.40f),
                  color, h * 0.10f, 0, t);
      dl->AddRect(at(centre, h, -0.29f, -0.40f), at(centre, h, 0.29f, -0.24f),
                  color, h * 0.04f, 0, t);
      stroke(dl, at(centre, h, -0.22f, 0.18f), at(centre, h, 0.22f, 0.18f),
             color, t * 0.8f);
      return;
    }
    case Icon::Thermometer: {
      dl->AddCircle(at(centre, h, 0.0f, 0.29f), h * 0.15f, color, 16, t);
      dl->AddRect(at(centre, h, -0.08f, -0.40f), at(centre, h, 0.08f, 0.24f),
                  color, h * 0.08f, 0, t);
      stroke(dl, at(centre, h, 0.08f, -0.20f), at(centre, h, 0.22f, -0.20f),
             color, t * 0.8f);
      stroke(dl, at(centre, h, 0.08f, 0.00f), at(centre, h, 0.20f, 0.00f),
             color, t * 0.8f);
      stroke(dl, at(centre, h, 0.0f, 0.25f), at(centre, h, 0.0f, -0.16f), color,
             t * 0.8f);
      return;
    }
    case Icon::Balance: {
      stroke(dl, at(centre, h, -0.40f, -0.20f), at(centre, h, 0.40f, -0.20f),
             color, t);
      const ImVec2 fulcrum[] = {{0.0f, -0.18f}, {-0.18f, 0.36f},
                                {0.18f, 0.36f}};
      polyline(dl, centre, h, fulcrum, 3, color, t, true);
      for (int side = -1; side <= 1; side += 2) {
        const float x = static_cast<float>(side) * 0.28f;
        stroke(dl, at(centre, h, x, -0.20f), at(centre, h, x, 0.08f), color,
               t * 0.8f);
        arc(dl, at(centre, h, x, 0.08f), h * 0.16f, 0.0f, 3.14159265f, color,
            t);
      }
      return;
    }
    case Icon::Molecule: {
      stroke(dl, at(centre, h, -0.23f, 0.20f), at(centre, h, 0.0f, -0.20f),
             color, t);
      stroke(dl, at(centre, h, 0.0f, -0.20f), at(centre, h, 0.29f, 0.17f),
             color, t);
      dl->AddCircle(at(centre, h, -0.27f, 0.25f), h * 0.12f, color, 14, t);
      dl->AddCircle(at(centre, h, 0.0f, -0.27f), h * 0.13f, color, 14, t);
      dl->AddCircle(at(centre, h, 0.34f, 0.23f), h * 0.11f, color, 14, t);
      return;
    }
    case Icon::Reaction:
      arrow(dl, centre, h, -0.40f, 0.0f, 0.40f, 0.0f, color, t);
      return;
    case Icon::Retro:
      arrow(dl, centre, h, 0.40f, -0.10f, -0.40f, -0.10f, color, t * 0.9f,
            0.16f);
      arrow(dl, centre, h, 0.40f, 0.10f, -0.40f, 0.10f, color, t * 0.9f,
            0.16f);
      return;
    case Icon::Ph: {
      dl->PathLineTo(at(centre, h, 0.0f, -0.42f));
      dl->PathBezierCubicCurveTo(at(centre, h, -0.08f, -0.20f),
                                 at(centre, h, -0.20f, -0.03f),
                                 at(centre, h, -0.20f, 0.09f));
      dl->PathBezierCubicCurveTo(at(centre, h, -0.20f, 0.30f),
                                 at(centre, h, 0.20f, 0.30f),
                                 at(centre, h, 0.20f, 0.09f));
      dl->PathBezierCubicCurveTo(at(centre, h, 0.20f, -0.03f),
                                 at(centre, h, 0.08f, -0.20f),
                                 at(centre, h, 0.0f, -0.42f));
      dl->PathStroke(color, ImDrawFlags_Closed, t);
      stroke(dl, at(centre, h, -0.38f, 0.38f), at(centre, h, 0.38f, 0.38f),
             color, t);
      stroke(dl, at(centre, h, 0.0f, 0.30f), at(centre, h, 0.0f, 0.44f), color,
             t);
      return;
    }
    case Icon::Energy: {
      const ImVec2 bolt[] = {{0.10f, -0.44f}, {-0.27f, 0.05f},
                             {-0.02f, 0.05f}, {-0.14f, 0.44f},
                             {0.30f, -0.10f}, {0.05f, -0.10f}};
      polyline(dl, centre, h, bolt, 6, color, t, true);
      return;
    }
    case Icon::Flame: {
      dl->PathLineTo(at(centre, h, 0.04f, -0.44f));
      dl->PathBezierCubicCurveTo(at(centre, h, 0.20f, -0.15f),
                                 at(centre, h, 0.38f, 0.02f),
                                 at(centre, h, 0.30f, 0.25f));
      dl->PathBezierCubicCurveTo(at(centre, h, 0.20f, 0.48f),
                                 at(centre, h, -0.24f, 0.48f),
                                 at(centre, h, -0.30f, 0.23f));
      dl->PathBezierCubicCurveTo(at(centre, h, -0.35f, 0.02f),
                                 at(centre, h, -0.08f, -0.08f),
                                 at(centre, h, 0.04f, -0.44f));
      dl->PathStroke(color, ImDrawFlags_Closed, t);
      const ImVec2 tongue[] = {{0.02f, -0.05f}, {-0.10f, 0.18f},
                               {0.0f, 0.34f},   {0.12f, 0.16f}};
      polyline(dl, centre, h, tongue, 4, color, t * 0.8f, false);
      return;
    }
    case Icon::Snowflake: {
      constexpr float pi = 3.14159265f;
      for (int axis = 0; axis < 3; ++axis) {
        const float angle = axis * pi / 3.0f;
        const float ux = std::cos(angle);
        const float uy = std::sin(angle);
        stroke(dl, at(centre, h, -ux * 0.40f, -uy * 0.40f),
               at(centre, h, ux * 0.40f, uy * 0.40f), color, t);
        for (int end = -1; end <= 1; end += 2) {
          const float ex = ux * 0.30f * static_cast<float>(end);
          const float ey = uy * 0.30f * static_cast<float>(end);
          for (int barb = -1; barb <= 1; barb += 2) {
            const float bx = ex - ux * 0.11f * static_cast<float>(end) -
                             uy * 0.10f * static_cast<float>(barb);
            const float by = ey - uy * 0.11f * static_cast<float>(end) +
                             ux * 0.10f * static_cast<float>(barb);
            stroke(dl, at(centre, h, ex, ey), at(centre, h, bx, by), color,
                   t * 0.75f);
          }
        }
      }
      return;
    }
    case Icon::Shake: {
      const ImVec2 vessel[] = {{-0.20f, -0.31f}, {-0.30f, 0.30f},
                               {0.30f, 0.30f},   {0.20f, -0.31f}};
      polyline(dl, centre, h, vessel, 4, color, t, false);
      stroke(dl, at(centre, h, -0.24f, -0.31f), at(centre, h, 0.24f, -0.31f),
             color, t);
      arc(dl, at(centre, h, -0.28f, 0.0f), h * 0.20f, 1.75f, 4.55f, color,
          t * 0.8f, 10);
      arc(dl, at(centre, h, 0.28f, 0.0f), h * 0.20f, -1.40f, 1.40f, color,
          t * 0.8f, 10);
      return;
    }
    case Icon::Timer:
      dl->AddCircle(at(centre, h, 0.0f, 0.06f), h * 0.34f, color, 24, t);
      stroke(dl, at(centre, h, 0.0f, 0.06f), at(centre, h, 0.0f, -0.18f), color,
             t);
      stroke(dl, at(centre, h, 0.0f, 0.06f), at(centre, h, 0.17f, 0.16f), color,
             t);
      stroke(dl, at(centre, h, 0.0f, -0.39f), at(centre, h, 0.0f, -0.28f), color,
             t);
      stroke(dl, at(centre, h, -0.10f, -0.39f), at(centre, h, 0.10f, -0.39f),
             color, t);
      return;
    case Icon::ChartBars:
      stroke(dl, at(centre, h, -0.38f, -0.34f), at(centre, h, -0.38f, 0.36f),
             color, t);
      stroke(dl, at(centre, h, -0.38f, 0.36f), at(centre, h, 0.40f, 0.36f),
             color, t);
      stroke(dl, at(centre, h, -0.20f, 0.32f), at(centre, h, -0.20f, 0.06f),
             color, t * 2.2f);
      stroke(dl, at(centre, h, 0.02f, 0.32f), at(centre, h, 0.02f, -0.18f), color,
             t * 2.2f);
      stroke(dl, at(centre, h, 0.25f, 0.32f), at(centre, h, 0.25f, -0.02f), color,
             t * 2.2f);
      return;
    case Icon::ChartLine: {
      stroke(dl, at(centre, h, -0.38f, -0.34f), at(centre, h, -0.38f, 0.36f),
             color, t);
      stroke(dl, at(centre, h, -0.38f, 0.36f), at(centre, h, 0.40f, 0.36f),
             color, t);
      const ImVec2 points[] = {{-0.30f, 0.20f}, {-0.10f, -0.04f},
                               {0.08f, 0.10f},  {0.35f, -0.26f}};
      polyline(dl, centre, h, points, 4, color, t, false);
      return;
    }
    case Icon::ChartScatter: {
      stroke(dl, at(centre, h, -0.38f, -0.34f), at(centre, h, -0.38f, 0.36f),
             color, t);
      stroke(dl, at(centre, h, -0.38f, 0.36f), at(centre, h, 0.40f, 0.36f),
             color, t);
      const ImVec2 points[] = {{-0.20f, 0.17f}, {-0.02f, -0.08f},
                               {0.18f, 0.08f},  {0.32f, -0.25f}};
      for (const ImVec2& p : points) {
        dl->AddCircleFilled(at(centre, h, p.x, p.y), h * 0.055f, color, 8);
      }
      return;
    }
    case Icon::Gauge:
      arc(dl, at(centre, h, 0.0f, 0.08f), h * 0.36f, 2.6179939f, 6.8067841f,
          color, t, 24);
      stroke(dl, at(centre, h, 0.0f, 0.08f), at(centre, h, 0.22f, -0.15f), color,
             t);
      dl->AddCircleFilled(at(centre, h, 0.0f, 0.08f), h * 0.06f, color, 10);
      return;
    case Icon::Layers: {
      const ImVec2 layer[] = {{-0.36f, -0.10f}, {-0.10f, -0.28f},
                              {0.36f, -0.10f},  {0.10f, 0.08f}};
      for (int i = 0; i < 3; ++i) {
        polyline(dl, at(centre, h, 0.0f, i * 0.17f - 0.12f), h, layer, 4, color,
                 t * 0.85f, true);
      }
      return;
    }
    case Icon::Grid:
      for (int y = 0; y < 2; ++y) {
        for (int x = 0; x < 2; ++x) {
          const float x0 = -0.38f + x * 0.42f;
          const float y0 = -0.38f + y * 0.42f;
          dl->AddRect(at(centre, h, x0, y0), at(centre, h, x0 + 0.32f, y0 + 0.32f),
                      color, h * 0.03f, 0, t);
        }
      }
      return;
    case Icon::List:
      for (int row = -1; row <= 1; ++row) {
        const float y = row * 0.25f;
        dl->AddCircleFilled(at(centre, h, -0.34f, y), h * 0.045f, color, 8);
        stroke(dl, at(centre, h, -0.20f, y), at(centre, h, 0.38f, y), color, t);
      }
      return;
    case Icon::Cube: {
      ImVec2 vertices[6];
      polygonVerts(centre, h * 0.42f, 6, vertices);
      dl->AddPolyline(vertices, 6, color, ImDrawFlags_Closed, t);
      stroke(dl, centre, vertices[0], color, t);
      stroke(dl, centre, vertices[2], color, t);
      stroke(dl, centre, vertices[4], color, t);
      return;
    }
    case Icon::Ruler:
      dl->AddRect(at(centre, h, -0.42f, -0.25f), at(centre, h, 0.42f, 0.25f),
                  color, h * 0.04f, 0, t);
      for (int i = 0; i < 5; ++i) {
        const float x = -0.28f + i * 0.14f;
        const float y = i % 2 == 0 ? -0.02f : -0.10f;
        stroke(dl, at(centre, h, x, -0.25f), at(centre, h, x, y), color,
               t * 0.8f);
      }
      return;
    case Icon::Table:
      dl->AddRect(at(centre, h, -0.40f, -0.34f), at(centre, h, 0.40f, 0.34f),
                  color, h * 0.03f, 0, t);
      stroke(dl, at(centre, h, 0.0f, -0.34f), at(centre, h, 0.0f, 0.34f), color,
             t);
      stroke(dl, at(centre, h, -0.40f, 0.0f), at(centre, h, 0.40f, 0.0f), color,
             t);
      return;
    case Icon::Play:
      dl->AddTriangleFilled(at(centre, h, -0.24f, -0.38f),
                            at(centre, h, -0.24f, 0.38f),
                            at(centre, h, 0.38f, 0.0f), color);
      return;
    case Icon::Pause:
      dl->AddRectFilled(at(centre, h, -0.30f, -0.38f),
                        at(centre, h, -0.08f, 0.38f), color, h * 0.03f);
      dl->AddRectFilled(at(centre, h, 0.08f, -0.38f),
                        at(centre, h, 0.30f, 0.38f), color, h * 0.03f);
      return;
    case Icon::Stop:
      dl->AddRectFilled(at(centre, h, -0.34f, -0.34f),
                        at(centre, h, 0.34f, 0.34f), color, h * 0.04f);
      return;
    case Icon::Rewind:
      dl->AddTriangleFilled(at(centre, h, -0.42f, 0.0f),
                            at(centre, h, -0.04f, -0.34f),
                            at(centre, h, -0.04f, 0.34f), color);
      dl->AddTriangleFilled(at(centre, h, -0.02f, 0.0f),
                            at(centre, h, 0.36f, -0.34f),
                            at(centre, h, 0.36f, 0.34f), color);
      return;
    case Icon::StepForward:
      dl->AddTriangleFilled(at(centre, h, -0.36f, -0.34f),
                            at(centre, h, -0.36f, 0.34f),
                            at(centre, h, 0.20f, 0.0f), color);
      stroke(dl, at(centre, h, 0.32f, -0.36f), at(centre, h, 0.32f, 0.36f),
             color, t * 1.6f);
      return;
    case Icon::Check:
      stroke(dl, at(centre, h, -0.36f, 0.02f), at(centre, h, -0.10f, 0.28f),
             color, t);
      stroke(dl, at(centre, h, -0.10f, 0.28f), at(centre, h, 0.40f, -0.30f),
             color, t);
      return;
    case Icon::Filter: {
      const ImVec2 funnel[] = {{-0.40f, -0.34f}, {0.40f, -0.34f},
                               {0.12f, -0.02f},  {0.12f, 0.34f},
                               {-0.08f, 0.42f},  {-0.08f, -0.02f}};
      polyline(dl, centre, h, funnel, 6, color, t, true);
      return;
    }
    case Icon::Trash:
      dl->AddRect(at(centre, h, -0.28f, -0.16f), at(centre, h, 0.28f, 0.40f),
                  color, h * 0.04f, 0, t);
      stroke(dl, at(centre, h, -0.36f, -0.24f), at(centre, h, 0.36f, -0.24f),
             color, t);
      stroke(dl, at(centre, h, -0.14f, -0.36f), at(centre, h, 0.14f, -0.36f),
             color, t);
      stroke(dl, at(centre, h, 0.0f, -0.36f), at(centre, h, 0.0f, -0.24f), color,
             t);
      stroke(dl, at(centre, h, -0.10f, -0.05f), at(centre, h, -0.10f, 0.28f),
             color, t * 0.75f);
      stroke(dl, at(centre, h, 0.10f, -0.05f), at(centre, h, 0.10f, 0.28f), color,
             t * 0.75f);
      return;
    case Icon::Save:
      dl->AddRect(at(centre, h, -0.38f, -0.40f), at(centre, h, 0.38f, 0.40f),
                  color, h * 0.04f, 0, t);
      dl->AddRect(at(centre, h, -0.22f, -0.40f), at(centre, h, 0.22f, -0.10f),
                  color, 0.0f, 0, t);
      dl->AddRect(at(centre, h, -0.22f, 0.10f), at(centre, h, 0.22f, 0.40f),
                  color, h * 0.04f, 0, t);
      return;
    case Icon::Folder: {
      const ImVec2 folder[] = {{-0.40f, -0.25f}, {-0.13f, -0.25f},
                               {-0.03f, -0.12f}, {0.40f, -0.12f},
                               {0.34f, 0.34f},   {-0.40f, 0.34f}};
      polyline(dl, centre, h, folder, 6, color, t, true);
      return;
    }
    case Icon::Undo:
      arc(dl, at(centre, h, 0.03f, 0.07f), h * 0.32f, 3.45f, 6.42f, color, t,
          18);
      stroke(dl, at(centre, h, -0.29f, -0.03f), at(centre, h, -0.38f, -0.27f),
             color, t);
      stroke(dl, at(centre, h, -0.29f, -0.03f), at(centre, h, -0.10f, -0.20f),
             color, t);
      return;
    case Icon::Redo:
      arc(dl, at(centre, h, -0.03f, 0.07f), h * 0.32f, 3.00f, 5.97f, color, t,
          18);
      stroke(dl, at(centre, h, 0.29f, -0.03f), at(centre, h, 0.38f, -0.27f),
             color, t);
      stroke(dl, at(centre, h, 0.29f, -0.03f), at(centre, h, 0.10f, -0.20f),
             color, t);
      return;
    case Icon::Gear: {
      dl->AddCircle(centre, h * 0.30f, color, 22, t);
      dl->AddCircle(centre, h * 0.12f, color, 16, t);
      constexpr float pi = 3.14159265f;
      for (int tooth = 0; tooth < 8; ++tooth) {
        const float angle = tooth * pi / 4.0f;
        const float ux = std::cos(angle);
        const float uy = std::sin(angle);
        stroke(dl, at(centre, h, ux * 0.31f, uy * 0.31f),
               at(centre, h, ux * 0.43f, uy * 0.43f), color, t * 1.2f);
      }
      return;
    }
    case Icon::Eye:
    case Icon::EyeOff: {
      dl->PathLineTo(at(centre, h, -0.42f, 0.0f));
      dl->PathBezierCubicCurveTo(at(centre, h, -0.20f, -0.30f),
                                 at(centre, h, 0.20f, -0.30f),
                                 at(centre, h, 0.42f, 0.0f));
      dl->PathBezierCubicCurveTo(at(centre, h, 0.20f, 0.30f),
                                 at(centre, h, -0.20f, 0.30f),
                                 at(centre, h, -0.42f, 0.0f));
      dl->PathStroke(color, ImDrawFlags_Closed, t);
      dl->AddCircle(centre, h * 0.13f, color, 16, t);
      if (icon == Icon::EyeOff) {
        stroke(dl, at(centre, h, -0.38f, -0.38f), at(centre, h, 0.38f, 0.38f),
               color, t * 1.1f);
      }
      return;
    }
    case Icon::Lock:
      arc(dl, at(centre, h, 0.0f, -0.12f), h * 0.23f, 3.14159265f, 6.2831853f,
          color, t, 14);
      stroke(dl, at(centre, h, -0.23f, -0.12f), at(centre, h, -0.23f, 0.02f),
             color, t);
      stroke(dl, at(centre, h, 0.23f, -0.12f), at(centre, h, 0.23f, 0.02f),
             color, t);
      dl->AddRect(at(centre, h, -0.34f, 0.0f), at(centre, h, 0.34f, 0.40f),
                  color, h * 0.06f, 0, t);
      return;
    case Icon::Link:
      dl->AddRect(at(centre, h, -0.42f, -0.20f), at(centre, h, 0.10f, 0.10f),
                  color, h * 0.15f, 0, t);
      dl->AddRect(at(centre, h, -0.10f, -0.10f), at(centre, h, 0.42f, 0.20f),
                  color, h * 0.15f, 0, t);
      stroke(dl, at(centre, h, -0.12f, 0.0f), at(centre, h, 0.12f, 0.0f), color,
             t);
      return;
    case Icon::Star: {
      ImVec2 points[10];
      for (int i = 0; i < 10; ++i) {
        const float angle = -1.5707963f + i * 0.6283185f;
        const float radius = h * (i % 2 == 0 ? 0.42f : 0.18f);
        points[i] = ImVec2(centre.x + std::cos(angle) * radius,
                           centre.y + std::sin(angle) * radius);
      }
      dl->AddPolyline(points, 10, color, ImDrawFlags_Closed, t);
      return;
    }
    case Icon::Warning:
      dl->AddTriangle(at(centre, h, 0.0f, -0.43f),
                      at(centre, h, -0.42f, 0.36f),
                      at(centre, h, 0.42f, 0.36f), color, t);
      stroke(dl, at(centre, h, 0.0f, -0.17f), at(centre, h, 0.0f, 0.10f), color,
             t);
      dl->AddCircleFilled(at(centre, h, 0.0f, 0.24f), h * 0.04f, color, 8);
      return;
    case Icon::Info:
      dl->AddCircle(centre, h * 0.41f, color, 24, t);
      dl->AddCircleFilled(at(centre, h, 0.0f, -0.20f), h * 0.045f, color, 8);
      stroke(dl, at(centre, h, 0.0f, -0.05f), at(centre, h, 0.0f, 0.24f), color,
             t);
      return;
    case Icon::Sparkle:
      fourPointStar(dl, at(centre, h, -0.08f, 0.02f), h, 0.82f, color, t);
      fourPointStar(dl, at(centre, h, 0.30f, -0.27f), h, 0.30f, color, t * 0.8f);
      return;
    case Icon::Book: {
      const ImVec2 left[] = {{-0.02f, -0.30f}, {-0.34f, -0.39f},
                             {-0.40f, 0.26f},  {-0.02f, 0.38f}};
      const ImVec2 right[] = {{0.02f, -0.30f}, {0.34f, -0.39f},
                              {0.40f, 0.26f},  {0.02f, 0.38f}};
      polyline(dl, centre, h, left, 4, color, t, false);
      polyline(dl, centre, h, right, 4, color, t, false);
      stroke(dl, at(centre, h, 0.0f, -0.30f), at(centre, h, 0.0f, 0.38f), color,
             t);
      return;
    }
    case Icon::Crosshair:
      dl->AddCircle(centre, h * 0.25f, color, 20, t);
      stroke(dl, at(centre, h, -0.42f, 0.0f), at(centre, h, -0.25f, 0.0f), color,
             t);
      stroke(dl, at(centre, h, 0.25f, 0.0f), at(centre, h, 0.42f, 0.0f), color,
             t);
      stroke(dl, at(centre, h, 0.0f, -0.42f), at(centre, h, 0.0f, -0.25f), color,
             t);
      stroke(dl, at(centre, h, 0.0f, 0.25f), at(centre, h, 0.0f, 0.42f), color,
             t);
      return;
    case Icon::ArrowLeft:
      arrow(dl, centre, h, 0.36f, 0.0f, -0.30f, 0.0f, color, t, 0.20f);
      return;
    case Icon::ChevronDown: {
      const ImVec2 points[] = {{-0.28f, -0.14f}, {0.0f, 0.16f},
                               {0.28f, -0.14f}};
      polyline(dl, centre, h, points, 3, color, t, false);
      return;
    }
    case Icon::ChevronRight: {
      const ImVec2 points[] = {{-0.14f, -0.28f}, {0.16f, 0.0f},
                               {-0.14f, 0.28f}};
      polyline(dl, centre, h, points, 3, color, t, false);
      return;
    }
    case Icon::ChevronLeft: {
      const ImVec2 points[] = {{0.14f, -0.28f}, {-0.16f, 0.0f},
                               {0.14f, 0.28f}};
      polyline(dl, centre, h, points, 3, color, t, false);
      return;
    }
    case Icon::ChevronUp: {
      const ImVec2 points[] = {{-0.28f, 0.14f}, {0.0f, -0.16f},
                               {0.28f, 0.14f}};
      polyline(dl, centre, h, points, 3, color, t, false);
      return;
    }
    case Icon::DragHandle:
      for (int row = -1; row <= 1; row += 2) {
        for (int column = -1; column <= 1; ++column) {
          dl->AddCircleFilled(
              at(centre, h, static_cast<float>(column) * 0.20f,
                 static_cast<float>(row) * 0.16f),
              h * 0.045f, color, 8);
        }
      }
      return;
    case Icon::Minus:
      stroke(dl, at(centre, h, -0.30f, 0.0f), at(centre, h, 0.30f, 0.0f), color, t);
      return;
    case Icon::ZoomFit: {
      // Four corner brackets: fit-to-window.
      const float c = 0.36f;
      const float a = 0.16f;
      for (int sx = -1; sx <= 1; sx += 2) {
        for (int sy = -1; sy <= 1; sy += 2) {
          const float fx = static_cast<float>(sx), fy = static_cast<float>(sy);
          stroke(dl, at(centre, h, fx * c, fy * (c - a)), at(centre, h, fx * c, fy * c), color,
                 t);
          stroke(dl, at(centre, h, fx * (c - a), fy * c), at(centre, h, fx * c, fy * c), color,
                 t);
        }
      }
      return;
    }
    case Icon::Plus:
      stroke(dl, at(centre, h, -0.30f, 0.0f), at(centre, h, 0.30f, 0.0f), color, t);
      stroke(dl, at(centre, h, 0.0f, -0.30f), at(centre, h, 0.0f, 0.30f), color, t);
      return;
    case Icon::Close:
      stroke(dl, at(centre, h, -0.26f, -0.26f), at(centre, h, 0.26f, 0.26f), color, t);
      stroke(dl, at(centre, h, -0.26f, 0.26f), at(centre, h, 0.26f, -0.26f), color, t);
      return;
    case Icon::Search: {
      const ImVec2 c = at(centre, h, -0.10f, -0.10f);
      dl->AddCircle(c, h * 0.28f, color, 20, t);
      stroke(dl, at(centre, h, 0.12f, 0.12f), at(centre, h, 0.40f, 0.40f), color,
             t * 1.1f);
      return;
    }
    case Icon::Copy: {
      const float r = h * 0.08f;
      dl->AddRect(at(centre, h, -0.34f, -0.40f), at(centre, h, 0.10f, 0.06f), color,
                  r, 0, t * 0.9f);
      dl->AddRect(at(centre, h, -0.10f, -0.06f), at(centre, h, 0.34f, 0.40f), color,
                  r, 0, t);
      return;
    }
    case Icon::ArrowRight:
      stroke(dl, at(centre, h, -0.36f, 0.0f), at(centre, h, 0.30f, 0.0f), color, t);
      stroke(dl, at(centre, h, 0.30f, 0.0f), at(centre, h, 0.10f, -0.18f), color, t);
      stroke(dl, at(centre, h, 0.30f, 0.0f), at(centre, h, 0.10f, 0.18f), color, t);
      return;
    case Icon::Logo:
      hexagon(dl, centre, h * 0.46f, color, t);
      dl->AddCircle(centre, h * 0.22f, color, 24, t * 0.9f);
      return;
  }
}

}  // namespace chemcad::ui::icons

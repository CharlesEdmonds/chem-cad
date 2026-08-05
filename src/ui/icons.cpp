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

void hexagon(ImDrawList* dl, ImVec2 c, float radius, ImU32 col, float t) {
  ImVec2 pts[6];
  for (int i = 0; i < 6; ++i) {
    const float angle = (-90.0f + 60.0f * static_cast<float>(i)) * 3.14159265f / 180.0f;
    pts[i] = ImVec2(c.x + radius * std::cos(angle), c.y + radius * std::sin(angle));
  }
  dl->AddPolyline(pts, 6, col, ImDrawFlags_Closed, t);
}

}  // namespace

void draw(ImDrawList* dl, Icon icon, ImVec2 centre, float size, ImU32 color,
          float thickness) {
  if (!dl || size <= 0.0f) return;
  // Glyph coordinates are normalized around ±0.4. Using `size` as the unit
  // lets the drawing occupy ~80% of the requested box instead of only ~40%.
  const float h = size;
  const float t = thickness > 0.0f ? thickness : std::max(1.4f, size * 0.078f);

  switch (icon) {
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
    case Icon::Ring:
      hexagon(dl, centre, h * 0.44f, color, t);
      return;
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

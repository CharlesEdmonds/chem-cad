#include "core/sprout.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace chemcad::core {
namespace {
constexpr float kPi = std::numbers::pi_v<float>;
constexpr float kDeg2Rad = kPi / 180.0f;

float wrapTwoPi(float a) {
  while (a < 0) a += 2 * kPi;
  while (a >= 2 * kPi) a -= 2 * kPi;
  return a;
}
}  // namespace

float length(Vec2 v) { return std::sqrt(v.x * v.x + v.y * v.y); }

Vec2 normalize(Vec2 v) {
  const float len = length(v);
  if (len < 1e-6f) return {1.0f, 0.0f};
  return {v.x / len, v.y / len};
}

float angleOf(Vec2 v) { return std::atan2(v.y, v.x); }

Vec2 fromAngle(float radians, float len) {
  return {std::cos(radians) * len, std::sin(radians) * len};
}

Vec2 snapAngle(Vec2 rawDir) {
  const float step = kSnapDegrees * kDeg2Rad;
  const float snapped = std::round(angleOf(rawDir) / step) * step;
  return fromAngle(snapped);
}

Vec2 sproutDirection(const Molecule& mol, AtomId from) {
  const Atom* self = mol.atom(from);
  if (!self) return fromAngle(30.0f * kDeg2Rad);

  const std::vector<AtomId> nbrs = mol.neighbors(from);

  // Free atom: the canonical starting direction, 30 degrees up from +x.
  if (nbrs.empty()) return fromAngle(30.0f * kDeg2Rad);

  if (nbrs.size() == 1) {
    const Atom* n = mol.atom(nbrs[0]);
    if (!n) return fromAngle(30.0f * kDeg2Rad);
    // Direction pointing back at the existing neighbour.
    const float toNbr = angleOf({n->pos.x - self->pos.x, n->pos.y - self->pos.y});
    const float optA = toNbr + 120.0f * kDeg2Rad;
    const float optB = toNbr - 120.0f * kDeg2Rad;

    // Continue a trans zig-zag: prefer the option pointing away from the
    // neighbour's own substituent, so chains alternate instead of curling.
    const std::vector<AtomId> nn = mol.neighbors(nbrs[0]);
    for (AtomId far : nn) {
      if (far == from) continue;
      const Atom* f = mol.atom(far);
      if (!f) continue;
      const Vec2 refDir{f->pos.x - n->pos.x, f->pos.y - n->pos.y};
      const Vec2 dirA = fromAngle(optA);
      const Vec2 dirB = fromAngle(optB);
      // Pick the option least aligned with the neighbour's other bond.
      const float dotA = dirA.x * refDir.x + dirA.y * refDir.y;
      const float dotB = dirB.x * refDir.x + dirB.y * refDir.y;
      return dotA < dotB ? dirA : dirB;
    }
    // Terminal neighbour: take the upward option.
    return std::sin(optA) >= std::sin(optB) ? fromAngle(optA) : fromAngle(optB);
  }

  // Crowded atom: bisect the widest angular gap between existing bonds.
  std::vector<float> angles;
  angles.reserve(nbrs.size());
  for (AtomId n : nbrs) {
    const Atom* a = mol.atom(n);
    if (!a) continue;
    angles.push_back(wrapTwoPi(angleOf({a->pos.x - self->pos.x, a->pos.y - self->pos.y})));
  }
  if (angles.empty()) return fromAngle(30.0f * kDeg2Rad);
  std::sort(angles.begin(), angles.end());

  float bestGap = -1.0f, bestMid = 30.0f * kDeg2Rad;
  for (size_t i = 0; i < angles.size(); ++i) {
    const float cur = angles[i];
    const float next = (i + 1 < angles.size()) ? angles[i + 1] : angles[0] + 2 * kPi;
    const float gap = next - cur;
    if (gap > bestGap) {
      bestGap = gap;
      bestMid = cur + gap * 0.5f;
    }
  }
  return fromAngle(bestMid);
}

Vec2 sproutPosition(const Molecule& mol, AtomId from) {
  const Atom* self = mol.atom(from);
  if (!self) return {0, 0};
  const Vec2 d = sproutDirection(mol, from);
  return {self->pos.x + d.x * kBondLength, self->pos.y + d.y * kBondLength};
}

}  // namespace chemcad::core

#pragma once
// Automatic bond-angle geometry -- the "auto bond angles" behaviour that makes
// a sketcher feel like ChemDraw. Pure math over core::Molecule, no RDKit.

#include "core/model.hpp"

namespace chemcad::core {

inline constexpr float kBondLength = 1.0f;  // world units
inline constexpr float kSnapDegrees = 15.0f;

// Unit direction for a NEW bond sprouting from `from`:
//   0 neighbours  -> 30 degrees above the +x axis
//   1 neighbour   -> +/-120 degrees off the existing bond, choosing the side
//                    that continues a trans zig-zag (uses the neighbour's own
//                    other bond when it has one, otherwise the upward option)
//   2+ neighbours -> bisector of the largest angular gap
Vec2 sproutDirection(const Molecule&, AtomId from);

// Snap an arbitrary drag direction to the nearest 15 degree increment.
Vec2 snapAngle(Vec2 rawDir);

// Small vector helpers shared by the geometry and the canvas.
Vec2 normalize(Vec2 v);
float length(Vec2 v);
float angleOf(Vec2 v);              // radians, atan2(y, x)
Vec2 fromAngle(float radians, float len = 1.0f);

// Position a new atom one standard bond length away from `from`.
Vec2 sproutPosition(const Molecule&, AtomId from);

}  // namespace chemcad::core

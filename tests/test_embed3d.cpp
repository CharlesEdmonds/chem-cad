// Tests for chem::embed3D: RDKit ETKDG + MMFF behind the RDKit-free
// boundary. These defend geometry the viewer depends on, not RDKit internals.
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <cmath>

#include "chem/bridge.hpp"
#include "chem/embed3d.hpp"

namespace chem = chemcad::chem;
namespace core = chemcad::core;

namespace {

float distance(const chem::Atom3D& a, const chem::Atom3D& b) {
  const float dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

}  // namespace

TEST_CASE("benzene embeds as a flat regular hexagon of aromatic bonds") {
  const chem::Embedded3D model = chem::embed3D(chem::fromSmiles("c1ccccc1"));
  REQUIRE(model.atoms.size() == 12);  // CPK models keep hydrogens: 6 C + 6 H
  REQUIRE(model.bonds.size() == 12);
  int ringBonds = 0;
  for (const chem::Bond3D& b : model.bonds) {
    const chem::Atom3D& a = model.atoms[static_cast<size_t>(b.a)];
    const chem::Atom3D& c = model.atoms[static_cast<size_t>(b.b)];
    if (a.atomicNumber != 6 || c.atomicNumber != 6) continue;  // skip C-H
    ++ringBonds;
    CHECK(b.order == 4);  // aromatic marker
    // Aromatic C-C sits between single and double bond lengths.
    const float d = distance(a, c);
    CHECK(d == doctest::Approx(1.39f).epsilon(0.08));
  }
  CHECK(ringBonds == 6);
  // Aromatic rings are planar.
  for (const chem::Atom3D& a : model.atoms) {
    CHECK(std::fabs(a.z) < 0.15f);
  }
  CHECK(model.radius > 1.0f);
}

TEST_CASE("acetic acid keeps its bond orders and sane lengths") {
  const chem::Embedded3D model = chem::embed3D(chem::fromSmiles("CC(=O)O"));
  REQUIRE(model.atoms.size() == 8);  // 2 C + 2 O + 4 H (hydrogens kept)
  REQUIRE(model.bonds.size() == 7);
  int singles = 0, doubles = 0;
  for (const chem::Bond3D& b : model.bonds) {
    const chem::Atom3D& a = model.atoms[static_cast<size_t>(b.a)];
    const chem::Atom3D& c = model.atoms[static_cast<size_t>(b.b)];
    if (a.atomicNumber == 1 || c.atomicNumber == 1) continue;  // skip X-H
    if (b.order == 1) ++singles;
    if (b.order == 2) ++doubles;
  }
  CHECK(singles == 2);   // C-C and C-OH
  CHECK(doubles == 1);   // C=O
  // The carbonyl is measurably shorter than the C-O single bond.
  for (const chem::Bond3D& b : model.bonds) {
    const float d = distance(model.atoms[static_cast<size_t>(b.a)],
                             model.atoms[static_cast<size_t>(b.b)]);
    if (b.order == 2) {
      CHECK(d == doctest::Approx(1.22f).epsilon(0.10));
    }
  }
}

TEST_CASE("embedding is deterministic and centred") {
  const chem::Embedded3D first = chem::embed3D(chem::fromSmiles("CCO"));
  const chem::Embedded3D second = chem::embed3D(chem::fromSmiles("CCO"));
  REQUIRE(first.atoms.size() == second.atoms.size());
  for (size_t i = 0; i < first.atoms.size(); ++i) {
    CHECK(first.atoms[i].x == doctest::Approx(second.atoms[i].x));
    CHECK(first.atoms[i].y == doctest::Approx(second.atoms[i].y));
    CHECK(first.atoms[i].z == doctest::Approx(second.atoms[i].z));
  }
  // Centroid sits at the origin, so the turntable rotates about the model.
  float cx = 0.0f, cy = 0.0f, cz = 0.0f;
  for (const chem::Atom3D& a : first.atoms) {
    cx += a.x;
    cy += a.y;
    cz += a.z;
  }
  const float n = static_cast<float>(first.atoms.size());
  CHECK(std::fabs(cx / n) < 1e-3f);
  CHECK(std::fabs(cy / n) < 1e-3f);
  CHECK(std::fabs(cz / n) < 1e-3f);
}

TEST_CASE("an empty molecule cannot be embedded") {
  const core::Molecule empty;
  CHECK_THROWS_AS(chem::embed3D(empty), chem::ChemError);
}

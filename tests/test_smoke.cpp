#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "core/model.hpp"
#include "core/sprout.hpp"

using namespace chemcad;

// Wave-0 sanity: the model links and the graph invariants hold. The full core
// suite lives in test_core.cpp.
TEST_CASE("molecule graph basics") {
  core::Molecule m;
  const core::AtomId a = m.addAtom({});
  const core::AtomId b = m.addAtom({});
  REQUIRE(a != core::kInvalidAtom);
  REQUIRE(b != a);
  const core::BondId bond = m.addBond(a, b, core::BondOrder::Single);
  CHECK(bond != core::kInvalidBond);
  CHECK(m.atomCount() == 2);
  CHECK(m.bondCount() == 1);

  m.removeAtom(a);
  CHECK(m.atomCount() == 1);
  CHECK(m.bondCount() == 0);  // incident bonds go with the atom
}

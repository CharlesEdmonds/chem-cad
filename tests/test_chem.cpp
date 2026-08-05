#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <cmath>
#include <string>

#include "chem/bridge.hpp"
#include "core/model.hpp"

namespace chem = chemcad::chem;
namespace core = chemcad::core;

TEST_CASE("SMILES roundtrip is canonical") {
  const std::string input = "c1ccccc1CC(=O)O";
  const core::Molecule molecule = chem::fromSmiles(input);
  CHECK(chem::toSmiles(molecule) == chem::canonicalize(input));
}

TEST_CASE("SMILES import creates normalized two-dimensional coordinates") {
  const core::Molecule molecule = chem::fromSmiles("c1ccccc1CC(=O)O");
  REQUIRE(molecule.atomCount() > 1);
  REQUIRE(molecule.bondCount() > 0);

  bool anyNonzero = false;
  bool anyDifferent = false;
  const core::Vec2 first = molecule.atoms().front().pos;
  for (const core::Atom& atom : molecule.atoms()) {
    anyNonzero = anyNonzero || std::abs(atom.pos.x) > 1e-5f || std::abs(atom.pos.y) > 1e-5f;
    anyDifferent = anyDifferent || std::abs(atom.pos.x - first.x) > 1e-5f ||
                                     std::abs(atom.pos.y - first.y) > 1e-5f;
  }
  CHECK(anyNonzero);
  CHECK(anyDifferent);

  double totalLength = 0.0;
  for (const core::Bond& bond : molecule.bonds()) {
    const core::Atom* begin = molecule.atom(bond.a);
    const core::Atom* end = molecule.atom(bond.b);
    REQUIRE(begin != nullptr);
    REQUIRE(end != nullptr);
    totalLength += std::hypot(static_cast<double>(end->pos.x - begin->pos.x),
                              static_cast<double>(end->pos.y - begin->pos.y));
  }
  const double meanLength = totalLength / static_cast<double>(molecule.bondCount());
  CHECK(meanLength == doctest::Approx(1.0).epsilon(0.2));
}

TEST_CASE("properties report formula mass and rings") {
  const chem::Properties ethanol = chem::computeProperties(chem::fromSmiles("CCO"));
  CHECK(ethanol.formula == "C2H6O");
  CHECK(ethanol.mw == doctest::Approx(46.0419).epsilon(0.001));
  CHECK(ethanol.rings == 0);

  const chem::Properties benzene = chem::computeProperties(chem::fromSmiles("c1ccccc1"));
  CHECK(benzene.rings == 1);
}

TEST_CASE("terminal single-bonded carbons are methyl groups") {
  const core::Molecule ethane = chem::fromSmiles("CC");
  REQUIRE(ethane.atomCount() == 2);
  for (const core::Atom& atom : ethane.atoms()) {
    CHECK(atom.atomicNumber == 6);
    CHECK(ethane.degree(atom.id) == 1);
    CHECK(chem::implicitHCount(ethane, atom.id) == 3);
  }
}

TEST_CASE("invalid valence is translated to ChemError") {
  core::Molecule molecule;
  const core::AtomId center = molecule.addAtom({});
  for (int i = 0; i < 5; ++i) {
    const core::AtomId neighbor = molecule.addAtom({});
    molecule.addBond(center, neighbor, core::BondOrder::Single);
  }

  CHECK_THROWS_AS(chem::computeProperties(molecule), chem::ChemError);
  CHECK_NOTHROW(chem::implicitHCount(molecule, center));
  CHECK(chem::implicitHCount(molecule, core::kInvalidAtom) == 0);
  CHECK(chem::implicitHCount(molecule, 999999) == 0);
}

TEST_CASE("invalid SMILES is rejected") {
  CHECK_THROWS_AS(chem::fromSmiles("not a smiles at all"), chem::ChemError);
}

TEST_CASE("SVG export contains a complete SVG element") {
  const std::string svg = chem::toSvg(chem::fromSmiles("CCO"), 320, 200);
  CHECK(svg.find("<svg") != std::string::npos);
  CHECK(svg.find("</svg>") != std::string::npos);
}

TEST_CASE("MolBlock roundtrip preserves chemical identity") {
  const core::Molecule ethanol = chem::fromSmiles("CCO");
  const core::Molecule imported = chem::fromMolBlock(chem::toMolBlock(ethanol));
  CHECK(chem::toSmiles(imported) == chem::canonicalize("CCO"));
}

TEST_CASE("periodic table helpers normalize symbols") {
  CHECK(std::string(chem::symbolFor(6)) == "C");
  CHECK(std::string(chem::symbolFor(17)) == "Cl");
  CHECK(chem::atomicNumberFor("Cl") == 17);
  CHECK(chem::atomicNumberFor("cl") == 17);
  CHECK(chem::atomicNumberFor("Xx") == 0);
}

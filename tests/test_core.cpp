#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "core/model.hpp"
#include "core/sprout.hpp"
#include "ui/camera.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <vector>

namespace {

using chemcad::core::Atom;
using chemcad::core::AtomId;
using chemcad::core::BondOrder;
using chemcad::core::Document;
using chemcad::core::Molecule;
using chemcad::core::Vec2;

constexpr float kPi = std::numbers::pi_v<float>;

AtomId addAtomAt(Molecule& molecule, float x, float y) {
  Atom atom;
  atom.pos = {x, y};
  return molecule.addAtom(atom);
}

float dot(Vec2 a, Vec2 b) { return a.x * b.x + a.y * b.y; }

float cross(Vec2 a, Vec2 b) { return a.x * b.y - a.y * b.x; }

Vec2 displacement(Vec2 from, Vec2 to) { return {to.x - from.x, to.y - from.y}; }

void checkUnit(Vec2 value, double epsilon = 1e-4) {
  CHECK(chemcad::core::length(value) == doctest::Approx(1.0).epsilon(epsilon));
}

void checkSameDocument(const Document& actual, const Document& expected) {
  REQUIRE(actual.molecules.size() == expected.molecules.size());
  for (size_t moleculeIndex = 0; moleculeIndex < actual.molecules.size(); ++moleculeIndex) {
    const Molecule& actualMolecule = actual.molecules[moleculeIndex];
    const Molecule& expectedMolecule = expected.molecules[moleculeIndex];
    REQUIRE(actualMolecule.atomCount() == expectedMolecule.atomCount());
    REQUIRE(actualMolecule.bondCount() == expectedMolecule.bondCount());
    for (size_t atomIndex = 0; atomIndex < actualMolecule.atomCount(); ++atomIndex) {
      CHECK(actualMolecule.atoms()[atomIndex].pos.x ==
            doctest::Approx(expectedMolecule.atoms()[atomIndex].pos.x));
      CHECK(actualMolecule.atoms()[atomIndex].pos.y ==
            doctest::Approx(expectedMolecule.atoms()[atomIndex].pos.y));
    }
  }
}

}  // namespace

TEST_CASE("Molecule maintains graph invariants") {
  Molecule molecule;
  const AtomId carbon = addAtomAt(molecule, 0.0f, 0.0f);
  const AtomId oxygen = addAtomAt(molecule, 1.0f, 0.0f);
  const AtomId nitrogen = addAtomAt(molecule, 0.0f, 1.0f);

  CHECK(carbon != chemcad::core::kInvalidAtom);
  CHECK(oxygen != chemcad::core::kInvalidAtom);
  CHECK(nitrogen != chemcad::core::kInvalidAtom);
  CHECK(carbon != oxygen);
  CHECK(carbon != nitrogen);
  CHECK(oxygen != nitrogen);

  const auto carbonOxygen = molecule.addBond(carbon, oxygen, BondOrder::Single);
  REQUIRE(carbonOxygen != chemcad::core::kInvalidBond);
  CHECK(molecule.bondCount() == 1);

  const auto retyped = molecule.addBond(oxygen, carbon, BondOrder::Double);
  CHECK(retyped == carbonOxygen);
  CHECK(molecule.bondCount() == 1);
  REQUIRE(molecule.bond(carbonOxygen) != nullptr);
  CHECK(molecule.bond(carbonOxygen)->order == BondOrder::Double);

  CHECK(molecule.addBond(carbon, carbon, BondOrder::Single) == chemcad::core::kInvalidBond);
  CHECK(molecule.addBond(carbon, 999999, BondOrder::Single) == chemcad::core::kInvalidBond);
  CHECK(molecule.addBond(999998, oxygen, BondOrder::Single) == chemcad::core::kInvalidBond);
  CHECK(molecule.addBond(chemcad::core::kInvalidAtom, oxygen, BondOrder::Single) ==
        chemcad::core::kInvalidBond);
  CHECK(molecule.bondBetween(oxygen, nitrogen) == chemcad::core::kInvalidBond);

  const auto carbonNitrogen = molecule.addBond(carbon, nitrogen, BondOrder::Single);
  REQUIRE(carbonNitrogen != chemcad::core::kInvalidBond);
  const std::vector<AtomId> neighbors = molecule.neighbors(carbon);
  const std::vector<chemcad::core::BondId> incident = molecule.incidentBonds(carbon);
  CHECK(neighbors.size() == 2);
  CHECK(incident.size() == neighbors.size());
  CHECK(molecule.degree(carbon) == static_cast<int>(neighbors.size()));
  CHECK(std::find(neighbors.begin(), neighbors.end(), oxygen) != neighbors.end());
  CHECK(std::find(neighbors.begin(), neighbors.end(), nitrogen) != neighbors.end());
  CHECK(std::find(incident.begin(), incident.end(), carbonOxygen) != incident.end());
  CHECK(std::find(incident.begin(), incident.end(), carbonNitrogen) != incident.end());

  molecule.removeBond(carbonNitrogen);
  CHECK(molecule.atomCount() == 3);
  CHECK(molecule.bondCount() == 1);
  CHECK(molecule.bondBetween(carbon, nitrogen) == chemcad::core::kInvalidBond);

  REQUIRE(molecule.addBond(carbon, nitrogen, BondOrder::Single) !=
          chemcad::core::kInvalidBond);
  molecule.removeAtom(carbon);
  CHECK(molecule.atomCount() == 2);
  CHECK(molecule.bondCount() == 0);
  CHECK(molecule.atom(oxygen) != nullptr);
  CHECK(molecule.atom(nitrogen) != nullptr);
  CHECK(molecule.neighbors(oxygen).empty());
  CHECK(molecule.incidentBonds(oxygen).empty());
  CHECK(molecule.degree(oxygen) == 0);
}

TEST_CASE("UndoStack restores snapshots and tracks branches") {
  Document document;
  document.molecules.emplace_back();
  Molecule& molecule = document.molecules.back();
  const AtomId first = addAtomAt(molecule, -1.0f, 2.0f);
  const AtomId second = addAtomAt(molecule, 1.0f, 2.0f);
  REQUIRE(molecule.addBond(first, second, BondOrder::Single) !=
          chemcad::core::kInvalidBond);

  const Document original = document;
  chemcad::core::UndoStack undo;
  CHECK_FALSE(undo.canUndo());
  CHECK_FALSE(undo.canRedo());

  undo.push(document);
  CHECK(undo.canUndo());
  CHECK_FALSE(undo.canRedo());
  molecule.atom(first)->pos = {7.0f, -3.0f};
  addAtomAt(molecule, 8.0f, -3.0f);
  const Document edited = document;

  CHECK(undo.undo(document));
  CHECK_FALSE(undo.canUndo());
  CHECK(undo.canRedo());
  checkSameDocument(document, original);

  CHECK(undo.redo(document));
  CHECK(undo.canUndo());
  CHECK_FALSE(undo.canRedo());
  checkSameDocument(document, edited);

  undo.clear();
  CHECK_FALSE(undo.canUndo());
  CHECK_FALSE(undo.canRedo());

  document = original;
  undo.push(document);
  document.molecules[0].atom(first)->pos.x = 3.0f;
  undo.push(document);
  document.molecules[0].atom(first)->pos.x = 4.0f;
  REQUIRE(undo.undo(document));
  CHECK(undo.canRedo());
  CHECK(document.molecules[0].atom(first)->pos.x == doctest::Approx(3.0));

  undo.push(document);
  CHECK_FALSE(undo.canRedo());
  document.molecules[0].atom(first)->pos.x = 5.0f;
  CHECK_FALSE(undo.redo(document));
  REQUIRE(undo.undo(document));
  CHECK(document.molecules[0].atom(first)->pos.x == doctest::Approx(3.0));
}

TEST_CASE("UndoStack bounds retained history") {
  Document document;
  document.molecules.emplace_back();
  const AtomId atom = addAtomAt(document.molecules[0], 0.0f, 0.0f);
  chemcad::core::UndoStack undo;

  for (size_t i = 0; i < chemcad::core::UndoStack::kCapacity + 10; ++i) {
    undo.push(document);
    document.molecules[0].atom(atom)->pos.x = static_cast<float>(i + 1);
  }

  size_t undoCount = 0;
  while (undo.undo(document)) ++undoCount;
  CHECK(undoCount == chemcad::core::UndoStack::kCapacity);
  CHECK_FALSE(undo.canUndo());
  CHECK(undo.canRedo());
  CHECK(document.molecules[0].atom(atom)->pos.x == doctest::Approx(10.0));
  CHECK_FALSE(undo.undo(document));
}

TEST_CASE("sproutDirection handles free, terminal, and crowded atoms") {
  Molecule freeMolecule;
  const AtomId freeAtom = addAtomAt(freeMolecule, 4.0f, -2.0f);
  const Vec2 freeDirection = chemcad::core::sproutDirection(freeMolecule, freeAtom);
  CHECK(freeDirection.x == doctest::Approx(std::cos(kPi / 6.0f)).epsilon(1e-4));
  CHECK(freeDirection.y == doctest::Approx(std::sin(kPi / 6.0f)).epsilon(1e-4));
  checkUnit(freeDirection);

  Molecule terminalMolecule;
  const AtomId neighbor = addAtomAt(terminalMolecule, 0.0f, 0.0f);
  const AtomId terminal = addAtomAt(terminalMolecule, 1.0f, 0.0f);
  REQUIRE(terminalMolecule.addBond(neighbor, terminal, BondOrder::Single) !=
          chemcad::core::kInvalidBond);
  const Vec2 terminalDirection = chemcad::core::sproutDirection(terminalMolecule, terminal);
  const Vec2 backToNeighbor = chemcad::core::normalize({-1.0f, 0.0f});
  const float cosine = std::clamp(dot(terminalDirection, backToNeighbor), -1.0f, 1.0f);
  const float angleDegrees = std::acos(cosine) * 180.0f / kPi;
  CHECK(angleDegrees == doctest::Approx(120.0).epsilon(1e-3));
  checkUnit(terminalDirection);

  Molecule crowdedMolecule;
  const AtomId center = addAtomAt(crowdedMolecule, 0.0f, 0.0f);
  const std::array<AtomId, 3> crowdedNeighbors = {
      addAtomAt(crowdedMolecule, 1.0f, 0.0f), addAtomAt(crowdedMolecule, 0.0f, 1.0f),
      addAtomAt(crowdedMolecule, -1.0f, 0.0f)};
  for (AtomId id : crowdedNeighbors) {
    REQUIRE(crowdedMolecule.addBond(center, id, BondOrder::Single) !=
            chemcad::core::kInvalidBond);
  }
  const Vec2 crowdedDirection = chemcad::core::sproutDirection(crowdedMolecule, center);
  CHECK(crowdedDirection.x == doctest::Approx(0.0).epsilon(1e-3));
  CHECK(crowdedDirection.y == doctest::Approx(-1.0).epsilon(1e-3));
  checkUnit(crowdedDirection);
}

TEST_CASE("repeated sprouts form an alternating zig-zag") {
  Molecule chain;
  AtomId tip = addAtomAt(chain, 0.0f, 0.0f);
  std::vector<Vec2> positions{{0.0f, 0.0f}};
  std::vector<Vec2> bondDirections;

  for (int i = 0; i < 3; ++i) {
    const Vec2 direction = chemcad::core::sproutDirection(chain, tip);
    checkUnit(direction);
    const Vec2 nextPosition = chemcad::core::sproutPosition(chain, tip);
    CHECK(displacement(positions.back(), nextPosition).x ==
          doctest::Approx(direction.x).epsilon(1e-4));
    CHECK(displacement(positions.back(), nextPosition).y ==
          doctest::Approx(direction.y).epsilon(1e-4));
    const AtomId next = addAtomAt(chain, nextPosition.x, nextPosition.y);
    REQUIRE(chain.addBond(tip, next, BondOrder::Single) != chemcad::core::kInvalidBond);
    positions.push_back(nextPosition);
    bondDirections.push_back(direction);
    tip = next;
  }

  REQUIRE(positions.size() == 4);
  REQUIRE(bondDirections.size() == 3);
  const float firstTurn = cross(bondDirections[0], bondDirections[1]);
  const float secondTurn = cross(bondDirections[1], bondDirections[2]);
  CHECK(std::abs(firstTurn) > 1e-3f);
  CHECK(std::abs(secondTurn) > 1e-3f);
  CHECK(firstTurn * secondTurn < 0.0f);
  CHECK(std::abs(cross(displacement(positions[0], positions[1]),
                       displacement(positions[0], positions[3]))) > 1e-3f);
}

TEST_CASE("snapAngle rounds to the nearest 15 degrees and normalizes") {
  const auto directionAt = [](float degrees, float magnitude = 1.0f) {
    return chemcad::core::fromAngle(degrees * kPi / 180.0f, magnitude);
  };
  const auto angleDegrees = [](Vec2 value) {
    return chemcad::core::angleOf(value) * 180.0f / kPi;
  };

  const Vec2 seventeen = chemcad::core::snapAngle(directionAt(17.0f));
  const Vec2 twentyTwo = chemcad::core::snapAngle(directionAt(22.0f));
  const Vec2 twentyThree = chemcad::core::snapAngle(directionAt(23.0f));
  const Vec2 negativeSeven = chemcad::core::snapAngle(directionAt(-7.0f));
  const Vec2 nonUnit = chemcad::core::snapAngle(directionAt(23.0f, 19.0f));

  CHECK(angleDegrees(seventeen) == doctest::Approx(15.0).epsilon(1e-4));
  CHECK(angleDegrees(twentyTwo) == doctest::Approx(15.0).epsilon(1e-4));
  CHECK(angleDegrees(twentyThree) == doctest::Approx(30.0).epsilon(1e-4));
  CHECK(angleDegrees(negativeSeven) == doctest::Approx(0.0).epsilon(1e-4));
  CHECK(angleDegrees(nonUnit) == doctest::Approx(30.0).epsilon(1e-4));
  checkUnit(seventeen);
  checkUnit(twentyTwo);
  checkUnit(twentyThree);
  checkUnit(negativeSeven);
  checkUnit(nonUnit);
}

TEST_CASE("Camera2D transforms are inverse and flip the canvas Y axis") {
  chemcad::ui::Camera2D camera;
  camera.pan = {-1.25f, 3.5f};
  const Vec2 origin{120.0f, 75.0f};
  const std::array<Vec2, 4> points = {
      Vec2{0.0f, 0.0f}, Vec2{-4.5f, 2.25f}, Vec2{100.0f, -50.0f}, Vec2{0.01f, 0.02f}};

  for (float zoom : {chemcad::ui::kMinZoom, 1.0f, 3.75f, chemcad::ui::kMaxZoom}) {
    camera.zoom = zoom;
    for (Vec2 point : points) {
      const Vec2 roundTrip = camera.screenToWorld(camera.worldToScreen(point, origin), origin);
      CHECK(roundTrip.x == doctest::Approx(point.x).epsilon(1e-3));
      CHECK(roundTrip.y == doctest::Approx(point.y).epsilon(1e-3));
    }
  }

  camera.zoom = 1.0f;
  const Vec2 lowerWorld = camera.worldToScreen({2.0f, -1.0f}, origin);
  const Vec2 higherWorld = camera.worldToScreen({2.0f, 1.0f}, origin);
  CHECK(higherWorld.y < lowerWorld.y);
}

TEST_CASE("Camera2D zoomAt preserves the point under the cursor") {
  chemcad::ui::Camera2D camera;
  camera.pan = {-3.0f, 2.0f};
  camera.zoom = 0.8f;
  const Vec2 origin{45.0f, 80.0f};
  const Vec2 cursor{627.0f, 419.0f};
  const Vec2 before = camera.screenToWorld(cursor, origin);

  camera.zoomAt(2.25f, cursor, origin);
  const Vec2 after = camera.screenToWorld(cursor, origin);
  CHECK(after.x == doctest::Approx(before.x).epsilon(1e-4));
  CHECK(after.y == doctest::Approx(before.y).epsilon(1e-4));
  CHECK(camera.zoom == doctest::Approx(1.8).epsilon(1e-5));
}

TEST_CASE("Camera2D fit contains atoms and handles an empty document") {
  Document document;
  document.molecules.emplace_back();
  addAtomAt(document.molecules[0], -3.0f, -2.0f);
  addAtomAt(document.molecules[0], 4.0f, 5.0f);
  addAtomAt(document.molecules[0], 1.5f, -1.0f);

  chemcad::ui::Camera2D camera;
  const Vec2 viewport{800.0f, 600.0f};
  camera.fit(document, viewport);
  CHECK(camera.zoom >= chemcad::ui::kMinZoom);
  CHECK(camera.zoom <= chemcad::ui::kMaxZoom);
  CHECK(std::isfinite(camera.zoom));
  for (const Molecule& molecule : document.molecules) {
    for (const Atom& atom : molecule.atoms()) {
      const Vec2 screen = camera.worldToScreen(atom.pos, {0.0f, 0.0f});
      CHECK(screen.x >= 0.0f);
      CHECK(screen.x <= viewport.x);
      CHECK(screen.y >= 0.0f);
      CHECK(screen.y <= viewport.y);
    }
  }

  Document empty;
  camera.pan = {99.0f, -99.0f};
  camera.zoom = 0.0f;
  camera.fit(empty, viewport);
  CHECK(std::isfinite(camera.zoom));
  CHECK(camera.zoom >= chemcad::ui::kMinZoom);
  CHECK(camera.zoom <= chemcad::ui::kMaxZoom);
  CHECK(std::isfinite(camera.pan.x));
  CHECK(std::isfinite(camera.pan.y));
}

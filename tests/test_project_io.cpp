#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <nlohmann/json.hpp>

#include "app/project_io.hpp"

namespace {

using chemcad::app::Project;
using chemcad::core::Atom;
using chemcad::core::AtomId;
using chemcad::core::BondOrder;
using chemcad::core::BondStereo;
using chemcad::core::Molecule;
using BondSignature = std::tuple<std::size_t, std::size_t, int, int>;

std::filesystem::path temporaryPath(const std::string& suffix) {
  static unsigned counter = 0;
  const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
  return std::filesystem::temp_directory_path() /
         ("chemcad-project-io-" + std::to_string(stamp) + "-" +
          std::to_string(counter++) + suffix);
}

struct RemoveOnExit {
  std::vector<std::filesystem::path> paths;

  ~RemoveOnExit() {
    for (const auto& path : paths) {
      std::error_code ignored;
      std::filesystem::remove(path, ignored);
      std::filesystem::remove(path.string() + ".tmp", ignored);
    }
  }
};

std::vector<BondSignature> bondSignatures(const Molecule& molecule) {
  std::unordered_map<AtomId, std::size_t> atomIndices;
  for (std::size_t index = 0; index < molecule.atoms().size(); ++index) {
    atomIndices.emplace(molecule.atoms()[index].id, index);
  }

  std::vector<BondSignature> result;
  for (const auto& bond : molecule.bonds()) {
    std::size_t a = atomIndices.at(bond.a);
    std::size_t b = atomIndices.at(bond.b);
    if (a > b) {
      std::swap(a, b);
    }
    result.emplace_back(a, b, static_cast<int>(bond.order), static_cast<int>(bond.stereo));
  }
  std::sort(result.begin(), result.end());
  return result;
}

void checkProjectsMatch(const Project& expected, const Project& actual) {
  REQUIRE(actual.doc.molecules.size() == expected.doc.molecules.size());
  for (std::size_t moleculeIndex = 0; moleculeIndex < expected.doc.molecules.size();
       ++moleculeIndex) {
    const Molecule& expectedMolecule = expected.doc.molecules[moleculeIndex];
    const Molecule& actualMolecule = actual.doc.molecules[moleculeIndex];
    REQUIRE(actualMolecule.atomCount() == expectedMolecule.atomCount());
    REQUIRE(actualMolecule.bondCount() == expectedMolecule.bondCount());

    for (std::size_t atomIndex = 0; atomIndex < expectedMolecule.atoms().size(); ++atomIndex) {
      const Atom& expectedAtom = expectedMolecule.atoms()[atomIndex];
      const Atom& actualAtom = actualMolecule.atoms()[atomIndex];
      CHECK(actualAtom.atomicNumber == expectedAtom.atomicNumber);
      CHECK(actualAtom.charge == expectedAtom.charge);
      CHECK(actualAtom.isotope == expectedAtom.isotope);
      CHECK(actualAtom.explicitH == expectedAtom.explicitH);
      CHECK(actualAtom.pos.x == expectedAtom.pos.x);
      CHECK(actualAtom.pos.y == expectedAtom.pos.y);
    }
    CHECK(bondSignatures(actualMolecule) == bondSignatures(expectedMolecule));
  }
  CHECK(actual.plannerStarts == expected.plannerStarts);
  CHECK(actual.plannerTarget == expected.plannerTarget);
}

Project exampleProject() {
  Project project;

  Molecule labelledFragment;
  Atom labelled;
  labelled.atomicNumber = 7;
  labelled.charge = 1;
  labelled.isotope = 15;
  labelled.explicitH = 2;
  labelled.pos = {-1.25F, 0.375F};
  const AtomId labelledId = labelledFragment.addAtom(labelled);

  Atom removed;
  removed.atomicNumber = 6;
  removed.pos = {50.0F, 50.0F};
  const AtomId removedId = labelledFragment.addAtom(removed);
  labelledFragment.removeAtom(removedId);

  Atom oxygen;
  oxygen.atomicNumber = 8;
  oxygen.pos = {0.125F, -2.75F};
  const AtomId oxygenId = labelledFragment.addAtom(oxygen);
  const auto wedgeId = labelledFragment.addBond(labelledId, oxygenId, BondOrder::Single);
  labelledFragment.bond(wedgeId)->stereo = BondStereo::Wedge;
  // A wavy bond (unknown stereochemistry) must survive the round-trip too.
  Atom nitrogen;
  nitrogen.atomicNumber = 7;
  nitrogen.pos = {1.5F, -2.75F};
  const AtomId nitrogenId = labelledFragment.addAtom(nitrogen);
  const auto wavyId = labelledFragment.addBond(oxygenId, nitrogenId, BondOrder::Single);
  labelledFragment.bond(wavyId)->stereo = BondStereo::Wavy;
  REQUIRE(labelledFragment.atoms()[1].id > labelledFragment.atoms()[0].id + 1);
  project.doc.molecules.push_back(std::move(labelledFragment));

  Molecule ringFragment;
  std::vector<AtomId> ringAtoms;
  for (const chemcad::core::Vec2 position :
       {chemcad::core::Vec2{0.0F, 1.0F}, chemcad::core::Vec2{1.0F, 0.0F},
        chemcad::core::Vec2{0.0F, -1.0F}, chemcad::core::Vec2{-1.0F, 0.0F}}) {
    Atom carbon;
    carbon.pos = position;
    ringAtoms.push_back(ringFragment.addAtom(carbon));
  }
  ringFragment.addBond(ringAtoms[0], ringAtoms[1], BondOrder::Single);
  ringFragment.addBond(ringAtoms[1], ringAtoms[2], BondOrder::Double);
  ringFragment.addBond(ringAtoms[2], ringAtoms[3], BondOrder::Single);
  ringFragment.addBond(ringAtoms[3], ringAtoms[0], BondOrder::Double);
  project.doc.molecules.push_back(std::move(ringFragment));

  project.plannerStarts = {"[15NH2+]O", "C1=CCC=C1"};
  project.plannerTarget = "CC(=O)O";
  return project;
}

}  // namespace

TEST_CASE("project JSON round-trip uses positional atom indices") {
  const Project expected = exampleProject();
  const std::string serialized = chemcad::app::serializeProject(expected);
  const Project actual = chemcad::app::deserializeProject(serialized);
  checkProjectsMatch(expected, actual);

  const nlohmann::json json = nlohmann::json::parse(serialized);
  CHECK(json.at("molecules").at(0).at("bonds").at(0).at("a") == 0);
  CHECK(json.at("molecules").at(0).at("bonds").at(0).at("b") == 1);
}

TEST_CASE("project file save is valid JSON and loads identically") {
  const Project expected = exampleProject();
  const std::filesystem::path path = temporaryPath(".chemcad");
  RemoveOnExit cleanup{{path}};

  chemcad::app::saveProject(expected, path.string());
  CHECK_FALSE(std::filesystem::exists(path.string() + ".tmp"));

  std::ifstream input(path, std::ios::binary);
  REQUIRE(input.good());
  const nlohmann::json saved = nlohmann::json::parse(input);
  CHECK(saved.at("version") == 1);

  const Project actual = chemcad::app::loadProject(path.string());
  checkProjectsMatch(expected, actual);
}

TEST_CASE("project format reports actionable errors") {
  const std::filesystem::path missing = temporaryPath("-missing.chemcad");
  const std::filesystem::path noVersion = temporaryPath("-no-version.chemcad");
  RemoveOnExit cleanup{{missing, noVersion}};
  {
    std::ofstream output(noVersion);
    output << "{}";
  }
  const std::string missingError = "cannot open " + missing.string();

  CHECK_THROWS_WITH_AS(chemcad::app::loadProject(missing.string()),
                       missingError.c_str(), std::runtime_error);
  CHECK_THROWS_WITH_AS(chemcad::app::loadProject(noVersion.string()),
                       "not a chemcad project (missing \"version\")", std::runtime_error);
  CHECK_THROWS_WITH_AS(chemcad::app::deserializeProject(R"({"version":999})"),
                       "unsupported project version 999", std::runtime_error);

  const std::string badBond = R"({
    "version": 1,
    "molecules": [{
      "atoms": [{"z": 6, "x": 0.0, "y": 0.0, "q": 0, "iso": 0, "h": -1}],
      "bonds": [{"a": 0, "b": 1, "o": 1, "s": 0}]
    }],
    "planner": {"starts": [], "target": ""}
  })";
  CHECK_THROWS_WITH_AS(chemcad::app::deserializeProject(badBond),
                       "malformed project: bond index out of range", std::runtime_error);
}

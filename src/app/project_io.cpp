#include "app/project_io.hpp"

#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <system_error>
#include <unordered_map>

#include <nlohmann/json.hpp>

namespace chemcad::app {
namespace {

using json = nlohmann::json;

template <class T>
T checkedInteger(const json& value, const char* field) {
  if (!value.is_number_integer() && !value.is_number_unsigned()) {
    throw std::runtime_error(std::string("malformed project: ") + field +
                             " must be an integer");
  }

  const auto number = value.get<std::int64_t>();
  if (number < static_cast<std::int64_t>(std::numeric_limits<T>::min()) ||
      number > static_cast<std::int64_t>(std::numeric_limits<T>::max())) {
    throw std::runtime_error(std::string("malformed project: ") + field +
                             " is out of range");
  }
  return static_cast<T>(number);
}

std::size_t bondIndex(const json& value, const char* field) {
  if (!value.is_number_unsigned()) {
    throw std::runtime_error(std::string("malformed project: bond ") + field +
                             " must be a non-negative integer");
  }
  return value.get<std::size_t>();
}

void removeTemporaryFile(const std::filesystem::path& path) {
  std::error_code ignored;
  std::filesystem::remove(path, ignored);
}

}  // namespace

void saveProject(const Project& project, const std::string& path) {
  const std::filesystem::path target(path);
  const std::filesystem::path temporary(path + ".tmp");
  const std::string contents = serializeProject(project);

  std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
  if (!output) {
    throw std::runtime_error("cannot open " + path);
  }

  output.write(contents.data(), static_cast<std::streamsize>(contents.size()));
  output.close();
  if (!output) {
    removeTemporaryFile(temporary);
    throw std::runtime_error("cannot write " + path);
  }

  std::error_code error;
  std::filesystem::rename(temporary, target, error);
  if (error) {
    removeTemporaryFile(temporary);
    throw std::runtime_error("cannot replace " + path + ": " + error.message());
  }
}

Project loadProject(const std::string& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::runtime_error("cannot open " + path);
  }

  std::string contents((std::istreambuf_iterator<char>(input)),
                       std::istreambuf_iterator<char>());
  if (input.bad()) {
    throw std::runtime_error("cannot read " + path);
  }
  return deserializeProject(contents);
}

std::string serializeProject(const Project& project) {
  json molecules = json::array();
  for (const core::Molecule& molecule : project.doc.molecules) {
    json atoms = json::array();
    std::unordered_map<core::AtomId, std::size_t> atomIndices;
    atomIndices.reserve(molecule.atomCount());

    for (const core::Atom& atom : molecule.atoms()) {
      atomIndices.emplace(atom.id, atomIndices.size());
      atoms.push_back({{"z", atom.atomicNumber},
                       {"x", atom.pos.x},
                       {"y", atom.pos.y},
                       {"q", atom.charge},
                       {"iso", atom.isotope},
                       {"h", atom.explicitH}});
    }

    json bonds = json::array();
    for (const core::Bond& bond : molecule.bonds()) {
      const auto a = atomIndices.find(bond.a);
      const auto b = atomIndices.find(bond.b);
      if (a == atomIndices.end() || b == atomIndices.end()) {
        throw std::runtime_error("cannot serialize project: bond references a missing atom");
      }
      bonds.push_back({{"a", a->second},
                       {"b", b->second},
                       {"o", static_cast<int>(bond.order)},
                       {"s", static_cast<int>(bond.stereo)}});
    }

    molecules.push_back({{"atoms", std::move(atoms)}, {"bonds", std::move(bonds)}});
  }

  const json root = {{"version", 1},
                     {"molecules", std::move(molecules)},
                     {"planner", {{"starts", project.plannerStarts},
                                  {"target", project.plannerTarget}}}};
  return root.dump(2);
}

Project deserializeProject(const std::string& text) {
  json root;
  try {
    root = json::parse(text);
  } catch (const json::exception& error) {
    throw std::runtime_error(std::string("malformed project: ") + error.what());
  }

  if (!root.is_object() || !root.contains("version")) {
    throw std::runtime_error("not a chemcad project (missing \"version\")");
  }

  int version = 0;
  try {
    version = root.at("version").get<int>();
  } catch (const json::exception& error) {
    throw std::runtime_error(std::string("malformed project: ") + error.what());
  }
  if (version != 1) {
    throw std::runtime_error("unsupported project version " + std::to_string(version));
  }

  Project project;
  try {
    const json& serializedMolecules = root.at("molecules");
    if (!serializedMolecules.is_array()) {
      throw std::runtime_error("malformed project: molecules must be an array");
    }

    project.doc.molecules.reserve(serializedMolecules.size());
    for (const json& serializedMolecule : serializedMolecules) {
      const json& serializedAtoms = serializedMolecule.at("atoms");
      const json& serializedBonds = serializedMolecule.at("bonds");
      if (!serializedAtoms.is_array() || !serializedBonds.is_array()) {
        throw std::runtime_error("malformed project: atoms and bonds must be arrays");
      }

      core::Molecule molecule;
      std::vector<core::AtomId> atomIds;
      atomIds.reserve(serializedAtoms.size());
      for (const json& serializedAtom : serializedAtoms) {
        core::Atom atom;
        atom.atomicNumber = checkedInteger<std::uint8_t>(serializedAtom.at("z"), "z");
        atom.pos.x = serializedAtom.at("x").get<float>();
        atom.pos.y = serializedAtom.at("y").get<float>();
        atom.charge = checkedInteger<std::int8_t>(serializedAtom.at("q"), "q");
        atom.isotope = checkedInteger<std::uint16_t>(serializedAtom.at("iso"), "iso");
        atom.explicitH = checkedInteger<std::int8_t>(serializedAtom.at("h"), "h");
        atomIds.push_back(molecule.addAtom(atom));
      }

      for (const json& serializedBond : serializedBonds) {
        const std::size_t a = bondIndex(serializedBond.at("a"), "a");
        const std::size_t b = bondIndex(serializedBond.at("b"), "b");
        if (a >= atomIds.size() || b >= atomIds.size()) {
          throw std::runtime_error("malformed project: bond index out of range");
        }

        const auto orderValue = checkedInteger<std::uint8_t>(serializedBond.at("o"), "o");
        const auto stereoValue = checkedInteger<std::uint8_t>(serializedBond.at("s"), "s");
        if (orderValue < static_cast<std::uint8_t>(core::BondOrder::Single) ||
            orderValue > static_cast<std::uint8_t>(core::BondOrder::Aromatic)) {
          throw std::runtime_error("malformed project: bond order is out of range");
        }
        if (stereoValue > static_cast<std::uint8_t>(core::BondStereo::Wavy)) {
          throw std::runtime_error("malformed project: bond stereo is out of range");
        }

        const core::BondId bondId = molecule.addBond(
            atomIds[a], atomIds[b], static_cast<core::BondOrder>(orderValue));
        molecule.bond(bondId)->stereo = static_cast<core::BondStereo>(stereoValue);
      }
      project.doc.molecules.push_back(std::move(molecule));
    }

    const json& planner = root.at("planner");
    project.plannerStarts = planner.at("starts").get<std::vector<std::string>>();
    project.plannerTarget = planner.at("target").get<std::string>();
  } catch (const json::exception& error) {
    throw std::runtime_error(std::string("malformed project: ") + error.what());
  }

  return project;
}

}  // namespace chemcad::app

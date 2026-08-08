#include "chem/bridge.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <exception>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Geometry/point.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/FileWriters.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/SanitException.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace chemcad::chem {
namespace {

using AtomIndexMap = std::unordered_map<core::AtomId, unsigned int>;

RDKit::Bond::BondType toRDKitBondType(core::BondOrder order) {
  switch (order) {
    case core::BondOrder::Single:
      return RDKit::Bond::SINGLE;
    case core::BondOrder::Double:
      return RDKit::Bond::DOUBLE;
    case core::BondOrder::Triple:
      return RDKit::Bond::TRIPLE;
    case core::BondOrder::Aromatic:
      return RDKit::Bond::AROMATIC;
  }
  return RDKit::Bond::SINGLE;
}

core::BondOrder fromRDKitBondType(const RDKit::Bond& bond) {
  if (bond.getIsAromatic() || bond.getBondType() == RDKit::Bond::AROMATIC) {
    return core::BondOrder::Aromatic;
  }
  switch (bond.getBondType()) {
    case RDKit::Bond::DOUBLE:
      return core::BondOrder::Double;
    case RDKit::Bond::TRIPLE:
      return core::BondOrder::Triple;
    default:
      return core::BondOrder::Single;
  }
}

RDKit::RWMol toRDKit(const core::Molecule& source, AtomIndexMap* indexMap = nullptr,
                     std::vector<core::AtomId>* atomIds = nullptr) {
  RDKit::RWMol result;
  AtomIndexMap localMap;
  localMap.reserve(source.atomCount());
  if (atomIds) {
    atomIds->clear();
    atomIds->reserve(source.atomCount());
  }

  for (const core::Atom& sourceAtom : source.atoms()) {
    auto atom = std::make_unique<RDKit::Atom>(sourceAtom.atomicNumber);
    atom->setFormalCharge(sourceAtom.charge);
    if (sourceAtom.isotope != 0) {
      atom->setIsotope(sourceAtom.isotope);
    }
    if (sourceAtom.explicitH >= 0) {
      atom->setNoImplicit(true);
      atom->setNumExplicitHs(static_cast<unsigned int>(sourceAtom.explicitH));
    }
    const unsigned int index = result.addAtom(atom.release(), true, true);
    localMap.emplace(sourceAtom.id, index);
    if (atomIds) {
      atomIds->push_back(sourceAtom.id);
    }
  }

  for (const core::Bond& sourceBond : source.bonds()) {
    const auto begin = localMap.find(sourceBond.a);
    const auto end = localMap.find(sourceBond.b);
    if (begin == localMap.end() || end == localMap.end()) {
      throw ChemError("bond refers to an unknown atom");
    }
    result.addBond(begin->second, end->second, toRDKitBondType(sourceBond.order));
    RDKit::Bond* bond = result.getBondBetweenAtoms(begin->second, end->second);
    if (sourceBond.order == core::BondOrder::Aromatic) {
      bond->setIsAromatic(true);
      result.getAtomWithIdx(begin->second)->setIsAromatic(true);
      result.getAtomWithIdx(end->second)->setIsAromatic(true);
    }
    if (sourceBond.stereo == core::BondStereo::Wedge) {
      bond->setBondDir(RDKit::Bond::BEGINWEDGE);
    } else if (sourceBond.stereo == core::BondStereo::Hash) {
      bond->setBondDir(RDKit::Bond::BEGINDASH);
    } else if (sourceBond.stereo == core::BondStereo::Wavy) {
      bond->setBondDir(RDKit::Bond::UNKNOWN);  // MDL "either": renders wavy
    }
  }

  if (result.getNumAtoms() != 0) {
    auto conformer = std::make_unique<RDKit::Conformer>(result.getNumAtoms());
    conformer->set3D(false);
    for (const core::Atom& sourceAtom : source.atoms()) {
      const unsigned int index = localMap.at(sourceAtom.id);
      conformer->setAtomPos(index, RDGeom::Point3D(sourceAtom.pos.x, sourceAtom.pos.y, 0.0));
    }
    result.addConformer(conformer.release(), true);
  }

  if (indexMap) {
    *indexMap = std::move(localMap);
  }
  return result;
}

core::Molecule fromRDKit(const RDKit::ROMol& source) {
  core::Molecule result;
  std::vector<core::AtomId> atomIds;
  atomIds.reserve(source.getNumAtoms());
  const RDKit::Conformer* conformer = source.getNumConformers() == 0 ? nullptr : &source.getConformer();

  for (const RDKit::Atom* sourceAtom : source.atoms()) {
    core::Atom atom;
    atom.atomicNumber = static_cast<uint8_t>(sourceAtom->getAtomicNum());
    atom.charge = static_cast<int8_t>(sourceAtom->getFormalCharge());
    atom.isotope = static_cast<uint16_t>(sourceAtom->getIsotope());
    if (sourceAtom->getNoImplicit() || sourceAtom->getNumExplicitHs() != 0) {
      atom.explicitH = static_cast<int8_t>(sourceAtom->getNumExplicitHs());
    }
    if (conformer) {
      const RDGeom::Point3D& position = conformer->getAtomPos(sourceAtom->getIdx());
      atom.pos = {static_cast<float>(position.x), static_cast<float>(position.y)};
    }
    atomIds.push_back(result.addAtom(atom));
  }

  for (const RDKit::Bond* sourceBond : source.bonds()) {
    const core::BondId id = result.addBond(atomIds.at(sourceBond->getBeginAtomIdx()),
                                           atomIds.at(sourceBond->getEndAtomIdx()),
                                           fromRDKitBondType(*sourceBond));
    core::Bond* bond = result.bond(id);
    if (sourceBond->getBondDir() == RDKit::Bond::BEGINWEDGE) {
      bond->stereo = core::BondStereo::Wedge;
    } else if (sourceBond->getBondDir() == RDKit::Bond::BEGINDASH) {
      bond->stereo = core::BondStereo::Hash;
    } else if (sourceBond->getBondDir() == RDKit::Bond::UNKNOWN) {
      bond->stereo = core::BondStereo::Wavy;
    }
  }
  return result;
}

std::string aromaticError() {
  return "cannot kekulize aromatic ring - check the ring is a valid aromatic system";
}

int explicitValenceForMessage(const RDKit::ROMol& molecule, unsigned int atomIndex) {
  const RDKit::Atom* atom = molecule.getAtomWithIdx(atomIndex);
  double valence = atom->getNumExplicitHs();
  for (const RDKit::Bond* bond : molecule.atomBonds(atom)) {
    valence += bond->getBondTypeAsDouble();
  }
  return static_cast<int>(std::lround(valence));
}

void sanitize(RDKit::RWMol& molecule) {
  try {
    RDKit::MolOps::sanitizeMol(molecule);
  } catch (const RDKit::KekulizeException&) {
    throw ChemError(aromaticError());
  } catch (const RDKit::AtomKekulizeException&) {
    throw ChemError(aromaticError());
  } catch (const RDKit::AtomValenceException& error) {
    const unsigned int atomIndex = error.getAtomIdx();
    if (atomIndex < molecule.getNumAtoms()) {
      const RDKit::Atom* atom = molecule.getAtomWithIdx(atomIndex);
      std::string symbol = RDKit::PeriodicTable::getTable()->getElementSymbol(atom->getAtomicNum());
      throw ChemError(symbol + " has too many bonds (valence " +
                      std::to_string(explicitValenceForMessage(molecule, atomIndex)) + ")");
    }
    throw ChemError("atom has too many bonds");
  } catch (const RDKit::MolSanitizeException& error) {
    throw ChemError(error.what());
  } catch (const std::exception& error) {
    throw ChemError(error.what());
  }
}

template <typename Function>
decltype(auto) translateErrors(Function&& function) {
  try {
    return std::forward<Function>(function)();
  } catch (const ChemError&) {
    throw;
  } catch (const RDKit::KekulizeException&) {
    throw ChemError(aromaticError());
  } catch (const RDKit::AtomKekulizeException&) {
    throw ChemError(aromaticError());
  } catch (const std::exception& error) {
    throw ChemError(error.what());
  }
}

std::unique_ptr<RDKit::RWMol> parseSmiles(const std::string& smiles) {
  try {
    std::unique_ptr<RDKit::RWMol> molecule(RDKit::SmilesToMol(smiles, 0, false));
    if (!molecule) {
      throw ChemError("could not parse SMILES");
    }
    sanitize(*molecule);
    return molecule;
  } catch (const ChemError&) {
    throw;
  } catch (const std::exception&) {
    throw ChemError("could not parse SMILES");
  }
}

}  // namespace

std::string toSmiles(const core::Molecule& molecule) {
  if (molecule.empty()) {
    return {};
  }
  return translateErrors([&] {
    RDKit::RWMol rdMolecule = toRDKit(molecule);
    sanitize(rdMolecule);
    return RDKit::MolToSmiles(rdMolecule, true, false, -1, true);
  });
}

core::Molecule fromSmiles(const std::string& smiles) {
  return translateErrors([&] {
    std::unique_ptr<RDKit::RWMol> rdMolecule = parseSmiles(smiles);
    core::Molecule result = fromRDKit(*rdMolecule);
    layout(result);
    return result;
  });
}

std::string toMolBlock(const core::Molecule& molecule) {
  return translateErrors([&] {
    RDKit::RWMol rdMolecule = toRDKit(molecule);
    sanitize(rdMolecule);
    return RDKit::MolToMolBlock(rdMolecule, true, -1, true, false);
  });
}

core::Molecule fromMolBlock(const std::string& molBlock) {
  if (molBlock.empty()) {
    throw ChemError("could not parse MolBlock");
  }
  return translateErrors([&] {
    std::unique_ptr<RDKit::RWMol> rdMolecule(
        RDKit::MolBlockToMol(molBlock, false, false, true));
    if (!rdMolecule) {
      throw ChemError("could not parse MolBlock");
    }
    sanitize(*rdMolecule);
    return fromRDKit(*rdMolecule);
  });
}

void layout(core::Molecule& molecule) {
  if (molecule.atomCount() <= 1) {
    return;
  }
  translateErrors([&] {
    AtomIndexMap indexMap;
    RDKit::RWMol rdMolecule = toRDKit(molecule, &indexMap);
    sanitize(rdMolecule);
    RDDepict::preferCoordGen = true;
    RDDepict::compute2DCoords(rdMolecule);

    const RDKit::Conformer& conformer = rdMolecule.getConformer();
    for (core::Atom& atom : molecule.mutableAtoms()) {
      const RDGeom::Point3D& position = conformer.getAtomPos(indexMap.at(atom.id));
      atom.pos = {static_cast<float>(position.x), static_cast<float>(position.y)};
    }

    if (molecule.bondCount() == 0) {
      return;
    }
    double lengthSum = 0.0;
    std::size_t lengthCount = 0;
    for (const core::Bond& bond : molecule.bonds()) {
      const core::Atom* begin = molecule.atom(bond.a);
      const core::Atom* end = molecule.atom(bond.b);
      if (!begin || !end) {
        continue;
      }
      lengthSum += std::hypot(static_cast<double>(end->pos.x - begin->pos.x),
                              static_cast<double>(end->pos.y - begin->pos.y));
      ++lengthCount;
    }
    if (lengthCount == 0 || lengthSum == 0.0) {
      return;
    }

    const double scale = static_cast<double>(lengthCount) / lengthSum;
    core::Vec2 center{};
    for (const core::Atom& atom : molecule.atoms()) {
      center.x += atom.pos.x;
      center.y += atom.pos.y;
    }
    center.x /= static_cast<float>(molecule.atomCount());
    center.y /= static_cast<float>(molecule.atomCount());
    for (core::Atom& atom : molecule.mutableAtoms()) {
      atom.pos.x = center.x + static_cast<float>((atom.pos.x - center.x) * scale);
      atom.pos.y = center.y + static_cast<float>((atom.pos.y - center.y) * scale);
    }
  });
}

Properties computeProperties(const core::Molecule& molecule) {
  if (molecule.empty()) {
    return {};
  }
  return translateErrors([&] {
    RDKit::RWMol rdMolecule = toRDKit(molecule);
    sanitize(rdMolecule);
    RDKit::MolOps::symmetrizeSSSR(rdMolecule);

    Properties result;
    result.formula = RDKit::Descriptors::calcMolFormula(rdMolecule);
    result.mw = RDKit::Descriptors::calcExactMW(rdMolecule);
    double molarRefractivity = 0.0;
    RDKit::Descriptors::calcCrippenDescriptors(rdMolecule, result.logP, molarRefractivity);
    result.rings = static_cast<int>(rdMolecule.getRingInfo()->numRings());
    return result;
  });
}

std::string canonicalize(const std::string& smiles) {
  if (smiles.empty()) {
    return {};
  }
  return translateErrors([&] {
    std::unique_ptr<RDKit::RWMol> molecule = parseSmiles(smiles);
    return RDKit::MolToSmiles(*molecule, true, false, -1, true);
  });
}

std::string toSvg(const core::Molecule& molecule, int width, int height) {
  return translateErrors([&] {
    RDKit::RWMol rdMolecule = toRDKit(molecule);
    sanitize(rdMolecule);
    RDKit::MolDraw2DUtils::prepareMolForDrawing(rdMolecule);
    RDKit::MolDraw2DSVG drawer(width, height);
    drawer.drawMolecule(rdMolecule);
    drawer.finishDrawing();
    return drawer.getDrawingText();
  });
}

int implicitHCount(const core::Molecule& molecule, core::AtomId atomId) {
  if (!molecule.atom(atomId)) {
    return 0;
  }
  try {
    AtomIndexMap indexMap;
    RDKit::RWMol rdMolecule = toRDKit(molecule, &indexMap);
    const auto found = indexMap.find(atomId);
    if (found == indexMap.end()) {
      return 0;
    }
    sanitize(rdMolecule);
    return static_cast<int>(rdMolecule.getAtomWithIdx(found->second)->getTotalNumHs());
  } catch (...) {
    return 0;
  }
}

const char* symbolFor(uint8_t atomicNumber) {
  static const std::array<std::string, 256> symbols = [] {
    std::array<std::string, 256> values;
    try {
      const RDKit::PeriodicTable* table = RDKit::PeriodicTable::getTable();
      for (unsigned int number = 0; number <= 118; ++number) {
        try {
          values[number] = table->getElementSymbol(number);
        } catch (...) {
          values[number].clear();
        }
      }
    } catch (...) {
    }
    return values;
  }();
  return symbols[atomicNumber].c_str();
}

uint8_t atomicNumberFor(const std::string& symbol) {
  if (symbol.empty()) {
    return 0;
  }
  std::string normalized = symbol;
  normalized[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(normalized[0])));
  std::transform(normalized.begin() + 1, normalized.end(), normalized.begin() + 1,
                 [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
  try {
    const int atomicNumber = RDKit::PeriodicTable::getTable()->getAtomicNumber(normalized);
    if (atomicNumber <= 0 || atomicNumber > 255) {
      return 0;
    }
    return static_cast<uint8_t>(atomicNumber);
  } catch (...) {
    return 0;
  }
}

}  // namespace chemcad::chem

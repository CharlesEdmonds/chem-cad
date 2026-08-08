// 3D conformer generation for the molecule viewer. Round-trips through
// canonical SMILES (deterministic, loses nothing the viewer renders), adds
// explicit hydrogens so the force field sees the full valence, embeds with
// fixed-seed ETKDG v3 and relaxes with MMFF94.

#include "chem/embed3d.hpp"

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <algorithm>
#include <cmath>
#include <memory>

namespace chemcad::chem {

Embedded3D embed3D(const core::Molecule& molecule) {
  if (molecule.empty()) throw ChemError("Nothing to display in 3D.");

  const std::string smiles = toSmiles(molecule);
  std::unique_ptr<RDKit::ROMol> parsed(RDKit::SmilesToMol(smiles));
  if (!parsed) throw ChemError("Could not rebuild molecule for 3D display.");
  RDKit::RWMol working = *RDKit::MolOps::addHs(*parsed);

  // Fixed seed: the viewer must not jump to a new conformer every time the
  // sketch changes elsewhere and back.
  RDKit::DGeomHelpers::EmbedParameters params = RDKit::DGeomHelpers::ETKDGv3;
  params.randomSeed = 0xC0FFEE;
  if (RDKit::DGeomHelpers::EmbedMolecule(working, params) < 0) {
    throw ChemError("3D embedding failed for this molecule.");
  }

  // MMFF can refuse exotic valences where ETKDG still produces sane
  // geometry; keep the embedded coordinates in that case rather than fail
  // the whole view.
  try {
    RDKit::MMFF::MMFFOptimizeMolecule(working, 200);
  } catch (...) {
  }

  const RDKit::Conformer& conformer = working.getConformer();

  // Readout keeps only non-hydrogen atoms, remapping indices as it goes.
  Embedded3D result;
  std::vector<int> indexMap(working.getNumAtoms(), -1);
  for (unsigned i = 0; i < working.getNumAtoms(); ++i) {
    const RDKit::Atom* atom = working.getAtomWithIdx(i);
    if (atom->getAtomicNum() <= 1) continue;
    const RDGeom::Point3D& p = conformer.getAtomPos(i);
    Atom3D out;
    out.atomicNumber = static_cast<uint8_t>(atom->getAtomicNum());
    out.x = static_cast<float>(p.x);
    out.y = static_cast<float>(p.y);
    out.z = static_cast<float>(p.z);
    indexMap[i] = static_cast<int>(result.atoms.size());
    result.atoms.push_back(out);
  }
  for (const RDKit::Bond* bond : working.bonds()) {
    const int a = indexMap[bond->getBeginAtomIdx()];
    const int b = indexMap[bond->getEndAtomIdx()];
    if (a < 0 || b < 0) continue;
    Bond3D out;
    out.a = a;
    out.b = b;
    switch (bond->getBondType()) {
      case RDKit::Bond::DOUBLE: out.order = 2; break;
      case RDKit::Bond::TRIPLE: out.order = 3; break;
      case RDKit::Bond::AROMATIC: out.order = 4; break;
      default: out.order = 1; break;
    }
    result.bonds.push_back(out);
  }

  // Centre on the origin and measure the bounding radius for framing.
  float cx = 0.0f, cy = 0.0f, cz = 0.0f;
  for (const Atom3D& a : result.atoms) {
    cx += a.x;
    cy += a.y;
    cz += a.z;
  }
  const float n = static_cast<float>(std::max<size_t>(result.atoms.size(), 1));
  cx /= n;
  cy /= n;
  cz /= n;
  float radius = 0.0f;
  for (Atom3D& a : result.atoms) {
    a.x -= cx;
    a.y -= cy;
    a.z -= cz;
    radius = std::max(radius, std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z));
  }
  result.radius = std::max(radius, 0.5f);
  return result;
}

}  // namespace chemcad::chem

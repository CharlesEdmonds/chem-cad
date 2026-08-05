#pragma once
// The single boundary between our editable graph and RDKit. RDKit headers are
// included ONLY in bridge.cpp so that the UI, tests and tooling stay cheap to
// compile and free of RDKit types.

#include <stdexcept>
#include <string>

#include "core/model.hpp"

namespace chemcad::chem {

struct ChemError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

std::string toSmiles(const core::Molecule&);      // canonical SMILES
core::Molecule fromSmiles(const std::string&);    // throws ChemError; lays out coords
std::string toMolBlock(const core::Molecule&);
core::Molecule fromMolBlock(const std::string&);  // throws ChemError

// Full CoordGen relayout -- the "Clean Up Structure" command.
void layout(core::Molecule&);

struct Properties {
  std::string formula;
  double mw = 0;
  double logP = 0;
  int rings = 0;
};

Properties computeProperties(const core::Molecule&);  // throws ChemError

// Canonical form used as the identity key for route matching and name caches.
std::string canonicalize(const std::string& smiles);  // throws ChemError

std::string toSvg(const core::Molecule&, int w, int h);

// Implicit hydrogen count for label rendering (CH3, NH2, OH ...).
// Returns 0 rather than throwing when the atom is invalid or the valence model
// cannot be applied -- rendering must never fail on a half-drawn structure.
int implicitHCount(const core::Molecule&, core::AtomId);

// Element helpers backed by RDKit's periodic table.
const char* symbolFor(uint8_t atomicNumber);
uint8_t atomicNumberFor(const std::string& symbol);  // 0 when unknown

}  // namespace chemcad::chem

#pragma once
// 3D coordinate generation for the molecule viewer: RDKit ETKDG conformer
// embedding + MMFF94 relaxation behind the same RDKit-free boundary as
// bridge.hpp. The result is self-contained (element + bond order +
// positions), so the UI never needs RDKit types or atom-id mapping.

#include <cstdint>
#include <vector>

#include "chem/bridge.hpp"
#include "core/model.hpp"

namespace chemcad::chem {

struct Atom3D {
  uint8_t atomicNumber = 0;
  float x = 0.0f;
  float y = 0.0f;
  float z = 0.0f;
};

struct Bond3D {
  int a = 0;
  int b = 0;
  int order = 1;  // 1, 2, 3; 4 = aromatic (drawn single + offset inner line)
};

struct Embedded3D {
  std::vector<Atom3D> atoms;
  std::vector<Bond3D> bonds;
  float radius = 1.0f;  // bounding radius of the centred model, in Angstrom
};

// Deterministic (fixed-seed ETKDG) embedding, so re-opening the view shows
// the same conformer. Implicit/explicit hydrogens are used during embedding
// (they steer the geometry) but are NOT returned -- the viewer follows the
// skeletal convention of hiding them.
// Throws ChemError for an empty molecule or a failed embedding.
Embedded3D embed3D(const core::Molecule&);

}  // namespace chemcad::chem

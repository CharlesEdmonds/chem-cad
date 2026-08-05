// STUB -- implemented by the rdkit-bridge workstream.
#include "chem/bridge.hpp"

namespace chemcad::chem {

std::string toSmiles(const core::Molecule&) { throw ChemError("unimplemented"); }
core::Molecule fromSmiles(const std::string&) { throw ChemError("unimplemented"); }
std::string toMolBlock(const core::Molecule&) { throw ChemError("unimplemented"); }
core::Molecule fromMolBlock(const std::string&) { throw ChemError("unimplemented"); }
void layout(core::Molecule&) { throw ChemError("unimplemented"); }
Properties computeProperties(const core::Molecule&) { throw ChemError("unimplemented"); }
std::string canonicalize(const std::string&) { throw ChemError("unimplemented"); }
std::string toSvg(const core::Molecule&, int, int) { throw ChemError("unimplemented"); }
int implicitHCount(const core::Molecule&, core::AtomId) { return 0; }
const char* symbolFor(uint8_t) { return "C"; }
uint8_t atomicNumberFor(const std::string&) { return 0; }

}  // namespace chemcad::chem

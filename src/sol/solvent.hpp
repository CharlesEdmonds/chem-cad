#pragma once
// Solvent property database. Pure data + lookup; no chemistry model here.

#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace chemcad::sol {

struct SolError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// Hansen solubility parameters, MPa^0.5.
struct Hansen {
  double dispersion = 0.0;    // deltaD
  double polar = 0.0;         // deltaP
  double hydrogenBond = 0.0;  // deltaH
};

struct Solvent {
  std::string id;             // stable key, e.g. "ethyl_acetate"
  std::string name;           // display name, e.g. "Ethyl acetate"
  std::string smiles;
  std::string family;         // "alkane" | "aromatic" | "alcohol" | "ether" | "ester" | "ketone" | "halogenated" | "amide" | "nitrile" | "sulfoxide" | "acid" | "water" | "amine"
  double molarMass = 0.0;     // g/mol
  double density = 0.0;       // g/mL at 25 C
  double molarVolume = 0.0;   // cm3/mol
  Hansen hansen;
  double dielectric = 0.0;
  double boilingPoint = 0.0;  // deg C
  double refractiveIndex = 0.0;
  bool waterMiscible = false;
};

// Loaded once from <data dir>/solvents.json (core::dataDir()). Throws SolError
// when the file is missing or malformed; returns the same vector afterwards.
const std::vector<Solvent>& solvents();

// nullptr when no solvent has that id.
const Solvent* findSolvent(std::string_view id);

}  // namespace chemcad::sol

#pragma once
// Solvent property database. Pure data + lookup; no chemistry model here.
//
// JSON schema and provenance:
// - id, name, SMILES and family are unitless identity/classification strings.
// - molar mass [g/mol], density [g/mL at 25 C], molar volume [cm3/mol],
//   dielectric constant [unitless], boiling point [deg C], refractive index
//   [unitless], Hansen parameters [MPa^0.5], and water miscibility [bool] are
//   required legacy fields compiled from standard solvent-property literature.
// - kappa_t [GPa^-1 at 25 C] and kappa_t_source are optional literature data.
// - melting_point and flash_point are optional [deg C] physical constants from
//   NIST or the standard handbook named by property_source. A flash point of
//   zero means not applicable for a non-flammable solvent.
// - chem21_safety, chem21_health and chem21_environment are optional unitless
//   1..10 CHEM21 scores (higher is worse); chem21_class is the final CHEM21
//   ranking. cost_tier is an internal unitless 1..4 bulk-cost band.
// - peroxide_former [bool] and hazard_note [text] are optional handling aids;
//   property_source [text] records provenance or explicitly says "unrated".
// All selector fields are optional on input and retain the defaults below.

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
  double kappaT = 0.0;        // isothermal compressibility, GPa^-1 at 25 C
  std::string kappaTSource;   // "NIST WTT" | "literature" | "" when unknown
  bool waterMiscible = false;
  double meltingPointC = 0.0;         // freezing point, deg C
  double flashPointC = 0.0;           // deg C, 0 when not applicable (e.g. water)
  int chem21Safety = 0;               // CHEM21 S score, 1..10, higher = worse
  int chem21Health = 0;               // CHEM21 H score, 1..10, higher = worse
  int chem21Environment = 0;          // CHEM21 E score, 1..10, higher = worse
  std::string chem21Class;            // "recommended" | "problematic" | "hazardous" | "highly hazardous"
  int costTier = 2;                   // 1 cheap bulk .. 4 expensive
  bool peroxideFormer = false;        // ethers and similar: peroxide accumulation risk
  std::string hazardNote;             // one short decision-relevant clause
  std::string propertySource;         // provenance, e.g. "CHEM21 2016; NIST"
};

// Loaded once from <data dir>/solvents.json (core::dataDir()). Throws SolError
// when the file is missing or malformed; returns the same vector afterwards.
const std::vector<Solvent>& solvents();

// nullptr when no solvent has that id.
const Solvent* findSolvent(std::string_view id);

// True when the pair forms a single homogeneous phase at room temperature.
// The database only records water miscibility, so the rule is deliberately
// simple: water against a water-immiscible partner separates into two
// phases; every other pairing is treated as mixable. Used to warn when a
// blend prediction assumes a homogeneous phase that does not exist.
bool miscibleWith(const Solvent& a, const Solvent& b);

}  // namespace chemcad::sol

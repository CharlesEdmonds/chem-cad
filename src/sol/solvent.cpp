#include "sol/solvent.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <mutex>

#include <nlohmann/json.hpp>

#include "core/paths.hpp"

namespace chemcad::sol {
namespace {

double requireNumber(const nlohmann::json& entry, const char* key, const std::string& id) {
  if (!entry.contains(key) || !entry.at(key).is_number()) {
    throw SolError("Solvent '" + id + "' is missing required numeric field '" + key + "'");
  }
  return entry.at(key).get<double>();
}

std::string requireString(const nlohmann::json& entry, const char* key, const std::string& id) {
  if (!entry.contains(key) || !entry.at(key).is_string()) {
    throw SolError("Solvent '" + id + "' is missing required string field '" + key + "'");
  }
  return entry.at(key).get<std::string>();
}

bool optionalNumber(const nlohmann::json& entry, const char* key, const std::string& id,
                    double& destination) {
  if (!entry.contains(key)) return false;
  if (!entry.at(key).is_number()) {
    throw SolError("Solvent '" + id + "' has a non-numeric optional field '" + key + "'");
  }
  destination = entry.at(key).get<double>();
  return true;
}

bool optionalInteger(const nlohmann::json& entry, const char* key, const std::string& id,
                     int& destination) {
  if (!entry.contains(key)) return false;
  if (!entry.at(key).is_number_integer()) {
    throw SolError("Solvent '" + id + "' has a non-integer optional field '" + key + "'");
  }
  destination = entry.at(key).get<int>();
  return true;
}

bool optionalBoolean(const nlohmann::json& entry, const char* key, const std::string& id,
                     bool& destination) {
  if (!entry.contains(key)) return false;
  if (!entry.at(key).is_boolean()) {
    throw SolError("Solvent '" + id + "' has a non-boolean optional field '" + key + "'");
  }
  destination = entry.at(key).get<bool>();
  return true;
}

bool optionalString(const nlohmann::json& entry, const char* key, const std::string& id,
                    std::string& destination) {
  if (!entry.contains(key)) return false;
  if (!entry.at(key).is_string()) {
    throw SolError("Solvent '" + id + "' has a non-string optional field '" + key + "'");
  }
  destination = entry.at(key).get<std::string>();
  return true;
}

Solvent parseSolvent(const nlohmann::json& entry) {
  if (!entry.is_object()) {
    throw SolError("solvents.json: each solvent entry must be a JSON object");
  }
  if (!entry.contains("id") || !entry.at("id").is_string()) {
    throw SolError("solvents.json: a solvent entry is missing its 'id' string field");
  }
  const std::string id = entry.at("id").get<std::string>();

  Solvent s;
  s.id = id;
  s.name = requireString(entry, "name", id);
  s.smiles = requireString(entry, "smiles", id);
  s.family = requireString(entry, "family", id);
  s.molarMass = requireNumber(entry, "molar_mass", id);
  s.density = requireNumber(entry, "density", id);
  s.molarVolume = requireNumber(entry, "molar_volume", id);
  s.hansen.dispersion = requireNumber(entry, "delta_d", id);
  s.hansen.polar = requireNumber(entry, "delta_p", id);
  s.hansen.hydrogenBond = requireNumber(entry, "delta_h", id);
  s.dielectric = requireNumber(entry, "dielectric", id);
  s.boilingPoint = requireNumber(entry, "boiling_point", id);
  s.refractiveIndex = requireNumber(entry, "refractive_index", id);
  if (!entry.contains("water_miscible") || !entry.at("water_miscible").is_boolean()) {
    throw SolError("Solvent '" + id + "' is missing required boolean field 'water_miscible'");
  }
  s.waterMiscible = entry.at("water_miscible").get<bool>();

  // Optional Kirkwood-Buff inputs: not every solvent has measured
  // compressibility; absence just disables the KB readout for that solvent.
  if (entry.contains("kappa_t")) {
    if (!entry.at("kappa_t").is_number()) {
      throw SolError("Solvent '" + id + "' has a non-numeric 'kappa_t'");
    }
    s.kappaT = entry.at("kappa_t").get<double>();
    s.kappaTSource = entry.value("kappa_t_source", std::string("literature"));
  }

  // Practical-selection data was added after the original database schema.
  // Keep every key optional so an older user-provided database remains valid.
  optionalNumber(entry, "melting_point", id, s.meltingPointC);
  optionalNumber(entry, "flash_point", id, s.flashPointC);
  bool hasChem21Rating = optionalInteger(entry, "chem21_safety", id, s.chem21Safety);
  hasChem21Rating |= optionalInteger(entry, "chem21_health", id, s.chem21Health);
  hasChem21Rating |=
      optionalInteger(entry, "chem21_environment", id, s.chem21Environment);
  hasChem21Rating |= optionalString(entry, "chem21_class", id, s.chem21Class);
  optionalInteger(entry, "cost_tier", id, s.costTier);
  optionalBoolean(entry, "peroxide_former", id, s.peroxideFormer);
  optionalString(entry, "hazard_note", id, s.hazardNote);
  if (!optionalString(entry, "property_source", id, s.propertySource)) {
    s.propertySource = hasChem21Rating
                           ? "CHEM21 rating present; physical-property source not recorded"
                           : "unrated (no CHEM21 entry)";
  }

  if ((s.chem21Safety != 0 && (s.chem21Safety < 1 || s.chem21Safety > 10)) ||
      (s.chem21Health != 0 && (s.chem21Health < 1 || s.chem21Health > 10)) ||
      (s.chem21Environment != 0 &&
       (s.chem21Environment < 1 || s.chem21Environment > 10))) {
    throw SolError("Solvent '" + id + "' has a CHEM21 score outside 1..10");
  }
  if (s.costTier < 1 || s.costTier > 4) {
    throw SolError("Solvent '" + id + "' has a cost_tier outside 1..4");
  }
  if (!s.chem21Class.empty() && s.chem21Class != "recommended" &&
      s.chem21Class != "problematic" && s.chem21Class != "hazardous" &&
      s.chem21Class != "highly hazardous") {
    throw SolError("Solvent '" + id + "' has an invalid chem21_class");
  }

  // Sanity-check the numbers that the solubility model divides by / relies on
  // being physically meaningful; a bad data entry should fail loudly here
  // rather than silently corrupting every prediction that touches it.
  if (!std::isfinite(s.molarMass) || s.molarMass <= 0.0) {
    throw SolError("Solvent '" + id + "' has a non-finite or non-positive molar_mass");
  }
  if (!std::isfinite(s.density) || s.density <= 0.0) {
    throw SolError("Solvent '" + id + "' has a non-finite or non-positive density");
  }
  if (!std::isfinite(s.molarVolume) || s.molarVolume <= 0.0) {
    throw SolError("Solvent '" + id + "' has a non-finite or non-positive molar_volume");
  }
  return s;
}

std::vector<Solvent> loadSolvents() {
  const auto path = core::dataDir() / "solvents.json";
  std::ifstream input(path);
  if (!input) {
    throw SolError("Could not open solvent database " + path.string());
  }

  nlohmann::json root;
  try {
    input >> root;
  } catch (const std::exception& e) {
    throw SolError("Could not parse solvent database " + path.string() + ": " + e.what());
  }

  if (!root.is_object() || !root.contains("solvents") || !root.at("solvents").is_array()) {
    throw SolError("solvent database " + path.string() + " must contain a 'solvents' array");
  }

  std::vector<Solvent> loaded;
  loaded.reserve(root.at("solvents").size());
  for (const nlohmann::json& entry : root.at("solvents")) {
    Solvent parsed = parseSolvent(entry);
    for (const Solvent& existing : loaded) {
      if (existing.id == parsed.id) {
        throw SolError("Duplicate solvent id '" + parsed.id + "' in " + path.string());
      }
    }
    loaded.push_back(parsed);
  }

  std::sort(loaded.begin(), loaded.end(),
            [](const Solvent& a, const Solvent& b) { return a.name < b.name; });
  return loaded;
}

}  // namespace

const std::vector<Solvent>& solvents() {
  // Loaded exactly once behind call_once; a parse failure is captured and
  // rethrown on every later call instead of retrying a broken file each
  // frame. After the first call this is a lock-free flag check, safe to call
  // from the render loop.
  static std::vector<Solvent> table;
  static std::once_flag once;
  static std::exception_ptr loadError;

  std::call_once(once, [] {
    try {
      table = loadSolvents();
    } catch (...) {
      loadError = std::current_exception();
    }
  });

  if (loadError) std::rethrow_exception(loadError);
  return table;
}

const Solvent* findSolvent(std::string_view id) {
  const std::vector<Solvent>& table = solvents();  // propagates a load failure
  for (const Solvent& s : table) {
    if (s.id == id) return &s;
  }
  return nullptr;
}

bool miscibleWith(const Solvent& a, const Solvent& b) {
  // The database records water miscibility only, so the rule is deliberately
  // simple: water against a water-immiscible partner separates; any other
  // pairing is treated as mixable. Water+water takes the first branch and is
  // true because water is waterMiscible.
  if (a.family == "water") return b.waterMiscible;
  if (b.family == "water") return a.waterMiscible;
  return true;
}

}  // namespace chemcad::sol

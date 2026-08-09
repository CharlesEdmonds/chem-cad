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

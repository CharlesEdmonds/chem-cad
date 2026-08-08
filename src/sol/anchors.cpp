// Literature anchors and aqueous ionic chemistry. See anchors.hpp for the
// model contract. Loading mirrors solvent.cpp: sticky call_once table with a
// persistent exception on load failure.

#include "sol/anchors.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <mutex>
#include <utility>

#include "chem/bridge.hpp"
#include "core/paths.hpp"

namespace chemcad::sol {
namespace {

struct AnchorTables {
  std::vector<SolubilityAnchor> anchors;
  std::vector<Salt> salts;
  std::vector<Electrolyte> electrolytes;
};

// Keys in the JSON are author-written SMILES; canonicalising both sides at
// load means a salt or solute matches no matter which valid SMILES the user
// drew. A key RDKit cannot read is a data bug and fails the load.
std::string canonicalKey(const nlohmann::json& entry, const char* file) {
  if (!entry.contains("solute") || !entry.at("solute").is_string()) {
    throw SolError(std::string(file) + ": an entry is missing its 'solute' SMILES");
  }
  const std::string smiles = entry.at("solute").get<std::string>();
  try {
    return chem::canonicalize(smiles);
  } catch (const std::exception& e) {
    throw SolError(std::string(file) + ": solute SMILES '" + smiles + "' does not parse: " +
                   e.what());
  }
}

AnchorTables loadTables() {
  const auto path = core::dataDir() / "solubility_anchors.json";
  std::ifstream input(path);
  if (!input) {
    throw SolError("Could not open solubility database " + path.string());
  }
  nlohmann::json root;
  try {
    input >> root;
  } catch (const std::exception& e) {
    throw SolError("Could not parse solubility database " + path.string() + ": " + e.what());
  }

  AnchorTables tables;
  for (const nlohmann::json& entry : root.value("anchors", nlohmann::json::array())) {
    SolubilityAnchor anchor;
    anchor.soluteSmiles = canonicalKey(entry, "solubility_anchors.json");
    anchor.solventId = entry.value("solvent", "");
    anchor.temperatureC = entry.value("temperature_c", 25.0);
    anchor.gramsPerMillilitre = entry.value("g_per_ml", 0.0);
    anchor.soluteMeltingPointC = entry.value("solute_tm_c", 25.0);
    anchor.note = entry.value("note", "");
    if (anchor.solventId.empty() || anchor.gramsPerMillilitre <= 0.0) {
      throw SolError("solubility_anchors.json: an anchor is missing solvent or g_per_ml");
    }
    if (!findSolvent(anchor.solventId)) {
      throw SolError("solubility_anchors.json: unknown solvent id '" + anchor.solventId + "'");
    }
    tables.anchors.push_back(std::move(anchor));
  }
  for (const nlohmann::json& entry : root.value("salts", nlohmann::json::array())) {
    Salt salt;
    salt.soluteSmiles = canonicalKey(entry, "solubility_anchors.json");
    salt.name = entry.value("name", "");
    salt.cation = entry.value("cation", "");
    salt.anion = entry.value("anion", "");
    salt.molarMass = entry.value("molar_mass", 0.0);
    salt.solubilityGPerMl25 = entry.value("g_per_ml_water25", 0.0);
    salt.ksp25 = entry.value("ksp25", 0.0);
    salt.dissolutionEnthalpyKj = entry.value("dh_solution_kj", 0.0);
    salt.note = entry.value("note", "");
    if (salt.name.empty() || salt.molarMass <= 0.0 || salt.solubilityGPerMl25 <= 0.0 ||
        salt.ksp25 <= 0.0) {
      throw SolError("solubility_anchors.json: salt '" + salt.name + "' has incomplete fields");
    }
    tables.salts.push_back(std::move(salt));
  }
  for (const nlohmann::json& entry : root.value("electrolytes", nlohmann::json::array())) {
    Electrolyte electrolyte;
    electrolyte.id = entry.value("id", "");
    electrolyte.name = entry.value("name", "");
    electrolyte.cation = entry.value("cation", "");
    electrolyte.anion = entry.value("anion", "");
    electrolyte.molarMass = entry.value("molar_mass", 0.0);
    if (electrolyte.id.empty() || electrolyte.cation.empty() || electrolyte.anion.empty()) {
      throw SolError("solubility_anchors.json: an electrolyte has incomplete fields");
    }
    tables.electrolytes.push_back(std::move(electrolyte));
  }
  std::sort(tables.electrolytes.begin(), tables.electrolytes.end(),
            [](const Electrolyte& a, const Electrolyte& b) { return a.name < b.name; });
  return tables;
}

const AnchorTables& tables() {
  static AnchorTables loaded;
  static std::exception_ptr loadError;
  static std::once_flag once;
  std::call_once(once, [] {
    try {
      loaded = loadTables();
    } catch (...) {
      loadError = std::current_exception();
    }
  });
  if (loadError) std::rethrow_exception(loadError);
  return loaded;
}

constexpr double kGasConstant = 8.314462618;  // J/(mol*K)

}  // namespace

// Davies equation: mean ionic activity of a 1:1 electrolyte. Valid to about
// I = 0.5; the term is clamped there so concentrated brines get the 0.5 M
// correction instead of an unphysical one.
double daviesGamma(double ionicStrength) {
  const double i = std::min(std::max(ionicStrength, 0.0), 0.5);
  const double root = std::sqrt(i);
  const double logGamma = -0.509 * (root / (1.0 + root) - 0.3 * i);
  return std::pow(10.0, logGamma);
}

const std::vector<SolubilityAnchor>& anchors() { return tables().anchors; }
const std::vector<Salt>& salts() { return tables().salts; }
const std::vector<Electrolyte>& electrolytes() { return tables().electrolytes; }

const SolubilityAnchor* findAnchor(std::string_view canonicalSmiles, std::string_view solventId,
                                   double temperatureC, double toleranceC) {
  for (const SolubilityAnchor& anchor : tables().anchors) {
    if (anchor.soluteSmiles == canonicalSmiles && anchor.solventId == solventId &&
        std::fabs(anchor.temperatureC - temperatureC) <= toleranceC) {
      return &anchor;
    }
  }
  return nullptr;
}

const Salt* findSalt(std::string_view canonicalSmiles) {
  for (const Salt& salt : tables().salts) {
    if (salt.soluteSmiles == canonicalSmiles) return &salt;
  }
  return nullptr;
}

const Electrolyte* findElectrolyte(std::string_view id) {
  for (const Electrolyte& electrolyte : tables().electrolytes) {
    if (electrolyte.id == id) return &electrolyte;
  }
  return nullptr;
}

double saltSolubilityMolar(const Salt& salt, const Electrolyte* background, double backgroundM,
                           double temperatureC) {
  // van't Hoff: endothermic salts dissolve better hot.
  const double t = temperatureC + 273.15;
  const double ksp = salt.ksp25 * std::exp(-(salt.dissolutionEnthalpyKj * 1000.0 / kGasConstant) *
                                           (1.0 / t - 1.0 / 298.15));

  const double commonM = (background && backgroundM > 0.0 &&
                          (background->cation == salt.cation || background->anion == salt.anion))
                             ? backgroundM
                             : 0.0;
  const double ionicBackground = (background && backgroundM > 0.0) ? backgroundM : 0.0;

  // gamma(I)^2 * (s + common) * s = Ksp, I = background + s. The left side is
  // monotone increasing in s, so bisection brackets it robustly -- but the
  // activity correction (gamma < 1) pushes the root ABOVE sqrt(Ksp), so the
  // bracket must be expanded until it actually contains the sign change.
  auto residual = [&](double s) {
    const double gamma = daviesGamma(ionicBackground + s);
    return gamma * gamma * (s + commonM) * s - ksp;
  };
  double lo = 0.0;
  double hi = std::sqrt(ksp) + commonM + 1.0;
  for (int expand = 0; expand < 24 && residual(hi) < 0.0; ++expand) {
    hi = hi * 2.0 + 1.0;
  }
  for (int i = 0; i < 80; ++i) {
    const double mid = 0.5 * (lo + hi);
    if (residual(mid) > 0.0) {
      hi = mid;
    } else {
      lo = mid;
    }
  }
  return 0.5 * (lo + hi);
}

}  // namespace chemcad::sol

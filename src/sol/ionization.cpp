#include "sol/ionization.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace chemcad::sol {

namespace {

constexpr double kGasConstant = 8.314;  // J/(mol K)
constexpr double kLn10 = 2.302585092994046;
// N_A * e^2 / (8 pi eps_0), J.m/mol: the Born prefactor, so the self-energy
// of a mole of unit charges of radius r is kBornPrefactor / r.
constexpr double kBornPrefactor = 6.947e-5;
// Latimer-style radius correction: bare crystallographic radii make the Born
// equation overestimate hydration/transfer free energies; adding ~0.08 nm to
// the ionic radius is the standard empirical repair (Latimer, Pitzer &
// Slansky, J. Chem. Phys. 1939, 7, 108).
constexpr double kLatimerShiftNm = 0.08;

bool isAromaticAtom(const core::Molecule& molecule, core::AtomId id) {
  for (core::BondId bondId : molecule.incidentBonds(id)) {
    const core::Bond* bond = molecule.bond(bondId);
    if (bond && bond->order == core::BondOrder::Aromatic) return true;
  }
  return false;
}

// True when `carbon` is a carbonyl (or thiocarbonyl) carbon: the electron
// sink that turns a neighbouring nitrogen into an amide instead of an amine.
bool isCarbonylCarbon(const core::Molecule& molecule, core::AtomId carbon) {
  const core::Atom* atom = molecule.atom(carbon);
  if (!atom || atom->atomicNumber != 6) return false;
  for (core::BondId bondId : molecule.incidentBonds(carbon)) {
    const core::Bond* bond = molecule.bond(bondId);
    if (!bond || bond->order != core::BondOrder::Double) continue;
    core::AtomId other = bond->a == carbon ? bond->b : bond->a;
    const core::Atom* neighbour = molecule.atom(other);
    if (neighbour && (neighbour->atomicNumber == 8 || neighbour->atomicNumber == 16)) return true;
  }
  return false;
}

bool hasDoubleBondedOxygen(const core::Molecule& molecule, core::AtomId id) {
  for (core::BondId bondId : molecule.incidentBonds(id)) {
    const core::Bond* bond = molecule.bond(bondId);
    if (!bond || bond->order != core::BondOrder::Double) continue;
    core::AtomId other = bond->a == id ? bond->b : bond->a;
    const core::Atom* neighbour = molecule.atom(other);
    if (neighbour && neighbour->atomicNumber == 8) return true;
  }
  return false;
}

int heavyDegree(const core::Molecule& molecule, core::AtomId id) {
  int degree = 0;
  for (core::AtomId nb : molecule.neighbors(id)) {
    const core::Atom* atom = molecule.atom(nb);
    if (atom && atom->atomicNumber != 1) ++degree;
  }
  return degree;
}

// Size of the smallest ring through `id`, or 0 when the atom is acyclic.
// Needed because a five-membered aromatic nitrogen (imidazole, pKa 7) and a
// six-membered one (pyridine, pKa 5.2) are different bases, and a bridge test
// alone cannot tell them apart. Bounded breadth-first search between each
// pair of neighbours, ignoring the centre atom.
int smallestRingSize(const core::Molecule& molecule, core::AtomId id) {
  const std::vector<core::AtomId> neighbours = molecule.neighbors(id);
  int best = 0;
  for (size_t i = 0; i < neighbours.size(); ++i) {
    for (size_t j = i + 1; j < neighbours.size(); ++j) {
      // Shortest path neighbours[i] -> neighbours[j] that avoids `id`.
      std::unordered_map<core::AtomId, int> depth;
      std::deque<core::AtomId> queue;
      depth[neighbours[i]] = 0;
      queue.push_back(neighbours[i]);
      int found = -1;
      while (!queue.empty() && found < 0) {
        const core::AtomId current = queue.front();
        queue.pop_front();
        const int d = depth[current];
        if (d > 6) continue;  // rings larger than 8 atoms are irrelevant here
        for (core::AtomId nb : molecule.neighbors(current)) {
          if (nb == id) continue;
          if (nb == neighbours[j]) {
            found = d + 1;
            break;
          }
          if (depth.find(nb) != depth.end()) continue;
          depth[nb] = d + 1;
          queue.push_back(nb);
        }
      }
      if (found > 0) {
        const int cycle = found + 2;  // + the two bonds through `id`
        if (best == 0 || cycle < best) best = cycle;
      }
    }
  }
  return best;
}

// Electron-withdrawing carbonyls attached to the ring system a nitrogen sits
// in. Each one drains basicity: caffeine's imidazole nitrogen would look like
// a pKa-7 base on ring type alone, but the two ring carbonyls of the xanthine
// core pull its real pKa down to 0.6. Two decades per carbonyl is the
// Hammett-scale order of magnitude and reproduces that ordering.
int ringCarbonylCount(const core::Molecule& molecule, core::AtomId nitrogen) {
  // Walk the aromatic/ring neighbourhood two bonds out and count carbonyls.
  std::unordered_set<core::AtomId> seen{nitrogen};
  std::deque<std::pair<core::AtomId, int>> queue{{nitrogen, 0}};
  int carbonyls = 0;
  while (!queue.empty()) {
    const auto [current, depth] = queue.front();
    queue.pop_front();
    if (depth >= 3) continue;
    for (core::AtomId nb : molecule.neighbors(current)) {
      if (!seen.insert(nb).second) continue;
      const core::Atom* atom = molecule.atom(nb);
      if (!atom) continue;
      if (atom->atomicNumber == 6 && hasDoubleBondedOxygen(molecule, nb)) ++carbonyls;
      queue.emplace_back(nb, depth + 1);
    }
  }
  return carbonyls;
}

// ---- group pKa table --------------------------------------------------
// Textbook aqueous values (Perrin's tables / any physical-organic text), not
// a fitted regression. Bases carry the pKa of their conjugate acid BH+.
namespace pKaTable {
constexpr double kGuanidine = 13.6;
constexpr double kPrimaryAmine = 10.6;
constexpr double kSecondaryAmine = 11.1;
constexpr double kTertiaryAmine = 9.8;
constexpr double kAniline = 4.6;
constexpr double kPyridine = 5.2;   // six-membered aromatic N
constexpr double kAzole = 7.0;      // five-membered aromatic N (imidazole)
constexpr double kRingCarbonylPenalty = 2.0;

constexpr double kAromaticAcid = 4.2;   // benzoic acid
constexpr double kAliphaticAcid = 4.8;  // acetic acid
constexpr double kSulfonicAcid = 1.0;
constexpr double kPhenol = 10.0;
constexpr double kThiol = 10.5;
}  // namespace pKaTable

struct Site {
  double pKa = 0.0;
  const char* label = "";
};

// Strongest base site (highest pKa) on the skeleton, if any.
bool strongestBase(const core::Molecule& molecule, Site& out) {
  bool found = false;
  for (const core::Atom& atom : molecule.atoms()) {
    if (atom.atomicNumber != 7) continue;

    // Amide, sulfonamide, nitrile and nitro nitrogens are not bases.
    bool amideLike = false;
    bool nitrileOrNitro = false;
    for (core::BondId bondId : molecule.incidentBonds(atom.id)) {
      const core::Bond* bond = molecule.bond(bondId);
      if (!bond) continue;
      const core::AtomId other = bond->a == atom.id ? bond->b : bond->a;
      const core::Atom* neighbour = molecule.atom(other);
      if (!neighbour) continue;
      if (bond->order == core::BondOrder::Triple) nitrileOrNitro = true;
      if (neighbour->atomicNumber == 8 && bond->order == core::BondOrder::Double) {
        nitrileOrNitro = true;  // nitro
      }
      if (neighbour->atomicNumber == 16 && hasDoubleBondedOxygen(molecule, other)) {
        amideLike = true;  // sulfonamide
      }
      if (isCarbonylCarbon(molecule, other)) amideLike = true;
    }
    if (nitrileOrNitro || amideLike) continue;

    Site site;
    if (isAromaticAtom(molecule, atom.id)) {
      // A substituted aromatic nitrogen (pyrrole-type, three heavy
      // neighbours) has its lone pair in the ring and is not basic.
      if (heavyDegree(molecule, atom.id) >= 3) continue;
      const int ringSize = smallestRingSize(molecule, atom.id);
      site.pKa = ringSize == 5 ? pKaTable::kAzole : pKaTable::kPyridine;
      site.label = ringSize == 5 ? "azole ring N" : "pyridine ring N";
      site.pKa -= pKaTable::kRingCarbonylPenalty *
                  static_cast<double>(ringCarbonylCount(molecule, atom.id));
    } else {
      // Guanidine: nitrogen on a carbon bearing two or more further
      // nitrogens.
      bool guanidine = false;
      bool onAromatic = false;
      for (core::AtomId nb : molecule.neighbors(atom.id)) {
        const core::Atom* neighbour = molecule.atom(nb);
        if (!neighbour || neighbour->atomicNumber != 6) continue;
        if (isAromaticAtom(molecule, nb)) onAromatic = true;
        int nitrogens = 0;
        for (core::AtomId sub : molecule.neighbors(nb)) {
          const core::Atom* subAtom = molecule.atom(sub);
          if (subAtom && subAtom->atomicNumber == 7) ++nitrogens;
        }
        if (nitrogens >= 3) guanidine = true;
      }
      const int degree = heavyDegree(molecule, atom.id);
      if (guanidine) {
        site.pKa = pKaTable::kGuanidine;
        site.label = "guanidine";
      } else if (onAromatic) {
        site.pKa = pKaTable::kAniline;
        site.label = "aromatic amine";
      } else if (degree <= 1) {
        site.pKa = pKaTable::kPrimaryAmine;
        site.label = "primary amine";
      } else if (degree == 2) {
        site.pKa = pKaTable::kSecondaryAmine;
        site.label = "secondary amine";
      } else {
        site.pKa = pKaTable::kTertiaryAmine;
        site.label = "tertiary amine";
      }
    }
    if (!found || site.pKa > out.pKa) {
      out = site;
      found = true;
    }
  }
  return found;
}

// Strongest acid site (lowest pKa) on the skeleton, if any.
bool strongestAcid(const core::Molecule& molecule, Site& out) {
  bool found = false;
  const auto consider = [&](const Site& site) {
    if (!found || site.pKa < out.pKa) {
      out = site;
      found = true;
    }
  };

  for (const core::Atom& atom : molecule.atoms()) {
    if (atom.atomicNumber == 16) {
      // Sulfonic acid: S with two double-bonded oxygens and a hydroxyl.
      int doubleO = 0;
      bool hydroxyl = false;
      for (core::BondId bondId : molecule.incidentBonds(atom.id)) {
        const core::Bond* bond = molecule.bond(bondId);
        if (!bond) continue;
        const core::AtomId other = bond->a == atom.id ? bond->b : bond->a;
        const core::Atom* neighbour = molecule.atom(other);
        if (!neighbour || neighbour->atomicNumber != 8) continue;
        if (bond->order == core::BondOrder::Double) ++doubleO;
        else if (heavyDegree(molecule, other) <= 1) hydroxyl = true;
      }
      if (doubleO >= 2 && hydroxyl) consider(Site{pKaTable::kSulfonicAcid, "sulfonic acid"});
      // Thiol: sulfur with one SINGLE bond to a heavy atom. A thiocarbonyl
      // (C=S) also has one heavy neighbour and is not acidic.
      bool singleBonded = false;
      for (core::BondId bondId : molecule.incidentBonds(atom.id)) {
        const core::Bond* bond = molecule.bond(bondId);
        if (bond && bond->order == core::BondOrder::Single) singleBonded = true;
      }
      if (heavyDegree(molecule, atom.id) <= 1 && doubleO == 0 && singleBonded) {
        consider(Site{pKaTable::kThiol, "thiol"});
      }
      continue;
    }
    if (atom.atomicNumber != 8) continue;
    if (heavyDegree(molecule, atom.id) != 1) continue;  // ether/ester bridge O

    // The single heavy neighbour decides what this hydroxyl is -- and the bond
    // to it must be SINGLE. A carbonyl oxygen (ester, ketone, amide, acid)
    // also has exactly one heavy neighbour, and reading it as a hydroxyl is
    // what made an ester-rich alkaloid look like a carboxylic acid, hence a
    // zwitterion, hence not a base at all.
    core::AtomId host = core::kInvalidAtom;
    for (core::BondId bondId : molecule.incidentBonds(atom.id)) {
      const core::Bond* bond = molecule.bond(bondId);
      if (!bond || bond->order != core::BondOrder::Single) continue;
      const core::AtomId other = bond->a == atom.id ? bond->b : bond->a;
      const core::Atom* neighbour = molecule.atom(other);
      if (neighbour && neighbour->atomicNumber != 1) host = other;
    }
    if (host == core::kInvalidAtom) continue;
    const core::Atom* hostAtom = molecule.atom(host);
    if (!hostAtom) continue;

    if (hostAtom->atomicNumber == 6 && hasDoubleBondedOxygen(molecule, host)) {
      // Carboxylic acid; aromatic attachment is the stronger acid.
      bool aromaticAttached = false;
      for (core::AtomId nb : molecule.neighbors(host)) {
        if (nb != atom.id && isAromaticAtom(molecule, nb)) aromaticAttached = true;
      }
      consider(Site{aromaticAttached ? pKaTable::kAromaticAcid : pKaTable::kAliphaticAcid,
                    "carboxylic acid"});
    } else if (hostAtom->atomicNumber == 6 && isAromaticAtom(molecule, host)) {
      consider(Site{pKaTable::kPhenol, "phenol"});
    }
    // Aliphatic alcohols (pKa ~16) never ionise in a usable window.
  }
  return found;
}

int heavyAtomCount(const core::Molecule& molecule) {
  int count = 0;
  for (const core::Atom& atom : molecule.atoms()) {
    if (atom.atomicNumber != 1) ++count;
  }
  return count;
}

const char* elementSymbol(uint8_t atomicNumber) {
  switch (atomicNumber) {
    case 1: return "H";
    case 6: return "C";
    case 7: return "N";
    case 8: return "O";
    case 9: return "F";
    case 11: return "Na";
    case 12: return "Mg";
    case 15: return "P";
    case 16: return "S";
    case 17: return "Cl";
    case 19: return "K";
    case 20: return "Ca";
    case 35: return "Br";
    case 53: return "I";
    default: return "X";
  }
}

std::string formulaOf(const core::Molecule& molecule) {
  std::unordered_map<uint8_t, int> counts;
  for (const core::Atom& atom : molecule.atoms()) ++counts[atom.atomicNumber];
  // Hill-ish order: C, H, then everything else by atomic number.
  std::vector<uint8_t> order;
  if (counts.count(6)) order.push_back(6);
  if (counts.count(1)) order.push_back(1);
  std::vector<uint8_t> rest;
  for (const auto& [z, n] : counts) {
    if (z != 1 && z != 6) rest.push_back(z);
  }
  std::sort(rest.begin(), rest.end());
  order.insert(order.end(), rest.begin(), rest.end());

  std::string formula;
  for (uint8_t z : order) {
    formula += elementSymbol(z);
    const int n = counts[z];
    if (n > 1) formula += std::to_string(n);
  }
  return formula;
}

// Crystallographic (Shannon/Pauling) radii of the counter-ions a bench
// chemist actually sees, nm. Anything unrecognised uses 0.20 nm, the middle
// of the range for small organic anions.
double counterIonRadius(const core::Molecule& ion) {
  double smallest = 0.0;
  for (const core::Atom& atom : ion.atoms()) {
    double r = 0.0;
    switch (atom.atomicNumber) {
      case 1: r = 0.10; break;   // proton (hydrated, effective)
      case 9: r = 0.133; break;  // F-
      case 11: r = 0.102; break; // Na+
      case 17: r = 0.181; break; // Cl-
      case 19: r = 0.138; break; // K+
      case 35: r = 0.196; break; // Br-
      case 53: r = 0.220; break; // I-
      case 16: r = 0.230; break; // sulfate / sulfonate
      case 7: r = 0.200; break;  // nitrate
      default: r = 0.0; break;
    }
    if (r > smallest) smallest = r;
  }
  return smallest > 0.0 ? smallest : 0.20;
}

// Does this fragment behave as a counter-ion rather than a second molecule?
bool looksLikeCounterIon(const core::Molecule& fragment) {
  int carbons = 0;
  int netCharge = 0;
  for (const core::Atom& atom : fragment.atoms()) {
    if (atom.atomicNumber == 6) ++carbons;
    netCharge += atom.charge;
  }
  const int heavy = heavyAtomCount(fragment);
  if (heavy == 0) return false;
  if (netCharge != 0) return true;             // explicitly drawn as an ion
  if (carbons == 0 && heavy <= 6) return true;  // HCl, HBr, H2SO4, NaOH, ...
  if (heavy > 12) return false;

  // Small organic acid used as a salt former (acetate, mesylate, tosylate,
  // maleate, tartrate): it must carry an acid group of its own.
  Site acid;
  return strongestAcid(fragment, acid) && acid.pKa <= 5.0;
}

}  // namespace

std::vector<core::Molecule> splitComponents(const core::Molecule& molecule) {
  std::vector<core::Molecule> components;
  if (molecule.empty()) return components;

  std::unordered_set<core::AtomId> visited;
  for (const core::Atom& seed : molecule.atoms()) {
    if (visited.count(seed.id)) continue;

    // Breadth-first sweep over one connected component.
    std::vector<core::AtomId> members;
    std::deque<core::AtomId> queue{seed.id};
    visited.insert(seed.id);
    while (!queue.empty()) {
      const core::AtomId current = queue.front();
      queue.pop_front();
      members.push_back(current);
      for (core::AtomId nb : molecule.neighbors(current)) {
        if (visited.insert(nb).second) queue.push_back(nb);
      }
    }

    core::Molecule component;
    std::unordered_map<core::AtomId, core::AtomId> remap;
    remap.reserve(members.size());
    for (core::AtomId id : members) {
      const core::Atom* atom = molecule.atom(id);
      if (!atom) continue;
      core::Atom copy = *atom;
      remap.emplace(id, component.addAtom(copy));
    }
    for (core::AtomId id : members) {
      for (core::BondId bondId : molecule.incidentBonds(id)) {
        const core::Bond* bond = molecule.bond(bondId);
        if (!bond || bond->a != id) continue;  // add each bond once, from its `a` end
        auto a = remap.find(bond->a);
        auto b = remap.find(bond->b);
        if (a == remap.end() || b == remap.end()) continue;
        const core::BondId added = component.addBond(a->second, b->second, bond->order);
        if (core::Bond* copy = component.bond(added)) copy->stereo = bond->stereo;
      }
    }
    components.push_back(std::move(component));
  }

  std::stable_sort(components.begin(), components.end(),
                   [](const core::Molecule& a, const core::Molecule& b) {
                     return heavyAtomCount(a) > heavyAtomCount(b);
                   });
  return components;
}

Ionization analyseIonization(const core::Molecule& structure) {
  Ionization result;
  const std::vector<core::Molecule> components = splitComponents(structure);
  if (components.empty()) return result;

  const core::Molecule& skeleton = components.front();
  result.fragmentCount = static_cast<int>(components.size());
  for (const core::Atom& atom : skeleton.atoms()) result.netCharge += atom.charge;

  Site base;
  Site acid;
  const bool hasBase = strongestBase(skeleton, base);
  const bool hasAcid = strongestAcid(skeleton, acid);

  if (hasBase && hasAcid && base.pKa >= 8.0 && acid.pKa <= 5.0) {
    // Both sites strong: the species is a zwitterion over a wide pH window,
    // net neutral, and the Henderson-Hasselbalch enhancement below does not
    // describe it (its aqueous solubility is already the zwitterion's).
    result.ionClass = IonClass::Zwitterion;
    result.pKa = base.pKa;
    result.siteLabel = std::string(acid.label) + " + " + base.label;
  } else if (hasBase && (!hasAcid || (base.pKa - 7.0) >= (7.0 - acid.pKa))) {
    result.ionClass = IonClass::Base;
    result.pKa = base.pKa;
    result.siteLabel = base.label;
  } else if (hasAcid) {
    result.ionClass = IonClass::Acid;
    result.pKa = acid.pKa;
    result.siteLabel = acid.label;
  }

  // Counter-ions: every component but the skeleton that reads as an ion.
  for (size_t i = 1; i < components.size(); ++i) {
    if (!looksLikeCounterIon(components[i])) continue;
    const std::string formula = formulaOf(components[i]);
    if (result.counterIon.empty()) {
      result.counterIon = formula;
      result.counterIonRadiusNm = counterIonRadius(components[i]);
    } else if (result.counterIon != formula) {
      result.counterIon += "." + formula;
    }
    result.saltForm = true;
  }
  // A drawn formal charge on the skeleton is a salt even without a separate
  // counter-ion fragment (the user simply did not draw it).
  if (result.netCharge != 0) result.saltForm = true;
  // A counter-ion next to a species with nothing to ionise is a spectator,
  // not a salt former.
  if (result.ionClass == IonClass::Neutral && result.netCharge == 0) result.saltForm = false;

  return result;
}

double ionisedRatio(const Ionization& ionization, double pH) {
  double exponent = 0.0;
  switch (ionization.ionClass) {
    case IonClass::Base: exponent = ionization.pKa - pH; break;
    case IonClass::Acid: exponent = pH - ionization.pKa; break;
    case IonClass::Neutral:
    case IonClass::Zwitterion: return 0.0;
  }
  return std::min(std::pow(10.0, std::clamp(exponent, -12.0, 12.0)), 1e12);
}

double selfBufferedPH(const Ionization& ionization, double molarity) {
  const double c = std::clamp(molarity, 1e-9, 30.0);
  const double logC = std::log10(c);
  double pH = 7.0;
  switch (ionization.ionClass) {
    case IonClass::Base:
      // Salt of a base (BH+ hydrolysis) is acidic; the free base is basic.
      pH = ionization.saltForm ? 0.5 * (ionization.pKa - logC)
                               : 14.0 - 0.5 * ((14.0 - ionization.pKa) - logC);
      break;
    case IonClass::Acid:
      // Salt of an acid (A- hydrolysis) is basic; the free acid is acidic.
      pH = ionization.saltForm ? 7.0 + 0.5 * (ionization.pKa + logC)
                               : 0.5 * (ionization.pKa - logC);
      break;
    case IonClass::Neutral:
    case IonClass::Zwitterion:
      pH = 7.0;
      break;
  }
  return std::clamp(pH, 0.0, 14.0);
}

double bornPenaltyDecades(double dielectric, double soluteRadiusNm, double counterIonRadiusNm,
                          double temperatureK, double waterDielectric) {
  const double epsSolvent = std::max(dielectric, 1.0);
  const double epsWater = std::max(waterDielectric, 1.0);
  const double inverseGap = 1.0 / epsSolvent - 1.0 / epsWater;
  if (inverseGap <= 0.0) return 0.0;  // as polar as water or more: no penalty

  const double rSolute = std::max(soluteRadiusNm, 0.05) + kLatimerShiftNm;
  const double rCounter = std::max(counterIonRadiusNm, 0.05) + kLatimerShiftNm;
  // Self-energy of both ions of a 1:1 salt, J/mol (radii nm -> m).
  const double selfEnergy = kBornPrefactor * (1.0 / (rSolute * 1e-9) + 1.0 / (rCounter * 1e-9));
  const double deltaG = selfEnergy * inverseGap;  // J/mol, positive = penalty
  const double decades = deltaG / (kLn10 * kGasConstant * std::max(temperatureK, 1.0));
  return std::max(decades, 0.0);
}

double bornRadiusFromMolarVolumeNm(double molarVolumeCm3PerMol) {
  const double volumePerMolecule =
      std::max(molarVolumeCm3PerMol, 1.0) * 1e-6 / 6.02214076e23;  // m^3
  const double radiusM = std::cbrt(3.0 * volumePerMolecule / (4.0 * 3.14159265358979323846));
  return radiusM * 1e9;
}

}  // namespace chemcad::sol

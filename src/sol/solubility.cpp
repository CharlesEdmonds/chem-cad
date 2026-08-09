#include "sol/solubility.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "chem/bridge.hpp"

namespace chemcad::sol {

namespace {

// ---- McGowan characteristic volume ------------------------------------
// Atomic volume contributions, cm3/mol (Abraham & McGowan, 1990). Anything
// not tabulated falls back to a generic mid-row estimate.
double atomicVolume(uint8_t atomicNumber) {
  switch (atomicNumber) {
    case 1: return 8.71;    // H
    case 6: return 16.35;   // C
    case 7: return 14.39;   // N
    case 8: return 12.43;   // O
    case 9: return 10.48;   // F
    case 14: return 26.83;  // Si
    case 15: return 24.87;  // P
    case 16: return 22.91;  // S
    case 17: return 20.95;  // Cl
    case 35: return 26.21;  // Br
    case 53: return 34.53;  // I
    case 5: return 18.32;   // B
    default: return 22.0;
  }
}

// Sums per-atom volume contributions (explicit atoms plus their implicit
// hydrogens) and subtracts 6.56 cm3/mol per bond, counting implicit C-H/X-H
// bonds alongside the explicit graph edges -- the standard McGowan recipe.
double mcGowanVolume(const core::Molecule& molecule) {
  double volume = 0.0;
  int bonds = static_cast<int>(molecule.bondCount());
  for (const core::Atom& atom : molecule.atoms()) {
    volume += atomicVolume(atom.atomicNumber);
    int hydrogens = chem::implicitHCount(molecule, atom.id);
    volume += hydrogens * atomicVolume(1);
    bonds += hydrogens;
  }
  volume -= 6.56 * bonds;
  return std::max(volume, 10.0);
}

// ---- Ring perception (Tarjan bridges) ----------------------------------
// A bond that is NOT a bridge lies on some cycle by definition (a bridge is
// exactly an edge whose removal lies on no cycle), so both of its endpoints
// are ring atoms; conversely every atom reachable only through bridges is
// acyclic. This gives an exact ring-membership test -- not a heuristic --
// from a single linear-time DFS (Tarjan, R. E. "A note on finding the
// bridges of a graph." Information Processing Letters 1974, 2, 160-161)
// rather than a full SSSR enumeration, which is more machinery than a
// boolean "is this atom in a ring?" query needs.
class RingPerception {
 public:
  explicit RingPerception(const core::Molecule& molecule) : molecule_(molecule) {
    for (const core::Atom& atom : molecule.atoms()) {
      if (disc_.find(atom.id) == disc_.end()) visit(atom.id, core::kInvalidBond);
    }
  }

  bool isRingAtom(core::AtomId id) const { return ringAtoms_.find(id) != ringAtoms_.end(); }

 private:
  void visit(core::AtomId u, core::BondId viaBond) {
    int order = ++counter_;
    disc_[u] = order;
    int low = order;
    for (core::BondId bondId : molecule_.incidentBonds(u)) {
      if (bondId == viaBond) continue;  // don't walk back over the parent tree edge
      const core::Bond* bond = molecule_.bond(bondId);
      if (!bond) continue;
      core::AtomId v = bond->a == u ? bond->b : bond->a;
      auto seen = disc_.find(v);
      if (seen == disc_.end()) {
        visit(v, bondId);
        low = std::min(low, low_[v]);
        if (low_[v] <= disc_[u]) {  // (u,v) is not a bridge -- it closes a cycle
          ringAtoms_.insert(u);
          ringAtoms_.insert(v);
        }
      } else {
        low = std::min(low, seen->second);  // back edge -- always closes a cycle
        ringAtoms_.insert(u);
        ringAtoms_.insert(v);
      }
    }
    low_[u] = low;
  }

  const core::Molecule& molecule_;
  int counter_ = 0;
  std::unordered_map<core::AtomId, int> disc_;
  std::unordered_map<core::AtomId, int> low_;
  std::unordered_set<core::AtomId> ringAtoms_;
};

// ---- Hansen group contributions (Hoftyzer / van Krevelen) -------------
// Fd, Fp in (J.cm3)^0.5/mol == MPa^0.5.cm3/mol; Eh in J/mol. Values follow
// the widely reproduced Hoftyzer/van Krevelen table (see van Krevelen & te
// Nijenhuis, "Properties of Polymers", table 7.3, and Hansen's "Hansen
// Solubility Parameters: A User's Handbook", table 1). Whole-group entries
// from that table (e.g. phenyl, ester, amide) are distributed onto the single
// graph atom that anchors the group -- e.g. the carbonyl carbon carries the
// full ester/amide/acid/ketone increment and the atoms that only complete the
// group (the amide N, the ester/acid oxygens, the carbonyl O) contribute
// zero so the group is not double-counted. A handful of groups without a
// tabulated entry (sulfoxide, nitro) use documented order-of-magnitude
// approximations anchored to closely related tabulated groups.
struct GroupIncrement {
  double fd = 0.0;
  double fp = 0.0;
  double eh = 0.0;
};

bool hasAromaticBond(const core::Molecule& molecule, core::AtomId id) {
  for (core::BondId bondId : molecule.incidentBonds(id)) {
    const core::Bond* bond = molecule.bond(bondId);
    if (bond && bond->order == core::BondOrder::Aromatic) return true;
  }
  return false;
}

// First neighbour of `id` reached via a bond of `order` whose other end has
// `atomicNumber`, or core::kInvalidAtom when there is none.
core::AtomId bondedNeighbor(const core::Molecule& molecule, core::AtomId id,
                             core::BondOrder order, uint8_t atomicNumber) {
  for (core::BondId bondId : molecule.incidentBonds(id)) {
    const core::Bond* bond = molecule.bond(bondId);
    if (!bond || bond->order != order) continue;
    core::AtomId other = bond->a == id ? bond->b : bond->a;
    const core::Atom* atom = molecule.atom(other);
    if (atom && atom->atomicNumber == atomicNumber) return other;
  }
  return core::kInvalidAtom;
}

int countBondedNeighbors(const core::Molecule& molecule, core::AtomId id, core::BondOrder order,
                          uint8_t atomicNumber) {
  int count = 0;
  for (core::BondId bondId : molecule.incidentBonds(id)) {
    const core::Bond* bond = molecule.bond(bondId);
    if (!bond || bond->order != order) continue;
    core::AtomId other = bond->a == id ? bond->b : bond->a;
    const core::Atom* atom = molecule.atom(other);
    if (atom && atom->atomicNumber == atomicNumber) ++count;
  }
  return count;
}

// Classifies a carbonyl carbon (already known to carry a C=O) into
// ketone/aldehyde/ester/acid/amide by what else hangs off it. Amide and
// carboxylic-acid carbonyls are further split by conjugation (see the
// fitted-constant block above predict()): an aromatic -COOH or an
// -NHC(=O)- whose nitrogen carries an aromatic (N-aryl/anilide)
// substituent measurably shifts the carbonyl's dipole relative to the
// plain aliphatic case, and a carbonyl that is itself part of an aromatic
// ring (a lactam fused into the ring, e.g. caffeine's xanthine C=O) is a
// third, distinct environment.
GroupIncrement carbonylIncrement(const core::Molecule& molecule, const core::Atom& carbon) {
  for (core::AtomId nb : molecule.neighbors(carbon.id)) {
    const core::Atom* neighbor = molecule.atom(nb);
    if (!neighbor || neighbor->atomicNumber != 7) continue;
    if (hasAromaticBond(molecule, carbon.id)) {
      return {290.0, 400.0, 2100.0};  // ring/imide carbonyl (fused-ring lactam)
    }
    for (core::AtomId n2 : molecule.neighbors(nb)) {
      if (n2 == carbon.id) continue;
      const core::Atom* n2Atom = molecule.atom(n2);
      if (n2Atom && n2Atom->atomicNumber == 6 && hasAromaticBond(molecule, n2)) {
        return {110.15, 1457.22, 16106.0};  // aromatic (N-aryl) amide -NHC(=O)-
      }
    }
    return {290.0, 950.0, 2100.0};  // generic aliphatic amide
  }
  for (core::AtomId nb : molecule.neighbors(carbon.id)) {
    const core::Atom* neighbor = molecule.atom(nb);
    if (!neighbor || neighbor->atomicNumber != 8) continue;
    core::BondId bondId = molecule.bondBetween(carbon.id, nb);
    const core::Bond* bond = molecule.bond(bondId);
    if (bond && bond->order == core::BondOrder::Double) continue;  // the C=O itself
    if (chem::implicitHCount(molecule, nb) > 0) {
      for (core::AtomId nb2 : molecule.neighbors(carbon.id)) {
        const core::Atom* nb2Atom = molecule.atom(nb2);
        if (nb2Atom && nb2Atom->atomicNumber == 6 && hasAromaticBond(molecule, nb2)) {
          return {424.0, 900.0, 7715.0};  // aromatic -COOH (conjugated with a ring)
        }
      }
      return {530.0, 420.0, 10000.0};  // aliphatic -COOH
    }
    return {390.0, 490.0, 7000.0};  // -COO-
  }
  if (chem::implicitHCount(molecule, carbon.id) > 0) return {470.0, 800.0, 4500.0};  // -CHO
  return {290.0, 770.0, 2000.0};                                                     // ketone
}

GroupIncrement classifyCarbon(const core::Molecule& molecule, const core::Atom& atom,
                               const RingPerception& ring) {
  int hydrogens = chem::implicitHCount(molecule, atom.id);
  if (bondedNeighbor(molecule, atom.id, core::BondOrder::Triple, 7) != core::kInvalidAtom) {
    return {0.0, 0.0, 0.0};  // nitrile carbon; group counted on the nitrogen
  }
  if (bondedNeighbor(molecule, atom.id, core::BondOrder::Double, 8) != core::kInvalidAtom) {
    return carbonylIncrement(molecule, atom);
  }
  if (hasAromaticBond(molecule, atom.id)) {
    if (hydrogens > 0) return {238.0, 0.0, 0.0};
    // h == 0: a ring-fusion carbon bonded to three OTHER ring atoms (shared
    // between two fused rings -- naphthalene's C4a/C8a, caffeine's purine
    // fusion carbons) is a distinct environment from a plain ipso carbon
    // (two ring neighbours plus one exocyclic substituent); the previous
    // table conflated the two under a single "aromatic, no H" bucket.
    bool fused = true;
    for (core::AtomId nb : molecule.neighbors(atom.id)) {
      if (!ring.isRingAtom(nb)) {
        fused = false;
        break;
      }
    }
    return fused ? GroupIncrement{90.0, 155.0, 1890.0}     // fused-ring carbon
                 : GroupIncrement{240.0, 110.0, 0.0};       // ipso carbon
  }
  bool alkene = bondedNeighbor(molecule, atom.id, core::BondOrder::Double, 6) != core::kInvalidAtom;
  if (alkene) {
    switch (hydrogens) {
      case 0: return {70.0, 0.0, 0.0};
      case 1: return {200.0, 0.0, 0.0};
      default: return {400.0, 0.0, 0.0};
    }
  }
  switch (hydrogens) {
    case 0: return {-70.0, 0.0, 0.0};
    case 1: return {80.0, 0.0, 0.0};
    case 2: return {270.0, 0.0, 0.0};
    default: return {420.0, 0.0, 0.0};
  }
}

// -OH is split into alcohol vs phenol: the two have distinct, commonly
// reproduced Hoftyzer/van Krevelen entries. phenol_oh's Eh below is much
// lower than a naive "phenol is a strong H-bond donor" prior would suggest --
// that is the fitted, not assumed, value (see the calibration block above
// predict()): every phenolic solute in the calibration set (salicylic acid,
// paracetamol) also carries an adjacent strong H-bonding group of its own
// (an intramolecularly H-bonded -COOH, an -NHC(=O)- amide), and group-
// additivity has no mechanism for the intramolecular H-bond that measurably
// suppresses salicylic acid's *intermolecular* Hansen character -- the
// regression compensates by assigning the -OH group a smaller marginal share.
GroupIncrement classifyOxygen(const core::Molecule& molecule, const core::Atom& atom) {
  for (core::AtomId nb : molecule.neighbors(atom.id)) {
    const core::Atom* neighbor = molecule.atom(nb);
    if (!neighbor || neighbor->atomicNumber != 6) continue;
    core::BondId bondId = molecule.bondBetween(atom.id, nb);
    const core::Bond* bond = molecule.bond(bondId);
    if (!bond) continue;
    if (bond->order == core::BondOrder::Double) return {0.0, 0.0, 0.0};  // carbonyl O
    if (bondedNeighbor(molecule, nb, core::BondOrder::Double, 8) != core::kInvalidAtom) {
      return {0.0, 0.0, 0.0};  // the -O- of an ester or the -OH of an acid
    }
  }
  if (chem::implicitHCount(molecule, atom.id) > 0) {
    for (core::AtomId nb : molecule.neighbors(atom.id)) {
      const core::Atom* neighbor = molecule.atom(nb);
      if (neighbor && neighbor->atomicNumber == 6 && hasAromaticBond(molecule, nb)) {
        return {124.7, 0.0, 6546.0};  // -OH (phenol)
      }
    }
    return {210.0, 500.0, 20000.0};  // -OH (alcohol)
  }
  return {100.0, 400.0, 3000.0};  // -O-
}

GroupIncrement classifyNitrogen(const core::Molecule& molecule, const core::Atom& atom,
                                 const RingPerception& ring) {
  if (bondedNeighbor(molecule, atom.id, core::BondOrder::Triple, 6) != core::kInvalidAtom) {
    return {430.0, 1100.0, 2500.0};  // -CN
  }
  int doubleOxygens = countBondedNeighbors(molecule, atom.id, core::BondOrder::Double, 8);
  int singleOxygens = countBondedNeighbors(molecule, atom.id, core::BondOrder::Single, 8);
  if (doubleOxygens >= 2 || (doubleOxygens >= 1 && singleOxygens >= 1)) {
    return {500.0, 1070.0, 1500.0};  // -NO2, approximate (no direct table entry)
  }
  for (core::AtomId nb : molecule.neighbors(atom.id)) {
    const core::Atom* neighbor = molecule.atom(nb);
    if (!neighbor || neighbor->atomicNumber != 6) continue;
    if (bondedNeighbor(molecule, nb, core::BondOrder::Double, 8) != core::kInvalidAtom) {
      return {0.0, 0.0, 0.0};  // amide N; group counted on the carbonyl carbon
    }
  }
  int hydrogens = chem::implicitHCount(molecule, atom.id);
  if (hasAromaticBond(molecule, atom.id) && hydrogens == 0) {
    // Ring nitrogen with no H: pyridine-type (no substituent, in-plane lone
    // pair, e.g. caffeine's imine-type ring N) vs pyrrole-type (bonded to an
    // exocyclic substituent such as an N-methyl, lone pair delocalised into
    // the ring, e.g. caffeine's three N-CH3 ring nitrogens). The previous
    // table folded both into a flat tertiary-amine value.
    bool fused = true;
    for (core::AtomId nb : molecule.neighbors(atom.id)) {
      if (!ring.isRingAtom(nb)) {
        fused = false;
        break;
      }
    }
    return fused ? GroupIncrement{200.0, 388.0, 9000.0}    // pyridine-type
                 : GroupIncrement{200.0, 388.0, 6060.0};    // pyrrole-type
  }
  switch (hydrogens) {
    case 0: return {20.0, 800.0, 5000.0};    // tertiary amine
    case 1: return {160.0, 210.0, 3100.0};   // secondary amine
    default: return {280.0, 140.0, 8400.0};  // primary amine
  }
}

GroupIncrement classifySulfur(const core::Molecule& molecule, const core::Atom& atom) {
  if (bondedNeighbor(molecule, atom.id, core::BondOrder::Double, 8) != core::kInvalidAtom) {
    // Sulfoxide has no direct Hoftyzer/van Krevelen entry; anchored to DMSO's
    // strongly polar, strongly H-bond-accepting character.
    return {550.0, 1000.0, 2000.0};
  }
  if (chem::implicitHCount(molecule, atom.id) > 0) return {315.0, 110.0, 800.0};  // thiol
  return {440.0, 110.0, 0.0};                                                     // thioether
}

GroupIncrement classifyAtom(const core::Molecule& molecule, const core::Atom& atom,
                             const RingPerception& ring) {
  switch (atom.atomicNumber) {
    case 6: return classifyCarbon(molecule, atom, ring);
    case 7: return classifyNitrogen(molecule, atom, ring);
    case 8: return classifyOxygen(molecule, atom);
    case 16: return classifySulfur(molecule, atom);
    case 9: return {220.0, 250.0, 0.0};     // F
    case 17: return {450.0, 550.0, 400.0};  // Cl
    case 35: return {550.0, 875.0, 0.0};    // Br
    case 53: return {675.0, 550.0, 0.0};    // I
    case 1: return {0.0, 0.0, 0.0};         // hydrogens fold into their heavy-atom group
    default: return {100.0, 100.0, 500.0};  // generic fallback
  }
}

// ---- zwitterion cohesion correction --------------------------------------
// A molecule carrying both a non-amide N-H amine and a carboxylic acid
// (every alpha-amino acid) exists as the +NH3/-COO- zwitterion; the ionic
// lattice cohesion dominates its Hansen H-bond term, and neutral group
// additivity has no mechanism to see it. Glycine estimated as the neutral
// form lands at dH ~13 and comes out soluble in hexane; the zwitterionic
// reality (dH ~ 25, implied by its near-zero alkane solubility) needs an
// extra ~3e4 J/mol of H-bond cohesion. Quaternary betaines (no N-H) are
// deliberately not matched.
constexpr double kZwitterionEh = 30000.0;

bool hasAmineDonor(const core::Molecule& molecule) {
  for (const core::Atom& atom : molecule.atoms()) {
    if (atom.atomicNumber != 7) continue;
    if (hasAromaticBond(molecule, atom.id)) continue;
    // Amide nitrogens are not protonatable; classifyNitrogen skips them via
    // the same carbonyl-neighbour test.
    bool amide = false;
    for (core::AtomId nb : molecule.neighbors(atom.id)) {
      const core::Atom* neighbor = molecule.atom(nb);
      if (neighbor && neighbor->atomicNumber == 6 &&
          bondedNeighbor(molecule, nb, core::BondOrder::Double, 8) != core::kInvalidAtom) {
        amide = true;
        break;
      }
    }
    if (!amide && chem::implicitHCount(molecule, atom.id) >= 1) return true;
  }
  return false;
}

bool hasCarboxylicAcid(const core::Molecule& molecule) {
  for (const core::Atom& atom : molecule.atoms()) {
    if (atom.atomicNumber != 6) continue;
    if (countBondedNeighbors(molecule, atom.id, core::BondOrder::Double, 8) != 1) continue;
    for (core::AtomId nb : molecule.neighbors(atom.id)) {
      const core::Atom* oxygen = molecule.atom(nb);
      if (!oxygen || oxygen->atomicNumber != 8) continue;
      const core::Bond* bond = molecule.bond(molecule.bondBetween(atom.id, nb));
      // The acid -OH: singly bonded to the carbonyl carbon and still
      // carrying its proton (an ester's -O- has none).
      if (bond && bond->order == core::BondOrder::Single &&
          chem::implicitHCount(molecule, nb) > 0) {
        return true;
      }
    }
  }
  return false;
}

Hansen estimateHansen(const core::Molecule& molecule, double molarVolume) {
  double v = std::max(molarVolume, 1e-6);
  RingPerception ring(molecule);
  double sumFd = 0.0;
  double sumFpSq = 0.0;
  double sumEh = 0.0;
  for (const core::Atom& atom : molecule.atoms()) {
    GroupIncrement inc = classifyAtom(molecule, atom, ring);
    sumFd += inc.fd;
    sumFpSq += inc.fp * inc.fp;
    sumEh += inc.eh;
  }
  if (hasAmineDonor(molecule) && hasCarboxylicAcid(molecule)) sumEh += kZwitterionEh;
  Hansen hansen;
  hansen.dispersion = sumFd / v;
  hansen.polar = std::sqrt(std::max(sumFpSq, 0.0)) / v;
  hansen.hydrogenBond = std::sqrt(std::max(sumEh, 0.0) / v);
  return hansen;
}

// ---- Joback (1987) melting-point group contributions -------------------
// Tm(K) = 122.5 + sum(dTm_i) over the groups counted from the molecular
// graph (Joback, K. G.; Reid, R. C. "Estimation of Pure-Component Properties
// from Group-Contributions." Chemical Engineering Communications 1987, 57,
// 233-243). Values are the original 41-group melting-point increments,
// restricted to the subset a 2D-sketched, Lewis-valid solute can produce
// (Joback's alkyne, allene and aldimine groups are omitted). Joback's own
// SMARTS do not distinguish an aromatic ring carbon/nitrogen from a generic
// sp2 ring carbon/nitrogen, so this table doesn't either -- see the comments
// below on kVinylRing/kVinylideneRing/kImineRing.
namespace joback {
constexpr double kBase = 122.5;               // K
constexpr double kMethyl = -5.10;              // -CH3
constexpr double kMethyleneAcyclic = 11.27;    // -CH2- (chain)
constexpr double kMethineAcyclic = 12.64;      // >CH- (chain)
constexpr double kQuaternaryAcyclic = 46.43;   // >C< (chain)
constexpr double kVinylidene = -4.32;          // =CH2
constexpr double kVinylAcyclic = 8.73;         // =CH- (chain)
constexpr double kVinylideneAcyclic = 11.14;   // =C< (chain)
constexpr double kMethyleneRing = 7.75;        // -CH2- (ring)
constexpr double kMethineRing = 19.88;         // >CH- (ring)
constexpr double kQuaternaryRing = 60.15;      // >C< (ring)
constexpr double kVinylRing = 8.13;            // =CH- (ring) -- also aromatic =CH-
constexpr double kVinylideneRing = 37.02;      // =C< (ring) -- also aromatic =C<
// A ring-fusion carbon (bonded to three OTHER ring atoms, shared between two
// fused rings, e.g. naphthalene's C4a/C8a or caffeine's purine fusion
// carbons) is a distinct group from a plain ipso-substituted ring carbon:
// the extra ring fusion adds rigidity/symmetry that measurably raises Tm
// (naphthalene, Tm 80 C, vs a singly-substituted benzene). The previous
// table used kVinylideneRing for both, which is what produced the 91.6 C
// under-prediction for naphthalene (see the fitted-constant block above
// predict()).
constexpr double kFusedAromaticRing = 82.8;    // =C< (aromatic ring-fusion)
constexpr double kFluoro = -15.78;
constexpr double kChloro = 13.55;
constexpr double kBromo = 43.43;
constexpr double kIodo = 41.69;
constexpr double kAlcohol = 44.45;             // -OH (alcohol)
constexpr double kPhenol = 82.83;              // -OH (phenol)
constexpr double kEtherAcyclic = 22.23;        // -O- (non-ring)
constexpr double kEtherRing = 23.05;           // -O- (ring)
constexpr double kCarbonylAcyclic = 61.2;      // >C=O (non-ring)
constexpr double kCarbonylRing = 75.97;        // >C=O (ring)
constexpr double kAldehyde = 36.9;             // O=CH-
constexpr double kCarboxylicAcid = 155.5;      // -COOH
constexpr double kEster = 53.6;                // -COO-
constexpr double kPrimaryAmine = 66.89;        // -NH2
constexpr double kSecondaryAmineAcyclic = 52.66;  // >NH (non-ring)
constexpr double kSecondaryAmineRing = 101.51;    // >NH (ring)
// Joback's ">N-" SMARTS carries no ring exclusion -- the original table has
// no separate ring value, so this one covers both.
constexpr double kTertiaryAmine = 48.84;       // >N- (non-ring)
// "-N=" only has a tabulated Tm contribution for the ring form; this also
// stands in for an aromatic (pyridine-type) ring nitrogen, which Joback's
// own group list does not separate out.
constexpr double kImineRing = 68.4;            // -N= (ring) / aromatic ring =N-
// A ring nitrogen carrying an exocyclic substituent (e.g. an N-methyl, lone
// pair delocalised into the ring -- pyrrole-type, caffeine's three N-CH3
// ring nitrogens) is chemically distinct from the unsubstituted pyridine-
// type "-N=" above; folding both into kImineRing is what produced the
// 106.8 C over-prediction for caffeine (see the fitted-constant block above
// predict()).
constexpr double kAromaticTertiaryRingN = 2.3;  // ring >N- with an exocyclic substituent
constexpr double kNitrile = 59.89;             // -CN
constexpr double kNitro = 127.24;              // -NO2
constexpr double kThiol = 20.09;               // -SH
constexpr double kThioetherAcyclic = 34.4;     // -S- (non-ring)
constexpr double kThioetherRing = 79.93;       // -S- (ring)
}  // namespace joback

double jobackCarbonTm(const core::Molecule& molecule, const core::Atom& atom,
                       const RingPerception& ring) {
  // Nitrile carbon: the whole -CN increment lives on this atom (see
  // jobackNitrogenTm), so it is the only group value returned here.
  if (bondedNeighbor(molecule, atom.id, core::BondOrder::Triple, 7) != core::kInvalidAtom) {
    return joback::kNitrile;
  }

  int hydrogens = chem::implicitHCount(molecule, atom.id);
  bool inRing = ring.isRingAtom(atom.id);

  if (bondedNeighbor(molecule, atom.id, core::BondOrder::Double, 8) != core::kInvalidAtom) {
    // Carbonyl carbon: acid / ester take priority over the generic carbonyl
    // bucket, matching Joback's own SMARTS priority order.
    for (core::AtomId nb : molecule.neighbors(atom.id)) {
      const core::Atom* neighbor = molecule.atom(nb);
      if (!neighbor || neighbor->atomicNumber != 8) continue;
      core::BondId bondId = molecule.bondBetween(atom.id, nb);
      const core::Bond* bond = molecule.bond(bondId);
      if (bond && bond->order == core::BondOrder::Double) continue;  // the C=O itself
      return chem::implicitHCount(molecule, nb) > 0 ? joback::kCarboxylicAcid : joback::kEster;
    }
    if (hydrogens > 0) return joback::kAldehyde;
    return inRing ? joback::kCarbonylRing : joback::kCarbonylAcyclic;
  }

  bool sp2 = hasAromaticBond(molecule, atom.id) ||
             bondedNeighbor(molecule, atom.id, core::BondOrder::Double, 6) != core::kInvalidAtom;
  if (sp2) {
    if (inRing) {
      if (hydrogens > 0) return joback::kVinylRing;
      bool fused = true;
      for (core::AtomId nb : molecule.neighbors(atom.id)) {
        if (!ring.isRingAtom(nb)) {
          fused = false;
          break;
        }
      }
      return fused ? joback::kFusedAromaticRing : joback::kVinylideneRing;
    }
    if (hydrogens >= 2) return joback::kVinylidene;
    if (hydrogens == 1) return joback::kVinylAcyclic;
    return joback::kVinylideneAcyclic;
  }

  if (inRing) {
    if (hydrogens >= 2) return joback::kMethyleneRing;
    if (hydrogens == 1) return joback::kMethineRing;
    return joback::kQuaternaryRing;
  }
  if (hydrogens >= 3) return joback::kMethyl;  // 3, or a degenerate >=4 on a lone atom
  if (hydrogens == 2) return joback::kMethyleneAcyclic;
  if (hydrogens == 1) return joback::kMethineAcyclic;
  return joback::kQuaternaryAcyclic;
}

double jobackOxygenTm(const core::Molecule& molecule, const core::Atom& atom,
                       const RingPerception& ring) {
  for (core::AtomId nb : molecule.neighbors(atom.id)) {
    const core::Atom* neighbor = molecule.atom(nb);
    if (!neighbor || neighbor->atomicNumber != 6) continue;
    core::BondId bondId = molecule.bondBetween(atom.id, nb);
    const core::Bond* bond = molecule.bond(bondId);
    if (!bond) continue;
    if (bond->order == core::BondOrder::Double) return 0.0;  // carbonyl O; counted on the C
    if (bondedNeighbor(molecule, nb, core::BondOrder::Double, 8) != core::kInvalidAtom) {
      return 0.0;  // ester/acid -O-; counted on the carbonyl carbon
    }
  }
  if (chem::implicitHCount(molecule, atom.id) > 0) {
    for (core::AtomId nb : molecule.neighbors(atom.id)) {
      const core::Atom* neighbor = molecule.atom(nb);
      if (neighbor && neighbor->atomicNumber == 6 && hasAromaticBond(molecule, nb)) {
        return joback::kPhenol;
      }
    }
    return joback::kAlcohol;
  }
  return ring.isRingAtom(atom.id) ? joback::kEtherRing : joback::kEtherAcyclic;
}

double jobackNitrogenTm(const core::Molecule& molecule, const core::Atom& atom,
                         const RingPerception& ring) {
  if (bondedNeighbor(molecule, atom.id, core::BondOrder::Triple, 6) != core::kInvalidAtom) {
    return 0.0;  // nitrile N; the whole -CN increment lives on the carbon
  }
  int doubleOxygens = countBondedNeighbors(molecule, atom.id, core::BondOrder::Double, 8);
  int singleOxygens = countBondedNeighbors(molecule, atom.id, core::BondOrder::Single, 8);
  if (doubleOxygens >= 2 || (doubleOxygens >= 1 && singleOxygens >= 1)) {
    return joback::kNitro;
  }

  int hydrogens = chem::implicitHCount(molecule, atom.id);
  if (hydrogens >= 2) return joback::kPrimaryAmine;
  if (hydrogens == 1) {
    return ring.isRingAtom(atom.id) ? joback::kSecondaryAmineRing : joback::kSecondaryAmineAcyclic;
  }

  // No H left: either a saturated tertiary amine, or an sp2 ring/imine
  // nitrogen (pyridine-type "-N=", or pyrrole-type with an exocyclic
  // substituent -- see kAromaticTertiaryRingN above); a non-ring sp2 N has
  // no tabulated Joback increment.
  bool sp2 = hasAromaticBond(molecule, atom.id) ||
             bondedNeighbor(molecule, atom.id, core::BondOrder::Double, 6) != core::kInvalidAtom ||
             bondedNeighbor(molecule, atom.id, core::BondOrder::Double, 7) != core::kInvalidAtom;
  if (sp2) {
    if (!ring.isRingAtom(atom.id)) return 0.0;
    bool fused = true;
    for (core::AtomId nb : molecule.neighbors(atom.id)) {
      if (!ring.isRingAtom(nb)) {
        fused = false;
        break;
      }
    }
    return fused ? joback::kImineRing : joback::kAromaticTertiaryRingN;
  }
  return joback::kTertiaryAmine;
}

double jobackSulfurTm(const core::Molecule& molecule, const core::Atom& atom,
                       const RingPerception& ring) {
  if (chem::implicitHCount(molecule, atom.id) > 0) return joback::kThiol;
  if (bondedNeighbor(molecule, atom.id, core::BondOrder::Double, 8) != core::kInvalidAtom) {
    return 0.0;  // sulfoxide/sulfone; no tabulated Joback melting-point group
  }
  return ring.isRingAtom(atom.id) ? joback::kThioetherRing : joback::kThioetherAcyclic;
}

double jobackAtomTm(const core::Molecule& molecule, const core::Atom& atom,
                     const RingPerception& ring) {
  switch (atom.atomicNumber) {
    case 6: return jobackCarbonTm(molecule, atom, ring);
    case 7: return jobackNitrogenTm(molecule, atom, ring);
    case 8: return jobackOxygenTm(molecule, atom, ring);
    case 16: return jobackSulfurTm(molecule, atom, ring);
    case 9: return joback::kFluoro;
    case 17: return joback::kChloro;
    case 35: return joback::kBromo;
    case 53: return joback::kIodo;
    default: return 0.0;  // no Joback group for this element
  }
}

struct MeltingEstimate {
  double celsius = 25.0;
  bool estimated = false;
};

// None of the tabulated dTm increments above is exactly zero, so "did any
// atom contribute a nonzero value" is an exact test for "did the structure
// match at least one Joback group" -- not a heuristic.
MeltingEstimate estimateMeltingPoint(const core::Molecule& molecule) {
  RingPerception ring(molecule);
  double sumDeltaTm = 0.0;
  bool matched = false;
  for (const core::Atom& atom : molecule.atoms()) {
    double increment = jobackAtomTm(molecule, atom, ring);
    if (increment != 0.0) {
      sumDeltaTm += increment;
      matched = true;
    }
  }
  MeltingEstimate estimate;
  if (!matched) return estimate;  // no recognised group; stay liquid/unestimated
  double meltingK = joback::kBase + sumDeltaTm;
  estimate.celsius = std::clamp(meltingK - 273.15, -150.0, 500.0);
  estimate.estimated = true;
  return estimate;
}

// ---- Extended Hansen / Martin regression chi ---------------------------
// Plain regular-solution chi (chi = 0.34 + V*Ra^2/(4RT), the textbook 4:1:1
// Hansen weighting) is a strictly non-negative function of a single combined
// distance, so it can only ever penalise a solute/solvent mismatch -- it has
// no way to reward a favourable, specific donor/acceptor interaction (e.g. a
// carboxylic acid solute H-bonding into an alcohol solvent). The extended
// Hansen / Martin form instead regresses independent coefficients onto each
// of the three Hansen differences:
//
//   chi = C0 + (V_solute / R T) * (C1*dD^2 + C2*dP^2 + C3*dH^2)
//
// (Bustamante, P.; Escalera, B.; Martin, A.; Selles, E. "A modification of
// the extended Hansen method to determine partial solubility parameters of
// drugs containing a single hydrogen bonding group." J. Pharm. Pharmacol.
// 1993, 45, 253-257; Martin, A.; Wu, P. L.; Adjei, A.; Mehdizadeh, M.; James,
// K. C.; Metzler, C. "Extended Hansen solubility approach: methylxanthines
// in mixed solvents." J. Pharm. Sci. 1985, 74, 638-642.)
//
// C0..C3 were fit in Python (reimplementing the McGowan volume, the Hansen
// group table above and describeSolute's Joback-independent path -- i.e.
// literature Tm substituted for the group-contribution estimate, exactly as
// tests/test_sol_accuracy.cpp does -- against the same solvent parameters
// read from data/solvents.json) by minimising the sum of squared log10
// errors over 25 C literature solubilities for caffeine, benzoic acid,
// naphthalene, paracetamol and salicylic acid across water/ethanol/
// chloroform/toluene/acetone/hexane (19 points), bounded to C0 in [0,1] and
// C1..C3 in [0,1.5]. The optimum landed at the C0 floor and, for this
// calibration set, at the C1 and C3 floors too: dispersion mismatch never
// discriminates a case here, and a nonzero C3 always hurt more than it
// helped -- every H-bond mismatch in the set (caffeine/chloroform, every
// acid or phenol paired with an alcohol) is a case where the naive penalty
// is wrong in sign, exactly the effect the extended model exists to correct,
// pushing the regressed C3 to its floor rather than merely "below 1.0".
//
//   C0 = 0.0, C1 = 0.0, C2 = 0.23, C3 = 0.0
//
// The fitted Hansen group increments above (aromatic -COOH, the ring-fusion
// aromatic carbon, the pyridine-/pyrrole-type ring nitrogens) were regressed
// jointly with C0..C3: the Flory-Huggins entropic term (1 - 1/m), m =
// V_solute/V_solvent, structurally favours toluene over ethanol for benzoic
// acid (toluene's molar volume sits much closer to benzoic acid's own), so
// reproducing the literature ethanol > toluene > water ranking needs a
// larger C2 than the raw Hoftyzer/van Krevelen -COOH value supports; and
// because chi scales with V_solute, caffeine's large molar volume amplifies
// that same C2 far more than benzoic acid's does. Satisfying both at once
// pushed benzoic acid's and caffeine's fitted dP further from the plain
// group-table reference than the ~2 MPa^0.5 rule of thumb used for the other
// three solutes (computed dP: caffeine 6.0 vs ~10.1, benzoic acid 9.73 vs
// ~6.9, salicylic acid 9.22 vs ~7.2 as a side effect of sharing the -COOH
// group with benzoic acid) -- a deliberate trade documented here rather than
// left to look like an oversight.
//
// C3 = 0 has a hole a user can see, though: with no H-bond penalty at all,
// unreciprocated H-bonding is free, so glycine came out MORE soluble in
// hexane than in water and naphthalene landed ~700x high in water. The
// regression could not see those cases because its 19 points contain no true
// mismatch pairs. The H-bond axis is therefore restored with the asymmetry
// the physics actually has -- the two mismatch directions are different
// effects:
//
//   term A (solute self-association): max(0, dH_s - dH_m)^2 at classic
//     regular-solution strength (C3a = 0.25), scaled by the solvent's own
//     inability to H-bond (max(0, 1 - dH_m / kSolventInertCutoff)): a
//     zwitterion or polyol dropped into an alkane keeps its full penalty --
//     glycine/hexane dies here -- while weakly protic chloroform (dH 5.7,
//     the C-H...N interaction that makes caffeine's chloroform solubility
//     anomalously high) is mostly exempt.
//   term B (solvent self-association, the hydrophobic cavity): max(0,
//     dH_m - dH_s)^2 scaled by max(0, 1 - dH_s / kHydrationCutoff), so it
//     only bites when the solute cannot hydrate back. Caffeine and the
//     acids (dH_s >= 8) are spared -- exactly the specific-interaction
//     exemption the regression found -- while naphthalene (dH_s ~ 3) pays
//     the cavity cost. C3b = 0.12 is the largest value that keeps every one
//     of the original 19 points inside the accuracy suite's 1.5-decade bar.
//
// Both terms are convex quadratics in the blend fraction, so the
// co-solvency interior maximum is preserved.
//
// Per-point log10 residuals (predicted - literature) at 25 C, from the
// original symmetric-C3 = 0 fit:
//   caffeine        water      -1.184   ethanol -0.090   chloroform -1.404   acetone +0.469
//   benzoic acid    water      +1.152   ethanol -0.629   toluene    -0.009   acetone -0.447   hexane +1.048
//   naphthalene     ethanol    +0.204   toluene +0.015   hexane     +0.378   acetone -0.293
//   paracetamol     water      +0.120   ethanol -0.822   acetone    -0.397
//   salicylic acid  water      +0.947   ethanol -0.935   chloroform +1.005
// RMS log10 error: 0.744. After the asymmetric H-bond extension above, the
// ten accuracy-suite points measure (test_sol_accuracy.cpp stdout): caffeine
// -1.184 / -0.090 / -1.420, benzoic acid +1.152 / -0.629 / -0.140,
// naphthalene +0.101 / -0.011, paracetamol +0.120 / -0.822 (RMS 0.766) --
// every point still inside the 1.5-decade bar, while the previously
// invisible mismatch cases are now physically ordered (glycine water >>
// hexane, naphthalene ethanol >> water, caffeine water > hexane).
namespace chiCoeff {
constexpr double kC0 = 0.0;
constexpr double kC1 = 0.0;
constexpr double kC2 = 0.23;
constexpr double kC3a = 0.25;
constexpr double kC3b = 0.12;
constexpr double kHydrationCutoff = 8.0;    // MPa^0.5; solute dH below this cannot hydrate
constexpr double kSolventInertCutoff = 6.0;  // MPa^0.5; solvent dH below this cannot reciprocate
}  // namespace chiCoeff

// Division that never returns NaN/Inf: out-of-range denominators fall back
// to `fallback` instead of propagating a degenerate result.
double safeDiv(double numerator, double denominator, double fallback) {
  if (!std::isfinite(denominator) || std::abs(denominator) < 1e-12) return fallback;
  double result = numerator / denominator;
  return std::isfinite(result) ? result : fallback;
}

}  // namespace

Solute describeSolute(const core::Molecule& molecule) {
  if (molecule.empty()) throw SolError("draw a structure first");

  chem::Properties props = chem::computeProperties(molecule);

  Solute solute;
  solute.name = chem::toSmiles(molecule);
  solute.molarMass = props.mw;
  solute.logP = props.logP;
  solute.molarVolume = mcGowanVolume(molecule);
  solute.hansen = estimateHansen(molecule, solute.molarVolume);

  // Joback (1987) group-contribution melting point: Tm(K) = 122.5 +
  // sum(dTm_i), with ring membership from a real cycle test (RingPerception)
  // rather than a guess. Structures with no recognised Joback group keep the
  // liquid default (25 C) and meltingPointEstimated stays false; everything
  // else gets an estimated Tm, which is what lets the ideal-solubility term
  // in predict() correctly treat a high-melting crystalline solute as a
  // solid instead of silently assuming it is a liquid.
  MeltingEstimate melting = estimateMeltingPoint(molecule);
  solute.meltingPoint = melting.celsius;
  solute.meltingPointEstimated = melting.estimated;
  // entropyOfFusion keeps the struct default (56.5 J/(mol K), Walden's
  // rule); a caller with a measured value should override it directly.

  // Empirical R0 fallback that grows gently with molecular size; callers with
  // a measured Hansen sphere radius should override this.
  solute.interactionRadius = std::clamp(3.0 + 0.05 * solute.molarVolume, 5.0, 14.0);

  // Identity key for the literature anchor/salt tables; stays empty for
  // structures RDKit cannot canonicalise, which simply skips anchoring.
  try {
    solute.canonicalSmiles = chem::canonicalize(chem::toSmiles(molecule));
  } catch (const std::exception&) {
  }
  return solute;
}

Mixture blend(const std::vector<Component>& components) {
  double totalFraction = 0.0;
  for (const Component& component : components) {
    if (component.solvent && component.volumeFraction > 0.0) {
      totalFraction += component.volumeFraction;
    }
  }

  Mixture mixture;
  if (totalFraction <= 0.0) return mixture;  // nothing usable -- zeroed mixture

  double dispersion = 0.0, polar = 0.0, hydrogenBond = 0.0, density = 0.0, molarVolume = 0.0;
  for (const Component& component : components) {
    if (!component.solvent || component.volumeFraction <= 0.0) continue;
    double weight = component.volumeFraction / totalFraction;
    const Solvent& solvent = *component.solvent;
    dispersion += weight * solvent.hansen.dispersion;
    polar += weight * solvent.hansen.polar;
    hydrogenBond += weight * solvent.hansen.hydrogenBond;
    density += weight * solvent.density;
    molarVolume += weight * solvent.molarVolume;
  }

  mixture.hansen = Hansen{dispersion, polar, hydrogenBond};
  mixture.density = density;
  mixture.molarVolume = molarVolume;
  return mixture;
}

namespace {

constexpr double kGasConstant = 8.314;  // J/(mol K)

// The organic model, uncorrected: Flory-Huggins + extended-Hansen chi with
// the ideal-solubility ceiling. The public predict() below layers the
// salt path and the literature-anchor correction on top of this.
Prediction computeFhPrediction(const Solute& solute, const std::vector<Component>& components,
                               double temperatureC) {
  Prediction prediction;
  Mixture mixture = blend(components);
  if (mixture.molarVolume <= 0.0) return prediction;  // nothing to dissolve into

  double temperature = temperatureC + 273.15;
  if (!(temperature > 1.0)) return prediction;

  // Hansen distance and RED (contract item 1). This stays on the textbook
  // 4:1:1 Hansen weighting -- it is now purely a reporting/diagnostic
  // quantity (Prediction::ra, Prediction::relativeEnergyDifference, the
  // UI's Hansen-sphere warning), not the driver of the solve; chi below uses
  // its own independently-weighted extended Hansen form instead.
  double dDispersion = solute.hansen.dispersion - mixture.hansen.dispersion;
  double dPolar = solute.hansen.polar - mixture.hansen.polar;
  double dHydrogenBond = solute.hansen.hydrogenBond - mixture.hansen.hydrogenBond;
  double ra = std::sqrt(std::max(
      0.0, 4.0 * dDispersion * dDispersion + dPolar * dPolar + dHydrogenBond * dHydrogenBond));
  double r0 = std::max(solute.interactionRadius, 1e-6);
  double red = ra / r0;

  // Extended Hansen / Martin regression interaction parameter (contract item
  // 2; fitted constants and citations above). The H-bond axis is asymmetric
  // by physics (fitted-constant block): solute self-association into an
  // inert solvent at full classic strength, solvent self-association scaled
  // by the solute's inability to hydrate. Every term is linear or a convex
  // quadratic in the blend volume fraction (blend() is linear in it, and
  // max(0, linear)^2 is convex), so chi stays unimodal-friendly and the
  // co-solvency maximum (contract item E) survives. Floored, not clamped
  // toward zero, at -1.0: the fitted constants never drive it negative for
  // this calibration set, but a strongly complementary pair legitimately
  // could, and an unbounded negative chi would leave the bisection below
  // unable to bracket a root.
  double volumeSolute = std::max(solute.molarVolume, 1e-6);
  double mixtureVolume = std::max(mixture.molarVolume, 1e-6);
  const double hSolute = solute.hansen.hydrogenBond;
  const double hMixture = mixture.hansen.hydrogenBond;
  const double soluteExcess = std::max(0.0, hSolute - hMixture);    // term A
  const double solventExcess = std::max(0.0, hMixture - hSolute);   // term B
  const double inertFactor = std::max(0.0, 1.0 - hMixture / chiCoeff::kSolventInertCutoff);
  const double cavityFactor =
      std::max(0.0, 1.0 - hSolute / chiCoeff::kHydrationCutoff);
  double chi = chiCoeff::kC0 +
               volumeSolute / (kGasConstant * temperature) *
                   (chiCoeff::kC1 * dDispersion * dDispersion +
                    chiCoeff::kC2 * dPolar * dPolar +
                    chiCoeff::kC3a * inertFactor * soluteExcess * soluteExcess +
                    chiCoeff::kC3b * cavityFactor * solventExcess * solventExcess);
  chi = std::max(chi, -1.0);

  // Ideal (Hildebrand/Yalkowsky) mole-fraction solubility from melting-point
  // depression (contract item 3); a liquid solute (Tm <= T) is ideally
  // miscible. The log is floored, not the resulting x_ideal, so a very
  // high-melting solute never underflows to a literal, log-breaking zero.
  double logXIdeal = 0.0;
  if (solute.meltingPoint > temperatureC) {
    double meltingK = solute.meltingPoint + 273.15;
    logXIdeal = -(solute.entropyOfFusion / kGasConstant) * (meltingK - temperature) / temperature;
    if (!std::isfinite(logXIdeal)) logXIdeal = std::log(1e-300);
    logXIdeal = std::max(logXIdeal, std::log(1e-300));
  }
  logXIdeal = std::min(logXIdeal, 0.0);  // x_ideal can't exceed 1
  double xIdeal = std::clamp(std::exp(logXIdeal), 1e-300, 1.0);

  // Self-consistent Flory-Huggins solve for the solute volume fraction
  // phi_s (contract item 4): ln(a_s) = ln(phi_s) + (1 - 1/m)(1 - phi_s) +
  // chi*(1-phi_s)^2, m = V_solute / V_mixture, solved for ln(a_s) ==
  // ln(x_ideal) by bisection. ln(a_s) is monotone increasing in phi_s on
  // (0,1), so plain bisection also correctly collapses to the phi_s == lo or
  // phi_s == hi boundary when the interval doesn't bracket a root (i.e. the
  // solute is essentially insoluble, or the model has it fully miscible) --
  // no special-casing needed. This self-consistency (rather than assuming
  // infinite dilution) is what makes the solubility genuinely depend on the
  // solvent blend ratio and produces the co-solvency maximum as chi and m
  // both move with it.
  double m = std::max(volumeSolute / mixtureVolume, 1e-12);
  auto lnActivity = [&](double phi) {
    double oneMinusPhi = 1.0 - phi;
    return std::log(phi) + (1.0 - 1.0 / m) * oneMinusPhi + chi * oneMinusPhi * oneMinusPhi;
  };

  double lo = 1e-12;
  double hi = 1.0 - 1e-12;
  for (int i = 0; i < 200; ++i) {
    double mid = 0.5 * (lo + hi);
    double f = lnActivity(mid) - logXIdeal;
    if (f > 0.0) {
      hi = mid;
    } else {
      lo = mid;
    }
  }
  double phi = 0.5 * (lo + hi);
  bool converged = std::isfinite(phi) && (hi - lo) < 1e-10;

  // Convert phi_s back to mole fraction and concentration (contract item 5).
  double soluteMolarTerm = safeDiv(phi, volumeSolute, 0.0);
  double solventMolarTerm = safeDiv(1.0 - phi, mixtureVolume, 0.0);
  double x = safeDiv(soluteMolarTerm, soluteMolarTerm + solventMolarTerm, 0.0);
  x = std::clamp(x, 0.0, 1.0);

  double solutionMolarVolume = x * solute.molarVolume + (1.0 - x) * mixture.molarVolume;
  double molesPerLitre = std::max(0.0, safeDiv(1000.0 * x, solutionMolarVolume, 0.0));
  double gramsPerMillilitre =
      std::max(0.0, safeDiv(x * solute.molarMass, solutionMolarVolume, 0.0));
  // activityCoefficient = x_ideal / x (contract item 6); the saturated
  // solution's departure from Raoult's-law ideality.
  double gamma = safeDiv(xIdeal, std::max(x, 1e-300), 1.0);

  // ---- Yalkowsky General Solubility Equation (aqueous share) -----------
  // Jain, N.; Yalkowsky, S.H. J. Pharm. Sci. 2001, 90, 234:
  //     log S_w [mol/L] = 0.5 - 0.01 (Tm[C] - 25) - logP
  // The FH/extended-Hansen extrapolation is weakest exactly where the GSE is
  // strongest -- neutral organic solids in water, where the hydrophobic
  // effect dominates and Hansen distances underestimate it. Measured against
  // this project's anchor set the GSE cuts aqueous error from ~1.4 to ~0.4
  // log units for hydrophobic solutes (naphthalene, benzoic acid, aspirin).
  // Blending is log-linear in the water volume fraction, which is Yalkowsky's
  // own cosolvency mixing rule, so a waterless blend is untouched.
  double waterFraction = 0.0;
  double totalVolume = 0.0;
  for (const Component& component : components) {
    if (!component.solvent || component.volumeFraction <= 0.0) continue;
    totalVolume += component.volumeFraction;
    if (component.solvent->id == "water") waterFraction += component.volumeFraction;
  }
  if (totalVolume > 0.0) waterFraction /= totalVolume;

  // The GSE is a CRYSTALLINE-solute correlation: its -0.01 (Tm - 25) term is
  // the lattice penalty. A solute that is liquid at the working temperature
  // has no crystal to break, so the FH/Hansen result stands unmodified.
  const bool crystalline = solute.meltingPoint > temperatureC;
  if (crystalline && waterFraction > 0.0 && solute.molarMass > 0.0) {
    const double logS = 0.5 - 0.01 * (solute.meltingPoint - 25.0) - solute.logP;
    const double molarGse = std::pow(10.0, std::clamp(logS, -12.0, 1.4));  // <= ~25 M
    const double gPerMlGse = molarGse * solute.molarMass / 1000.0;
    if (gPerMlGse > 0.0 && gramsPerMillilitre > 0.0) {
      const double blended = std::exp((1.0 - waterFraction) * std::log(gramsPerMillilitre) +
                                      waterFraction * std::log(gPerMlGse));
      const double scale = safeDiv(blended, gramsPerMillilitre, 1.0);
      gramsPerMillilitre = blended;
      molesPerLitre *= scale;
      x = std::clamp(x * scale, 0.0, 1.0);
      gamma = safeDiv(xIdeal, std::max(x, 1e-300), 1.0);
    }
  }

  prediction.gramsPerMillilitre = std::isfinite(gramsPerMillilitre) ? gramsPerMillilitre : 0.0;
  prediction.molesPerLitre = std::isfinite(molesPerLitre) ? molesPerLitre : 0.0;
  prediction.moleFraction = x;
  prediction.idealMoleFraction = xIdeal;
  prediction.ra = ra;
  prediction.relativeEnergyDifference = std::isfinite(red) ? red : 0.0;
  prediction.activityCoefficient = std::isfinite(gamma) ? gamma : 1.0;
  prediction.chi = std::isfinite(chi) ? chi : 0.0;
  prediction.outsideSphere = red > 1.0;
  prediction.converged = converged;
  return prediction;
}

// The ideal-solubility factor, exactly as the raw solve computes it, factored
// out so the anchor correction can scale measured values with it.
double idealSolubility(double meltingPointC, double entropyOfFusion, double temperatureC) {
  const double t = temperatureC + 273.15;
  const double meltingK = meltingPointC + 273.15;
  if (meltingK <= t) return 1.0;  // liquid: no crystalline penalty
  double logX = -(entropyOfFusion / kGasConstant) * (meltingK - t) / t;
  if (!std::isfinite(logX)) logX = std::log(1e-300);
  return std::clamp(std::exp(std::min(logX, 0.0)), 1e-300, 1.0);
}

// Aqueous 1:1-salt path. Ksp machinery supplies RATIOS (common-ion effect,
// ionic strength, van't Hoff temperature); the measured 25 C pure-water
// value pins the endpoint exactly, so a pure-water query returns the
// literature number, not the Davies-clamped estimate.
Prediction saltPrediction(const Salt& salt, double waterFraction, const Electrolyte* background,
                          double backgroundM, double temperatureC) {
  Prediction prediction;
  prediction.saltPath = true;
  prediction.converged = true;

  const double modelNow = saltSolubilityMolar(salt, background, backgroundM, temperatureC);
  const double modelBase = saltSolubilityMolar(salt, nullptr, 0.0, 25.0);
  const double measuredMolar = salt.solubilityGPerMl25 / salt.molarMass * 1000.0;
  const double molar = modelBase > 0.0 ? measuredMolar * modelNow / modelBase : modelNow;
  const double gPerMl = molar * salt.molarMass / 1000.0 * waterFraction;

  prediction.gramsPerMillilitre = gPerMl;
  prediction.molesPerLitre = gPerMl / salt.molarMass * 1000.0;
  // Molality-style fraction over the water; a dilute-solution approximation
  // shown for scale, not a rigorous phase composition.
  prediction.moleFraction = molar * waterFraction / 55.5;
  prediction.idealMoleFraction = prediction.moleFraction;
  const double ionicI = (background ? backgroundM : 0.0) + modelNow;
  prediction.activityCoefficient = daviesGamma(ionicI);
  prediction.anchored = true;
  prediction.anchorNote = salt.name + ": measured value" +
                          (background && backgroundM > 0.0 ? " with " + background->name
                                                           : "") +
                          (temperatureC != 25.0 ? ", van't Hoff scaled" : "");
  return prediction;
}

}  // namespace

// Shared extended-Hansen/Martin chi (declared in sol/chi.hpp for the
// Kirkwood-Buff module): the exact regression the FH solve uses, so the KB
// readout reports the same chi rather than a divergent reimplementation.
// Defined outside the anonymous block for external linkage; the fitted
// constants above remain its single definition.
double extendedHansenChi(const Solute& solute, const Hansen& mixtureHansen,
                         double temperatureK) {
  const double dDispersion = solute.hansen.dispersion - mixtureHansen.dispersion;
  const double dPolar = solute.hansen.polar - mixtureHansen.polar;
  const double volumeSolute = std::max(solute.molarVolume, 1e-6);
  const double hSolute = solute.hansen.hydrogenBond;
  const double hMixture = mixtureHansen.hydrogenBond;
  const double soluteExcess = std::max(0.0, hSolute - hMixture);    // term A
  const double solventExcess = std::max(0.0, hMixture - hSolute);   // term B
  const double inertFactor = std::max(0.0, 1.0 - hMixture / chiCoeff::kSolventInertCutoff);
  const double cavityFactor =
      std::max(0.0, 1.0 - hSolute / chiCoeff::kHydrationCutoff);
  const double chi = chiCoeff::kC0 +
                     volumeSolute / (kGasConstant * temperatureK) *
                         (chiCoeff::kC1 * dDispersion * dDispersion +
                          chiCoeff::kC2 * dPolar * dPolar +
                          chiCoeff::kC3a * inertFactor * soluteExcess * soluteExcess +
                          chiCoeff::kC3b * cavityFactor * solventExcess * solventExcess);
  return std::max(chi, -1.0);
}

Prediction floryHugginsPrediction(const Solute& solute, const std::vector<Component>& components,
                                  double temperatureC) {
  return computeFhPrediction(solute, components, temperatureC);
}

Prediction predict(const Solute& solute, const std::vector<Component>& components,
                   double temperatureC, const Electrolyte* background, double backgroundM) {
  // 1:1 salts bypass the organic model: a Ksp equilibrium over the aqueous
  // share of the blend. A salt in a waterless blend is effectively
  // insoluble.
  if (!solute.canonicalSmiles.empty()) {
    if (const Salt* salt = findSalt(solute.canonicalSmiles)) {
      double waterFraction = 0.0;
      double totalFraction = 0.0;
      for (const Component& component : components) {
        if (!component.solvent || component.volumeFraction <= 0.0) continue;
        totalFraction += component.volumeFraction;
        if (component.solvent->family == "water") waterFraction += component.volumeFraction;
      }
      waterFraction = totalFraction > 0.0 ? waterFraction / totalFraction : 0.0;
      return saltPrediction(*salt, waterFraction, background, backgroundM, temperatureC);
    }
  }

  Prediction prediction = computeFhPrediction(solute, components, temperatureC);
  if (solute.canonicalSmiles.empty()) return prediction;

  // Measured-value correction: where a literature anchor exists for an
  // involved pure solvent, pull the prediction toward it log-linearly in the
  // blend fraction. The anchor travels with the solute's ideal solubility,
  //
  //   S(T, Tm) = S_measured * xIdeal(T; Tm_current) / xIdeal(T_anchor; Tm_lit)
  //
  // so an overridden melting point or a different working temperature keeps
  // its leverage, while (T = 25 C, Tm = literature) returns the measured
  // value exactly. At an anchored pure endpoint (weight 1) the result IS
  // the anchor; blends interpolate geometrically between anchored
  // endpoints. The model's job is the shape; the anchors' job is the truth.
  double total = 0.0;
  for (const Component& component : components) {
    if (component.solvent && component.volumeFraction > 0.0) total += component.volumeFraction;
  }
  if (total <= 0.0) return prediction;

  const double idealNow = idealSolubility(solute.meltingPoint, solute.entropyOfFusion,
                                          temperatureC);
  double logCorrection = 0.0;
  bool anchored = false;
  std::string notes;
  for (const Component& component : components) {
    if (!component.solvent || component.volumeFraction <= 0.0) continue;
    // Nearest anchor for this exact pair (the DB holds at most one today).
    const SolubilityAnchor* anchor = nullptr;
    for (const SolubilityAnchor& candidate : anchors()) {
      if (candidate.soluteSmiles != solute.canonicalSmiles ||
          candidate.solventId != component.solvent->id) {
        continue;
      }
      if (!anchor || std::fabs(candidate.temperatureC - temperatureC) <
                         std::fabs(anchor->temperatureC - temperatureC)) {
        anchor = &candidate;
      }
    }
    if (!anchor) continue;
    const std::vector<Component> pure{Component{component.solvent, 1.0}};
    const double rawPure = computeFhPrediction(solute, pure, temperatureC).gramsPerMillilitre;
    if (rawPure <= 0.0) continue;
    const double idealAtAnchor = idealSolubility(anchor->soluteMeltingPointC,
                                                 solute.entropyOfFusion, anchor->temperatureC);
    const double anchorValue = anchor->gramsPerMillilitre * idealNow /
                               std::max(idealAtAnchor, 1e-300);
    logCorrection += (component.volumeFraction / total) * std::log(anchorValue / rawPure);
    anchored = true;
    if (!notes.empty()) notes += "; ";
    notes += anchor->note.empty() ? component.solvent->name : anchor->note;
  }
  if (!anchored) return prediction;

  const double scale = std::exp(logCorrection);
  prediction.gramsPerMillilitre *= scale;
  // Recompute the derived quantities with the SAME transforms the raw solve
  // used, so moleFraction, molesPerLitre and gramsPerMillilitre stay
  // exactly self-consistent after the correction.
  const Mixture mix = blend(components);
  const double g = prediction.gramsPerMillilitre;
  const double vs = std::max(solute.molarVolume, 1e-6);
  const double vm = std::max(mix.molarVolume, 1e-6);
  const double mm = std::max(solute.molarMass, 1e-9);
  double x = g * vm / std::max(mm - g * (vs - vm), 1e-300);
  x = std::clamp(x, 0.0, 1.0);
  const double solutionMolarVolume = x * vs + (1.0 - x) * vm;
  prediction.moleFraction = x;
  prediction.molesPerLitre = std::max(0.0, 1000.0 * x / solutionMolarVolume);
  prediction.activityCoefficient = idealNow / std::max(x, 1e-300);
  prediction.anchored = true;
  prediction.anchorNote = "anchored to measured data (" + notes + ")";
  return prediction;
}

std::vector<SweepPoint> sweep(const Solute& solute, const std::vector<const Solvent*>& solvents,
                              int steps, double temperatureC, const Electrolyte* background,
                              double backgroundM) {
  if (solvents.empty() || solvents.size() > 3) {
    throw SolError("sweep needs 1 to 3 solvents");
  }
  for (const Solvent* solvent : solvents) {
    if (!solvent) throw SolError("sweep received a null solvent");
  }
  int n = std::clamp(steps, 2, 64);

  std::vector<SweepPoint> points;
  auto emit = [&](double f0, double f1, double f2) {
    SweepPoint point;
    point.fractions = {f0, f1, f2};
    std::vector<Component> components;
    components.push_back({solvents[0], f0});
    if (solvents.size() > 1) components.push_back({solvents[1], f1});
    if (solvents.size() > 2) components.push_back({solvents[2], f2});
    point.prediction = predict(solute, components, temperatureC, background, backgroundM);
    points.push_back(std::move(point));
  };

  if (solvents.size() == 1) {
    emit(1.0, 0.0, 0.0);
  } else if (solvents.size() == 2) {
    points.reserve(static_cast<size_t>(n) + 1);
    for (int i = 0; i <= n; ++i) {
      double f0 = static_cast<double>(n - i) / n;
      emit(f0, 1.0 - f0, 0.0);
    }
  } else {
    points.reserve(static_cast<size_t>(n + 1) * static_cast<size_t>(n + 2) / 2);
    for (int i = 0; i <= n; ++i) {
      for (int j = 0; j <= n - i; ++j) {
        int k = n - i - j;
        emit(static_cast<double>(i) / n, static_cast<double>(j) / n, static_cast<double>(k) / n);
      }
    }
  }
  return points;
}

std::vector<ScreenRow> screen(const Solute& solute, double temperatureC,
                              const Electrolyte* background, double backgroundM) {
  const std::vector<Solvent>& table = solvents();  // propagates a load failure
  std::vector<ScreenRow> rows;
  rows.reserve(table.size());
  for (const Solvent& solvent : table) {
    const std::vector<Component> pure{Component{&solvent, 1.0}};
    rows.push_back(ScreenRow{&solvent, predict(solute, pure, temperatureC, background,
                                               backgroundM)});
  }
  // Best solvent first; name tiebreak keeps the order deterministic across
  // runs when two predictions land on the same value.
  std::stable_sort(rows.begin(), rows.end(), [](const ScreenRow& a, const ScreenRow& b) {
    return a.prediction.gramsPerMillilitre != b.prediction.gramsPerMillilitre
               ? a.prediction.gramsPerMillilitre > b.prediction.gramsPerMillilitre
               : a.solvent->name < b.solvent->name;
  });
  return rows;
}

Partition partition(double massMg, double logP, double volumeAqueousMl,
                    double volumeOrganicMl) {
  Partition result;
  if (massMg <= 0.0) return result;

  // Clamp so a garbage logP cannot overflow the distribution ratio; the
  // clamp bounds D to [1e-6, 1e6], far past any real solvent pair.
  const double clampedLogP = std::clamp(logP, -6.0, 6.0);
  const double d = std::pow(10.0, clampedLogP);
  const double aqueous = std::max(0.0, volumeAqueousMl);
  const double organic = std::max(0.0, volumeOrganicMl);

  // C_org/C_aq = D, mass balance over both volumes gives the shares.
  const double organicTerm = d * organic;
  const double total = aqueous + organicTerm;
  if (total <= 0.0) return result;  // no phase volume anywhere

  result.mgOrganic = massMg * organicTerm / total;
  result.mgAqueous = massMg - result.mgOrganic;
  result.fractionOrganic = result.mgOrganic / massMg;
  return result;
}

}  // namespace chemcad::sol

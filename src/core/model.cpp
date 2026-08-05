#include "core/model.hpp"

#include <algorithm>
#include <utility>

namespace chemcad::core {

AtomId Molecule::addAtom(Atom a) {
  a.id = nextAtomId_++;
  atoms_.push_back(a);
  return a.id;
}

BondId Molecule::addBond(AtomId a, AtomId b, BondOrder o) {
  if (a == kInvalidAtom || b == kInvalidAtom || a == b) return kInvalidBond;
  if (!atom(a) || !atom(b)) return kInvalidBond;
  if (BondId existing = bondBetween(a, b); existing != kInvalidBond) {
    bond(existing)->order = o;  // re-bonding the same pair retypes it
    return existing;
  }
  Bond bnd;
  bnd.id = nextBondId_++;
  bnd.a = a;
  bnd.b = b;
  bnd.order = o;
  bonds_.push_back(bnd);
  return bnd.id;
}

void Molecule::removeAtom(AtomId id) {
  std::erase_if(bonds_, [&](const Bond& b) { return b.a == id || b.b == id; });
  std::erase_if(atoms_, [&](const Atom& a) { return a.id == id; });
}

void Molecule::removeBond(BondId id) {
  std::erase_if(bonds_, [&](const Bond& b) { return b.id == id; });
}

Atom* Molecule::atom(AtomId id) {
  return const_cast<Atom*>(std::as_const(*this).atom(id));
}

const Atom* Molecule::atom(AtomId id) const {
  auto it = std::find_if(atoms_.begin(), atoms_.end(),
                         [&](const Atom& a) { return a.id == id; });
  return it == atoms_.end() ? nullptr : &*it;
}

Bond* Molecule::bond(BondId id) {
  return const_cast<Bond*>(std::as_const(*this).bond(id));
}

const Bond* Molecule::bond(BondId id) const {
  auto it = std::find_if(bonds_.begin(), bonds_.end(),
                         [&](const Bond& b) { return b.id == id; });
  return it == bonds_.end() ? nullptr : &*it;
}

std::vector<AtomId> Molecule::neighbors(AtomId id) const {
  std::vector<AtomId> out;
  for (const Bond& b : bonds_) {
    if (b.a == id) out.push_back(b.b);
    else if (b.b == id) out.push_back(b.a);
  }
  return out;
}

std::vector<BondId> Molecule::incidentBonds(AtomId id) const {
  std::vector<BondId> out;
  for (const Bond& b : bonds_)
    if (b.a == id || b.b == id) out.push_back(b.id);
  return out;
}

BondId Molecule::bondBetween(AtomId a, AtomId b) const {
  for (const Bond& bnd : bonds_)
    if ((bnd.a == a && bnd.b == b) || (bnd.a == b && bnd.b == a)) return bnd.id;
  return kInvalidBond;
}

int Molecule::degree(AtomId id) const {
  int n = 0;
  for (const Bond& b : bonds_)
    if (b.a == id || b.b == id) ++n;
  return n;
}

void Molecule::clear() {
  atoms_.clear();
  bonds_.clear();
  nextAtomId_ = 1;
  nextBondId_ = 1;
}

bool Document::empty() const {
  return std::all_of(molecules.begin(), molecules.end(),
                     [](const Molecule& m) { return m.empty(); });
}

int Document::largestMoleculeIndex() const {
  int best = -1;
  size_t bestCount = 0;
  for (size_t i = 0; i < molecules.size(); ++i) {
    if (molecules[i].atomCount() > bestCount) {
      bestCount = molecules[i].atomCount();
      best = static_cast<int>(i);
    }
  }
  return best;
}

void UndoStack::push(const Document& doc) {
  past_.push_back(doc);
  if (past_.size() > kCapacity) past_.pop_front();
  future_.clear();
}

bool UndoStack::undo(Document& doc) {
  if (past_.empty()) return false;
  future_.push_back(doc);
  if (future_.size() > kCapacity) future_.pop_front();
  doc = past_.back();
  past_.pop_back();
  return true;
}

bool UndoStack::redo(Document& doc) {
  if (future_.empty()) return false;
  past_.push_back(doc);
  if (past_.size() > kCapacity) past_.pop_front();
  doc = future_.back();
  future_.pop_back();
  return true;
}

void UndoStack::clear() {
  past_.clear();
  future_.clear();
}

}  // namespace chemcad::core

#pragma once
// Editable molecule graph used by the canvas. Deliberately independent of
// RDKit: this is the fast, undo-friendly structure the UI mutates every frame.
// Conversion to/from RDKit happens at the boundary in chemcad::chem.

#include <cstdint>
#include <deque>
#include <span>
#include <vector>

namespace chemcad::core {

struct Vec2 {
  float x = 0, y = 0;
};

using AtomId = uint32_t;
using BondId = uint32_t;
inline constexpr AtomId kInvalidAtom = 0;  // 0 is never a valid id
inline constexpr BondId kInvalidBond = 0;

enum class BondOrder : uint8_t { Single = 1, Double = 2, Triple = 3, Aromatic = 4 };
enum class BondStereo : uint8_t { None = 0, Wedge = 1, Hash = 2 };

struct Atom {
  AtomId id{};
  uint8_t atomicNumber = 6;
  int8_t charge = 0;
  uint16_t isotope = 0;  // 0 = natural abundance
  Vec2 pos;
  int8_t explicitH = -1;  // -1 = let the toolkit decide
};

struct Bond {
  BondId id{};
  AtomId a{}, b{};
  BondOrder order = BondOrder::Single;
  BondStereo stereo = BondStereo::None;
};

class Molecule {
 public:
  AtomId addAtom(Atom a);                        // assigns and returns a fresh id
  BondId addBond(AtomId a, AtomId b, BondOrder o);
  void removeAtom(AtomId);                       // also removes incident bonds
  void removeBond(BondId);

  Atom* atom(AtomId);
  const Atom* atom(AtomId) const;
  Bond* bond(BondId);
  const Bond* bond(BondId) const;

  std::span<const Atom> atoms() const { return atoms_; }
  std::span<const Bond> bonds() const { return bonds_; }
  std::span<Atom> mutableAtoms() { return atoms_; }

  std::vector<AtomId> neighbors(AtomId) const;
  std::vector<BondId> incidentBonds(AtomId) const;
  BondId bondBetween(AtomId, AtomId) const;      // kInvalidBond when unbonded
  int degree(AtomId) const;

  bool empty() const { return atoms_.empty(); }
  size_t atomCount() const { return atoms_.size(); }
  size_t bondCount() const { return bonds_.size(); }
  void clear();

 private:
  std::vector<Atom> atoms_;
  std::vector<Bond> bonds_;
  AtomId nextAtomId_ = 1;
  BondId nextBondId_ = 1;
};

struct Document {
  std::vector<Molecule> molecules;

  bool empty() const;
  void clear() { molecules.clear(); }
  // Index of the molecule with the most atoms, or -1 when the document is empty.
  int largestMoleculeIndex() const;
};

// Snapshot undo. Documents are tiny (hundreds of atoms at most), so copying the
// whole document per gesture is far simpler than a command log and never drifts.
class UndoStack {
 public:
  static constexpr size_t kCapacity = 256;

  void push(const Document&);  // call BEFORE mutating; clears the redo branch
  bool undo(Document&);
  bool redo(Document&);
  bool canUndo() const { return !past_.empty(); }
  bool canRedo() const { return !future_.empty(); }
  void clear();

 private:
  std::deque<Document> past_;
  std::deque<Document> future_;
};

}  // namespace chemcad::core

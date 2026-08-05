#include "ui/app_state.hpp"

#include <algorithm>

namespace chemcad::ui {

bool Selection::contains(const AtomRef& r) const {
  return std::find(atoms.begin(), atoms.end(), r) != atoms.end();
}

bool Selection::contains(const BondRef& r) const {
  return std::find(bonds.begin(), bonds.end(), r) != bonds.end();
}

core::Molecule* AppState::molecule(int index) {
  if (index < 0 || index >= static_cast<int>(doc.molecules.size())) return nullptr;
  return &doc.molecules[static_cast<size_t>(index)];
}

core::Atom* AppState::atomAt(const AtomRef& r) {
  core::Molecule* m = molecule(r.mol);
  return m ? m->atom(r.id) : nullptr;
}

core::Bond* AppState::bondAt(const BondRef& r) {
  core::Molecule* m = molecule(r.mol);
  return m ? m->bond(r.id) : nullptr;
}

}  // namespace chemcad::ui

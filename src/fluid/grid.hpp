#pragma once
// Uniform-cell neighbour search for the SPH solver.
//
// Cell edge equals the kernel support radius H, so every neighbour within
// support lies in the 27-cell block around a particle's own cell. The grid is
// rebuilt every substep (all particles move every substep, which is why a
// KD-tree is the wrong structure here).
//
// Determinism is a hard requirement: the solver has reproducibility tests, so
// cell contents are produced by a radix sort on a Morton key with a fixed
// coordinate bias, never by iterating a hash map. Particles are reordered into
// key order so that neighbour loops are cache-local, and the immutable
// `Particles::id` is what tests and diagnostics identify particles by.

#include <cstdint>
#include <vector>

#include "fluid/types.hpp"

namespace chemcad::fluid {

class NeighbourGrid {
 public:
  // Rebuilds the grid for the current particle positions. `support` is the
  // kernel support radius in metres and becomes the cell edge. Reorders the
  // particle arrays into Morton order as a side effect.
  void build(Particles& particles, double support);

  // Number of particles the last build saw.
  std::size_t size() const { return keys_.size(); }

  // Calls `fn(j)` for every particle j (including i itself) whose cell is in
  // the 27-cell block around particle i. Callers must still test r < H.
  // Iteration order is ascending j so pairwise accumulation is reproducible.
  template <typename Fn>
  void forEachCandidate(std::size_t i, Fn&& fn) const;

  // Same block query for an arbitrary point, used by the renderer's surface
  // sampling and by diagnostics.
  template <typename Fn>
  void forEachCandidateAt(float x, float y, float z, Fn&& fn) const;

  double support() const { return support_; }

 private:
  struct CellRange {
    uint64_t key = 0;
    uint32_t begin = 0;
    uint32_t end = 0;
  };

  // Cell coordinate of a position, with the fixed +2^20 bias that keeps keys
  // positive for the whole reachable vessel volume.
  void cellOf(float x, float y, float z, int32_t out[3]) const;
  uint64_t keyOf(const int32_t cell[3]) const;
  const CellRange* findCell(uint64_t key) const;

  double support_ = 8.0e-3;
  double inverseSupport_ = 125.0;
  std::vector<uint64_t> keys_;        // per particle, in sorted order
  std::vector<CellRange> cells_;      // ascending key, one entry per occupied cell
  std::vector<uint32_t> scratchKeys_; // radix-sort scratch, kept to avoid churn
  std::vector<uint32_t> order_;
};

}  // namespace chemcad::fluid

namespace chemcad::fluid {

template <typename Fn>
void NeighbourGrid::forEachCandidate(std::size_t i, Fn&& fn) const {
  if (i >= keys_.size()) return;

  // Reverse the 21-bit Morton interleave. Keeping the sorted key per particle
  // avoids retaining a second copy of the positions solely for grid queries.
  const auto compactByThree = [](uint64_t value) {
    value &= 0x1249249249249249ULL;
    value = (value ^ (value >> 2)) & 0x10c30c30c30c30c3ULL;
    value = (value ^ (value >> 4)) & 0x100f00f00f00f00fULL;
    value = (value ^ (value >> 8)) & 0x001f0000ff0000ffULL;
    value = (value ^ (value >> 16)) & 0x001f00000000ffffULL;
    value = (value ^ (value >> 32)) & 0x00000000001fffffULL;
    return static_cast<int32_t>(value) - (1 << 20);
  };

  const uint64_t key = keys_[i];
  const int32_t base[3] = {
      compactByThree(key), compactByThree(key >> 1), compactByThree(key >> 2)};
  const CellRange* ranges[27]{};
  std::size_t rangeCount = 0;
  for (int32_t dz = -1; dz <= 1; ++dz) {
    for (int32_t dy = -1; dy <= 1; ++dy) {
      for (int32_t dx = -1; dx <= 1; ++dx) {
        const int32_t cell[3] = {base[0] + dx, base[1] + dy, base[2] + dz};
        const CellRange* range = findCell(keyOf(cell));
        if (range == nullptr) continue;
        bool duplicate = false;
        for (std::size_t n = 0; n < rangeCount; ++n) {
          if (ranges[n] == range) duplicate = true;
        }
        if (!duplicate) ranges[rangeCount++] = range;
      }
    }
  }

  // Cell lookup above is fixed (dz,dy,dx). Reordering its at-most-27 ranges by
  // begin index then makes the actual gather globally ascending in particle j,
  // satisfying the reduction-order contract without allocating a candidate list.
  for (std::size_t n = 1; n < rangeCount; ++n) {
    const CellRange* value = ranges[n];
    std::size_t insertion = n;
    while (insertion > 0 && ranges[insertion - 1]->begin > value->begin) {
      ranges[insertion] = ranges[insertion - 1];
      --insertion;
    }
    ranges[insertion] = value;
  }
  for (std::size_t n = 0; n < rangeCount; ++n) {
    for (uint32_t j = ranges[n]->begin; j < ranges[n]->end; ++j) fn(j);
  }
}

template <typename Fn>
void NeighbourGrid::forEachCandidateAt(float x, float y, float z, Fn&& fn) const {
  if (!(support_ > 0.0)) return;
  int32_t base[3];
  cellOf(x, y, z, base);
  const CellRange* ranges[27]{};
  std::size_t rangeCount = 0;
  for (int32_t dz = -1; dz <= 1; ++dz) {
    for (int32_t dy = -1; dy <= 1; ++dy) {
      for (int32_t dx = -1; dx <= 1; ++dx) {
        const int32_t cell[3] = {base[0] + dx, base[1] + dy, base[2] + dz};
        const CellRange* range = findCell(keyOf(cell));
        if (range == nullptr) continue;
        bool duplicate = false;
        for (std::size_t n = 0; n < rangeCount; ++n) {
          if (ranges[n] == range) duplicate = true;
        }
        if (!duplicate) ranges[rangeCount++] = range;
      }
    }
  }
  for (std::size_t n = 1; n < rangeCount; ++n) {
    const CellRange* value = ranges[n];
    std::size_t insertion = n;
    while (insertion > 0 && ranges[insertion - 1]->begin > value->begin) {
      ranges[insertion] = ranges[insertion - 1];
      --insertion;
    }
    ranges[insertion] = value;
  }
  for (std::size_t n = 0; n < rangeCount; ++n) {
    for (uint32_t j = ranges[n]->begin; j < ranges[n]->end; ++j) fn(j);
  }
}

}  // namespace chemcad::fluid

#include "fluid/grid.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>

namespace chemcad::fluid {
namespace {

constexpr int32_t kCellBias = 1 << 20;
constexpr uint64_t kCoordinateMask = (uint64_t{1} << 21) - 1;

uint64_t spreadByThree(uint32_t value) {
  uint64_t x = static_cast<uint64_t>(value) & kCoordinateMask;
  x = (x | (x << 32)) & 0x001f00000000ffffULL;
  x = (x | (x << 16)) & 0x001f0000ff0000ffULL;
  x = (x | (x << 8)) & 0x100f00f00f00f00fULL;
  x = (x | (x << 4)) & 0x10c30c30c30c30c3ULL;
  x = (x | (x << 2)) & 0x1249249249249249ULL;
  return x;
}

template <typename T>
void swapParticleEntries(std::vector<T>& values, std::size_t a, std::size_t b,
                         std::size_t particleCount) {
  if (values.size() == particleCount) std::swap(values[a], values[b]);
}

}  // namespace

void NeighbourGrid::cellOf(float x, float y, float z, int32_t out[3]) const {
  out[0] = static_cast<int32_t>(std::floor(static_cast<double>(x) * inverseSupport_));
  out[1] = static_cast<int32_t>(std::floor(static_cast<double>(y) * inverseSupport_));
  out[2] = static_cast<int32_t>(std::floor(static_cast<double>(z) * inverseSupport_));
}

uint64_t NeighbourGrid::keyOf(const int32_t cell[3]) const {
  const auto biased = [](int32_t coordinate) {
    const int64_t value = static_cast<int64_t>(coordinate) + kCellBias;
    return static_cast<uint32_t>(std::clamp<int64_t>(value, 0, kCoordinateMask));
  };
  const uint64_t x = spreadByThree(biased(cell[0]));
  const uint64_t y = spreadByThree(biased(cell[1]));
  const uint64_t z = spreadByThree(biased(cell[2]));
  return x | (y << 1) | (z << 2);
}

const NeighbourGrid::CellRange* NeighbourGrid::findCell(uint64_t key) const {
  const auto it = std::lower_bound(
      cells_.begin(), cells_.end(), key,
      [](const CellRange& cell, uint64_t sought) { return cell.key < sought; });
  return it != cells_.end() && it->key == key ? &*it : nullptr;
}

void NeighbourGrid::build(Particles& particles, double support) {
  support_ = support;
  if (!(support > 0.0) || !std::isfinite(support)) {
    inverseSupport_ = 0.0;
    keys_.clear();
    cells_.clear();
    order_.clear();
    scratchKeys_.clear();
    return;
  }
  inverseSupport_ = 1.0 / support;

  const std::size_t count = particles.size();
  keys_.resize(count);
  order_.resize(count);
  scratchKeys_.resize(count);
  cells_.clear();
  cells_.reserve(count);

  for (std::size_t i = 0; i < count; ++i) {
    int32_t cell[3];
    cellOf(particles.px[i], particles.py[i], particles.pz[i], cell);
    keys_[i] = keyOf(cell);
    order_[i] = static_cast<uint32_t>(i);
  }

  // Stable least-significant-digit radix sort. Eight fixed byte passes avoid
  // comparison-sort implementation differences and preserve the prior index
  // order among particles sharing a cell.
  for (unsigned pass = 0; pass < 8; ++pass) {
    std::array<uint32_t, 256> histogram{};
    const unsigned shift = pass * 8;
    for (uint32_t index : order_) {
      ++histogram[static_cast<std::size_t>((keys_[index] >> shift) & 0xffULL)];
    }

    std::array<uint32_t, 256> next{};
    uint32_t offset = 0;
    for (std::size_t digit = 0; digit < histogram.size(); ++digit) {
      next[digit] = offset;
      offset += histogram[digit];
    }
    for (uint32_t index : order_) {
      const std::size_t digit = static_cast<std::size_t>((keys_[index] >> shift) & 0xffULL);
      scratchKeys_[next[digit]++] = index;
    }
    order_.swap(scratchKeys_);
  }

  // Convert new-position -> old-position into old-position -> new-position,
  // then apply one permutation cycle to every structure-of-arrays member at
  // once. This needs no temporary particle buffers and, critically, keeps the
  // immutable id travelling with the physical particle.
  for (std::size_t newIndex = 0; newIndex < count; ++newIndex) {
    scratchKeys_[order_[newIndex]] = static_cast<uint32_t>(newIndex);
  }
  for (std::size_t i = 0; i < count; ++i) {
    while (scratchKeys_[i] != i) {
      const std::size_t j = scratchKeys_[i];
      std::swap(keys_[i], keys_[j]);
      swapParticleEntries(particles.px, i, j, count);
      swapParticleEntries(particles.py, i, j, count);
      swapParticleEntries(particles.pz, i, j, count);
      swapParticleEntries(particles.vx, i, j, count);
      swapParticleEntries(particles.vy, i, j, count);
      swapParticleEntries(particles.vz, i, j, count);
      swapParticleEntries(particles.ax, i, j, count);
      swapParticleEntries(particles.ay, i, j, count);
      swapParticleEntries(particles.az, i, j, count);
      swapParticleEntries(particles.delta, i, j, count);
      swapParticleEntries(particles.pressure, i, j, count);
      swapParticleEntries(particles.colour, i, j, count);
      swapParticleEntries(particles.nx, i, j, count);
      swapParticleEntries(particles.ny, i, j, count);
      swapParticleEntries(particles.nz, i, j, count);
      swapParticleEntries(particles.phase, i, j, count);
      swapParticleEntries(particles.id, i, j, count);
      std::swap(scratchKeys_[i], scratchKeys_[j]);
    }
  }

  if (count == 0) return;
  uint32_t begin = 0;
  while (begin < count) {
    uint32_t end = begin + 1;
    while (end < count && keys_[end] == keys_[begin]) ++end;
    cells_.push_back(CellRange{keys_[begin], begin, end});
    begin = end;
  }
}

}  // namespace chemcad::fluid

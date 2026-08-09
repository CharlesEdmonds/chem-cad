#include "fluid/grid.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <limits>
#include <utility>
#include <thread>

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


unsigned configuredWorkerCount(std::size_t workItems) {
  constexpr unsigned kMaxWorkers = 8;
  unsigned requested = std::thread::hardware_concurrency();
  unsigned overrideWorkers = 0;
#ifdef _WIN32
  char* text = nullptr;
  std::size_t textLength = 0;
  if (_dupenv_s(&text, &textLength, "CHEMCAD_FLUID_WORKERS") == 0 && text != nullptr) {
    overrideWorkers = static_cast<unsigned>(std::strtoul(text, nullptr, 10));
    std::free(text);
  }
#else
  if (const char* text = std::getenv("CHEMCAD_FLUID_WORKERS")) {
    overrideWorkers = static_cast<unsigned>(std::strtoul(text, nullptr, 10));
  }
#endif
  if (overrideWorkers > 0) {
    const unsigned availableItems =
        static_cast<unsigned>(std::max<std::size_t>(1, workItems));
    return std::clamp(overrideWorkers, 1U,
                      std::min(kMaxWorkers, availableItems));
  }
  if (requested == 0) requested = 1;
  const unsigned useful =
      static_cast<unsigned>(std::max<std::size_t>(1, (workItems + 255) / 256));
  return std::clamp(requested, 1U, std::min(kMaxWorkers, useful));
}

template <typename Fn>
void parallelFor(std::size_t count, Fn&& fn) {
  const unsigned workers = configuredWorkerCount(count);
  std::array<std::thread, 7> threads;
  for (unsigned worker = 1; worker < workers; ++worker) {
    threads[worker - 1] = std::thread([&, worker] {
      fn(count * worker / workers, count * (worker + 1) / workers);
    });
  }
  fn(0, count / workers);
  for (unsigned worker = 1; worker < workers; ++worker) threads[worker - 1].join();
}

template <typename T>
void permuteByCycles(std::vector<T>& values, std::size_t particleCount,
                     const std::vector<uint32_t>& oldToNew,
                     const std::vector<uint32_t>& cycleLeaders, std::size_t leaderCount) {
  if (values.size() != particleCount) return;
  for (std::size_t cycle = 0; cycle < leaderCount; ++cycle) {
    const std::size_t leader = cycleLeaders[cycle];
    T carried = std::move(values[leader]);
    std::size_t current = leader;
    for (;;) {
      const std::size_t next = oldToNew[current];
      if (next == leader) {
        values[leader] = std::move(carried);
        break;
      }
      std::swap(carried, values[next]);
      current = next;
    }
  }
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

  parallelFor(count, [&](std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      int32_t cell[3];
      cellOf(particles.px[i], particles.py[i], particles.pz[i], cell);
      keys_[i] = keyOf(cell);
      order_[i] = static_cast<uint32_t>(i);
    }
  });

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

  // Convert new-position -> old-position into old-position -> new-position.
  // The high bit is available because a practical particle set is far below
  // 2^31 entries; use it briefly to find disjoint permutation cycles without
  // allocating another particle-sized buffer.
  for (std::size_t newIndex = 0; newIndex < count; ++newIndex) {
    scratchKeys_[order_[newIndex]] = static_cast<uint32_t>(newIndex);
  }
  constexpr uint32_t kVisited = uint32_t{1} << 31;
  constexpr uint32_t kIndexMask = ~kVisited;
  std::size_t leaderCount = 0;
  for (std::size_t i = 0; i < count; ++i) {
    if ((scratchKeys_[i] & kVisited) != 0) continue;
    order_[leaderCount++] = static_cast<uint32_t>(i);
    std::size_t current = i;
    do {
      const std::size_t next = scratchKeys_[current] & kIndexMask;
      scratchKeys_[current] |= kVisited;
      current = next;
    } while (current != i);
  }
  for (uint32_t& index : scratchKeys_) index &= kIndexMask;

  // Each SoA member is independent, so workers own complete arrays while
  // applying the same immutable cycle table. This parallelises permutation
  // without concurrent writes, atomics, or any change to stable Morton order.
  constexpr std::size_t kFieldCount = 18;
  parallelFor(kFieldCount * 256, [&](std::size_t scaledBegin, std::size_t scaledEnd) {
    const std::size_t fieldBegin = scaledBegin / 256;
    const std::size_t fieldEnd = scaledEnd / 256;
    for (std::size_t field = fieldBegin; field < fieldEnd; ++field) {
      switch (field) {
        case 0: permuteByCycles(keys_, count, scratchKeys_, order_, leaderCount); break;
        case 1: permuteByCycles(particles.px, count, scratchKeys_, order_, leaderCount); break;
        case 2: permuteByCycles(particles.py, count, scratchKeys_, order_, leaderCount); break;
        case 3: permuteByCycles(particles.pz, count, scratchKeys_, order_, leaderCount); break;
        case 4: permuteByCycles(particles.vx, count, scratchKeys_, order_, leaderCount); break;
        case 5: permuteByCycles(particles.vy, count, scratchKeys_, order_, leaderCount); break;
        case 6: permuteByCycles(particles.vz, count, scratchKeys_, order_, leaderCount); break;
        case 7: permuteByCycles(particles.ax, count, scratchKeys_, order_, leaderCount); break;
        case 8: permuteByCycles(particles.ay, count, scratchKeys_, order_, leaderCount); break;
        case 9: permuteByCycles(particles.az, count, scratchKeys_, order_, leaderCount); break;
        case 10: permuteByCycles(particles.delta, count, scratchKeys_, order_, leaderCount); break;
        case 11: permuteByCycles(particles.pressure, count, scratchKeys_, order_, leaderCount); break;
        case 12: permuteByCycles(particles.colour, count, scratchKeys_, order_, leaderCount); break;
        case 13: permuteByCycles(particles.nx, count, scratchKeys_, order_, leaderCount); break;
        case 14: permuteByCycles(particles.ny, count, scratchKeys_, order_, leaderCount); break;
        case 15: permuteByCycles(particles.nz, count, scratchKeys_, order_, leaderCount); break;
        case 16: permuteByCycles(particles.phase, count, scratchKeys_, order_, leaderCount); break;
        case 17: permuteByCycles(particles.id, count, scratchKeys_, order_, leaderCount); break;
      }
    }
  });

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

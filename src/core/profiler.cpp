#include "core/profiler.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <thread>
#include <utility>

namespace chemcad::core {
namespace {

using Clock = std::chrono::steady_clock;

constexpr std::size_t kMaxThreads = 64;
constexpr std::size_t kMaxDepth = 64;
constexpr std::size_t kMaxZonesPerThread = 512;
constexpr std::size_t kMaxCountersPerThread = 128;
constexpr std::size_t kNameCapacity = 96;
constexpr std::uint64_t kHashOffset = 1469598103934665603ull;
constexpr std::uint64_t kHashPrime = 1099511628211ull;

struct ZoneAccum {
  std::uint64_t id = 0;
  std::uint64_t parentId = 0;
  char name[kNameCapacity]{};
  std::uint64_t inclusiveNs = 0;
  std::uint64_t exclusiveNs = 0;
  std::uint32_t calls = 0;
};

struct CounterAccum {
  char name[kNameCapacity]{};
  double value = 0.0;
};

struct ThreadBuffer {
  std::atomic<std::uint32_t> writers{0};
  std::array<ZoneAccum, kMaxZonesPerThread> zones{};
  std::array<CounterAccum, kMaxCountersPerThread> counters{};
  std::uint16_t zoneCount = 0;
  std::uint16_t counterCount = 0;
};

struct StackEntry {
  Clock::time_point started{};
  std::uint64_t childNs = 0;
  std::uint64_t id = 0;
  std::uint64_t parentId = 0;
  char name[kNameCapacity]{};
};

struct ThreadSlot {
  std::atomic<std::uintptr_t> owner{0};
  std::array<ThreadBuffer, 2> buffers{};
};

std::array<ThreadSlot, kMaxThreads> gThreadSlots;
std::atomic<bool> gEnabled{false};
std::atomic<std::uint64_t> gEpoch{0};
std::atomic<bool> gFrameOpen{false};
Clock::time_point gFrameStarted{};
std::uint64_t gFrameSequence = 0;

std::array<ProfileFrame, Profiler::frameCapacity> gFrames;
std::size_t gFrameHead = 0;
std::size_t gFrameCount = 0;
std::mutex gFramesMutex;

struct LocalSlot {
  ThreadSlot* slot = nullptr;
  std::array<StackEntry, kMaxDepth> stack{};
  std::uint16_t depth = 0;
  std::uint16_t overflowDepth = 0;

  ~LocalSlot() {
    if (slot) slot->owner.store(0, std::memory_order_release);
  }
};

thread_local LocalSlot gLocalSlot;

void copyName(char (&destination)[kNameCapacity], const char* source) {
  if (!source) source = "(unnamed)";
  std::strncpy(destination, source, kNameCapacity - 1);
  destination[kNameCapacity - 1] = '\0';
}

std::uint64_t hashPath(std::uint64_t parent, const char* name) {
  std::uint64_t hash = parent == 0 ? kHashOffset : parent;
  const auto* cursor = reinterpret_cast<const unsigned char*>(name ? name : "(unnamed)");
  while (*cursor) {
    hash ^= *cursor++;
    hash *= kHashPrime;
  }
  hash ^= 0xffu;
  hash *= kHashPrime;
  return hash;
}

ThreadSlot* acquireThreadSlot() {
  if (gLocalSlot.slot) return gLocalSlot.slot;

  const std::uintptr_t token = reinterpret_cast<std::uintptr_t>(&gLocalSlot);
  for (ThreadSlot& slot : gThreadSlots) {
    std::uintptr_t expected = 0;
    if (slot.owner.compare_exchange_strong(expected, token, std::memory_order_acq_rel)) {
      gLocalSlot.depth = 0;
      gLocalSlot.overflowDepth = 0;
      gLocalSlot.slot = &slot;
      return &slot;
    }
  }
  return nullptr;
}

ThreadBuffer& acquireCurrentBuffer(std::uint64_t& epoch) {
  ThreadSlot& slot = *gLocalSlot.slot;
  for (;;) {
    epoch = gEpoch.load(std::memory_order_acquire);
    ThreadBuffer& buffer = slot.buffers[epoch & 1u];
    buffer.writers.fetch_add(1, std::memory_order_acquire);
    if (gEpoch.load(std::memory_order_acquire) == epoch) return buffer;
    buffer.writers.fetch_sub(1, std::memory_order_release);
  }
}

ZoneAccum* ensureZone(ThreadBuffer& buffer, std::uint64_t id, std::uint64_t parentId,
                      const char* name) {
  for (std::uint16_t i = 0; i < buffer.zoneCount; ++i) {
    ZoneAccum& zone = buffer.zones[i];
    if (zone.id == id && zone.parentId == parentId && std::strcmp(zone.name, name) == 0) {
      return &zone;
    }
  }
  if (buffer.zoneCount >= buffer.zones.size()) return nullptr;

  ZoneAccum& zone = buffer.zones[buffer.zoneCount++];
  zone = {};
  zone.id = id;
  zone.parentId = parentId;
  copyName(zone.name, name);
  return &zone;
}

void mergeBuffer(ThreadBuffer& buffer, ProfileFrame& frame) {
  constexpr double kNsToMs = 1.0 / 1'000'000.0;
  for (std::uint16_t i = 0; i < buffer.zoneCount; ++i) {
    const ZoneAccum& source = buffer.zones[i];
    int parent = -1;
    if (source.parentId != 0) {
      const auto parentIt = std::find_if(
          frame.zones.begin(), frame.zones.end(),
          [&](const ProfileZone& candidate) { return candidate.id == source.parentId; });
      if (parentIt != frame.zones.end()) {
        parent = static_cast<int>(parentIt - frame.zones.begin());
      }
    }

    auto found = std::find_if(frame.zones.begin(), frame.zones.end(),
                              [&](const ProfileZone& candidate) {
                                return candidate.id == source.id && candidate.parent == parent &&
                                       candidate.name == source.name;
                              });
    if (found == frame.zones.end()) {
      ProfileZone zone;
      zone.id = source.id;
      zone.name = source.name;
      zone.parent = parent;
      frame.zones.push_back(std::move(zone));
      found = frame.zones.end() - 1;
    }
    found->calls += source.calls;
    found->inclusiveMs += static_cast<double>(source.inclusiveNs) * kNsToMs;
    found->exclusiveMs += static_cast<double>(source.exclusiveNs) * kNsToMs;
  }

  for (std::uint16_t i = 0; i < buffer.counterCount; ++i) {
    const CounterAccum& source = buffer.counters[i];
    auto found = std::find_if(frame.counters.begin(), frame.counters.end(),
                              [&](const ProfileCounter& candidate) {
                                return candidate.name == source.name;
                              });
    if (found == frame.counters.end()) {
      frame.counters.push_back(ProfileCounter{source.name, source.value});
    } else {
      found->value = source.value;
    }
  }

  buffer.zoneCount = 0;
  buffer.counterCount = 0;
}

std::string csvField(const std::string& value) {
  if (value.find_first_of(",\"\r\n") == std::string::npos) return value;
  std::string escaped;
  escaped.reserve(value.size() + 2);
  escaped.push_back('"');
  for (const char character : value) {
    if (character == '"') escaped.push_back('"');
    escaped.push_back(character);
  }
  escaped.push_back('"');
  return escaped;
}

}  // namespace

Profiler& profiler() {
  static Profiler instance;
  return instance;
}

void Profiler::beginZone(const char* name) {
  if (!gEnabled.load(std::memory_order_relaxed)) return;

  ThreadSlot* slot = acquireThreadSlot();
  if (!slot) return;
  if (gLocalSlot.overflowDepth > 0 || gLocalSlot.depth >= gLocalSlot.stack.size()) {
    ++gLocalSlot.overflowDepth;
    return;
  }

  StackEntry& entry = gLocalSlot.stack[gLocalSlot.depth++];
  entry.started = Clock::now();
  entry.childNs = 0;
  entry.parentId =
      gLocalSlot.depth > 1 ? gLocalSlot.stack[gLocalSlot.depth - 2].id : 0;
  copyName(entry.name, name);
  entry.id = hashPath(entry.parentId, entry.name);
}

void Profiler::endZone() {
  ThreadSlot* slot = gLocalSlot.slot;
  if (!slot) return;
  if (gLocalSlot.overflowDepth > 0) {
    --gLocalSlot.overflowDepth;
    return;
  }
  if (gLocalSlot.depth == 0) return;

  const Clock::time_point finished = Clock::now();
  StackEntry entry = gLocalSlot.stack[--gLocalSlot.depth];
  const std::uint64_t inclusiveNs = static_cast<std::uint64_t>(
      std::chrono::duration_cast<std::chrono::nanoseconds>(finished - entry.started).count());
  const std::uint64_t exclusiveNs = inclusiveNs >= entry.childNs ? inclusiveNs - entry.childNs : 0;
  if (gLocalSlot.depth > 0) gLocalSlot.stack[gLocalSlot.depth - 1].childNs += inclusiveNs;

  std::uint64_t epoch = 0;
  ThreadBuffer& buffer = acquireCurrentBuffer(epoch);
  for (std::uint16_t i = 0; i < gLocalSlot.depth; ++i) {
    const StackEntry& ancestor = gLocalSlot.stack[i];
    ensureZone(buffer, ancestor.id, ancestor.parentId, ancestor.name);
  }
  if (ZoneAccum* zone = ensureZone(buffer, entry.id, entry.parentId, entry.name)) {
    ++zone->calls;
    zone->inclusiveNs += inclusiveNs;
    zone->exclusiveNs += exclusiveNs;
  }
  buffer.writers.fetch_sub(1, std::memory_order_release);
}

void Profiler::beginFrame() {
  if (!gEnabled.load(std::memory_order_relaxed)) {
    gFrameOpen.store(false, std::memory_order_release);
    return;
  }
  gFrameStarted = Clock::now();
  gFrameOpen.store(true, std::memory_order_release);
}

void Profiler::endFrame() {
  if (!gFrameOpen.exchange(false, std::memory_order_acq_rel)) return;

  const Clock::time_point finished = Clock::now();
  const std::uint64_t completedEpoch = gEpoch.fetch_add(1, std::memory_order_acq_rel);
  const std::size_t completedBuffer = completedEpoch & 1u;

  ProfileFrame frame;
  frame.sequence = gFrameSequence++;
  frame.frameMs = std::chrono::duration<double, std::milli>(finished - gFrameStarted).count();

  for (ThreadSlot& slot : gThreadSlots) {
    ThreadBuffer& buffer = slot.buffers[completedBuffer];
    while (buffer.writers.load(std::memory_order_acquire) != 0) std::this_thread::yield();
    mergeBuffer(buffer, frame);
  }

  for (const ProfileZone& zone : frame.zones) frame.cpuMs += zone.exclusiveMs;

  std::lock_guard lock(gFramesMutex);
  if (gFrameCount < gFrames.size()) {
    const std::size_t index = (gFrameHead + gFrameCount) % gFrames.size();
    gFrames[index] = std::move(frame);
    ++gFrameCount;
  } else {
    gFrames[gFrameHead] = std::move(frame);
    gFrameHead = (gFrameHead + 1) % gFrames.size();
  }
}

void Profiler::setEnabled(bool enabledValue) {
  gEnabled.store(enabledValue, std::memory_order_release);
}

bool Profiler::enabled() const {
  return gEnabled.load(std::memory_order_relaxed);
}

void Profiler::counter(const char* name, double value) {
  if (!gEnabled.load(std::memory_order_relaxed)) return;
  ThreadSlot* slot = acquireThreadSlot();
  if (!slot) return;

  std::uint64_t epoch = 0;
  ThreadBuffer& buffer = acquireCurrentBuffer(epoch);
  char safeName[kNameCapacity];
  copyName(safeName, name);
  CounterAccum* destination = nullptr;
  for (std::uint16_t i = 0; i < buffer.counterCount; ++i) {
    if (std::strcmp(buffer.counters[i].name, safeName) == 0) {
      destination = &buffer.counters[i];
      break;
    }
  }
  if (!destination && buffer.counterCount < buffer.counters.size()) {
    destination = &buffer.counters[buffer.counterCount++];
    copyName(destination->name, safeName);
  }
  if (destination) destination->value = value;
  buffer.writers.fetch_sub(1, std::memory_order_release);
}

std::vector<ProfileFrame> Profiler::frames() const {
  std::lock_guard lock(gFramesMutex);
  std::vector<ProfileFrame> snapshot;
  snapshot.reserve(gFrameCount);
  for (std::size_t i = 0; i < gFrameCount; ++i) {
    snapshot.push_back(gFrames[(gFrameHead + i) % gFrames.size()]);
  }
  return snapshot;
}

void Profiler::reset() {
  std::lock_guard lock(gFramesMutex);
  gFrameHead = 0;
  gFrameCount = 0;
  for (ProfileFrame& frame : gFrames) frame = {};
}

std::string Profiler::exportCsv() const {
  std::lock_guard lock(gFramesMutex);
  std::ostringstream csv;
  csv << "frame,frame_ms,cpu_ms,zone,parent,calls,inclusive_ms,exclusive_ms,share_percent\n";
  csv << std::setprecision(9);
  for (std::size_t frameOffset = 0; frameOffset < gFrameCount; ++frameOffset) {
    const ProfileFrame& frame = gFrames[(gFrameHead + frameOffset) % gFrames.size()];
    for (const ProfileZone& zone : frame.zones) {
      const std::string parent =
          zone.parent >= 0 && zone.parent < static_cast<int>(frame.zones.size())
              ? frame.zones[static_cast<std::size_t>(zone.parent)].name
              : std::string{};
      const double share = frame.cpuMs > 0.0 ? zone.exclusiveMs * 100.0 / frame.cpuMs : 0.0;
      csv << frame.sequence << ',' << frame.frameMs << ',' << frame.cpuMs << ','
          << csvField(zone.name) << ',' << csvField(parent) << ',' << zone.calls << ','
          << zone.inclusiveMs << ',' << zone.exclusiveMs << ',' << share << '\n';
    }
  }
  return csv.str();
}

ScopedZone::ScopedZone(const char* name) : active_(profiler().enabled()) {
  if (active_) profiler().beginZone(name);
}

ScopedZone::~ScopedZone() {
  if (active_) profiler().endZone();
}

}  // namespace chemcad::core

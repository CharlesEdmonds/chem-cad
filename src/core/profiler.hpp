#pragma once
// Scoped CPU instrumentation keeps nesting stacks in thread_local storage and
// writes completed samples to fixed-capacity per-thread buffers. endFrame()
// merges those buffers, so worker threads never allocate or take a mutex on the
// zone hot path.

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#ifndef CHEMCAD_ENABLE_PROFILING
#define CHEMCAD_ENABLE_PROFILING 1
#endif

namespace chemcad::core {

struct ProfileZone {
  std::uint64_t id = 0;
  std::string name;
  std::uint32_t calls = 0;
  double inclusiveMs = 0.0;
  double exclusiveMs = 0.0;
  int parent = -1;
};

struct ProfileCounter {
  std::string name;
  double value = 0.0;
};

struct ProfileFrame {
  std::uint64_t sequence = 0;
  double frameMs = 0.0;
  double cpuMs = 0.0;
  std::vector<ProfileZone> zones;
  std::vector<ProfileCounter> counters;
};

class Profiler {
 public:
  static constexpr std::size_t frameCapacity = 240;

  void beginZone(const char* name);
  void endZone();
  void beginFrame();
  void endFrame();

  void setEnabled(bool enabled);
  bool enabled() const;
  void counter(const char* name, double value);

  std::vector<ProfileFrame> frames() const;
  void reset();
  std::string exportCsv() const;

 private:
  Profiler() = default;
  friend Profiler& profiler();
};

Profiler& profiler();

class ScopedZone {
 public:
  explicit ScopedZone(const char* name);
  ~ScopedZone();

  ScopedZone(const ScopedZone&) = delete;
  ScopedZone& operator=(const ScopedZone&) = delete;

 private:
  bool active_ = false;
};

}  // namespace chemcad::core

#define CHEMCAD_PROFILE_JOIN_INNER(a, b) a##b
#define CHEMCAD_PROFILE_JOIN(a, b) CHEMCAD_PROFILE_JOIN_INNER(a, b)

#if CHEMCAD_ENABLE_PROFILING
#define CHEMCAD_PROFILE_ZONE(name) \
  ::chemcad::core::ScopedZone CHEMCAD_PROFILE_JOIN(chemcadProfileZone_, __COUNTER__)(name)
#else
#define CHEMCAD_PROFILE_ZONE(name) ((void)0)
#endif

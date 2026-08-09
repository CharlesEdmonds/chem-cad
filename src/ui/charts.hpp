#pragma once
// Live instrumentation across ChemCAD panels shares one visual language here.
// The module deliberately uses ImDrawList instead of a charting dependency
// because the Solubility Suite already defines the application's charting idiom.

#include "imgui.h"

#include <cstddef>
#include <vector>

#include "ui/theme.hpp"

namespace chemcad::ui::charts {

// Fixed-capacity ring of recent samples for a live instrument trace.
class Trace {
 public:
  explicit Trace(std::size_t capacity = 240);
  void push(double value);
  void clear();
  std::size_t size() const;
  bool empty() const;
  double operator[](std::size_t index) const;
  double latest() const;
  double minimum() const;
  double maximum() const;

 private:
  void recomputeExtrema();

  std::vector<double> values_;
  std::size_t head_ = 0;
  std::size_t size_ = 0;
  double minimum_ = 0.0;
  double maximum_ = 0.0;
};

struct SparklineStyle {
  ImVec4 accent = style::col::Accent;
  bool fill = true;
  bool showLatest = true;
  bool autoFloor = false;
  double floorValue = 0.0;
  double ceilingValue = 0.0;
};
void sparkline(const char* id, const Trace& trace, ImVec2 size, const SparklineStyle& style);

struct MeterStyle {
  ImVec4 accent = style::col::Teal;
  double warnAt = 0.0;
  double dangerAt = 0.0;
  bool showTrack = true;
};
void meter(const char* id, double fraction, ImVec2 size, const MeterStyle& style);

// One instrument tile: caption, big monospace value with unit, and a sparkline.
// `value` is already formatted by the caller; `unit` may be nullptr.
void instrument(const char* id, const char* caption, const char* value, const char* unit,
                const Trace& trace, ImVec2 size, const SparklineStyle& style);

// Ranked horizontal bars on one shared scale. Returns the clicked row index, else -1.
struct BarRow {
  const char* label = nullptr;
  double value = 0.0;
  const char* annotation = nullptr;
  ImVec4 accent = style::col::Teal;
  bool selected = false;
};
int rankedBars(const char* id, const BarRow* rows, int count, ImVec2 size);

// Proportional stacked bar with inline labels, e.g. phase volumes.
struct StackSegment {
  const char* label = nullptr;
  double value = 0.0;
  ImVec4 colour = style::col::Teal;
};
void stackedBar(const char* id, const StackSegment* segments, int count, ImVec2 size);

}  // namespace chemcad::ui::charts

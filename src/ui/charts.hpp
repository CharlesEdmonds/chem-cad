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

// Radial gauge: a 240-degree arc, a needle, and the value in the well. Reads at
// a glance from across a bench, which a numeric label does not.
struct GaugeStyle {
  ImVec4 accent = style::col::Accent;
  double warnAt = 0.0;    // fraction of full scale; 0 disables the band
  double dangerAt = 0.0;
  const char* minLabel = nullptr;
  const char* maxLabel = nullptr;
};
void gauge(const char* id, double fraction, const char* value, const char* caption,
           ImVec2 size, const GaugeStyle& style);

// Ring breakdown with a total in the hole. Same data as stackedBar, used where
// the question is "what is this made of" rather than "how does it stratify".
void donut(const char* id, const StackSegment* segments, int count, const char* centreValue,
           const char* centreCaption, ImVec2 size);

// One curve of an XY plot. `x` may be null, in which case samples are placed at
// their index. Ownership stays with the caller for the duration of the call.
struct Series {
  const char* label = nullptr;
  const double* x = nullptr;
  const double* y = nullptr;
  int count = 0;
  ImVec4 colour = style::col::Accent;
  bool fill = false;
  bool markers = false;
  bool dashed = false;
};

struct PlotStyle {
  const char* xLabel = nullptr;
  const char* yLabel = nullptr;
  bool logY = false;
  bool legend = true;
  bool grid = true;
  // Equal min/max means autoscale that axis over every series.
  double xMin = 0.0;
  double xMax = 0.0;
  double yMin = 0.0;
  double yMax = 0.0;
  // Vertical rule at a data-space x, e.g. the current operating point.
  bool hasCursor = false;
  double cursorX = 0.0;
};

// Axed line/area plot with a hover readout. Returns the data-space x under the
// pointer, or NaN when the pointer is elsewhere, so callers can drive a linked
// selection without owning any hit-testing of their own.
double linePlot(const char* id, const Series* series, int count, ImVec2 size,
                const PlotStyle& style);

// Row-major scalar field as a colour ramp. Used for composition surfaces and
// property matrices, where a table of numbers hides the shape of the data.
struct HeatmapStyle {
  ImVec4 low = style::col::BgSurface;
  ImVec4 mid = style::col::Teal;
  ImVec4 high = style::col::Accent;
  const char* xLabel = nullptr;
  const char* yLabel = nullptr;
  bool showScale = true;
  bool interpolate = true;
};
// Returns the linear index of the hovered cell, or -1.
int heatmap(const char* id, const double* values, int columns, int rows, ImVec2 size,
            const HeatmapStyle& style);

// Spider chart over shared normalised axes: the shape of a solvent's character
// (polarity, H-bonding, boiling point, safety) is recognisable where four
// separate bars are not.
struct RadarSeries {
  const char* label = nullptr;
  const double* values = nullptr;  // axisCount entries, each already 0..1
  ImVec4 colour = style::col::Teal;
  bool filled = true;
};
void radar(const char* id, const char* const* axisLabels, int axisCount,
           const RadarSeries* series, int seriesCount, ImVec2 size);

// A measured value against its target and its acceptable band. The bullet chart
// answers "is this in spec, and by how much" in one row.
void bullet(const char* id, double value, double target, double rangeMin, double rangeMax,
            ImVec2 size, ImVec4 accent);

}  // namespace chemcad::ui::charts

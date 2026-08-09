#include "ui/charts.hpp"
#include "ui/layout.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

namespace chemcad::ui::charts {

namespace {

// imgui.h does not export IM_PI -- that macro lives in imgui_internal.h, and
// this module stays off ImGui internals on purpose.
constexpr float kPi = 3.14159265358979323846f;

struct DrawRect {
  ImVec2 min;
  ImVec2 max;
  ImVec2 size;
};

ImVec2 validSize(ImVec2 size) {
  if (!std::isfinite(size.x) || size.x < 0.0f) size.x = 0.0f;
  if (!std::isfinite(size.y) || size.y < 0.0f) size.y = 0.0f;
  return size;
}

DrawRect reserveRect(const char* id, ImVec2 requestedSize) {
  const ImVec2 size = validSize(requestedSize);
  const ImVec2 min = ImGui::GetCursorScreenPos();
  if (id && size.x > 0.0f && size.y > 0.0f) {
    ImGui::InvisibleButton(id, size);
  } else {
    ImGui::Dummy(size);
  }
  return {min, ImVec2(min.x + size.x, min.y + size.y), size};
}

bool hasArea(const DrawRect& rect) {
  return rect.size.x > 0.0f && rect.size.y > 0.0f;
}

bool nearlyEqual(double left, double right) {
  const double scale = std::max({1.0, std::abs(left), std::abs(right)});
  return std::abs(left - right) <= std::numeric_limits<double>::epsilon() * scale * 8.0;
}

struct AxisBuilder {
  bool hasValue = false;
  double low = 0.0;
  double high = 0.0;

  void add(double value) {
    if (!std::isfinite(value)) return;
    if (!hasValue) {
      low = value;
      high = value;
      hasValue = true;
      return;
    }
    low = std::min(low, value);
    high = std::max(high, value);
  }

  Axis linear(bool includeZero, int targetTicks) const {
    if (!hasValue) return {};
    return niceAxis(includeZero ? std::min(0.0, low) : low,
                    includeZero ? std::max(0.0, high) : high, targetTicks);
  }
};

float sampleY(double value, const Axis& axis, float top, float bottom) {
  return bottom - static_cast<float>(axis.normalise(value)) * (bottom - top);
}

void drawSparkline(ImDrawList* drawList, const DrawRect& rect, const Trace& trace,
                   const SparklineStyle& sparkStyle) {
  if (!hasArea(rect)) return;

  const style::Metrics& metrics = style::metrics();
  const float fontSize = ImGui::GetFontSize();
  const float inset = std::max(metrics.hairline, fontSize * 0.04f);
  const float left = std::min(rect.max.x, rect.min.x + inset);
  const float right = std::max(left, rect.max.x - inset);
  const float top = std::min(rect.max.y, rect.min.y + inset);
  const float bottom = std::max(top, rect.max.y - inset);

  drawList->AddLine(ImVec2(left, bottom), ImVec2(right, bottom),
                    style::u32(style::col::GridLine, 0.65f), metrics.hairline);
  if (trace.empty()) return;

  const bool flatTrace = nearlyEqual(trace.minimum(), trace.maximum());
  AxisBuilder bounds;
  bounds.add(trace.minimum());
  bounds.add(trace.maximum());
  const double requestedFloor =
      std::isfinite(sparkStyle.floorValue) ? sparkStyle.floorValue : 0.0;
  if (!flatTrace) {
    if (!sparkStyle.autoFloor) bounds.add(requestedFloor);
    if (std::isfinite(sparkStyle.ceilingValue) &&
        sparkStyle.ceilingValue > requestedFloor) {
      bounds.add(sparkStyle.ceilingValue);
    }
  }
  const Axis axis = bounds.linear(false, 4);

  const std::size_t count = trace.size();
  const float span = right - left;
  const auto pointAt = [&](std::size_t index) {
    const float x = count > 1
                        ? left + span * static_cast<float>(index) / static_cast<float>(count - 1)
                        : (left + right) * 0.5f;
    return ImVec2(x, sampleY(trace[index], axis, top, bottom));
  };

  if (sparkStyle.fill && !flatTrace && count > 1) {
    ImVec2 previous = pointAt(0);
    for (std::size_t index = 1; index < count; ++index) {
      const ImVec2 current = pointAt(index);
      const ImVec2 quad[] = {
          previous, current, ImVec2(current.x, bottom), ImVec2(previous.x, bottom)};
      drawList->AddConvexPolyFilled(quad, 4, style::u32(sparkStyle.accent, 0.13f));
      previous = current;
    }
  }

  const float strokeWidth = metrics.hairline * 2.0f;
  if (count == 1) {
    const ImVec2 point = pointAt(0);
    const float halfLine = std::min(span * 0.5f, fontSize * 0.22f);
    drawList->AddLine(ImVec2(point.x - halfLine, point.y),
                      ImVec2(point.x + halfLine, point.y),
                      style::u32(sparkStyle.accent), strokeWidth);
  } else {
    for (std::size_t index = 0; index < count; ++index) {
      drawList->PathLineTo(pointAt(index));
    }
    drawList->PathStroke(style::u32(sparkStyle.accent), ImDrawFlags_None, strokeWidth);
  }

  if (sparkStyle.showLatest) {
    const ImVec2 latest = pointAt(count - 1);
    const float haloRadius = fontSize * 0.29f;
    const float pointRadius = fontSize * 0.13f;
    drawList->AddCircleFilled(latest, haloRadius, style::u32(sparkStyle.accent, 0.18f));
    drawList->AddCircleFilled(latest, pointRadius, style::u32(sparkStyle.accent));
  }
}

std::string uppercase(const char* text) {
  std::string result = text ? text : "";
  for (char& character : result) {
    character = static_cast<char>(std::toupper(static_cast<unsigned char>(character)));
  }
  return result;
}

float chartLabelFont(float requested) {
  return layout::labelFont(requested);
}

ImVec2 chartTextSize(const char* text, float requested) {
  const float fontSize = chartLabelFont(requested);
  return ImGui::GetFont()->CalcTextSizeA(
      fontSize, std::numeric_limits<float>::max(), 0.0f, text ? text : "");
}

void drawChartText(ImDrawList* drawList, ImVec2 position, ImU32 colour,
                   const char* text, float requested) {
  const float fontSize = chartLabelFont(requested);
  drawList->AddText(nullptr, fontSize, position, colour, text ? text : "");
}

bool chartTextFits(const char* text, float width, float requested) {
  return width >= 0.0f && chartTextSize(text, requested).x <= width;
}

}  // namespace

double Axis::normalise(double value) const {
  const double span = max - min;
  if (!std::isfinite(value)) return 0.0;
  if (!std::isfinite(span) || span <= 0.0) return 0.5;
  return std::clamp((value - min) / span, 0.0, 1.0);
}

Axis niceAxis(double low, double high, int targetTicks) {
  if (std::isfinite(low) && std::isfinite(high) && low > high) {
    std::swap(low, high);
  }

  if (!std::isfinite(low) || !std::isfinite(high) || nearlyEqual(low, high)) {
    if (!std::isfinite(low) && !std::isfinite(high)) return {};
    const double centre = std::isfinite(low) ? low : high;
    low = centre - 0.5;
    high = centre + 0.5;
  }

  targetTicks = std::max(1, targetTicks);
  const double rawStep = (high - low) / static_cast<double>(targetTicks);
  if (!std::isfinite(rawStep) || rawStep <= 0.0) return {};

  double exponent = std::floor(std::log10(rawStep));
  double power = std::pow(10.0, exponent);
  if (!std::isfinite(power) || power <= 0.0) return {};
  const double mantissa = rawStep / power;
  double snapped = 1.0;
  if (mantissa > 5.0) {
    snapped = 1.0;
    power *= 10.0;
    exponent += 1.0;
  } else if (mantissa > 2.0) {
    snapped = 5.0;
  } else if (mantissa > 1.0) {
    snapped = 2.0;
  }
  const double step = snapped * power;
  if (!std::isfinite(step) || step <= 0.0) return {};

  double axisMin = std::floor(low / step) * step;
  double axisMax = std::ceil(high / step) * step;
  if (!(axisMin < low)) axisMin -= step;
  if (!(axisMax > high)) axisMax += step;
  if (!std::isfinite(axisMin) || !std::isfinite(axisMax) ||
      !(axisMax > axisMin)) {
    return {};
  }

  const double rawTicks = std::round((axisMax - axisMin) / step) + 1.0;
  Axis axis;
  axis.min = axisMin;
  axis.max = axisMax;
  axis.step = step;
  axis.decimals =
      std::clamp(static_cast<int>(std::max(0.0, -exponent)), 0, 6);
  axis.ticks = std::clamp(
      static_cast<int>(std::min(rawTicks, static_cast<double>(std::numeric_limits<int>::max()))),
      2, std::numeric_limits<int>::max());
  return axis;
}

Axis niceLogAxis(double low, double high, double fallbackLow) {
  if (!std::isfinite(fallbackLow) || fallbackLow <= 0.0) fallbackLow = 1.0e-9;
  if (!std::isfinite(low) || low <= 0.0) low = fallbackLow;
  if (!std::isfinite(high) || high <= 0.0) high = fallbackLow;
  if (low > high) std::swap(low, high);

  double minDecade = std::floor(std::log10(low));
  double maxDecade = std::ceil(std::log10(high));
  if (!std::isfinite(minDecade) || !std::isfinite(maxDecade)) return {};
  if (maxDecade - minDecade < 1.0) {
    if (nearlyEqual(low, high)) {
      minDecade = std::floor(std::log10(low));
      maxDecade = minDecade + 1.0;
    } else {
      maxDecade = minDecade + 1.0;
    }
  }

  Axis axis;
  axis.min = minDecade;
  axis.max = maxDecade;
  axis.step = 1.0;
  axis.decimals = 0;
  axis.ticks = std::max(2, static_cast<int>(maxDecade - minDecade) + 1);
  return axis;
}

Axis axisFor(const double* values, int count, bool includeZero) {
  if (!values || count <= 0) return {};
  AxisBuilder bounds;
  for (int index = 0; index < count; ++index) bounds.add(values[index]);
  return bounds.linear(includeZero, 5);
}

Trace::Trace(std::size_t capacity) : values_(capacity) {}

void Trace::push(double value) {
  if (!std::isfinite(value) || values_.empty()) return;

  if (size_ == 0) {
    values_[0] = value;
    size_ = 1;
    minimum_ = value;
    maximum_ = value;
    return;
  }

  if (size_ < values_.size()) {
    const std::size_t index = (head_ + size_) % values_.size();
    values_[index] = value;
    ++size_;
    minimum_ = std::min(minimum_, value);
    maximum_ = std::max(maximum_, value);
    return;
  }

  const double overwritten = values_[head_];
  values_[head_] = value;
  head_ = (head_ + 1) % values_.size();

  // A scan is only needed when eviction may have removed an extremum; ordinary
  // live updates retain O(1) insertion without making every query pay for it.
  if (overwritten == minimum_ || overwritten == maximum_) {
    recomputeExtrema();
  } else {
    minimum_ = std::min(minimum_, value);
    maximum_ = std::max(maximum_, value);
  }
}

void Trace::clear() {
  head_ = 0;
  size_ = 0;
  minimum_ = 0.0;
  maximum_ = 0.0;
}

std::size_t Trace::size() const {
  return size_;
}

bool Trace::empty() const {
  return size_ == 0;
}

double Trace::operator[](std::size_t index) const {
  if (index >= size_ || values_.empty()) return 0.0;
  return values_[(head_ + index) % values_.size()];
}

double Trace::latest() const {
  return empty() ? 0.0 : (*this)[size_ - 1];
}

double Trace::minimum() const {
  return empty() ? 0.0 : minimum_;
}

double Trace::maximum() const {
  return empty() ? 0.0 : maximum_;
}

void Trace::recomputeExtrema() {
  if (empty()) {
    minimum_ = 0.0;
    maximum_ = 0.0;
    return;
  }

  minimum_ = (*this)[0];
  maximum_ = minimum_;
  for (std::size_t index = 1; index < size_; ++index) {
    minimum_ = std::min(minimum_, (*this)[index]);
    maximum_ = std::max(maximum_, (*this)[index]);
  }
}

void sparkline(const char* id, const Trace& trace, ImVec2 size,
               const SparklineStyle& sparkStyle) {
  const DrawRect rect = reserveRect(id, size);
  if (!hasArea(rect)) return;

  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawList->PushClipRect(rect.min, rect.max, true);
  drawSparkline(drawList, rect, trace, sparkStyle);
  drawList->PopClipRect();
}

void meter(const char* id, double fraction, ImVec2 size, const MeterStyle& meterStyle) {
  const DrawRect rect = reserveRect(id, size);
  if (!hasArea(rect)) return;

  if (!std::isfinite(fraction)) fraction = 0.0;
  const double clamped = std::clamp(fraction, 0.0, 1.0);
  ImVec4 fillColour = meterStyle.accent;
  if (std::isfinite(meterStyle.warnAt) && meterStyle.warnAt > 0.0 &&
      clamped >= meterStyle.warnAt) {
    fillColour = style::col::DataBright;
  }
  if (std::isfinite(meterStyle.dangerAt) && meterStyle.dangerAt > 0.0 &&
      clamped >= meterStyle.dangerAt) {
    fillColour = style::col::Danger;
  }

  const style::Metrics& metrics = style::metrics();
  const float radius = std::min(metrics.radiusSm, rect.size.y * 0.5f);
  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawList->PushClipRect(rect.min, rect.max, true);
  if (meterStyle.showTrack) {
    drawList->AddRectFilled(rect.min, rect.max, style::u32(style::col::DataDim, 0.18f), radius);
  }

  const float fillWidth = rect.size.x * static_cast<float>(clamped);
  if (fillWidth > 0.0f) {
    const float fillRadius = std::min(radius, fillWidth * 0.5f);
    drawList->AddRectFilled(rect.min, ImVec2(rect.min.x + fillWidth, rect.max.y),
                            style::u32(fillColour), fillRadius);
  }

  const auto drawThreshold = [&](double threshold) {
    if (!std::isfinite(threshold) || threshold <= 0.0) return;
    const float x = rect.min.x + rect.size.x * static_cast<float>(std::clamp(threshold, 0.0, 1.0));
    drawList->AddLine(ImVec2(x, rect.min.y), ImVec2(x, rect.max.y),
                      style::u32(style::col::GridLine), metrics.hairline);
  };
  drawThreshold(meterStyle.warnAt);
  drawThreshold(meterStyle.dangerAt);
  drawList->AddRect(rect.min, rect.max, style::u32(style::col::DataDim), radius,
                    metrics.hairline, ImDrawFlags_None);
  drawList->PopClipRect();
}

void instrument(const char* id, const char* caption, const char* value, const char* unit,
                const Trace& trace, ImVec2 size, const SparklineStyle& sparkStyle) {
  const DrawRect rect = reserveRect(id, size);
  if (!hasArea(rect)) return;

  const style::Metrics& metrics = style::metrics();
  const float fontSize = ImGui::GetFontSize();
  const float padding = metrics.gap * 0.9f;
  const float captionSize = chartLabelFont(fontSize * 0.72f);
  const std::string upperCaption = uppercase(caption);
  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawList->PushClipRect(rect.min, rect.max, true);
  drawList->AddRectFilled(rect.min, rect.max, style::u32(style::col::BgSurface), metrics.radiusMd);
  drawList->AddRect(rect.min, rect.max, style::u32(style::col::DataDim), metrics.radiusMd,
                    metrics.hairline, ImDrawFlags_None);

  const ImVec2 captionPosition(rect.min.x + padding, rect.min.y + padding);
  const float textWidth = std::max(0.0f, rect.size.x - padding * 2.0f);
  if (chartTextFits(upperCaption.c_str(), textWidth, captionSize)) {
    drawChartText(drawList, captionPosition, style::u32(style::col::DataDim),
                  upperCaption.c_str(), captionSize);
  }

  const float valueY = captionPosition.y + captionSize + metrics.gap * 0.35f;
  const char* shownValue = value ? value : "";
  const float valueFont = chartLabelFont(fontSize);
  const bool monoPushed = style::pushFont(style::fonts::mono());
  const ImVec2 valueSize = chartTextSize(shownValue, valueFont);
  const ImVec4 valueColour = trace.empty() ? style::col::DataDim : style::col::DataBright;
  if (valueSize.x <= textWidth) {
    drawChartText(drawList, ImVec2(rect.min.x + padding, valueY),
                  style::u32(valueColour), shownValue, valueFont);
  }
  style::popFont(monoPushed);

  if (unit && unit[0] != '\0') {
    const float unitFont = chartLabelFont(fontSize * 0.82f);
    const ImVec2 unitSize = chartTextSize(unit, unitFont);
    const float unitX = rect.min.x + padding + valueSize.x + fontSize * 0.28f;
    const float unitY = valueY + (valueFont - unitSize.y) * 0.5f;
    if (unitX + unitSize.x <= rect.max.x - padding) {
      drawChartText(drawList, ImVec2(unitX, unitY), style::u32(style::col::DataDim),
                    unit, unitFont);
    }
  }

  const float traceTop = valueY + fontSize + metrics.gap * 0.55f;
  const float traceBottom = rect.max.y - padding;
  if (traceBottom > traceTop) {
    const DrawRect traceRect = {
        ImVec2(rect.min.x + padding, traceTop),
        ImVec2(rect.max.x - padding, traceBottom),
        ImVec2(std::max(0.0f, rect.size.x - padding * 2.0f), traceBottom - traceTop)};
    drawSparkline(drawList, traceRect, trace, sparkStyle);
  }
  drawList->PopClipRect();
}

int rankedBars(const char* id, const BarRow* rows, int count, ImVec2 size) {
  const DrawRect rect = reserveRect(id, size);
  if (!hasArea(rect) || !rows || count <= 0) return -1;

  const style::Metrics& metrics = style::metrics();
  const float fontSize = ImGui::GetFontSize();
  const float rowHeight = ImGui::GetTextLineHeightWithSpacing();
  const float padding = metrics.gap * 0.75f;
  const float columnGap = metrics.gap * 0.75f;
  const float innerWidth = std::max(0.0f, rect.size.x - padding * 2.0f);

  std::vector<double> magnitudes(static_cast<std::size_t>(count));
  const float labelFont = chartLabelFont(fontSize * 0.88f);
  float widestLabel = 0.0f;
  for (int index = 0; index < count; ++index) {
    magnitudes[static_cast<std::size_t>(index)] =
        std::isfinite(rows[index].value) ? std::abs(rows[index].value)
                                         : std::numeric_limits<double>::quiet_NaN();
    if (rows[index].label) {
      widestLabel =
          std::max(widestLabel, chartTextSize(rows[index].label, labelFont).x);
    }
  }
  const Axis magnitudeAxis = axisFor(magnitudes.data(), count, true);

  float widestAnnotation = 0.0f;
  const bool monoPushed = style::pushFont(style::fonts::mono());
  for (int index = 0; index < count; ++index) {
    if (rows[index].annotation) {
      widestAnnotation = std::max(
          widestAnnotation, chartTextSize(rows[index].annotation, labelFont).x);
    }
  }
  style::popFont(monoPushed);

  const float labelWidth = std::min(widestLabel, innerWidth * 0.34f);
  const float annotationWidth = std::min(widestAnnotation, innerWidth * 0.28f);
  const float barLeft = rect.min.x + padding + labelWidth + columnGap;
  const float annotationLeft = rect.max.x - padding - annotationWidth;
  const float barRight = annotationWidth > 0.0f ? annotationLeft - columnGap : rect.max.x - padding;
  const float barSpan = std::max(0.0f, barRight - barLeft);

  int hoveredRow = -1;
  if (ImGui::IsItemHovered()) {
    const float relativeY = ImGui::GetIO().MousePos.y - rect.min.y;
    if (relativeY >= 0.0f) {
      const int candidate = static_cast<int>(relativeY / rowHeight);
      if (candidate >= 0 && candidate < count) hoveredRow = candidate;
    }
  }
  const int clicked = hoveredRow >= 0 && ImGui::IsItemClicked(ImGuiMouseButton_Left)
                          ? hoveredRow
                          : -1;

  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawList->PushClipRect(rect.min, rect.max, true);
  for (int index = 0; index < count; ++index) {
    const BarRow& row = rows[index];
    const float rowTop = rect.min.y + rowHeight * static_cast<float>(index);
    const float rowBottom = rowTop + rowHeight;
    if (rowTop >= rect.max.y) break;

    if (row.selected || index == hoveredRow) {
      const ImU32 background = index == hoveredRow
                                   ? style::mix(style::col::BgRaised, style::col::BorderStrong,
                                                row.selected ? 0.24f : 0.14f)
                                   : style::u32(style::col::BgRaised);
      drawList->AddRectFilled(ImVec2(rect.min.x, rowTop), ImVec2(rect.max.x, rowBottom),
                              background, metrics.radiusSm);
    }
    if (row.selected) {
      const float edgeWidth = metrics.hairline * 2.5f;
      drawList->AddRectFilled(ImVec2(rect.min.x, rowTop),
                              ImVec2(rect.min.x + edgeWidth, rowBottom),
                              style::u32(row.accent), metrics.radiusSm,
                              ImDrawFlags_RoundCornersLeft);
    }

    const float textY = rowTop + (rowHeight - labelFont) * 0.5f;
    if (row.label) {
      const float labelRight = std::max(rect.min.x + padding, barLeft - columnGap * 0.5f);
      const float available = labelRight - rect.min.x - padding;
      if (chartTextFits(row.label, available, labelFont)) {
        drawChartText(drawList, ImVec2(rect.min.x + padding, textY),
                      style::u32(style::col::DataBright), row.label, labelFont);
      }
    }

    if (barSpan > 0.0f) {
      const float barHeight = rowHeight * 0.28f;
      const float barTop = rowTop + (rowHeight - barHeight) * 0.5f;
      const float barBottom = barTop + barHeight;
      drawList->AddRectFilled(ImVec2(barLeft, barTop), ImVec2(barRight, barBottom),
                              style::u32(style::col::DataDim, 0.18f), barHeight * 0.5f);
      const double zero = magnitudeAxis.normalise(0.0);
      const double available = std::max(1.0 - zero, std::numeric_limits<double>::epsilon());
      const double normalized =
          (magnitudeAxis.normalise(magnitudes[static_cast<std::size_t>(index)]) - zero) /
          available;
      const float width = barSpan * static_cast<float>(std::clamp(normalized, 0.0, 1.0));
      if (width > 0.0f) {
        const ImU32 barColour = index == hoveredRow
                                   ? style::mix(row.accent, style::col::Text, 0.18f)
                                   : style::u32(row.accent);
        drawList->AddRectFilled(ImVec2(barLeft, barTop),
                                ImVec2(barLeft + width, barBottom), barColour,
                                std::min(barHeight * 0.5f, width * 0.5f));
      }
    }

    if (row.annotation) {
      const bool annotationMono = style::pushFont(style::fonts::mono());
      const ImVec2 annotationSize = chartTextSize(row.annotation, labelFont);
      if (annotationSize.x <= annotationWidth) {
        drawChartText(drawList,
                      ImVec2(rect.max.x - padding - annotationSize.x, textY),
                      style::u32(style::col::DataDim), row.annotation, labelFont);
      }
      style::popFont(annotationMono);
    }
  }
  drawList->PopClipRect();
  return clicked;
}

void stackedBar(const char* id, const StackSegment* segments, int count, ImVec2 size) {
  const DrawRect rect = reserveRect(id, size);
  if (!hasArea(rect)) return;

  const style::Metrics& metrics = style::metrics();
  const float fontSize = ImGui::GetFontSize();
  const float labelFont = chartLabelFont(fontSize * 0.78f);
  const float radius = std::min(metrics.radiusSm, rect.size.y * 0.5f);
  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawList->PushClipRect(rect.min, rect.max, true);
  drawList->AddRectFilled(rect.min, rect.max,
                          style::u32(style::col::DataDim, 0.14f), radius);

  const int sampleCount = segments && count > 0 ? count : 0;
  std::vector<double> values(static_cast<std::size_t>(sampleCount));
  for (int index = 0; index < sampleCount; ++index) {
    values[static_cast<std::size_t>(index)] =
        std::isfinite(segments[index].value) && segments[index].value > 0.0
            ? segments[index].value
            : std::numeric_limits<double>::quiet_NaN();
  }
  const Axis valueAxis = axisFor(values.data(), sampleCount, true);
  const double scale = valueAxis.max > 0.0 ? valueAxis.max : 1.0;
  double scaledTotal = 0.0;
  int first = -1;
  int last = -1;
  for (int index = 0; index < sampleCount; ++index) {
    const double sample = values[static_cast<std::size_t>(index)];
    if (!std::isfinite(sample)) continue;
    scaledTotal += sample / scale;
    if (first < 0) first = index;
    last = index;
  }

  if (scaledTotal > 0.0) {
    const Axis totalAxis = niceAxis(0.0, scaledTotal);
    const double axisZero = totalAxis.normalise(0.0);
    const double axisEnd = totalAxis.normalise(scaledTotal);
    const double axisSpan =
        std::max(axisEnd - axisZero, std::numeric_limits<double>::epsilon());
    double cumulative = 0.0;
    float x = rect.min.x;
    for (int index = 0; index < sampleCount; ++index) {
      const double sample = values[static_cast<std::size_t>(index)];
      if (!std::isfinite(sample)) continue;
      const double scaled = sample / scale;
      cumulative += scaled;
      const double proportion = scaled / scaledTotal;
      const double normalized =
          (totalAxis.normalise(cumulative) - axisZero) / axisSpan;
      const float nextX =
          index == last
              ? rect.max.x
              : rect.min.x + rect.size.x *
                                 static_cast<float>(std::clamp(normalized, 0.0, 1.0));

      ImDrawFlags corners = ImDrawFlags_RoundCornersNone;
      if (index == first && index == last) {
        corners = ImDrawFlags_RoundCornersAll;
      } else if (index == first) {
        corners = ImDrawFlags_RoundCornersLeft;
      } else if (index == last) {
        corners = ImDrawFlags_RoundCornersRight;
      }
      drawList->AddRectFilled(ImVec2(x, rect.min.y), ImVec2(nextX, rect.max.y),
                              style::u32(segments[index].colour), radius, corners);

      char label[256];
      if (segments[index].label && segments[index].label[0] != '\0') {
        std::snprintf(label, sizeof(label), "%s  %.0f%%", segments[index].label,
                      proportion * 100.0);
      } else {
        std::snprintf(label, sizeof(label), "%.0f%%", proportion * 100.0);
      }
      const ImVec2 textSize = chartTextSize(label, labelFont);
      const float textPadding = metrics.gap * 0.55f;
      if (nextX - x >= textSize.x + textPadding * 2.0f &&
          rect.size.y >= textSize.y) {
        drawChartText(
            drawList,
            ImVec2(x + (nextX - x - textSize.x) * 0.5f,
                   rect.min.y + (rect.size.y - textSize.y) * 0.5f),
            style::u32(style::col::OnAccent), label, labelFont);
      }
      if (index != last) {
        drawList->AddLine(ImVec2(nextX, rect.min.y), ImVec2(nextX, rect.max.y),
                          style::u32(style::col::GridLine), metrics.hairline);
      }
      x = nextX;
    }
  }

  drawList->AddRect(rect.min, rect.max, style::u32(style::col::DataDim), radius,
                    metrics.hairline, ImDrawFlags_None);
  drawList->PopClipRect();
}


void gauge(const char* id, double fraction, const char* value, const char* caption,
           ImVec2 size, const GaugeStyle& gaugeStyle) {
  const DrawRect rect = reserveRect(id, size);
  if (!hasArea(rect)) return;

  if (!std::isfinite(fraction)) fraction = 0.0;
  const float normalized = static_cast<float>(std::clamp(fraction, 0.0, 1.0));
  const style::Metrics& metrics = style::metrics();
  const float fontSize = ImGui::GetFontSize();
  const float padding = metrics.gap * 0.55f;
  const float captionSpace = fontSize * 1.35f;
  const float radius = std::max(
      0.0f, std::min(rect.size.x * 0.5f - padding,
                     (rect.size.y - captionSpace - padding) * 0.64f));
  const ImVec2 centre((rect.min.x + rect.max.x) * 0.5f,
                      rect.min.y + padding + radius);
  const float startAngle = kPi * (5.0f / 6.0f);
  const float endAngle = kPi * (13.0f / 6.0f);
  const float sweep = endAngle - startAngle;
  const float trackWidth = std::max(metrics.hairline * 2.0f, fontSize * 0.24f);
  const int arcSegments = std::max(
      12, static_cast<int>(radius / std::max(metrics.hairline, fontSize * 0.08f)));

  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawList->PushClipRect(rect.min, rect.max, true);
  drawList->AddRectFilled(rect.min, rect.max, style::u32(style::col::BgSurface),
                          metrics.radiusMd);
  drawList->AddRect(rect.min, rect.max, style::u32(style::col::DataDim),
                    metrics.radiusMd, metrics.hairline, ImDrawFlags_None);
  if (radius <= 0.0f) {
    drawList->PopClipRect();
    return;
  }

  const auto strokeArc = [&](float arcRadius, float from, float to, ImU32 colour,
                             float thickness) {
    if (to <= from || arcRadius <= 0.0f) return;
    drawList->PathArcTo(centre, arcRadius, from, to, arcSegments);
    drawList->PathStroke(colour, ImDrawFlags_None, thickness);
  };

  strokeArc(radius, startAngle, endAngle, style::u32(style::col::DataDim, 0.45f), trackWidth);

  const auto thresholdAngle = [&](double threshold) {
    if (!std::isfinite(threshold)) return startAngle;
    return startAngle + sweep *
                            static_cast<float>(std::clamp(threshold, 0.0, 1.0));
  };
  const bool hasWarn = std::isfinite(gaugeStyle.warnAt) && gaugeStyle.warnAt > 0.0;
  const bool hasDanger =
      std::isfinite(gaugeStyle.dangerAt) && gaugeStyle.dangerAt > 0.0;
  const float bandRadius = radius + trackWidth * 0.78f;
  const float bandWidth = std::max(metrics.hairline, trackWidth * 0.28f);
  if (hasWarn) {
    const float warnStart = thresholdAngle(gaugeStyle.warnAt);
    const float warnEnd =
        hasDanger && gaugeStyle.dangerAt > gaugeStyle.warnAt
            ? thresholdAngle(gaugeStyle.dangerAt)
            : endAngle;
    strokeArc(bandRadius, warnStart, warnEnd, style::u32(style::col::DataBright),
              bandWidth);
  }
  if (hasDanger) {
    strokeArc(bandRadius, thresholdAngle(gaugeStyle.dangerAt), endAngle,
              style::u32(style::col::Danger), bandWidth);
  }

  const float valueAngle = startAngle + sweep * normalized;
  strokeArc(radius, startAngle, valueAngle, style::u32(gaugeStyle.accent),
            trackWidth);

  const float needleLength = radius * 0.73f;
  const ImVec2 needleTip(centre.x + std::cos(valueAngle) * needleLength,
                         centre.y + std::sin(valueAngle) * needleLength);
  drawList->AddLine(centre, needleTip, style::u32(style::col::DataBright),
                    metrics.hairline * 1.5f);
  drawList->AddCircleFilled(centre, fontSize * 0.17f,
                            style::u32(gaugeStyle.accent));
  drawList->AddCircle(centre, fontSize * 0.17f, style::u32(style::col::DataDim),
                      0, metrics.hairline);

  const char* shownValue = value ? value : "";
  const float valueFont = chartLabelFont(fontSize);
  const bool monoPushed = style::pushFont(style::fonts::mono());
  const ImVec2 valueSize = chartTextSize(shownValue, valueFont);
  if (valueSize.x <= rect.size.x - padding * 2.0f) {
    drawChartText(
        drawList,
        ImVec2(centre.x - valueSize.x * 0.5f, centre.y + fontSize * 0.28f),
        style::u32(style::col::DataBright), shownValue, valueFont);
  }
  style::popFont(monoPushed);

  const char* shownCaption = caption ? caption : "";
  const float captionFont = chartLabelFont(fontSize * 0.78f);
  const ImVec2 captionSize = chartTextSize(shownCaption, captionFont);
  if (captionSize.x <= rect.size.x - padding * 2.0f) {
    drawChartText(
        drawList,
        ImVec2(centre.x - captionSize.x * 0.5f, centre.y + fontSize * 1.18f),
        style::u32(style::col::DataDim), shownCaption, captionFont);
  }

  const auto drawEndLabel = [&](const char* label, float angle, bool rightAligned) {
    if (!label || label[0] == '\0') return;
    const float endLabelFont = chartLabelFont(fontSize * 0.72f);
    const ImVec2 labelSize = chartTextSize(label, endLabelFont);
    if (labelSize.x > radius) return;
    const float labelRadius = std::max(0.0f, radius - trackWidth * 0.25f);
    const ImVec2 anchor(centre.x + std::cos(angle) * labelRadius,
                        centre.y + std::sin(angle) * labelRadius);
    const float x = rightAligned ? anchor.x - labelSize.x : anchor.x;
    drawChartText(drawList, ImVec2(x, anchor.y + fontSize * 0.16f),
                  style::u32(style::col::DataDim), label, endLabelFont);
  };
  drawEndLabel(gaugeStyle.minLabel, startAngle, true);
  drawEndLabel(gaugeStyle.maxLabel, endAngle, false);
  drawList->PopClipRect();
}

void donut(const char* id, const StackSegment* segments, int count,
           const char* centreValue, const char* centreCaption, ImVec2 size) {
  const DrawRect rect = reserveRect(id, size);
  if (!hasArea(rect)) return;

  const style::Metrics& metrics = style::metrics();
  const float fontSize = ImGui::GetFontSize();
  const ImVec2 centre((rect.min.x + rect.max.x) * 0.5f,
                      (rect.min.y + rect.max.y) * 0.5f);
  const float radius =
      std::max(0.0f, std::min(rect.size.x, rect.size.y) * 0.5f - metrics.gap * 0.45f);
  const float ringWidth = radius * 0.28f;
  const float ringRadius = std::max(0.0f, radius - ringWidth * 0.5f);
  const int arcSegments = std::max(
      12, static_cast<int>(radius / std::max(metrics.hairline, fontSize * 0.08f)));

  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawList->PushClipRect(rect.min, rect.max, true);
  if (ringRadius > 0.0f && ringWidth > 0.0f) {
    drawList->PathArcTo(centre, ringRadius, 0.0f, kPi * 2.0f, arcSegments);
    drawList->PathStroke(style::u32(style::col::DataDim, 0.24f),
                         ImDrawFlags_Closed, ringWidth);
  }

  const int sampleCount = segments && count > 0 ? count : 0;
  std::vector<double> values(static_cast<std::size_t>(sampleCount));
  for (int index = 0; index < sampleCount; ++index) {
    values[static_cast<std::size_t>(index)] =
        std::isfinite(segments[index].value) && segments[index].value > 0.0
            ? segments[index].value
            : std::numeric_limits<double>::quiet_NaN();
  }
  const Axis valueAxis = axisFor(values.data(), sampleCount, true);
  const double scale = valueAxis.max > 0.0 ? valueAxis.max : 1.0;
  double scaledTotal = 0.0;
  int validSegments = 0;
  for (double sample : values) {
    if (!std::isfinite(sample)) continue;
    scaledTotal += sample / scale;
    ++validSegments;
  }

  if (scaledTotal > 0.0 && ringRadius > 0.0f && ringWidth > 0.0f) {
    const Axis totalAxis = niceAxis(0.0, scaledTotal);
    const double axisZero = totalAxis.normalise(0.0);
    const double axisEnd = totalAxis.normalise(scaledTotal);
    const double axisSpan =
        std::max(axisEnd - axisZero, std::numeric_limits<double>::epsilon());
    const float gapAngle =
        validSegments > 1
            ? std::min(kPi * 0.025f,
                       metrics.hairline / std::max(ringRadius, metrics.hairline))
            : 0.0f;
    double cumulative = 0.0;
    float angle = -kPi * 0.5f;
    for (int index = 0; index < sampleCount; ++index) {
      const double sample = values[static_cast<std::size_t>(index)];
      if (!std::isfinite(sample)) continue;
      const double nextCumulative = cumulative + sample / scale;
      const float nextAngle =
          nextCumulative >= scaledTotal
              ? -kPi * 0.5f + kPi * 2.0f
              : -kPi * 0.5f +
                    kPi * 2.0f *
                        static_cast<float>(
                            (totalAxis.normalise(nextCumulative) - axisZero) / axisSpan);
      const float from = angle + gapAngle * 0.5f;
      const float to = nextAngle - gapAngle * 0.5f;
      if (to > from) {
        const int segmentCount = std::max(
            2, static_cast<int>(arcSegments * (to - from) / (kPi * 2.0f)));
        drawList->PathArcTo(centre, ringRadius, from, to, segmentCount);
        drawList->PathStroke(style::u32(segments[index].colour),
                             ImDrawFlags_None, ringWidth);
      }
      cumulative = nextCumulative;
      angle = nextAngle;
    }
  }

  const float holeWidth = std::max(0.0f, (ringRadius - ringWidth * 0.5f) * 2.0f);
  const char* shownValue = centreValue ? centreValue : "";
  const float valueFont = chartLabelFont(fontSize);
  const bool monoPushed = style::pushFont(style::fonts::mono());
  const ImVec2 valueSize = chartTextSize(shownValue, valueFont);
  if (valueSize.x <= holeWidth) {
    drawChartText(
        drawList,
        ImVec2(centre.x - valueSize.x * 0.5f, centre.y - fontSize * 0.78f),
        style::u32(style::col::DataBright), shownValue, valueFont);
  }
  style::popFont(monoPushed);

  const char* shownCaption = centreCaption ? centreCaption : "";
  const float captionFont = chartLabelFont(fontSize * 0.72f);
  const ImVec2 captionSize = chartTextSize(shownCaption, captionFont);
  if (captionSize.x <= holeWidth) {
    drawChartText(
        drawList,
        ImVec2(centre.x - captionSize.x * 0.5f, centre.y + fontSize * 0.18f),
        style::u32(style::col::DataDim), shownCaption, captionFont);
  }
  drawList->PopClipRect();
}

double linePlot(const char* id, const Series* series, int count, ImVec2 size,
                const PlotStyle& plotStyle) {
  const DrawRect rect = reserveRect(id, size);
  const double noHover = std::numeric_limits<double>::quiet_NaN();
  if (!hasArea(rect)) return noHover;

  const style::Metrics& metrics = style::metrics();
  const float fontSize = ImGui::GetFontSize();
  const float tickFont = chartLabelFont(fontSize * 0.72f);
  const float axisLabelFont = chartLabelFont(fontSize * 0.78f);
  const float leftGutter = fontSize * 3.25f;
  const float rightGutter = metrics.gap * 0.45f;
  const float topGutter = metrics.gap * 0.45f;
  const float bottomGutter =
      fontSize * (plotStyle.xLabel && plotStyle.xLabel[0] != '\0' ? 2.35f : 1.45f);
  const DrawRect plot = {
      ImVec2(std::min(rect.max.x, rect.min.x + leftGutter),
             std::min(rect.max.y, rect.min.y + topGutter)),
      ImVec2(std::max(rect.min.x, rect.max.x - rightGutter),
             std::max(rect.min.y, rect.max.y - bottomGutter)),
      ImVec2(std::max(0.0f, rect.size.x - leftGutter - rightGutter),
             std::max(0.0f, rect.size.y - topGutter - bottomGutter))};

  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawList->AddRectFilled(rect.min, rect.max, style::u32(style::col::BgPanel),
                          metrics.radiusMd);
  drawList->AddRect(rect.min, rect.max, style::u32(style::col::DataDim),
                    metrics.radiusMd, metrics.hairline, ImDrawFlags_None);
  if (!hasArea(plot)) return noHover;
  drawList->AddRectFilled(plot.min, plot.max, style::u32(style::col::BgSurface),
                          metrics.radiusSm);

  AxisBuilder xBounds;
  AxisBuilder yBounds;
  const bool explicitX = std::isfinite(plotStyle.xMin) &&
                         std::isfinite(plotStyle.xMax) &&
                         !nearlyEqual(plotStyle.xMin, plotStyle.xMax);
  const bool explicitY = std::isfinite(plotStyle.yMin) &&
                         std::isfinite(plotStyle.yMax) &&
                         !nearlyEqual(plotStyle.yMin, plotStyle.yMax) &&
                         (!plotStyle.logY ||
                          (plotStyle.yMin > 0.0 && plotStyle.yMax > 0.0));
  if (explicitX) {
    xBounds.add(plotStyle.xMin);
    xBounds.add(plotStyle.xMax);
  }
  if (explicitY) {
    yBounds.add(plotStyle.yMin);
    yBounds.add(plotStyle.yMax);
  }
  if (series && count > 0) {
    for (int seriesIndex = 0; seriesIndex < count; ++seriesIndex) {
      const Series& item = series[seriesIndex];
      if (!item.y || item.count <= 0) continue;
      for (int sample = 0; sample < item.count; ++sample) {
        const double x = item.x ? item.x[sample] : static_cast<double>(sample);
        const double y = item.y[sample];
        if (!std::isfinite(x) || !std::isfinite(y) ||
            (plotStyle.logY && y <= 0.0)) {
          continue;
        }
        xBounds.add(x);
        yBounds.add(y);
      }
    }
  }

  const Axis xAxis = xBounds.linear(false, 5);
  const Axis yAxis =
      plotStyle.logY
          ? (yBounds.hasValue
                 ? niceLogAxis(yBounds.low, yBounds.high, yBounds.low)
                 : Axis{})
          : yBounds.linear(false, 4);
  const auto mapX = [&](double x) {
    return plot.min.x + static_cast<float>(xAxis.normalise(x)) * plot.size.x;
  };
  const auto mapY = [&](double y) {
    const double transformed = plotStyle.logY ? std::log10(y) : y;
    return plot.max.y - static_cast<float>(yAxis.normalise(transformed)) * plot.size.y;
  };
  const auto formatTick = [](char* buffer, std::size_t capacity,
                             const Axis& axis, double tick, bool logarithmic) {
    if (logarithmic) {
      std::snprintf(buffer, capacity, "1e%.0f", tick);
      return;
    }
    const double magnitude = std::abs(tick);
    if (magnitude >= 1.0e6 || (magnitude > 0.0 && magnitude < 1.0e-4)) {
      std::snprintf(buffer, capacity, "%.3g", tick);
    } else {
      std::snprintf(buffer, capacity, "%.*f", axis.decimals, tick);
    }
  };

  const float yLabelWidth = std::max(0.0f, leftGutter - metrics.gap * 0.65f);
  if (plotStyle.grid) {
    float widestXTick = 0.0f;
    char label[64];
    for (int tick = 0; tick < xAxis.ticks; ++tick) {
      const double value = xAxis.min + xAxis.step * static_cast<double>(tick);
      formatTick(label, sizeof(label), xAxis, value, false);
      widestXTick = std::max(widestXTick, chartTextSize(label, tickFont).x);
    }
    const float xSlot =
        plot.size.x / static_cast<float>(std::max(1, xAxis.ticks - 1));
    const int xStride =
        std::max(1, static_cast<int>(std::ceil(
                        (widestXTick + metrics.gap * 0.35f) /
                        std::max(xSlot, metrics.hairline))));
    for (int tick = 0; tick < xAxis.ticks; ++tick) {
      const double value = xAxis.min + xAxis.step * static_cast<double>(tick);
      const float screenX = mapX(value);
      drawList->AddLine(ImVec2(screenX, plot.min.y), ImVec2(screenX, plot.max.y),
                        style::u32(style::col::GridLine), metrics.hairline);
      if (tick % xStride != 0 && tick != xAxis.ticks - 1) continue;
      formatTick(label, sizeof(label), xAxis, value, false);
      const ImVec2 labelSize = chartTextSize(label, tickFont);
      if (labelSize.x > xSlot * static_cast<float>(xStride)) continue;
      drawChartText(drawList,
                    ImVec2(screenX - labelSize.x * 0.5f,
                           plot.max.y + fontSize * 0.22f),
                    style::u32(style::col::DataDim), label, tickFont);
    }

    const float ySlot =
        plot.size.y / static_cast<float>(std::max(1, yAxis.ticks - 1));
    const int yStride =
        std::max(1, static_cast<int>(std::ceil(
                        (tickFont + metrics.gap * 0.2f) /
                        std::max(ySlot, metrics.hairline))));

    for (int tick = 0; tick < yAxis.ticks; ++tick) {
      const double value = yAxis.min + yAxis.step * static_cast<double>(tick);
      const float screenY =
          plot.max.y - static_cast<float>(yAxis.normalise(value)) * plot.size.y;
      drawList->AddLine(ImVec2(plot.min.x, screenY), ImVec2(plot.max.x, screenY),
                        style::u32(style::col::GridLine), metrics.hairline);
      if (tick % yStride != 0 && tick != yAxis.ticks - 1) continue;
      formatTick(label, sizeof(label), yAxis, value, plotStyle.logY);
      const ImVec2 labelSize = chartTextSize(label, tickFont);
      if (labelSize.x > yLabelWidth) continue;
      drawChartText(drawList,
                    ImVec2(plot.min.x - metrics.gap * 0.35f - labelSize.x,
                           screenY - labelSize.y * 0.5f),
                    style::u32(style::col::DataDim), label, tickFont);
    }
  }

  if (plotStyle.xLabel && plotStyle.xLabel[0] != '\0' &&
      chartTextFits(plotStyle.xLabel, plot.size.x, axisLabelFont)) {
    const ImVec2 labelSize = chartTextSize(plotStyle.xLabel, axisLabelFont);
    drawChartText(
        drawList,
        ImVec2(plot.min.x + (plot.size.x - labelSize.x) * 0.5f,
               rect.max.y - labelSize.y - metrics.gap * 0.22f),
        style::u32(style::col::DataDim), plotStyle.xLabel, axisLabelFont);
  }
  if (plotStyle.yLabel && plotStyle.yLabel[0] != '\0' &&
      chartTextFits(plotStyle.yLabel, yLabelWidth, axisLabelFont)) {
    drawChartText(drawList,
                  ImVec2(rect.min.x + metrics.gap * 0.25f,
                         rect.min.y + metrics.gap * 0.25f),
                  style::u32(style::col::DataDim), plotStyle.yLabel,
                  axisLabelFont);
  }

  drawList->PushClipRect(plot.min, plot.max, true);
  if (series && count > 0) {
    const float strokeWidth = metrics.hairline * 2.0f;
    const float markerRadius = fontSize * 0.12f;
    const float dashLength = fontSize * 0.46f;
    const float dashGap = fontSize * 0.27f;
    for (int seriesIndex = 0; seriesIndex < count; ++seriesIndex) {
      const Series& item = series[seriesIndex];
      if (!item.y || item.count <= 0) continue;
      const auto pointAt = [&](int sample, ImVec2& point) {
        const double x = item.x ? item.x[sample] : static_cast<double>(sample);
        const double y = item.y[sample];
        if (!std::isfinite(x) || !std::isfinite(y) ||
            (plotStyle.logY && y <= 0.0)) {
          return false;
        }
        point = ImVec2(mapX(x), mapY(y));
        return true;
      };

      int validPoints = 0;
      ImVec2 onlyPoint;
      for (int sample = 0; sample < item.count; ++sample) {
        ImVec2 point;
        if (pointAt(sample, point)) {
          ++validPoints;
          onlyPoint = point;
        }
      }
      for (int sample = 1; sample < item.count; ++sample) {
        ImVec2 previous;
        ImVec2 current;
        if (!pointAt(sample - 1, previous) || !pointAt(sample, current)) continue;
        if (item.fill) {
          const ImVec2 quad[] = {
              previous, current, ImVec2(current.x, plot.max.y),
              ImVec2(previous.x, plot.max.y)};
          drawList->AddConvexPolyFilled(quad, 4, style::u32(item.colour, 0.13f));
        }
        if (!item.dashed) {
          drawList->AddLine(previous, current, style::u32(item.colour), strokeWidth);
          continue;
        }

        const float deltaX = current.x - previous.x;
        const float deltaY = current.y - previous.y;
        const float length = std::sqrt(deltaX * deltaX + deltaY * deltaY);
        if (length <= 0.0f) continue;
        const float step = dashLength + dashGap;
        for (float offset = 0.0f; offset < length; offset += step) {
          const float finish = std::min(length, offset + dashLength);
          const float fromT = offset / length;
          const float toT = finish / length;
          drawList->AddLine(
              ImVec2(previous.x + deltaX * fromT, previous.y + deltaY * fromT),
              ImVec2(previous.x + deltaX * toT, previous.y + deltaY * toT),
              style::u32(item.colour), strokeWidth);
        }
      }
      if (item.markers) {
        for (int sample = 0; sample < item.count; ++sample) {
          ImVec2 point;
          if (pointAt(sample, point)) {
            drawList->AddCircleFilled(point, markerRadius, style::u32(item.colour));
          }
        }
      }
      if (validPoints == 1 && !item.markers) {
        drawList->AddCircleFilled(onlyPoint, markerRadius, style::u32(item.colour));
      }
    }
  }

  if (plotStyle.hasCursor && std::isfinite(plotStyle.cursorX)) {
    const float cursorX = mapX(plotStyle.cursorX);
    drawList->AddLine(ImVec2(cursorX, plot.min.y), ImVec2(cursorX, plot.max.y),
                      style::u32(style::col::Accent, 0.55f), metrics.hairline);
  }

  double hoveredX = noHover;
  if (ImGui::IsItemHovered() && series && count > 0 && series[0].count > 0) {
    const Series& first = series[0];
    const double mouseNormal =
        std::clamp(static_cast<double>((ImGui::GetIO().MousePos.x - plot.min.x) /
                                       std::max(plot.size.x, metrics.hairline)),
                   0.0, 1.0);
    const double mouseDataX =
        xAxis.min + mouseNormal * (xAxis.max - xAxis.min);
    double nearestDistance = std::numeric_limits<double>::infinity();
    for (int sample = 0; sample < first.count; ++sample) {
      const double x = first.x ? first.x[sample] : static_cast<double>(sample);
      const double y = first.y ? first.y[sample] : noHover;
      if (!std::isfinite(x) || !std::isfinite(y) ||
          (plotStyle.logY && y <= 0.0)) {
        continue;
      }
      const double distance = std::abs(x - mouseDataX);
      if (distance < nearestDistance) {
        nearestDistance = distance;
        hoveredX = x;
      }
    }

    if (std::isfinite(hoveredX)) {
      const float hoverScreenX = mapX(hoveredX);
      drawList->AddLine(ImVec2(hoverScreenX, plot.min.y),
                        ImVec2(hoverScreenX, plot.max.y),
                        style::u32(style::col::DataBright, 0.72f), metrics.hairline);
      char line[160];
      std::snprintf(line, sizeof(line), "x  %.3g", hoveredX);
      std::string tooltip = line;
      for (int seriesIndex = 0; seriesIndex < count; ++seriesIndex) {
        const Series& item = series[seriesIndex];
        int nearest = -1;
        double distance = std::numeric_limits<double>::infinity();
        if (item.y) {
          for (int sample = 0; sample < item.count; ++sample) {
            const double x = item.x ? item.x[sample] : static_cast<double>(sample);
            const double y = item.y[sample];
            if (!std::isfinite(x) || !std::isfinite(y) ||
                (plotStyle.logY && y <= 0.0)) {
              continue;
            }
            const double candidateDistance = std::abs(x - hoveredX);
            if (candidateDistance < distance) {
              distance = candidateDistance;
              nearest = sample;
            }
          }
        }
        const char* label =
            item.label && item.label[0] != '\0' ? item.label : "series";
        if (nearest >= 0) {
          const double x =
              item.x ? item.x[nearest] : static_cast<double>(nearest);
          const double y = item.y[nearest];
          drawList->AddCircleFilled(ImVec2(mapX(x), mapY(y)),
                                    fontSize * 0.18f,
                                    style::u32(item.colour));
          std::snprintf(line, sizeof(line), "\n%s  %.3g", label, y);
        } else {
          std::snprintf(line, sizeof(line), "\n%s  --", label);
        }
        tooltip += line;
      }
      ImGui::SetTooltip("%s", tooltip.c_str());
    }
  }
  drawList->PopClipRect();

  if (plotStyle.legend && series && count > 0) {
    const float legendFont = chartLabelFont(fontSize * 0.72f);
    const float swatch = legendFont * 0.58f;
    const float legendGap = metrics.gap * 0.45f;
    float legendX = plot.min.x + metrics.gap * 0.35f;
    const float legendY = plot.min.y + metrics.gap * 0.3f;
    for (int seriesIndex = 0; seriesIndex < count; ++seriesIndex) {
      const char* label = series[seriesIndex].label;
      if (!label || label[0] == '\0') continue;
      const ImVec2 textSize = chartTextSize(label, legendFont);
      const float entryWidth = swatch + metrics.gap * 0.28f + textSize.x;
      if (legendX + entryWidth > plot.max.x - metrics.gap * 0.35f) break;
      drawList->AddRectFilled(
          ImVec2(legendX, legendY + (legendFont - swatch) * 0.5f),
          ImVec2(legendX + swatch, legendY + (legendFont + swatch) * 0.5f),
          style::u32(series[seriesIndex].colour), metrics.radiusSm);
      legendX += swatch + metrics.gap * 0.28f;
      drawChartText(drawList, ImVec2(legendX, legendY),
                    style::u32(style::col::DataBright), label, legendFont);
      legendX += textSize.x + legendGap;
    }
  }

  drawList->AddRect(plot.min, plot.max, style::u32(style::col::DataDim),
                    metrics.radiusSm, metrics.hairline, ImDrawFlags_None);
  return hoveredX;
}

int heatmap(const char* id, const double* values, int columns, int rows, ImVec2 size,
            const HeatmapStyle& heatmapStyle) {
  const DrawRect rect = reserveRect(id, size);
  if (!hasArea(rect)) return -1;

  const style::Metrics& metrics = style::metrics();
  const float fontSize = ImGui::GetFontSize();
  const float labelFont = chartLabelFont(fontSize * 0.72f);
  const bool hasXLabel = heatmapStyle.xLabel && heatmapStyle.xLabel[0] != '\0';
  const bool hasYLabel = heatmapStyle.yLabel && heatmapStyle.yLabel[0] != '\0';
  const float leftGutter = hasYLabel ? fontSize * 1.1f : metrics.gap * 0.35f;
  const float bottomGutter = hasXLabel ? fontSize * 1.55f : metrics.gap * 0.35f;
  const float scaleGutter = heatmapStyle.showScale ? fontSize * 3.15f :
                                                      metrics.gap * 0.35f;
  const float topGutter = metrics.gap * 0.35f;
  const DrawRect field = {
      ImVec2(std::min(rect.max.x, rect.min.x + leftGutter),
             std::min(rect.max.y, rect.min.y + topGutter)),
      ImVec2(std::max(rect.min.x, rect.max.x - scaleGutter),
             std::max(rect.min.y, rect.max.y - bottomGutter)),
      ImVec2(std::max(0.0f, rect.size.x - leftGutter - scaleGutter),
             std::max(0.0f, rect.size.y - topGutter - bottomGutter))};

  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawList->PushClipRect(rect.min, rect.max, true);
  drawList->AddRectFilled(rect.min, rect.max, style::u32(style::col::BgPanel),
                          metrics.radiusMd);
  drawList->AddRect(rect.min, rect.max, style::u32(style::col::DataDim),
                    metrics.radiusMd, metrics.hairline, ImDrawFlags_None);
  if (columns <= 0 || rows <= 0 || !hasArea(field)) {
    drawList->PopClipRect();
    return -1;
  }

  const int sampleCount =
      columns <= std::numeric_limits<int>::max() / rows ? columns * rows : 0;
  const Axis valueAxis = axisFor(values, sampleCount, false);
  const auto colourFor = [&](double value) {
    if (!std::isfinite(value)) return style::u32(style::col::BgSurface);
    const double normalized = valueAxis.normalise(value);
    if (normalized <= 0.5) {
      return style::mix(heatmapStyle.low, heatmapStyle.mid,
                        static_cast<float>(normalized * 2.0));
    }
    return style::mix(heatmapStyle.mid, heatmapStyle.high,
                      static_cast<float>((normalized - 0.5) * 2.0));
  };
  const auto averageFinite = [](double a, double b, double c, double d) {
    double total = 0.0;
    int finiteCount = 0;
    const double candidates[] = {a, b, c, d};
    for (double candidate : candidates) {
      if (!std::isfinite(candidate)) continue;
      total += candidate;
      ++finiteCount;
    }
    return finiteCount > 0 ? total / static_cast<double>(finiteCount)
                           : std::numeric_limits<double>::quiet_NaN();
  };
  const auto cellValue = [&](int row, int column) {
    if (!values || row < 0 || row >= rows || column < 0 || column >= columns) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return values[row * columns + column];
  };
  const auto cornerValue = [&](int rowVertex, int columnVertex) {
    return averageFinite(cellValue(rowVertex - 1, columnVertex - 1),
                         cellValue(rowVertex - 1, columnVertex),
                         cellValue(rowVertex, columnVertex - 1),
                         cellValue(rowVertex, columnVertex));
  };

  const float cellWidth = field.size.x / static_cast<float>(columns);
  const float cellHeight = field.size.y / static_cast<float>(rows);
  for (int row = 0; row < rows; ++row) {
    for (int column = 0; column < columns; ++column) {
      const ImVec2 cellMin(field.min.x + cellWidth * static_cast<float>(column),
                           field.min.y + cellHeight * static_cast<float>(row));
      const ImVec2 cellMax(field.min.x + cellWidth * static_cast<float>(column + 1),
                           field.min.y + cellHeight * static_cast<float>(row + 1));
      if (!heatmapStyle.interpolate) {
        drawList->AddRectFilled(cellMin, cellMax,
                                colourFor(cellValue(row, column)));
        continue;
      }

      const double topLeft = cornerValue(row, column);
      const double topRight = cornerValue(row, column + 1);
      const double bottomLeft = cornerValue(row + 1, column);
      const double bottomRight = cornerValue(row + 1, column + 1);
      const float halfX = (cellMin.x + cellMax.x) * 0.5f;
      const float halfY = (cellMin.y + cellMax.y) * 0.5f;
      for (int subRow = 0; subRow < 2; ++subRow) {
        for (int subColumn = 0; subColumn < 2; ++subColumn) {
          const double u = (static_cast<double>(subColumn) + 0.5) * 0.5;
          const double v = (static_cast<double>(subRow) + 0.5) * 0.5;
          const double top = topLeft + (topRight - topLeft) * u;
          const double bottom = bottomLeft + (bottomRight - bottomLeft) * u;
          const double interpolated =
              std::isfinite(top) && std::isfinite(bottom)
                  ? top + (bottom - top) * v
                  : averageFinite(topLeft, topRight, bottomLeft, bottomRight);
          const ImVec2 subMin(subColumn == 0 ? cellMin.x : halfX,
                              subRow == 0 ? cellMin.y : halfY);
          const ImVec2 subMax(subColumn == 0 ? halfX : cellMax.x,
                              subRow == 0 ? halfY : cellMax.y);
          drawList->AddRectFilled(subMin, subMax, colourFor(interpolated));
        }
      }
    }
  }

  int hovered = -1;
  if (ImGui::IsItemHovered()) {
    const ImVec2 mouse = ImGui::GetIO().MousePos;
    if (mouse.x >= field.min.x && mouse.x < field.max.x &&
        mouse.y >= field.min.y && mouse.y < field.max.y) {
      const int column = std::clamp(
          static_cast<int>((mouse.x - field.min.x) / cellWidth), 0, columns - 1);
      const int row = std::clamp(
          static_cast<int>((mouse.y - field.min.y) / cellHeight), 0, rows - 1);
      hovered = row * columns + column;
      const ImVec2 hoverMin(field.min.x + cellWidth * static_cast<float>(column),
                            field.min.y + cellHeight * static_cast<float>(row));
      const ImVec2 hoverMax(
          field.min.x + cellWidth * static_cast<float>(column + 1),
          field.min.y + cellHeight * static_cast<float>(row + 1));
      drawList->AddRect(hoverMin, hoverMax, style::u32(style::col::Accent),
                        0.0f, metrics.hairline * 2.0f, ImDrawFlags_None);
      const double hoveredValue = cellValue(row, column);
      ImGui::SetTooltip("row %d  column %d\nvalue  %.3g", row, column,
                        hoveredValue);
    }
  }

  drawList->AddRect(field.min, field.max, style::u32(style::col::DataDim),
                    metrics.radiusSm, metrics.hairline, ImDrawFlags_None);
  if (heatmapStyle.showScale) {
    const float rampWidth = fontSize * 0.5f;
    const float rampX = field.max.x + metrics.gap * 0.45f;
    const int rampSteps = 32;
    for (int step = 0; step < rampSteps; ++step) {
      const float t0 = static_cast<float>(step) / static_cast<float>(rampSteps);
      const float t1 =
          static_cast<float>(step + 1) / static_cast<float>(rampSteps);
      const float y0 = field.max.y - field.size.y * t1;
      const float y1 = field.max.y - field.size.y * t0;
      const double value =
          valueAxis.min + (valueAxis.max - valueAxis.min) *
                              (static_cast<double>(t0 + t1) * 0.5);
      drawList->AddRectFilled(ImVec2(rampX, y0),
                              ImVec2(rampX + rampWidth, y1),
                              colourFor(value));
    }
    drawList->AddRect(ImVec2(rampX, field.min.y),
                      ImVec2(rampX + rampWidth, field.max.y),
                      style::u32(style::col::DataDim), metrics.radiusSm,
                      metrics.hairline, ImDrawFlags_None);
    const float labelX = rampX + rampWidth + metrics.gap * 0.25f;
    const float labelWidth = std::max(0.0f, rect.max.x - labelX);
    if (field.size.y >= labelFont * 2.0f) {
      char label[64];
      std::snprintf(label, sizeof(label), "%.*f", valueAxis.decimals, valueAxis.max);
      if (chartTextFits(label, labelWidth, labelFont)) {
        drawChartText(drawList, ImVec2(labelX, field.min.y),
                      style::u32(style::col::DataDim), label, labelFont);
      }
      std::snprintf(label, sizeof(label), "%.*f", valueAxis.decimals, valueAxis.min);
      const ImVec2 minLabelSize = chartTextSize(label, labelFont);
      if (minLabelSize.x <= labelWidth) {
        drawChartText(drawList, ImVec2(labelX, field.max.y - minLabelSize.y),
                      style::u32(style::col::DataDim), label, labelFont);
      }
    }
  }

  if (hasXLabel && chartTextFits(heatmapStyle.xLabel, field.size.x, labelFont)) {
    const ImVec2 labelSize = chartTextSize(heatmapStyle.xLabel, labelFont);
    drawChartText(
        drawList,
        ImVec2(field.min.x + (field.size.x - labelSize.x) * 0.5f,
               rect.max.y - labelSize.y - metrics.gap * 0.18f),
        style::u32(style::col::DataDim), heatmapStyle.xLabel, labelFont);
  }
  if (hasYLabel &&
      chartTextFits(heatmapStyle.yLabel,
                    std::max(0.0f, leftGutter - metrics.gap * 0.3f),
                    labelFont)) {
    const ImVec2 labelSize = chartTextSize(heatmapStyle.yLabel, labelFont);
    drawChartText(
        drawList,
        ImVec2(rect.min.x + metrics.gap * 0.15f,
               field.min.y + (field.size.y - labelSize.y) * 0.5f),
        style::u32(style::col::DataDim), heatmapStyle.yLabel, labelFont);
  }
  drawList->PopClipRect();
  return hovered;
}

void radar(const char* id, const char* const* axisLabels, int axisCount,
           const RadarSeries* series, int seriesCount, ImVec2 size) {
  const DrawRect rect = reserveRect(id, size);
  if (!hasArea(rect)) return;

  const style::Metrics& metrics = style::metrics();
  const float fontSize = ImGui::GetFontSize();
  const float labelFont = chartLabelFont(fontSize * 0.72f);
  const ImVec2 centre((rect.min.x + rect.max.x) * 0.5f,
                      (rect.min.y + rect.max.y) * 0.5f);
  const float radius =
      std::max(0.0f, std::min(rect.size.x, rect.size.y) * 0.5f - fontSize * 1.7f);

  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawList->PushClipRect(rect.min, rect.max, true);
  drawList->AddRectFilled(rect.min, rect.max, style::u32(style::col::BgPanel),
                          metrics.radiusMd);
  drawList->AddRect(rect.min, rect.max, style::u32(style::col::DataDim),
                    metrics.radiusMd, metrics.hairline, ImDrawFlags_None);
  if (axisCount < 3 || radius <= 0.0f) {
    drawList->PopClipRect();
    return;
  }

  const float angleStep = kPi * 2.0f / static_cast<float>(axisCount);
  const auto axisPoint = [&](int axis, float magnitude) {
    const float angle = -kPi * 0.5f + angleStep * static_cast<float>(axis);
    return ImVec2(centre.x + std::cos(angle) * radius * magnitude,
                  centre.y + std::sin(angle) * radius * magnitude);
  };
  for (int ring = 1; ring <= 4; ++ring) {
    const float magnitude = static_cast<float>(ring) / 4.0f;
    for (int axis = 0; axis < axisCount; ++axis) {
      drawList->PathLineTo(axisPoint(axis, magnitude));
    }
    drawList->PathStroke(style::u32(style::col::GridLine),
                         ImDrawFlags_Closed, metrics.hairline);
  }
  for (int axis = 0; axis < axisCount; ++axis) {
    drawList->AddLine(centre, axisPoint(axis, 1.0f),
                      style::u32(style::col::GridLine), metrics.hairline);
  }

  if (axisLabels) {
    const float labelSlot =
        std::max(labelFont, 2.0f * radius * std::sin(kPi / axisCount));
    for (int axis = 0; axis < axisCount; ++axis) {
      const char* label = axisLabels[axis];
      if (!label || label[0] == '\0' ||
          !chartTextFits(label, labelSlot, labelFont)) {
        continue;
      }
      const float angle = -kPi * 0.5f + angleStep * static_cast<float>(axis);
      const ImVec2 labelSize = chartTextSize(label, labelFont);
      const float labelRadius = radius + fontSize * 0.58f;
      ImVec2 position(centre.x + std::cos(angle) * labelRadius,
                      centre.y + std::sin(angle) * labelRadius -
                          labelSize.y * 0.5f);
      const float horizontal = std::cos(angle);
      if (horizontal < -0.18f) {
        position.x -= labelSize.x;
      } else if (horizontal <= 0.18f) {
        position.x -= labelSize.x * 0.5f;
      }
      if (position.x < rect.min.x || position.x + labelSize.x > rect.max.x ||
          position.y < rect.min.y || position.y + labelSize.y > rect.max.y) {
        continue;
      }
      drawChartText(drawList, position, style::u32(style::col::DataDim), label,
                    labelFont);
    }
  }

  if (series && seriesCount > 0) {
    const float pointRadius = fontSize * 0.11f;
    for (int seriesIndex = 0; seriesIndex < seriesCount; ++seriesIndex) {
      const RadarSeries& item = series[seriesIndex];
      if (!item.values) continue;
      const auto valid = [&](int axis) {
        return std::isfinite(item.values[axis]);
      };
      const auto seriesPoint = [&](int axis) {
        return axisPoint(
            axis, static_cast<float>(std::clamp(item.values[axis], 0.0, 1.0)));
      };
      for (int axis = 0; axis < axisCount; ++axis) {
        const int next = (axis + 1) % axisCount;
        if (!valid(axis) || !valid(next)) continue;
        if (item.filled) {
          drawList->AddTriangleFilled(centre, seriesPoint(axis),
                                      seriesPoint(next),
                                      style::u32(item.colour, 0.18f));
        }
        drawList->AddLine(seriesPoint(axis), seriesPoint(next),
                          style::u32(item.colour), metrics.hairline * 2.0f);
      }
      for (int axis = 0; axis < axisCount; ++axis) {
        if (!valid(axis)) continue;
        drawList->AddCircleFilled(seriesPoint(axis), pointRadius,
                                  style::u32(item.colour));
      }
    }
  }
  drawList->PopClipRect();
}

void bullet(const char* id, double value, double target, double rangeMin,
            double rangeMax, ImVec2 size, ImVec4 accent) {
  const DrawRect rect = reserveRect(id, size);
  if (!hasArea(rect)) return;

  const style::Metrics& metrics = style::metrics();
  const float fontSize = ImGui::GetFontSize();
  const float padding = metrics.gap * 0.55f;
  const float left = std::min(rect.max.x, rect.min.x + padding);
  const float right = std::max(left, rect.max.x - padding);
  const float centreY = (rect.min.y + rect.max.y) * 0.5f;
  const float trackHeight = fontSize * 0.34f;
  const float barHeight = fontSize * 0.58f;
  const float tickHeight = fontSize * 0.95f;
  const float radius = std::min(metrics.radiusSm, trackHeight * 0.5f);

  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawList->PushClipRect(rect.min, rect.max, true);
  drawList->AddRectFilled(ImVec2(left, centreY - trackHeight * 0.5f),
                          ImVec2(right, centreY + trackHeight * 0.5f),
                          style::u32(style::col::DataDim, 0.18f), radius);

  const double samples[] = {value, target, rangeMin, rangeMax};
  const Axis axis = axisFor(samples, 4, false);
  const auto mapValue = [&](double sample) {
    return left + (right - left) * static_cast<float>(axis.normalise(sample));
  };

  const bool hasRange = std::isfinite(rangeMin) && std::isfinite(rangeMax);
  if (hasRange && rangeMin > rangeMax) std::swap(rangeMin, rangeMax);
  if (hasRange && std::isfinite(target)) {
    const double clampedTarget = std::clamp(target, rangeMin, rangeMax);
    const double nearerEnd =
        std::abs(clampedTarget - rangeMin) <= std::abs(rangeMax - clampedTarget)
            ? rangeMin
            : rangeMax;
    const float targetX = mapValue(clampedTarget);
    const float bandEndX = mapValue(nearerEnd);
    drawList->AddRectFilled(
        ImVec2(std::min(targetX, bandEndX), centreY - trackHeight * 0.5f),
        ImVec2(std::max(targetX, bandEndX), centreY + trackHeight * 0.5f),
        style::u32(style::col::DataDim, 0.28f), radius);
  }

  if (std::isfinite(value)) {
    const bool outside = hasRange && (value < rangeMin || value > rangeMax);
    const float valueX = mapValue(value);
    const ImVec4 barColour = outside ? style::col::Danger : accent;
    if (valueX > left) {
      drawList->AddRectFilled(ImVec2(left, centreY - barHeight * 0.5f),
                              ImVec2(valueX, centreY + barHeight * 0.5f),
                              style::u32(barColour),
                              std::min(metrics.radiusSm, barHeight * 0.5f));
    } else {
      drawList->AddCircleFilled(ImVec2(left, centreY), barHeight * 0.5f,
                                style::u32(barColour));
    }
  }
  if (std::isfinite(target)) {
    const float targetX = mapValue(target);
    drawList->AddLine(ImVec2(targetX, centreY - tickHeight * 0.5f),
                      ImVec2(targetX, centreY + tickHeight * 0.5f),
                      style::u32(style::col::DataBright),
                      metrics.hairline * 2.0f);
  }
  drawList->PopClipRect();
}

}  // namespace chemcad::ui::charts

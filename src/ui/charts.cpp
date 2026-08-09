#include "ui/charts.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>

namespace chemcad::ui::charts {

namespace {

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

float sampleY(double value, double floor, double ceiling, bool degenerate,
              float top, float bottom) {
  if (degenerate) return (top + bottom) * 0.5f;
  const double normalized = std::clamp((value - floor) / (ceiling - floor), 0.0, 1.0);
  return bottom - static_cast<float>(normalized) * (bottom - top);
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
                    style::u32(style::col::Border, 0.65f), metrics.hairline);
  if (trace.empty()) return;

  const double traceMinimum = trace.minimum();
  const double traceMaximum = trace.maximum();
  double floor = sparkStyle.autoFloor ? traceMinimum : sparkStyle.floorValue;
  if (!std::isfinite(floor)) floor = 0.0;

  double ceiling = sparkStyle.ceilingValue;
  if (!std::isfinite(ceiling) || ceiling <= floor) ceiling = traceMaximum;

  const bool flatTrace = nearlyEqual(traceMinimum, traceMaximum);
  const bool invalidRange = !std::isfinite(ceiling) || ceiling <= floor || nearlyEqual(floor, ceiling);
  const bool degenerate = flatTrace || invalidRange;
  if (invalidRange) ceiling = floor;

  const std::size_t count = trace.size();
  const float span = right - left;
  const auto pointAt = [&](std::size_t index) {
    const float x = count > 1
                        ? left + span * static_cast<float>(index) / static_cast<float>(count - 1)
                        : (left + right) * 0.5f;
    return ImVec2(x, sampleY(trace[index], floor, ceiling, degenerate, top, bottom));
  };

  if (sparkStyle.fill && count > 1) {
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

double finiteMagnitude(double value) {
  return std::isfinite(value) ? std::abs(value) : 0.0;
}

}  // namespace

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
    fillColour = style::col::Accent;
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
    drawList->AddRectFilled(rect.min, rect.max, style::u32(style::col::BgSurface), radius);
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
                      style::u32(style::col::BorderStrong), metrics.hairline);
  };
  drawThreshold(meterStyle.warnAt);
  drawThreshold(meterStyle.dangerAt);
  drawList->AddRect(rect.min, rect.max, style::u32(style::col::Border), radius,
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
  const float captionSize = fontSize * 0.72f;
  const std::string upperCaption = uppercase(caption);
  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawList->PushClipRect(rect.min, rect.max, true);
  drawList->AddRectFilled(rect.min, rect.max, style::u32(style::col::BgSurface), metrics.radiusMd);
  drawList->AddRect(rect.min, rect.max, style::u32(style::col::Border), metrics.radiusMd,
                    metrics.hairline, ImDrawFlags_None);

  const ImVec2 captionPosition(rect.min.x + padding, rect.min.y + padding);
  drawList->AddText(nullptr, captionSize, captionPosition, style::u32(style::col::TextDim),
                    upperCaption.c_str());

  const float valueY = captionPosition.y + captionSize + metrics.gap * 0.35f;
  const char* shownValue = value ? value : "";
  const bool monoPushed = style::pushFont(style::fonts::mono());
  const ImVec2 valueSize = ImGui::CalcTextSize(shownValue);
  const ImVec4 valueColour = trace.empty() ? style::col::TextFaint : style::col::Text;
  drawList->AddText(ImVec2(rect.min.x + padding, valueY), style::u32(valueColour), shownValue);
  style::popFont(monoPushed);

  if (unit && unit[0] != '\0') {
    const ImVec2 unitSize = ImGui::CalcTextSize(unit);
    const float unitX = rect.min.x + padding + valueSize.x + fontSize * 0.28f;
    const float unitY = valueY + (fontSize - unitSize.y) * 0.5f;
    drawList->AddText(ImVec2(unitX, unitY), style::u32(style::col::TextFaint), unit);
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

  double maximumMagnitude = 0.0;
  float widestLabel = 0.0f;
  for (int index = 0; index < count; ++index) {
    maximumMagnitude = std::max(maximumMagnitude, finiteMagnitude(rows[index].value));
    if (rows[index].label) {
      widestLabel = std::max(widestLabel, ImGui::CalcTextSize(rows[index].label).x);
    }
  }

  float widestAnnotation = 0.0f;
  const bool monoPushed = style::pushFont(style::fonts::mono());
  for (int index = 0; index < count; ++index) {
    if (rows[index].annotation) {
      widestAnnotation =
          std::max(widestAnnotation, ImGui::CalcTextSize(rows[index].annotation).x);
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

    const float textY = rowTop + (rowHeight - fontSize) * 0.5f;
    if (row.label) {
      const float labelRight = std::max(rect.min.x + padding, barLeft - columnGap * 0.5f);
      drawList->PushClipRect(ImVec2(rect.min.x, rowTop), ImVec2(labelRight, rowBottom), true);
      drawList->AddText(ImVec2(rect.min.x + padding, textY), style::u32(style::col::Text),
                        row.label);
      drawList->PopClipRect();
    }

    if (barSpan > 0.0f) {
      const float barHeight = rowHeight * 0.28f;
      const float barTop = rowTop + (rowHeight - barHeight) * 0.5f;
      const float barBottom = barTop + barHeight;
      drawList->AddRectFilled(ImVec2(barLeft, barTop), ImVec2(barRight, barBottom),
                              style::u32(style::col::BgPanel), barHeight * 0.5f);
      const double normalized = maximumMagnitude > 0.0
                                    ? finiteMagnitude(row.value) / maximumMagnitude
                                    : 0.0;
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
      const ImVec2 annotationSize = ImGui::CalcTextSize(row.annotation);
      drawList->AddText(ImVec2(rect.max.x - padding - annotationSize.x, textY),
                        style::u32(style::col::TextDim), row.annotation);
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
  const float radius = std::min(metrics.radiusSm, rect.size.y * 0.5f);
  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawList->PushClipRect(rect.min, rect.max, true);
  drawList->AddRectFilled(rect.min, rect.max, style::u32(style::col::BgSurface), radius);

  double largest = 0.0;
  if (segments && count > 0) {
    for (int index = 0; index < count; ++index) {
      if (std::isfinite(segments[index].value) && segments[index].value > 0.0) {
        largest = std::max(largest, segments[index].value);
      }
    }
  }

  double scaledTotal = 0.0;
  int first = -1;
  int last = -1;
  if (largest > 0.0) {
    for (int index = 0; index < count; ++index) {
      const double value = segments[index].value;
      if (!std::isfinite(value) || value <= 0.0) continue;
      scaledTotal += value / largest;
      if (first < 0) first = index;
      last = index;
    }
  }

  float x = rect.min.x;
  if (scaledTotal > 0.0) {
    for (int index = 0; index < count; ++index) {
      const StackSegment& segment = segments[index];
      if (!std::isfinite(segment.value) || segment.value <= 0.0) continue;

      const double proportion = (segment.value / largest) / scaledTotal;
      const float segmentWidth = index == last
                                     ? rect.max.x - x
                                     : rect.size.x * static_cast<float>(proportion);
      const float nextX = std::min(rect.max.x, x + std::max(0.0f, segmentWidth));
      ImDrawFlags corners = ImDrawFlags_RoundCornersNone;
      if (index == first && index == last) {
        corners = ImDrawFlags_RoundCornersAll;
      } else if (index == first) {
        corners = ImDrawFlags_RoundCornersLeft;
      } else if (index == last) {
        corners = ImDrawFlags_RoundCornersRight;
      }
      drawList->AddRectFilled(ImVec2(x, rect.min.y), ImVec2(nextX, rect.max.y),
                              style::u32(segment.colour), radius, corners);

      char label[256];
      if (segment.label && segment.label[0] != '\0') {
        std::snprintf(label, sizeof(label), "%s  %.0f%%", segment.label, proportion * 100.0);
      } else {
        std::snprintf(label, sizeof(label), "%.0f%%", proportion * 100.0);
      }
      const ImVec2 textSize = ImGui::CalcTextSize(label);
      const float textPadding = metrics.gap * 0.55f;
      if (nextX - x >= textSize.x + textPadding * 2.0f &&
          rect.size.y >= textSize.y) {
        drawList->AddText(ImVec2(x + (nextX - x - textSize.x) * 0.5f,
                                 rect.min.y + (rect.size.y - textSize.y) * 0.5f),
                          style::u32(style::col::OnAccent), label);
      }

      if (index != last) {
        drawList->AddLine(ImVec2(nextX, rect.min.y), ImVec2(nextX, rect.max.y),
                          style::u32(style::col::Border), metrics.hairline);
      }
      x = nextX;
    }
  }

  drawList->AddRect(rect.min, rect.max, style::u32(style::col::Border), radius,
                    metrics.hairline, ImDrawFlags_None);
  drawList->PopClipRect();
}

}  // namespace chemcad::ui::charts

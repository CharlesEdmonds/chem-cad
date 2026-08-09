#include "ui/layout.hpp"

#include <algorithm>
#include <cmath>

#include "ui/theme.hpp"

namespace chemcad::ui::layout {

namespace {

float positiveFinite(const float value) {
  return std::isfinite(value) && value > 0.0f ? value : 0.0f;
}

}  // namespace

Frame measure() {
  return measure(ImGui::GetContentRegionAvail());
}

Frame measure(ImVec2 available) {
  Frame frame;
  frame.size = available;
  frame.em = ImGui::GetFontSize();
  frame.aspect = available.x / std::max(available.y, 1.0f);
  frame.wide = frame.aspect >= 1.4f;
  frame.tall = frame.aspect < 0.9f;

  if (frame.ems() < 46.0f || available.y < frame.em * 26.0f) {
    frame.density = Density::Compact;
  } else if (frame.ems() >= 96.0f && available.y >= frame.em * 40.0f) {
    frame.density = Density::Roomy;
  } else {
    frame.density = Density::Regular;
  }

  float spacingScale = 1.0f;
  if (frame.density == Density::Compact) {
    spacingScale = 0.65f;
  } else if (frame.density == Density::Roomy) {
    spacingScale = 1.35f;
  }
  frame.gap = style::metrics().gap * spacingScale;
  frame.pad = style::metrics().gap * spacingScale;
  frame.row = ImGui::GetTextLineHeightWithSpacing();
  frame.control = ImGui::GetFrameHeight();
  return frame;
}

float columnWidth(const Frame& frame, int columns, int span) {
  if (columns <= 0) return frame.size.x;

  const float gaps = static_cast<float>(columns - 1) * frame.gap;
  const float oneColumn = (frame.size.x - gaps) / static_cast<float>(columns);
  const float width = oneColumn * static_cast<float>(span) +
                      static_cast<float>(span - 1) * frame.gap;
  return std::max(width, 1.0f);
}

int columnsThatFit(const Frame& frame, float minEm) {
  if (!std::isfinite(frame.size.x) || frame.size.x <= 0.0f ||
      !std::isfinite(frame.em) || frame.em <= 0.0f ||
      !std::isfinite(minEm) || minEm <= 0.0f) {
    return 1;
  }

  const float minimum = minEm * frame.em;
  if (!std::isfinite(minimum) || minimum <= 1.0f) return 1;

  const float gap = positiveFinite(frame.gap);
  const float count = std::floor((frame.size.x + gap) / (minimum + gap));
  if (!std::isfinite(count) || count < 1.0f) return 1;
  return static_cast<int>(count);
}

void distribute(float budget, const float* weights, const float* minimums, int count,
                float gap, float* out) {
  if (count <= 0 || out == nullptr) return;

  const float gapTotal = static_cast<float>(count - 1) * gap;
  const float available = std::max(budget - gapTotal, 0.0f);

  double minimumTotal = 0.0;
  if (minimums != nullptr) {
    for (int i = 0; i < count; ++i) {
      minimumTotal += static_cast<double>(positiveFinite(minimums[i]));
    }
  }

  float assigned = 0.0f;
  if (minimums != nullptr && minimumTotal > static_cast<double>(available)) {
    const double scale = minimumTotal > 0.0
                             ? static_cast<double>(available) / minimumTotal
                             : 0.0;
    for (int i = 0; i < count - 1; ++i) {
      out[i] = std::floor(static_cast<float>(
          static_cast<double>(positiveFinite(minimums[i])) * scale));
      assigned += out[i];
    }
    // The last row absorbs fractional pixels so independent rounding cannot
    // overflow the content region and create a scrollbar.
    out[count - 1] = available - assigned;
    return;
  }

  double weightTotal = 0.0;
  for (int i = 0; i < count; ++i) {
    weightTotal += weights != nullptr
                       ? static_cast<double>(positiveFinite(weights[i]))
                       : 0.0;
  }
  const bool equalShares = !std::isfinite(weightTotal) || weightTotal <= 0.0;
  const float remainder = available - static_cast<float>(minimumTotal);

  for (int i = 0; i < count - 1; ++i) {
    const float base = minimums != nullptr ? positiveFinite(minimums[i]) : 0.0f;
    const double fraction = equalShares
                                ? 1.0 / static_cast<double>(count)
                                : static_cast<double>(positiveFinite(weights[i])) /
                                      weightTotal;
    out[i] = std::floor(base + remainder * static_cast<float>(fraction));
    assigned += out[i];
  }
  // Keeping the exact remainder in one place prevents sub-pixel losses from
  // accumulating across rows and unexpectedly enabling the window scrollbar.
  out[count - 1] = available - assigned;
}

void nextRow(float yPosition) {
  // A cursor move on its own never grows the content extent, so ImGui reports
  // "code uses SetCursorPos to extend window/parent boundaries" once per frame
  // whenever a band advance is the last thing in a window. Submitting a
  // zero-size item claims the position; positioning one ItemSpacing short lets
  // that item's own advance land the cursor exactly on `yPosition`.
  ImGui::SetCursorPosY(yPosition - ImGui::GetStyle().ItemSpacing.y);
  ImGui::Dummy(ImVec2(0.0f, 0.0f));
}

float pageHeight() {
  return ImGui::GetContentRegionAvail().y;
}

float minReadablePx() {
  const float em = ImGui::GetFontSize();
  // Below roughly 72% of body text, dense numeric labels stop being readable;
  // the second term keeps the legibility floor scaling with the user's font.
  return std::max(em * 0.72f, 11.0f * (em / 16.0f));
}

float labelFont(float requested) {
  const float em = ImGui::GetFontSize();
  if (!std::isfinite(requested)) return em;
  return std::clamp(requested, minReadablePx(), em);
}

bool fits(const char* text, float width) {
  if (text == nullptr || text[0] == '\0') return true;
  return ImGui::CalcTextSize(text).x <= width;
}

const char* bestLabel(float width, const char* const* candidates, int count) {
  if (candidates == nullptr || count <= 0) return "";

  const char* last = nullptr;
  for (int i = 0; i < count; ++i) {
    const char* candidate = candidates[i];
    if (candidate == nullptr) continue;
    last = candidate;
    if (fits(candidate, width)) return candidate;
  }
  return last != nullptr ? last : "";
}

}  // namespace chemcad::ui::layout

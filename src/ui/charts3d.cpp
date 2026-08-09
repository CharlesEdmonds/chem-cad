#include "ui/charts3d.hpp"

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <initializer_list>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "ui/charts.hpp"
#include "ui/layout.hpp"

namespace chemcad::ui::charts3d {

namespace {

namespace uiStyle = chemcad::ui::style;

constexpr float kPi = 3.14159265358979323846f;
constexpr float kDegreesToRadians = kPi / 180.0f;

struct Vec3 {
  float x = 0.0f;
  float y = 0.0f;
  float z = 0.0f;
};

Vec3 operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
Vec3 operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec3 operator*(Vec3 v, float scale) { return {v.x * scale, v.y * scale, v.z * scale}; }

float dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

Vec3 cross(Vec3 a, Vec3 b) {
  return {a.y * b.z - a.z * b.y,
          a.z * b.x - a.x * b.z,
          a.x * b.y - a.y * b.x};
}

float length(Vec3 v) { return std::sqrt(dot(v, v)); }

Vec3 normalise(Vec3 v) {
  const float magnitude = length(v);
  return magnitude > std::numeric_limits<float>::epsilon() ? v * (1.0f / magnitude)
                                                            : Vec3{0.0f, 0.0f, 1.0f};
}

Vec3 rotate(Vec3 p, const Orbit& orbit) {
  const float yaw = orbit.yawDeg * kDegreesToRadians;
  const float pitch = orbit.pitchDeg * kDegreesToRadians;
  const float cy = std::cos(yaw);
  const float sy = std::sin(yaw);
  const float cp = std::cos(pitch);
  const float sp = std::sin(pitch);
  const Vec3 yawed{cy * p.x + sy * p.z, p.y, -sy * p.x + cy * p.z};
  return {yawed.x, cp * yawed.y - sp * yawed.z, sp * yawed.y + cp * yawed.z};
}

struct DrawRect {
  ImVec2 min;
  ImVec2 max;
  ImVec2 size;
};

DrawRect reserveRect(const char* id, ImVec2 requested) {
  if (!std::isfinite(requested.x) || requested.x < 0.0f) requested.x = 0.0f;
  if (!std::isfinite(requested.y) || requested.y < 0.0f) requested.y = 0.0f;
  const ImVec2 min = ImGui::GetCursorScreenPos();
  if (id && requested.x > 0.0f && requested.y > 0.0f) {
    ImGui::InvisibleButton(id, requested);
  } else {
    ImGui::Dummy(requested);
  }
  return {min, ImVec2(min.x + requested.x, min.y + requested.y), requested};
}

bool hasArea(const DrawRect& rect) {
  return rect.size.x > std::numeric_limits<float>::epsilon() &&
         rect.size.y > std::numeric_limits<float>::epsilon();
}

void drawFrame(ImDrawList* drawList, const DrawRect& rect) {
  if (!hasArea(rect)) return;
  const auto& metrics = uiStyle::metrics();
  drawList->AddRectFilled(rect.min, rect.max, uiStyle::u32(uiStyle::col::BgSurface, 0.62f),
                          metrics.radiusSm);
  drawList->AddRect(rect.min, rect.max, uiStyle::u32(uiStyle::col::Border),
                    metrics.radiusSm, ImDrawFlags_None, metrics.hairline);
}

float chartLabelSize(float requestedScale = 0.76f) {
  return layout::labelFont(ImGui::GetFontSize() * requestedScale);
}

ImVec2 textSize(const char* text, float fontSize) {
  if (!text) return {};
  return ImGui::GetFont()->CalcTextSizeA(fontSize, FLT_MAX, 0.0f, text);
}

void addLabel(ImDrawList* drawList, ImVec2 position, ImVec4 colour, const char* text,
              float requestedScale = 0.76f) {
  if (!text || text[0] == '\0') return;
  const float fontSize = chartLabelSize(requestedScale);
  drawList->AddText(ImGui::GetFont(), fontSize, position, uiStyle::u32(colour), text);
}

struct Projection {
  ImVec2 centre{};
  float scale = 1.0f;
  float perspectiveDistance = 4.5f;
  Orbit orbit{};

  Vec3 camera(Vec3 p) const { return rotate(p, orbit); }

  ImVec2 projectCamera(Vec3 p) const {
    const float divisor = std::max(perspectiveDistance - p.z,
                                   perspectiveDistance * 0.16f);
    const float perspective = perspectiveDistance / divisor;
    return {centre.x + p.x * perspective * scale,
            centre.y - p.y * perspective * scale};
  }

  ImVec2 project(Vec3 p) const { return projectCamera(camera(p)); }
};

Projection fittedProjection(const DrawRect& rect, const Orbit& orbit,
                            const std::vector<Vec3>& points) {
  Projection projection;
  projection.orbit = orbit;
  projection.centre = ImVec2((rect.min.x + rect.max.x) * 0.5f,
                             (rect.min.y + rect.max.y) * 0.5f);
  if (!hasArea(rect) || points.empty()) return projection;

  float minX = std::numeric_limits<float>::max();
  float maxX = -std::numeric_limits<float>::max();
  float minY = std::numeric_limits<float>::max();
  float maxY = -std::numeric_limits<float>::max();
  for (Vec3 point : points) {
    const Vec3 camera = projection.camera(point);
    const float divisor = std::max(projection.perspectiveDistance - camera.z,
                                   projection.perspectiveDistance * 0.16f);
    const float perspective = projection.perspectiveDistance / divisor;
    minX = std::min(minX, camera.x * perspective);
    maxX = std::max(maxX, camera.x * perspective);
    minY = std::min(minY, camera.y * perspective);
    maxY = std::max(maxY, camera.y * perspective);
  }

  const float fontSize = ImGui::GetFontSize();
  const float inset = std::min(rect.size.x, rect.size.y) * 0.08f + fontSize * 0.35f;
  const float availableX = std::max(fontSize, rect.size.x - inset * 2.0f);
  const float availableY = std::max(fontSize, rect.size.y - inset * 2.0f);
  const float spanX = std::max(maxX - minX, std::numeric_limits<float>::epsilon());
  const float spanY = std::max(maxY - minY, std::numeric_limits<float>::epsilon());
  projection.scale = std::min(availableX / spanX, availableY / spanY) *
                     std::clamp(orbit.zoom, 0.4f, 4.0f);
  const float centreX = (minX + maxX) * 0.5f;
  const float centreY = (minY + maxY) * 0.5f;
  projection.centre.x -= centreX * projection.scale;
  projection.centre.y += centreY * projection.scale;
  return projection;
}

enum class PrimitiveKind { Fill, Line, Circle };

struct Primitive {
  PrimitiveKind kind = PrimitiveKind::Line;
  std::array<ImVec2, 4> points{};
  int count = 0;
  float depth = 0.0f;
  ImU32 colour = 0;
  ImU32 strokeColour = 0;
  float thickness = 0.0f;
  float radius = 0.0f;
};

void sortAndEmit(ImDrawList* drawList, std::vector<Primitive>& primitives) {
  std::stable_sort(primitives.begin(), primitives.end(),
                   [](const Primitive& left, const Primitive& right) {
                     return left.depth < right.depth;
                   });
  for (const Primitive& primitive : primitives) {
    if (primitive.kind == PrimitiveKind::Fill && primitive.count >= 3) {
      drawList->AddConvexPolyFilled(primitive.points.data(), primitive.count,
                                    primitive.colour);
      if (primitive.strokeColour != 0) {
        drawList->AddPolyline(primitive.points.data(), primitive.count,
                              primitive.strokeColour, ImDrawFlags_Closed,
                              primitive.thickness);
      }
    } else if (primitive.kind == PrimitiveKind::Line && primitive.count >= 2) {
      drawList->AddPolyline(primitive.points.data(), primitive.count,
                            primitive.colour, ImDrawFlags_None,
                            primitive.thickness);
    } else if (primitive.kind == PrimitiveKind::Circle) {
      drawList->AddCircleFilled(primitive.points[0], primitive.radius,
                                primitive.colour);
      if (primitive.strokeColour != 0) {
        drawList->AddCircle(primitive.points[0], primitive.radius,
                            primitive.strokeColour, 0, primitive.thickness);
      }
    }
  }
}

float representativeDepth(const Projection& projection, const Vec3* points, int count) {
  float depth = 0.0f;
  for (int index = 0; index < count; ++index) depth += projection.camera(points[index]).z;
  return count > 0 ? depth / static_cast<float>(count) : 0.0f;
}

float lambert(Vec3 normal, const Projection& projection) {
  const Vec3 cameraNormal = normalise(rotate(normal, projection.orbit));
  const Vec3 key = normalise(Vec3{-0.42f, 0.72f, 0.55f});
  return 0.28f + 0.72f * std::max(0.0f, dot(cameraNormal, key));
}

float distanceSquared(ImVec2 a, ImVec2 b) {
  const float dx = a.x - b.x;
  const float dy = a.y - b.y;
  return dx * dx + dy * dy;
}

float pointSegmentDistanceSquared(ImVec2 p, ImVec2 a, ImVec2 b) {
  const float dx = b.x - a.x;
  const float dy = b.y - a.y;
  const float denominator = dx * dx + dy * dy;
  if (denominator <= std::numeric_limits<float>::epsilon()) return distanceSquared(p, a);
  const float t = std::clamp(((p.x - a.x) * dx + (p.y - a.y) * dy) / denominator,
                             0.0f, 1.0f);
  return distanceSquared(p, ImVec2(a.x + dx * t, a.y + dy * t));
}

bool finite(double value) { return std::isfinite(value); }

bool nearlyEqual(double left, double right) {
  const double scale = std::max({1.0, std::abs(left), std::abs(right)});
  return std::abs(left - right) <= std::numeric_limits<double>::epsilon() * scale * 16.0;
}

float signedLog(double value) {
  if (!finite(value)) return 0.0f;
  return static_cast<float>(std::copysign(std::log1p(std::abs(value)), value));
}

ImVec4 colourAt(ImVec4 low, ImVec4 high, float t) {
  t = std::clamp(t, 0.0f, 1.0f);
  return {low.x + (high.x - low.x) * t,
          low.y + (high.y - low.y) * t,
          low.z + (high.z - low.z) * t,
          low.w + (high.w - low.w) * t};
}

bool rectsOverlap(ImVec2 aMin, ImVec2 aMax, ImVec2 bMin, ImVec2 bMax) {
  return aMin.x < bMax.x && aMax.x > bMin.x && aMin.y < bMax.y && aMax.y > bMin.y;
}

Vec3 spherical(float radius, float theta, float phi) {
  const float sinTheta = std::sin(theta);
  return {radius * sinTheta * std::cos(phi),
          radius * sinTheta * std::sin(phi),
          radius * std::cos(theta)};
}

long double binomial(int n, int k) {
  if (k < 0 || k > n) return 0.0L;
  k = std::min(k, n - k);
  long double result = 1.0L;
  for (int index = 1; index <= k; ++index) {
    result *= static_cast<long double>(n - k + index) / static_cast<long double>(index);
  }
  return result;
}

long double associatedLaguerre(int degree, int alpha, long double x) {
  long double sum = 0.0L;
  long double power = 1.0L;
  long double factorial = 1.0L;
  for (int index = 0; index <= degree; ++index) {
    if (index > 0) {
      power *= x;
      factorial *= static_cast<long double>(index);
    }
    const long double sign = (index & 1) ? -1.0L : 1.0L;
    sum += sign * binomial(degree + alpha, degree - index) * power / factorial;
  }
  return sum;
}

double radialHydrogenic(int n, int l, double radius) {
  const long double rho = 2.0L * static_cast<long double>(radius) /
                          static_cast<long double>(n);
  const int degree = n - l - 1;
  const long double envelope = std::exp(-rho * 0.5L) *
                               std::pow(std::max(0.0L, rho), l);
  return static_cast<double>(envelope * associatedLaguerre(degree, 2 * l + 1, rho));
}

double realSphericalHarmonic(int l, int m, double theta, double phi) {
  const double c = std::cos(theta);
  const double s = std::sin(theta);
  switch (l) {
    case 0:
      return std::sqrt(1.0 / (4.0 * kPi));
    case 1:
      if (m == 0) return std::sqrt(3.0 / (4.0 * kPi)) * c;
      if (m == 1) return std::sqrt(3.0 / (4.0 * kPi)) * s * std::cos(phi);
      if (m == -1) return std::sqrt(3.0 / (4.0 * kPi)) * s * std::sin(phi);
      break;
    case 2:
      if (m == 0) return std::sqrt(5.0 / (16.0 * kPi)) * (3.0 * c * c - 1.0);
      if (m == 1) return std::sqrt(15.0 / (4.0 * kPi)) * s * c * std::cos(phi);
      if (m == -1) return std::sqrt(15.0 / (4.0 * kPi)) * s * c * std::sin(phi);
      if (m == 2) return std::sqrt(15.0 / (16.0 * kPi)) * s * s * std::cos(2.0 * phi);
      if (m == -2) return std::sqrt(15.0 / (16.0 * kPi)) * s * s * std::sin(2.0 * phi);
      break;
    case 3:
      if (m == 0) return std::sqrt(7.0 / (16.0 * kPi)) * c * (5.0 * c * c - 3.0);
      if (m == 1) return std::sqrt(21.0 / (32.0 * kPi)) * s *
                         (5.0 * c * c - 1.0) * std::cos(phi);
      if (m == -1) return std::sqrt(21.0 / (32.0 * kPi)) * s *
                          (5.0 * c * c - 1.0) * std::sin(phi);
      if (m == 2) return std::sqrt(105.0 / (16.0 * kPi)) * s * s * c *
                         std::cos(2.0 * phi);
      if (m == -2) return std::sqrt(105.0 / (16.0 * kPi)) * s * s * c *
                          std::sin(2.0 * phi);
      if (m == 3) return std::sqrt(35.0 / (32.0 * kPi)) * s * s * s *
                         std::cos(3.0 * phi);
      if (m == -3) return std::sqrt(35.0 / (32.0 * kPi)) * s * s * s *
                          std::sin(3.0 * phi);
      break;
    default:
      break;
  }
  return 0.0;
}

std::vector<double> radialNodes(int n, int l, double maximumRadius) {
  std::vector<double> nodes;
  const int expected = n - l - 1;
  if (expected <= 0) return nodes;
  const int samples = std::max(512, expected * 192);
  auto nodePolynomial = [&](double radius) {
    const long double rho = 2.0L * static_cast<long double>(radius) /
                            static_cast<long double>(n);
    return static_cast<double>(associatedLaguerre(expected, 2 * l + 1, rho));
  };
  double previousRadius = 0.0;
  double previous = nodePolynomial(previousRadius);
  for (int index = 1; index <= samples && static_cast<int>(nodes.size()) < expected; ++index) {
    const double radius = maximumRadius * static_cast<double>(index) /
                          static_cast<double>(samples);
    const double current = nodePolynomial(radius);
    if (current == 0.0) {
      nodes.push_back(radius);
    } else if ((previous < 0.0 && current > 0.0) ||
               (previous > 0.0 && current < 0.0)) {
      double low = previousRadius;
      double high = radius;
      for (int iteration = 0; iteration < 36; ++iteration) {
        const double middle = (low + high) * 0.5;
        const double value = nodePolynomial(middle);
        if ((previous < 0.0 && value < 0.0) || (previous > 0.0 && value > 0.0)) {
          low = middle;
        } else {
          high = middle;
        }
      }
      nodes.push_back((low + high) * 0.5);
    }
    previousRadius = radius;
    previous = current;
  }
  return nodes;
}

std::vector<double> isoCrossings(int n, int l, double angularMagnitude,
                                 double target, double maximumRadius) {
  std::vector<double> crossings;
  if (angularMagnitude <= std::numeric_limits<double>::epsilon()) return crossings;
  const int radialNodeCount = n - l - 1;
  const int samples = std::max(320, (radialNodeCount + 1) * 128);
  auto field = [&](double radius) {
    return std::abs(radialHydrogenic(n, l, radius)) * angularMagnitude - target;
  };
  double previousRadius = 0.0;
  double previous = field(previousRadius);
  for (int index = 1; index <= samples; ++index) {
    const double fraction = static_cast<double>(index) / static_cast<double>(samples);
    const double radius = maximumRadius * fraction * fraction;
    const double current = field(radius);
    if ((previous < 0.0 && current >= 0.0) || (previous >= 0.0 && current < 0.0)) {
      double low = previousRadius;
      double high = radius;
      double lowValue = previous;
      for (int iteration = 0; iteration < 32; ++iteration) {
        const double middle = (low + high) * 0.5;
        const double middleValue = field(middle);
        if ((lowValue < 0.0) == (middleValue < 0.0)) {
          low = middle;
          lowValue = middleValue;
        } else {
          high = middle;
        }
      }
      crossings.push_back((low + high) * 0.5);
    }
    previousRadius = radius;
    previous = current;
  }
  return crossings;
}

void beginChartTooltip() {
  ImGui::BeginTooltip();
  const float base = std::max(ImGui::GetFontSize(), std::numeric_limits<float>::epsilon());
  ImGui::SetWindowFontScale(layout::labelFont(base * 0.82f) / base);
}

void tooltipValue(const char* label, double value, const char* unit = nullptr) {
  if (label && label[0] != '\0') {
    if (unit && unit[0] != '\0') ImGui::Text("%s: %.5g %s", label, value, unit);
    else ImGui::Text("%s: %.5g", label, value);
  }
}

}  // namespace

void handleOrbitInput(Orbit& orbit, bool hovered, bool active) {
  ImGuiIO& io = ImGui::GetIO();
  const bool dragging = active && ImGui::IsMouseDown(ImGuiMouseButton_Left);
  if (dragging) {
    const float sensitivity = 180.0f / std::max(ImGui::GetFontSize() * 18.0f,
                                                std::numeric_limits<float>::epsilon());
    orbit.yawDeg += io.MouseDelta.x * sensitivity;
    orbit.pitchDeg -= io.MouseDelta.y * sensitivity;
    orbit.pitchDeg = std::clamp(orbit.pitchDeg, -88.9f, 88.9f);
  }
  if (hovered && io.MouseWheel != 0.0f) {
    orbit.zoom *= std::pow(1.12f, io.MouseWheel);
    orbit.zoom = std::clamp(orbit.zoom, 0.4f, 4.0f);
  }
  if (orbit.spin && !dragging && std::isfinite(io.DeltaTime)) {
    orbit.yawDeg += orbit.spinRateDegS * io.DeltaTime;
  }
  orbit.yawDeg = std::fmod(orbit.yawDeg, 360.0f);
  if (orbit.yawDeg < 0.0f) orbit.yawDeg += 360.0f;
}

int surface(const char* id, const double* values, int columns, int rows, ImVec2 size,
            Orbit& orbit, const SurfaceStyle& style) {
  const DrawRect rect = reserveRect(id, size);
  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawFrame(drawList, rect);
  const bool hovered = ImGui::IsItemHovered();
  const bool active = ImGui::IsItemActive();
  handleOrbitInput(orbit, hovered, active);
  if (!hasArea(rect) || !values || columns <= 0 || rows <= 0) return -1;

  const int count = columns * rows;
  std::vector<float> heights(static_cast<std::size_t>(count));
  double minimum = std::numeric_limits<double>::max();
  double maximum = -std::numeric_limits<double>::max();
  int peakIndex = -1;
  for (int index = 0; index < count; ++index) {
    const double value = finite(values[index]) ? values[index] : 0.0;
    const float mapped = style.logHeight ? signedLog(value) : static_cast<float>(value);
    heights[static_cast<std::size_t>(index)] = mapped;
    if (mapped < minimum) minimum = mapped;
    if (mapped > maximum) {
      maximum = mapped;
      peakIndex = index;
    }
  }
  if (nearlyEqual(minimum, maximum)) return -1;
  const double span = std::max(maximum - minimum, std::numeric_limits<double>::epsilon());
  auto vertex = [&](int column, int row, bool base = false) {
    const int index = row * columns + column;
    const float x = columns > 1 ? -1.0f + 2.0f * static_cast<float>(column) /
                                               static_cast<float>(columns - 1) : 0.0f;
    const float y = rows > 1 ? -1.0f + 2.0f * static_cast<float>(row) /
                                          static_cast<float>(rows - 1) : 0.0f;
    const float z = base ? -0.72f
                         : -0.72f + 1.44f * static_cast<float>(
                             (heights[static_cast<std::size_t>(index)] - minimum) / span);
    return Vec3{x, y, z};
  };

  std::vector<Vec3> fitPoints;
  fitPoints.reserve(static_cast<std::size_t>(count + 8));
  for (int row = 0; row < rows; ++row) {
    for (int column = 0; column < columns; ++column) fitPoints.push_back(vertex(column, row));
  }
  for (float x : {-1.0f, 1.0f}) {
    for (float y : {-1.0f, 1.0f}) {
      fitPoints.push_back({x, y, -0.72f});
      fitPoints.push_back({x, y, 0.72f});
    }
  }
  const Projection projection = fittedProjection(rect, orbit, fitPoints);

  std::vector<Primitive> primitives;
  primitives.reserve(static_cast<std::size_t>(std::max(1, (rows - 1) * (columns - 1) * 2)));
  const float hairline = uiStyle::metrics().hairline;
  for (int row = 0; row + 1 < rows; ++row) {
    for (int column = 0; column + 1 < columns; ++column) {
      const Vec3 corners[] = {vertex(column, row), vertex(column + 1, row),
                              vertex(column + 1, row + 1), vertex(column, row + 1)};
      const int indices[] = {row * columns + column, row * columns + column + 1,
                             (row + 1) * columns + column + 1,
                             (row + 1) * columns + column};
      float normalizedHeight = 0.0f;
      for (int corner = 0; corner < 4; ++corner) {
        normalizedHeight += static_cast<float>(
            (heights[static_cast<std::size_t>(indices[corner])] - minimum) / span);
      }
      normalizedHeight *= 0.25f;
      const ImVec4 ramp = colourAt(style.low, style.high, normalizedHeight);
      const Vec3 normal = cross(corners[1] - corners[0], corners[3] - corners[0]);
      Primitive quad;
      quad.kind = PrimitiveKind::Fill;
      quad.count = 4;
      quad.depth = representativeDepth(projection, corners, 4);
      for (int corner = 0; corner < 4; ++corner) quad.points[corner] = projection.project(corners[corner]);
      quad.colour = uiStyle::mix(uiStyle::col::BgSurface, ramp, lambert(normal, projection), 0.94f);
      if (style.showGrid) {
        quad.strokeColour = uiStyle::u32(uiStyle::col::GridLine, 0.92f);
        quad.thickness = hairline;
      }
      primitives.push_back(quad);

      if (style.showFloor) {
        Vec3 floorCorners[] = {vertex(column, row, true), vertex(column + 1, row, true),
                               vertex(column + 1, row + 1, true), vertex(column, row + 1, true)};
        Primitive floor;
        floor.kind = PrimitiveKind::Fill;
        floor.count = 4;
        floor.depth = representativeDepth(projection, floorCorners, 4) - 0.02f;
        for (int corner = 0; corner < 4; ++corner) floor.points[corner] = projection.project(floorCorners[corner]);
        floor.colour = uiStyle::u32(ramp, 0.30f);
        floor.strokeColour = uiStyle::u32(uiStyle::col::GridLine, 0.42f);
        floor.thickness = hairline;
        primitives.push_back(floor);
      }
    }
  }
  Vec3 axisEndpoints[3][2]{};
  for (int dimension = 0; dimension < 3; ++dimension) {
    float nearestDepth = -std::numeric_limits<float>::max();
    for (int first = 0; first < 2; ++first) {
      for (int second = 0; second < 2; ++second) {
        const float a = first ? 1.0f : -1.0f;
        const float b = second ? 0.72f : -0.72f;
        Vec3 endpoints[2];
        if (dimension == 0) {
          endpoints[0] = {-1.0f, a, b};
          endpoints[1] = {1.0f, a, b};
        } else if (dimension == 1) {
          endpoints[0] = {a, -1.0f, b};
          endpoints[1] = {a, 1.0f, b};
        } else {
          const float y = second ? 1.0f : -1.0f;
          endpoints[0] = {a, y, -0.72f};
          endpoints[1] = {a, y, 0.72f};
        }
        const float depth = representativeDepth(projection, endpoints, 2);
        if (depth > nearestDepth) {
          nearestDepth = depth;
          axisEndpoints[dimension][0] = endpoints[0];
          axisEndpoints[dimension][1] = endpoints[1];
        }
      }
    }
    Primitive axisLine;
    axisLine.kind = PrimitiveKind::Line;
    axisLine.count = 2;
    axisLine.points[0] = projection.project(axisEndpoints[dimension][0]);
    axisLine.points[1] = projection.project(axisEndpoints[dimension][1]);
    axisLine.depth = nearestDepth;
    axisLine.colour = uiStyle::u32(uiStyle::col::DataDim);
    axisLine.thickness = hairline;
    primitives.push_back(axisLine);
  }
  drawList->PushClipRect(rect.min, rect.max, true);
  sortAndEmit(drawList, primitives);

  int hoveredIndex = -1;
  if (hovered) {
    const ImVec2 mouse = ImGui::GetIO().MousePos;
    float bestDistance = std::numeric_limits<float>::max();
    for (int row = 0; row < rows; ++row) {
      for (int column = 0; column < columns; ++column) {
        const int index = row * columns + column;
        const float distance = distanceSquared(mouse, projection.project(vertex(column, row)));
        if (distance < bestDistance) {
          bestDistance = distance;
          hoveredIndex = index;
        }
      }
    }
    const float reach = ImGui::GetFontSize() * 0.82f;
    if (bestDistance > reach * reach) hoveredIndex = -1;
  }

  if (hoveredIndex >= 0) {
    const int column = hoveredIndex % columns;
    const int row = hoveredIndex / columns;
    const ImVec2 point = projection.project(vertex(column, row));
    const float radius = ImGui::GetFontSize() * 0.22f;
    drawList->AddCircleFilled(point, radius, uiStyle::u32(uiStyle::col::Accent));
    drawList->AddCircle(point, radius * 1.75f, uiStyle::u32(uiStyle::col::Accent, 0.36f),
                        0, hairline);
    beginChartTooltip();
    ImGui::Text("u: %.5g", columns > 1 ? static_cast<double>(column) / (columns - 1) : 0.0);
    ImGui::Text("v: %.5g", rows > 1 ? static_cast<double>(row) / (rows - 1) : 0.0);
    ImGui::Text("w: %.5g", values[hoveredIndex]);
    ImGui::EndTooltip();
  }

  if (style.markPeak && peakIndex >= 0) {
    const int column = peakIndex % columns;
    const int row = peakIndex / columns;
    const ImVec2 peak = projection.project(vertex(column, row));
    const float radius = ImGui::GetFontSize() * 0.25f;
    const ImVec2 leaderEnd(peak.x + ImGui::GetFontSize() * 0.72f,
                           peak.y - ImGui::GetFontSize() * 0.72f);
    drawList->AddLine(peak, leaderEnd, uiStyle::u32(style.peak), hairline);
    drawList->AddCircleFilled(peak, radius, uiStyle::u32(style.peak));
  }

  const ImVec2 uEnd = projection.project(axisEndpoints[0][1]);
  const ImVec2 vEnd = projection.project(axisEndpoints[1][1]);
  const ImVec2 wEnd = projection.project(axisEndpoints[2][1]);
  addLabel(drawList, uEnd, uiStyle::col::TextDim, style.uLabel);
  addLabel(drawList, vEnd, uiStyle::col::TextDim, style.vLabel);
  addLabel(drawList, wEnd, uiStyle::col::TextDim, style.wLabel);
  drawList->PopClipRect();
  return hoveredIndex;
}

int ternary(const char* id, const double* points, const double* values, int count,
            ImVec2 size, const TernaryStyle& style, double* outA, double* outB,
            double* outC) {
  const DrawRect rect = reserveRect(id, size);
  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawFrame(drawList, rect);
  if (!hasArea(rect) || !points || !values || count <= 0) return -1;

  const float fontSize = ImGui::GetFontSize();
  const auto& metrics = uiStyle::metrics();
  const float labelSize = chartLabelSize();
  const float inset = fontSize * 1.35f + metrics.gap * 0.35f;
  const float availableWidth = std::max(0.0f, rect.size.x - inset * 2.0f);
  const float availableHeight = std::max(0.0f, rect.size.y - inset * 2.0f);
  const float side = std::min(availableWidth, availableHeight * 2.0f / std::sqrt(3.0f));
  const float triangleHeight = side * std::sqrt(3.0f) * 0.5f;
  const ImVec2 a((rect.min.x + rect.max.x) * 0.5f, rect.min.y + inset);
  const ImVec2 b(a.x - side * 0.5f, a.y + triangleHeight);
  const ImVec2 c(a.x + side * 0.5f, a.y + triangleHeight);
  if (side <= std::numeric_limits<float>::epsilon()) return -1;

  auto baryPoint = [&](double weightA, double weightB, double weightC) {
    const double sum = weightA + weightB + weightC;
    if (!finite(sum) || std::abs(sum) <= std::numeric_limits<double>::epsilon()) {
      return ImVec2((a.x + b.x + c.x) / 3.0f, (a.y + b.y + c.y) / 3.0f);
    }
    return ImVec2(static_cast<float>((weightA * a.x + weightB * b.x + weightC * c.x) / sum),
                  static_cast<float>((weightA * a.y + weightB * b.y + weightC * c.y) / sum));
  };

  double minimum = std::numeric_limits<double>::max();
  double maximum = -std::numeric_limits<double>::max();
  int validSampleCount = 0;
  for (int index = 0; index < count; ++index) {
    const double sum = points[index * 3] + points[index * 3 + 1] + points[index * 3 + 2];
    if (!finite(points[index * 3]) || !finite(points[index * 3 + 1]) ||
        !finite(points[index * 3 + 2]) || !finite(sum) ||
        std::abs(sum) <= std::numeric_limits<double>::epsilon() ||
        !finite(values[index])) continue;
    ++validSampleCount;
    minimum = std::min(minimum, values[index]);
    maximum = std::max(maximum, values[index]);
  }
  if (validSampleCount == 0 || nearlyEqual(minimum, maximum)) return -1;
  const double span = std::max(maximum - minimum, std::numeric_limits<double>::epsilon());
  auto interpolate = [&](double wa, double wb, double wc) {
    if (!points || !values || count <= 0) return 0.5;
    double weighted = 0.0;
    double weights = 0.0;
    for (int index = 0; index < count; ++index) {
      const double sum = points[index * 3] + points[index * 3 + 1] + points[index * 3 + 2];
      if (!finite(sum) || std::abs(sum) <= std::numeric_limits<double>::epsilon() ||
          !finite(values[index])) continue;
      const double pa = points[index * 3] / sum;
      const double pb = points[index * 3 + 1] / sum;
      const double pc = points[index * 3 + 2] / sum;
      const double distance = (pa - wa) * (pa - wa) + (pb - wb) * (pb - wb) +
                              (pc - wc) * (pc - wc);
      const double weight = 1.0 / std::max(distance, 1.0e-6);
      weighted += values[index] * weight;
      weights += weight;
    }
    return weights > 0.0 ? (weighted / weights - minimum) / span : 0.5;
  };

  drawList->PushClipRect(rect.min, rect.max, true);
  const int divisions = 18;
  for (int row = 0; row < divisions; ++row) {
    for (int column = 0; column < divisions - row; ++column) {
      const double inv = 1.0 / divisions;
      const double wa0 = 1.0 - (row + column) * inv;
      const double wb0 = column * inv;
      const double wc0 = row * inv;
      const std::array<std::array<double, 3>, 3> first{{
          {wa0, wb0, wc0}, {wa0 - inv, wb0 + inv, wc0},
          {wa0 - inv, wb0, wc0 + inv}}};
      ImVec2 triangle[3];
      double shade = 0.0;
      for (int vertexIndex = 0; vertexIndex < 3; ++vertexIndex) {
        const auto& weights = first[vertexIndex];
        triangle[vertexIndex] = baryPoint(weights[0], weights[1], weights[2]);
        shade += interpolate(weights[0], weights[1], weights[2]);
      }
      shade /= 3.0;
      drawList->AddConvexPolyFilled(triangle, 3,
          uiStyle::u32(colourAt(style.low, style.high,
                                static_cast<float>(shade)), 0.86f));

      if (column + row + 1 < divisions) {
        const std::array<std::array<double, 3>, 3> second{{
            {wa0 - inv, wb0 + inv, wc0},
            {wa0 - 2.0 * inv, wb0 + inv, wc0 + inv},
            {wa0 - inv, wb0, wc0 + inv}}};
        shade = 0.0;
        for (int vertexIndex = 0; vertexIndex < 3; ++vertexIndex) {
          const auto& weights = second[vertexIndex];
          triangle[vertexIndex] = baryPoint(weights[0], weights[1], weights[2]);
          shade += interpolate(weights[0], weights[1], weights[2]);
        }
        shade /= 3.0;
        drawList->AddConvexPolyFilled(triangle, 3,
            uiStyle::u32(colourAt(style.low, style.high,
                                  static_cast<float>(shade)), 0.86f));
      }
    }
  }

  if (style.isolines && values && count > 0) {
    const int contourCount = 5;
    const int grid = 24;
    auto drawContourTriangle = [&](const std::array<std::array<double, 3>, 3>& weights,
                                   double iso) {
      std::array<double, 3> samples{};
      for (int vertexIndex = 0; vertexIndex < 3; ++vertexIndex) {
        samples[vertexIndex] = interpolate(weights[vertexIndex][0],
                                            weights[vertexIndex][1],
                                            weights[vertexIndex][2]);
      }
      ImVec2 intersections[2];
      int intersectionCount = 0;
      for (int edge = 0; edge < 3 && intersectionCount < 2; ++edge) {
        const int next = (edge + 1) % 3;
        if ((samples[edge] < iso) == (samples[next] < iso) ||
            std::abs(samples[next] - samples[edge]) <=
                std::numeric_limits<double>::epsilon()) continue;
        const double t = (iso - samples[edge]) / (samples[next] - samples[edge]);
        const ImVec2 p0 = baryPoint(weights[edge][0], weights[edge][1], weights[edge][2]);
        const ImVec2 p1 = baryPoint(weights[next][0], weights[next][1], weights[next][2]);
        intersections[intersectionCount++] = ImVec2(
            p0.x + static_cast<float>(t) * (p1.x - p0.x),
            p0.y + static_cast<float>(t) * (p1.y - p0.y));
      }
      if (intersectionCount == 2) {
        drawList->AddLine(intersections[0], intersections[1],
                          uiStyle::u32(uiStyle::col::GridLine, 0.88f),
                          metrics.hairline);
      }
    };
    for (int contourIndex = 1; contourIndex < contourCount; ++contourIndex) {
      const double iso = static_cast<double>(contourIndex) / contourCount;
      for (int row = 0; row < grid; ++row) {
        for (int column = 0; column < grid - row; ++column) {
          const double inv = 1.0 / grid;
          const double wa0 = 1.0 - (row + column) * inv;
          const double wb0 = column * inv;
          const double wc0 = row * inv;
          drawContourTriangle({{{wa0, wb0, wc0},
                                {wa0 - inv, wb0 + inv, wc0},
                                {wa0 - inv, wb0, wc0 + inv}}}, iso);
          if (row + column + 1 < grid) {
            drawContourTriangle({{{wa0 - inv, wb0 + inv, wc0},
                                  {wa0 - 2.0 * inv, wb0 + inv, wc0 + inv},
                                  {wa0 - inv, wb0, wc0 + inv}}}, iso);
          }
        }
      }
    }
  }

  drawList->AddTriangle(a, b, c, uiStyle::u32(uiStyle::col::DataDim),
                        metrics.hairline * 1.5f);
  ImVec2 aText = textSize(style.aLabel, labelSize);
  ImVec2 bText = textSize(style.bLabel, labelSize);
  addLabel(drawList, ImVec2(a.x - aText.x * 0.5f, a.y - aText.y - metrics.gap * 0.25f),
           uiStyle::col::TextDim, style.aLabel);
  addLabel(drawList, ImVec2(b.x - bText.x - metrics.gap * 0.2f, b.y),
           uiStyle::col::TextDim, style.bLabel);
  addLabel(drawList, ImVec2(c.x + metrics.gap * 0.2f, c.y),
           uiStyle::col::TextDim, style.cLabel);

  if (style.hasMarker) {
    const ImVec2 marker = baryPoint(style.markerA, style.markerB,
                                    style.markerC);
    auto foot = [](ImVec2 p, ImVec2 start, ImVec2 end) {
      const float dx = end.x - start.x;
      const float dy = end.y - start.y;
      const float denominator = dx * dx + dy * dy;
      const float t = denominator > 0.0f
                          ? std::clamp(((p.x - start.x) * dx + (p.y - start.y) * dy) /
                                           denominator, 0.0f, 1.0f)
                          : 0.0f;
      return ImVec2(start.x + dx * t, start.y + dy * t);
    };
    drawList->AddLine(marker, foot(marker, a, b), uiStyle::u32(uiStyle::col::Accent, 0.42f),
                      metrics.hairline);
    drawList->AddLine(marker, foot(marker, b, c), uiStyle::u32(uiStyle::col::Accent, 0.42f),
                      metrics.hairline);
    drawList->AddLine(marker, foot(marker, c, a), uiStyle::u32(uiStyle::col::Accent, 0.42f),
                      metrics.hairline);
    drawList->AddCircleFilled(marker, fontSize * 0.24f, uiStyle::u32(uiStyle::col::Accent));
  }
  drawList->PopClipRect();

  const bool hovered = ImGui::IsItemHovered();
  const ImVec2 mouse = ImGui::GetIO().MousePos;
  const float denominator = (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y);
  const double mouseA = ((b.y - c.y) * (mouse.x - c.x) +
                         (c.x - b.x) * (mouse.y - c.y)) / denominator;
  const double mouseB = ((c.y - a.y) * (mouse.x - c.x) +
                         (a.x - c.x) * (mouse.y - c.y)) / denominator;
  const double mouseC = 1.0 - mouseA - mouseB;
  const bool inside = mouseA >= 0.0 && mouseB >= 0.0 && mouseC >= 0.0;
  if (hovered && inside && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
    if (outA) *outA = mouseA;
    if (outB) *outB = mouseB;
    if (outC) *outC = mouseC;
    return -2;
  }

  int hoveredIndex = -1;
  if (hovered && points && count > 0) {
    float best = std::numeric_limits<float>::max();
    for (int index = 0; index < count; ++index) {
      const double sum = points[index * 3] + points[index * 3 + 1] +
                         points[index * 3 + 2];
      if (!finite(sum) || std::abs(sum) <= std::numeric_limits<double>::epsilon() ||
          !finite(values[index])) continue;
      const ImVec2 point = baryPoint(points[index * 3], points[index * 3 + 1],
                                     points[index * 3 + 2]);
      const float distance = distanceSquared(mouse, point);
      if (distance < best) {
        best = distance;
        hoveredIndex = index;
      }
    }
    const float reach = fontSize * 0.72f;
    if (best > reach * reach) hoveredIndex = -1;
  }
  if (hoveredIndex >= 0) {
    const double sum = points[hoveredIndex * 3] + points[hoveredIndex * 3 + 1] +
                       points[hoveredIndex * 3 + 2];
    beginChartTooltip();
    if (std::abs(sum) > std::numeric_limits<double>::epsilon()) {
      ImGui::Text("%s: %.4g", style.aLabel ? style.aLabel : "A",
                  points[hoveredIndex * 3] / sum);
      ImGui::Text("%s: %.4g", style.bLabel ? style.bLabel : "B",
                  points[hoveredIndex * 3 + 1] / sum);
      ImGui::Text("%s: %.4g", style.cLabel ? style.cLabel : "C",
                  points[hoveredIndex * 3 + 2] / sum);
    }
    if (values) ImGui::Text("value: %.5g", values[hoveredIndex]);
    ImGui::EndTooltip();
  }
  return hoveredIndex;
}

int cloud(const char* id, const CloudPoint* points, int count, ImVec2 size, Orbit& orbit,
          const CloudStyle& style) {
  const DrawRect rect = reserveRect(id, size);
  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawFrame(drawList, rect);
  const bool hovered = ImGui::IsItemHovered();
  handleOrbitInput(orbit, hovered, ImGui::IsItemActive());
  if (!hasArea(rect) || !points || count <= 0) return -1;

  double xLow = std::numeric_limits<double>::max();
  double xHigh = -std::numeric_limits<double>::max();
  double yLow = xLow;
  double yHigh = xHigh;
  double zLow = xLow;
  double zHigh = xHigh;
  int validPointCount = 0;
  for (int index = 0; index < count; ++index) {
    if (!finite(points[index].x) || !finite(points[index].y) || !finite(points[index].z)) continue;
    ++validPointCount;
    xLow = std::min(xLow, points[index].x);
    xHigh = std::max(xHigh, points[index].x);
    yLow = std::min(yLow, points[index].y);
    yHigh = std::max(yHigh, points[index].y);
    zLow = std::min(zLow, points[index].z);
    zHigh = std::max(zHigh, points[index].z);
  }
  const bool validSphere = style.hasSphere && finite(style.sphereX) &&
      finite(style.sphereY) && finite(style.sphereZ) &&
      finite(style.sphereRadius) && style.sphereRadius > 0.0;
  if (validSphere) {
    xLow = std::min(xLow, style.sphereX - style.sphereRadius);
    xHigh = std::max(xHigh, style.sphereX + style.sphereRadius);
    yLow = std::min(yLow, style.sphereY - style.sphereRadius);
    yHigh = std::max(yHigh, style.sphereY + style.sphereRadius);
    zLow = std::min(zLow, style.sphereZ - style.sphereRadius);
    zHigh = std::max(zHigh, style.sphereZ + style.sphereRadius);
  }
  if (validPointCount == 0 && !validSphere) return -1;
  if (!validSphere && nearlyEqual(xLow, xHigh) && nearlyEqual(yLow, yHigh) &&
      nearlyEqual(zLow, zHigh)) return -1;
  const charts::Axis xAxis = charts::niceAxis(xLow, xHigh);
  const charts::Axis yAxis = charts::niceAxis(yLow, yHigh);
  const charts::Axis zAxis = charts::niceAxis(zLow, zHigh);
  const double xCentre = (xAxis.min + xAxis.max) * 0.5;
  const double yCentre = (yAxis.min + yAxis.max) * 0.5;
  const double zCentre = (zAxis.min + zAxis.max) * 0.5;
  const double commonSpan = std::max({xAxis.max - xAxis.min, yAxis.max - yAxis.min,
                                      zAxis.max - zAxis.min,
                                      std::numeric_limits<double>::epsilon()});
  auto model = [&](double x, double y, double z) {
    return Vec3{2.0f * static_cast<float>((x - xCentre) / commonSpan),
                2.0f * static_cast<float>((y - yCentre) / commonSpan),
                2.0f * static_cast<float>((z - zCentre) / commonSpan)};
  };
  const std::array<Vec3, 8> boxCorners{{
      model(xAxis.min, yAxis.min, zAxis.min), model(xAxis.max, yAxis.min, zAxis.min),
      model(xAxis.max, yAxis.max, zAxis.min), model(xAxis.min, yAxis.max, zAxis.min),
      model(xAxis.min, yAxis.min, zAxis.max), model(xAxis.max, yAxis.min, zAxis.max),
      model(xAxis.max, yAxis.max, zAxis.max), model(xAxis.min, yAxis.max, zAxis.max)}};

  std::vector<Vec3> fitPoints;
  fitPoints.reserve(static_cast<std::size_t>(validPointCount + 8));
  for (int index = 0; index < count; ++index) {
    if (!finite(points[index].x) || !finite(points[index].y) || !finite(points[index].z)) continue;
    fitPoints.push_back(model(points[index].x, points[index].y, points[index].z));
  }
  fitPoints.insert(fitPoints.end(), boxCorners.begin(), boxCorners.end());
  const Projection projection = fittedProjection(rect, orbit, fitPoints);
  const auto& metrics = uiStyle::metrics();
  const float fontSize = ImGui::GetFontSize();
  std::vector<Primitive> primitives;

  if (style.showAxes) {
    const auto& corners = boxCorners;
    const int edges[][2] = {{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},
                            {6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
    for (const auto& edge : edges) {
      Primitive line;
      line.kind = PrimitiveKind::Line;
      line.count = 2;
      const Vec3 endpoints[] = {corners[edge[0]], corners[edge[1]]};
      line.points[0] = projection.project(endpoints[0]);
      line.points[1] = projection.project(endpoints[1]);
      line.depth = representativeDepth(projection, endpoints, 2);
      line.colour = uiStyle::u32(uiStyle::col::GridLine, 0.90f);
      line.thickness = metrics.hairline;
      primitives.push_back(line);
    }
  }

  if (validSphere) {
    const int segments = 36;
    const int ringCount = 5;
    auto spherePoint = [&](Vec3 unit) {
      return model(style.sphereX + unit.x * style.sphereRadius,
                   style.sphereY + unit.y * style.sphereRadius,
                   style.sphereZ + unit.z * style.sphereRadius);
    };
    for (int ring = 1; ring < ringCount; ++ring) {
      const float theta = kPi * static_cast<float>(ring) / ringCount;
      for (int segment = 0; segment < segments; ++segment) {
        const float p0 = 2.0f * kPi * segment / segments;
        const float p1 = 2.0f * kPi * (segment + 1) / segments;
        const Vec3 endpoints[] = {
            spherePoint(spherical(1.0f, theta, p0)),
            spherePoint(spherical(1.0f, theta, p1))};
        Primitive line;
        line.kind = PrimitiveKind::Line;
        line.count = 2;
        line.points[0] = projection.project(endpoints[0]);
        line.points[1] = projection.project(endpoints[1]);
        line.depth = representativeDepth(projection, endpoints, 2);
        line.colour = uiStyle::u32(style.sphereColour, 0.30f);
        line.thickness = metrics.hairline;
        primitives.push_back(line);
      }
    }
    for (int longitude = 0; longitude < ringCount * 2; ++longitude) {
      const float phi = kPi * longitude / ringCount;
      for (int segment = 0; segment < segments / 2; ++segment) {
        const float t0 = kPi * segment / (segments / 2);
        const float t1 = kPi * (segment + 1) / (segments / 2);
        const Vec3 endpoints[] = {spherePoint(spherical(1.0f, t0, phi)),
                                  spherePoint(spherical(1.0f, t1, phi))};
        Primitive line;
        line.kind = PrimitiveKind::Line;
        line.count = 2;
        line.points[0] = projection.project(endpoints[0]);
        line.points[1] = projection.project(endpoints[1]);
        line.depth = representativeDepth(projection, endpoints, 2);
        line.colour = uiStyle::u32(style.sphereColour, 0.24f);
        line.thickness = metrics.hairline;
        primitives.push_back(line);
      }
    }
  }

  struct ScreenPoint {
    int index = -1;
    ImVec2 position{};
    float depth = 0.0f;
    float radius = 0.0f;
  };
  std::vector<ScreenPoint> screenPoints;
  screenPoints.reserve(static_cast<std::size_t>(count));
  for (int index = 0; index < count; ++index) {
    if (!finite(points[index].x) || !finite(points[index].y) || !finite(points[index].z)) continue;
    const Vec3 point = model(points[index].x, points[index].y, points[index].z);
    const Vec3 camera = projection.camera(point);
    const float perspective = projection.perspectiveDistance /
                              std::max(projection.perspectiveDistance - camera.z,
                                       projection.perspectiveDistance * 0.16f);
    const float radius = fontSize * 0.18f * std::max(0.25f, points[index].radius) *
                         std::clamp(perspective, 0.55f, 1.65f);
    screenPoints.push_back({index, projection.projectCamera(camera), camera.z, radius});
    Primitive marker;
    marker.kind = PrimitiveKind::Circle;
    marker.points[0] = screenPoints.back().position;
    marker.depth = camera.z;
    marker.radius = radius;
    marker.colour = uiStyle::u32(points[index].highlighted ? uiStyle::col::Accent
                                                           : points[index].colour,
                                 points[index].highlighted ? 1.0f : 0.88f);
    marker.strokeColour = uiStyle::u32(points[index].highlighted ? uiStyle::col::AccentHover
                                                                 : uiStyle::col::BgDeep,
                                       0.86f);
    marker.thickness = metrics.hairline;
    primitives.push_back(marker);
  }

  drawList->PushClipRect(rect.min, rect.max, true);
  sortAndEmit(drawList, primitives);

  if (style.showAxes) {
    auto drawTicks = [&](const charts::Axis& axis, int dimension, const char* label) {
      const int tickCount = std::max(2, axis.ticks);
      for (int tick = 0; tick < tickCount; ++tick) {
        const double value = axis.min + axis.step * tick;
        if (value > axis.max + axis.step * 0.25) break;
        Vec3 point = boxCorners[0];
        if (dimension == 0) point = model(value, yAxis.min, zAxis.min);
        if (dimension == 1) point = model(xAxis.min, value, zAxis.min);
        if (dimension == 2) point = model(xAxis.min, yAxis.min, value);
        const ImVec2 screen = projection.project(point);
        char buffer[48];
        std::snprintf(buffer, sizeof(buffer), "%.*f", axis.decimals, value);
        addLabel(drawList, ImVec2(screen.x + metrics.gap * 0.18f,
                                  screen.y + metrics.gap * 0.12f),
                 uiStyle::col::TextFaint, buffer, 0.66f);
      }
      Vec3 end = boxCorners[0];
      if (dimension == 0) end = model(xAxis.max, yAxis.min, zAxis.min);
      if (dimension == 1) end = model(xAxis.min, yAxis.max, zAxis.min);
      if (dimension == 2) end = model(xAxis.min, yAxis.min, zAxis.max);
      addLabel(drawList, projection.project(end), uiStyle::col::TextDim, label);
    };
    drawTicks(xAxis, 0, style.xLabel);
    drawTicks(yAxis, 1, style.yLabel);
    drawTicks(zAxis, 2, style.zLabel);
  }

  if (style.showLabels) {
    std::sort(screenPoints.begin(), screenPoints.end(),
              [](const ScreenPoint& left, const ScreenPoint& right) {
                return left.depth > right.depth;
              });
    std::vector<std::pair<ImVec2, ImVec2>> occupied;
    for (const ScreenPoint& screenPoint : screenPoints) {
      const char* label = points[screenPoint.index].label;
      if (!label || label[0] == '\0') continue;
      const float sizeLabel = chartLabelSize(0.70f);
      const ImVec2 extent = textSize(label, sizeLabel);
      const ImVec2 minimum(screenPoint.position.x + screenPoint.radius + metrics.gap * 0.2f,
                           screenPoint.position.y - extent.y * 0.5f);
      const ImVec2 maximum(minimum.x + extent.x, minimum.y + extent.y);
      bool collision = false;
      for (const auto& used : occupied) {
        if (rectsOverlap(minimum, maximum, used.first, used.second)) {
          collision = true;
          break;
        }
      }
      if (!collision && maximum.x <= rect.max.x && minimum.y >= rect.min.y &&
          maximum.y <= rect.max.y) {
        addLabel(drawList, minimum, uiStyle::col::TextDim, label, 0.70f);
        occupied.emplace_back(minimum, maximum);
      }
    }
  }
  drawList->PopClipRect();

  int hoveredIndex = -1;
  if (hovered) {
    const ImVec2 mouse = ImGui::GetIO().MousePos;
    float best = std::numeric_limits<float>::max();
    for (const ScreenPoint& screenPoint : screenPoints) {
      const float distance = distanceSquared(mouse, screenPoint.position);
      const float reach = std::max(screenPoint.radius, fontSize * 0.42f);
      if (distance <= reach * reach && distance < best) {
        best = distance;
        hoveredIndex = screenPoint.index;
      }
    }
  }
  if (hoveredIndex >= 0) {
    beginChartTooltip();
    if (points[hoveredIndex].label) ImGui::TextUnformatted(points[hoveredIndex].label);
    ImGui::Text("x: %.5g", points[hoveredIndex].x);
    ImGui::Text("y: %.5g", points[hoveredIndex].y);
    ImGui::Text("z: %.5g", points[hoveredIndex].z);
    ImGui::EndTooltip();
  }
  return hoveredIndex;
}

void orbital(const char* id, int n, int l, int m, ImVec2 size, Orbit& orbit,
             const OrbitalStyle& style) {
  const DrawRect rect = reserveRect(id, size);
  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawFrame(drawList, rect);
  const bool hovered = ImGui::IsItemHovered();
  handleOrbitInput(orbit, hovered, ImGui::IsItemActive());
  if (!hasArea(rect)) return;
  if (n < 1 || l < 0 || l >= n || l > 3 || std::abs(m) > l) {
    const char* message = "Invalid orbital quantum numbers";
    const float fontSize = chartLabelSize();
    const ImVec2 extent = textSize(message, fontSize);
    addLabel(drawList,
             ImVec2((rect.min.x + rect.max.x - extent.x) * 0.5f,
                    (rect.min.y + rect.max.y - extent.y) * 0.5f),
             uiStyle::col::TextFaint, message);
    return;
  }

  const int thetaSteps = 28;
  const int phiSteps = 56;
  const double maximumRadius = std::max(8.0, 5.5 * n * n);
  double radialMaximum = 0.0;
  for (int sample = 0; sample <= 1024; ++sample) {
    const double fraction = static_cast<double>(sample) / 1024.0;
    const double radius = maximumRadius * fraction * fraction;
    radialMaximum = std::max(radialMaximum, std::abs(radialHydrogenic(n, l, radius)));
  }
  double angularMaximum = 0.0;
  for (int theta = 0; theta <= thetaSteps; ++theta) {
    for (int phi = 0; phi < phiSteps; ++phi) {
      angularMaximum = std::max(angularMaximum,
          std::abs(realSphericalHarmonic(l, m, kPi * theta / thetaSteps,
                                        2.0 * kPi * phi / phiSteps)));
    }
  }
  const double iso = std::clamp(static_cast<double>(style.isoLevel), 0.02, 0.95);
  const double target = std::max(radialMaximum * angularMaximum * iso,
                                 std::numeric_limits<double>::epsilon());

  struct AngularVertex {
    double harmonic = 0.0;
    std::vector<double> radii;
  };
  std::vector<AngularVertex> angular(
      static_cast<std::size_t>((thetaSteps + 1) * phiSteps));
  int branches = 0;
  for (int thetaIndex = 0; thetaIndex <= thetaSteps; ++thetaIndex) {
    const double theta = kPi * thetaIndex / thetaSteps;
    for (int phiIndex = 0; phiIndex < phiSteps; ++phiIndex) {
      const double phi = 2.0 * kPi * phiIndex / phiSteps;
      AngularVertex& vertex = angular[static_cast<std::size_t>(thetaIndex * phiSteps + phiIndex)];
      vertex.harmonic = realSphericalHarmonic(l, m, theta, phi);
      vertex.radii = isoCrossings(n, l, std::abs(vertex.harmonic), target, maximumRadius);
      branches = std::max(branches, static_cast<int>(vertex.radii.size()));
    }
  }

  std::vector<Vec3> fitPoints;
  for (int thetaIndex = 0; thetaIndex <= thetaSteps; ++thetaIndex) {
    const float theta = kPi * thetaIndex / thetaSteps;
    for (int phiIndex = 0; phiIndex < phiSteps; ++phiIndex) {
      const float phi = 2.0f * kPi * phiIndex / phiSteps;
      const AngularVertex& vertex = angular[static_cast<std::size_t>(thetaIndex * phiSteps + phiIndex)];
      for (double radius : vertex.radii) {
        fitPoints.push_back(spherical(static_cast<float>(radius / maximumRadius), theta, phi));
      }
    }
  }
  if (fitPoints.empty()) {
    const char* message = "Iso level has no visible surface";
    const float fontSize = chartLabelSize();
    const ImVec2 extent = textSize(message, fontSize);
    addLabel(drawList,
             ImVec2((rect.min.x + rect.max.x - extent.x) * 0.5f,
                    (rect.min.y + rect.max.y - extent.y) * 0.5f),
             uiStyle::col::TextFaint, message);
    return;
  }
  float visibleRadius = 0.0f;
  for (Vec3 point : fitPoints) visibleRadius = std::max(visibleRadius, length(point));
  const float angularNodeRadius = visibleRadius * 1.03f;
  const float axisLength = visibleRadius * 1.14f;
  if (style.showNodes) {
    fitPoints.push_back({angularNodeRadius, 0.0f, 0.0f});
    fitPoints.push_back({-angularNodeRadius, 0.0f, 0.0f});
    fitPoints.push_back({0.0f, angularNodeRadius, 0.0f});
    fitPoints.push_back({0.0f, -angularNodeRadius, 0.0f});
    fitPoints.push_back({0.0f, 0.0f, angularNodeRadius});
    fitPoints.push_back({0.0f, 0.0f, -angularNodeRadius});
  }
  if (style.showAxes) {
    fitPoints.push_back({axisLength, 0.0f, 0.0f});
    fitPoints.push_back({0.0f, axisLength, 0.0f});
    fitPoints.push_back({0.0f, 0.0f, axisLength});
  }
  const Projection projection = fittedProjection(rect, orbit, fitPoints);
  const auto& metrics = uiStyle::metrics();
  std::vector<Primitive> primitives;
  primitives.reserve(static_cast<std::size_t>(branches * thetaSteps * phiSteps));

  for (int branch = 0; branch < branches; ++branch) {
    for (int thetaIndex = 0; thetaIndex < thetaSteps; ++thetaIndex) {
      const float theta0 = kPi * thetaIndex / thetaSteps;
      const float theta1 = kPi * (thetaIndex + 1) / thetaSteps;
      for (int phiIndex = 0; phiIndex < phiSteps; ++phiIndex) {
        const int phiNext = (phiIndex + 1) % phiSteps;
        const float phi0 = 2.0f * kPi * phiIndex / phiSteps;
        const float phi1 = 2.0f * kPi * (phiIndex + 1) / phiSteps;
        const AngularVertex* vertices[] = {
            &angular[static_cast<std::size_t>(thetaIndex * phiSteps + phiIndex)],
            &angular[static_cast<std::size_t>(thetaIndex * phiSteps + phiNext)],
            &angular[static_cast<std::size_t>((thetaIndex + 1) * phiSteps + phiNext)],
            &angular[static_cast<std::size_t>((thetaIndex + 1) * phiSteps + phiIndex)]};
        bool complete = true;
        for (const AngularVertex* vertex : vertices) {
          if (branch >= static_cast<int>(vertex->radii.size())) complete = false;
        }
        if (!complete) continue;
        const float radii[] = {
            static_cast<float>(vertices[0]->radii[branch] / maximumRadius),
            static_cast<float>(vertices[1]->radii[branch] / maximumRadius),
            static_cast<float>(vertices[2]->radii[branch] / maximumRadius),
            static_cast<float>(vertices[3]->radii[branch] / maximumRadius)};
        const Vec3 corners[] = {spherical(radii[0], theta0, phi0),
                                spherical(radii[1], theta0, phi1),
                                spherical(radii[2], theta1, phi1),
                                spherical(radii[3], theta1, phi0)};
        if (style.cutaway && representativeDepth(projection, corners, 4) > 0.0f) continue;
        const double averageRadius = (vertices[0]->radii[branch] + vertices[1]->radii[branch] +
                                      vertices[2]->radii[branch] + vertices[3]->radii[branch]) * 0.25;
        const double harmonic = (vertices[0]->harmonic + vertices[1]->harmonic +
                                 vertices[2]->harmonic + vertices[3]->harmonic) * 0.25;
        const bool positive = radialHydrogenic(n, l, averageRadius) * harmonic >= 0.0;
        const Vec3 normal = cross(corners[1] - corners[0], corners[3] - corners[0]);
        Primitive quad;
        quad.kind = PrimitiveKind::Fill;
        quad.count = 4;
        quad.depth = representativeDepth(projection, corners, 4);
        for (int corner = 0; corner < 4; ++corner) quad.points[corner] = projection.project(corners[corner]);
        const ImVec4 phase = positive ? style.positive : style.negative;
        quad.colour = uiStyle::mix(uiStyle::col::BgSurface, phase,
                                   lambert(normal, projection), 0.93f);
        quad.strokeColour = uiStyle::u32(uiStyle::col::GridLine, 0.18f);
        quad.thickness = metrics.hairline;
        primitives.push_back(quad);
      }
    }
  }

  if (style.showNodes) {
    const std::vector<double> nodes = radialNodes(n, l, maximumRadius);
    const int ringSegments = 48;
    for (double node : nodes) {
      const float radius = static_cast<float>(node / maximumRadius);
      for (int plane = 0; plane < 3; ++plane) {
        for (int segment = 0; segment < ringSegments; ++segment) {
          const float a0 = 2.0f * kPi * segment / ringSegments;
          const float a1 = 2.0f * kPi * (segment + 1) / ringSegments;
          Vec3 p0;
          Vec3 p1;
          if (plane == 0) {
            p0 = {radius * std::cos(a0), radius * std::sin(a0), 0.0f};
            p1 = {radius * std::cos(a1), radius * std::sin(a1), 0.0f};
          } else if (plane == 1) {
            p0 = {radius * std::cos(a0), 0.0f, radius * std::sin(a0)};
            p1 = {radius * std::cos(a1), 0.0f, radius * std::sin(a1)};
          } else {
            p0 = {0.0f, radius * std::cos(a0), radius * std::sin(a0)};
            p1 = {0.0f, radius * std::cos(a1), radius * std::sin(a1)};
          }
          const Vec3 endpoints[] = {p0, p1};
          Primitive line;
          line.kind = PrimitiveKind::Line;
          line.count = 2;
          line.points[0] = projection.project(p0);
          line.points[1] = projection.project(p1);
          line.depth = representativeDepth(projection, endpoints, 2);
          line.colour = uiStyle::u32(uiStyle::col::DataBright, 0.20f);
          line.thickness = metrics.hairline;
          primitives.push_back(line);
        }
      }
    }

    const int nodeTheta = 36;
    const int nodePhi = 72;
    for (int thetaIndex = 0; thetaIndex < nodeTheta; ++thetaIndex) {
      const double theta0 = kPi * thetaIndex / nodeTheta;
      const double theta1 = kPi * (thetaIndex + 1) / nodeTheta;
      for (int phiIndex = 0; phiIndex < nodePhi; ++phiIndex) {
        const double phi0 = 2.0 * kPi * phiIndex / nodePhi;
        const double phi1 = 2.0 * kPi * (phiIndex + 1) / nodePhi;
        const double samples[] = {
            realSphericalHarmonic(l, m, theta0, phi0),
            realSphericalHarmonic(l, m, theta0, phi1),
            realSphericalHarmonic(l, m, theta1, phi1),
            realSphericalHarmonic(l, m, theta1, phi0)};
        const double thetas[] = {theta0, theta0, theta1, theta1};
        const double phis[] = {phi0, phi1, phi1, phi0};
        Vec3 intersections[2];
        int intersectionCount = 0;
        for (int edge = 0; edge < 4 && intersectionCount < 2; ++edge) {
          const int next = (edge + 1) % 4;
          if ((samples[edge] < 0.0) == (samples[next] < 0.0) ||
              std::abs(samples[next] - samples[edge]) <=
                  std::numeric_limits<double>::epsilon()) continue;
          const double t = -samples[edge] / (samples[next] - samples[edge]);
          const float theta = static_cast<float>(thetas[edge] +
              t * (thetas[next] - thetas[edge]));
          const float phi = static_cast<float>(phis[edge] + t * (phis[next] - phis[edge]));
          intersections[intersectionCount++] = spherical(angularNodeRadius, theta, phi);
        }
        if (intersectionCount == 2) {
          Primitive line;
          line.kind = PrimitiveKind::Line;
          line.count = 2;
          line.points[0] = projection.project(intersections[0]);
          line.points[1] = projection.project(intersections[1]);
          line.depth = representativeDepth(projection, intersections, 2);
          line.colour = uiStyle::u32(uiStyle::col::DataBright, 0.24f);
          line.thickness = metrics.hairline;
          primitives.push_back(line);
        }
      }
    }
  }

  if (style.showAxes) {
    const Vec3 axes[][2] = {{{0.0f,0.0f,0.0f},{axisLength,0.0f,0.0f}},
                            {{0.0f,0.0f,0.0f},{0.0f,axisLength,0.0f}},
                            {{0.0f,0.0f,0.0f},{0.0f,0.0f,axisLength}}};
    for (const auto& axis : axes) {
      Primitive line;
      line.kind = PrimitiveKind::Line;
      line.count = 2;
      line.points[0] = projection.project(axis[0]);
      line.points[1] = projection.project(axis[1]);
      line.depth = representativeDepth(projection, axis, 2);
      line.colour = uiStyle::u32(uiStyle::col::DataDim, 0.78f);
      line.thickness = metrics.hairline;
      primitives.push_back(line);
    }
  }

  drawList->PushClipRect(rect.min, rect.max, true);
  sortAndEmit(drawList, primitives);
  if (style.showAxes) {
    addLabel(drawList, projection.project({axisLength, 0.0f, 0.0f}), uiStyle::col::TextDim, "x");
    addLabel(drawList, projection.project({0.0f, axisLength, 0.0f}), uiStyle::col::TextDim, "y");
    addLabel(drawList, projection.project({0.0f, 0.0f, axisLength}), uiStyle::col::TextDim, "z");
  }
  addLabel(drawList, ImVec2(rect.min.x + metrics.gap * 0.45f,
                            rect.min.y + metrics.gap * 0.35f),
           uiStyle::col::TextDim, orbitalName(n, l, m), 0.82f);
  drawList->PopClipRect();
}

const char* orbitalName(int n, int l, int m) {
  if (n < 1 || l < 0 || l >= n || l > 3 || std::abs(m) > l) return "invalid";
  static thread_local char name[32];
  const char* suffix = "";
  if (l == 0) suffix = "s";
  if (l == 1) {
    if (m == 0) suffix = "p_z";
    if (m == 1) suffix = "p_x";
    if (m == -1) suffix = "p_y";
  }
  if (l == 2) {
    if (m == 0) suffix = "d_z2";
    if (m == 1) suffix = "d_xz";
    if (m == -1) suffix = "d_yz";
    if (m == 2) suffix = "d_x2-y2";
    if (m == -2) suffix = "d_xy";
  }
  if (l == 3) {
    if (m == 0) suffix = "f_z3";
    if (m == 1) suffix = "f_xz2";
    if (m == -1) suffix = "f_yz2";
    if (m == 2) suffix = "f_z(x2-y2)";
    if (m == -2) suffix = "f_xyz";
    if (m == 3) suffix = "f_x(x2-3y2)";
    if (m == -3) suffix = "f_y(3x2-y2)";
  }
  std::snprintf(name, sizeof(name), "%d%s", n, suffix);
  return name;
}

int parallelCoordinates(const char* id, const ParallelAxis* axes, int axisCount,
                        const ParallelSeries* series, int seriesCount, ImVec2 size) {
  const DrawRect rect = reserveRect(id, size);
  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawFrame(drawList, rect);
  if (!hasArea(rect) || !axes || axisCount < 2 || !series || seriesCount <= 0) return -1;

  bool haveValue = false;
  bool allValuesEqual = true;
  double firstValue = 0.0;
  for (int seriesIndex = 0; seriesIndex < seriesCount; ++seriesIndex) {
    if (!series[seriesIndex].values) continue;
    for (int axisIndex = 0; axisIndex < axisCount; ++axisIndex) {
      const double value = series[seriesIndex].values[axisIndex];
      if (!finite(value)) continue;
      if (!haveValue) {
        firstValue = value;
        haveValue = true;
      } else if (!nearlyEqual(firstValue, value)) {
        allValuesEqual = false;
      }
    }
  }
  if (!haveValue || allValuesEqual) return -1;

  const auto& metrics = uiStyle::metrics();
  const float fontSize = ImGui::GetFontSize();
  const float labelSize = chartLabelSize(0.70f);
  const float top = rect.min.y + fontSize * 1.65f;
  const float bottom = rect.max.y - fontSize * 1.85f;
  const float left = rect.min.x + fontSize * 1.10f;
  const float right = rect.max.x - fontSize * 1.10f;
  if (right <= left || bottom <= top) return -1;
  auto axisX = [&](int index) {
    return left + (right - left) * static_cast<float>(index) /
                      static_cast<float>(axisCount - 1);
  };
  auto valueY = [&](int axisIndex, double value) {
    const double low = axes[axisIndex].min;
    const double high = axes[axisIndex].max;
    const double span = high - low;
    const double normalized = finite(value) && std::abs(span) >
        std::numeric_limits<double>::epsilon()
        ? std::clamp((value - low) / span, 0.0, 1.0) : 0.5;
    return bottom - static_cast<float>(normalized) * (bottom - top);
  };

  drawList->PushClipRect(rect.min, rect.max, true);
  for (int axisIndex = 0; axisIndex < axisCount; ++axisIndex) {
    const float x = axisX(axisIndex);
    const ImVec4 direction = axes[axisIndex].higherIsBetter
                                 ? uiStyle::col::Success : uiStyle::col::Danger;
    drawList->AddLine(ImVec2(x, top), ImVec2(x, bottom),
                      uiStyle::u32(uiStyle::col::DataDim), metrics.hairline);
    const float arrow = fontSize * 0.18f;
    if (axes[axisIndex].higherIsBetter) {
      drawList->AddTriangleFilled(ImVec2(x, top - arrow),
                                  ImVec2(x - arrow, top + arrow),
                                  ImVec2(x + arrow, top + arrow),
                                  uiStyle::u32(direction, 0.82f));
    } else {
      drawList->AddTriangleFilled(ImVec2(x, bottom + arrow),
                                  ImVec2(x - arrow, bottom - arrow),
                                  ImVec2(x + arrow, bottom - arrow),
                                  uiStyle::u32(direction, 0.82f));
    }
    std::string title = axes[axisIndex].label ? axes[axisIndex].label : "";
    if (axes[axisIndex].unit && axes[axisIndex].unit[0] != '\0') {
      title += " (";
      title += axes[axisIndex].unit;
      title += ")";
    }
    const ImVec2 titleExtent = textSize(title.c_str(), labelSize);
    addLabel(drawList, ImVec2(x - titleExtent.x * 0.5f, rect.min.y + metrics.gap * 0.22f),
             uiStyle::col::TextDim, title.c_str(), 0.70f);
    char maximum[48];
    char minimum[48];
    std::snprintf(maximum, sizeof(maximum), "%.4g", axes[axisIndex].max);
    std::snprintf(minimum, sizeof(minimum), "%.4g", axes[axisIndex].min);
    const ImVec2 maxExtent = textSize(maximum, labelSize);
    const ImVec2 minExtent = textSize(minimum, labelSize);
    addLabel(drawList, ImVec2(x - maxExtent.x * 0.5f, top - maxExtent.y - metrics.gap * 0.12f),
             uiStyle::col::TextFaint, maximum, 0.70f);
    addLabel(drawList, ImVec2(x - minExtent.x * 0.5f, bottom + metrics.gap * 0.16f),
             uiStyle::col::TextFaint, minimum, 0.70f);
  }

  const ImVec2 mouse = ImGui::GetIO().MousePos;
  int hoveredSeries = -1;
  float bestDistance = std::numeric_limits<float>::max();
  for (int pass = 0; pass < 2; ++pass) {
    for (int seriesIndex = 0; seriesIndex < seriesCount; ++seriesIndex) {
      if (!series[seriesIndex].values || static_cast<int>(series[seriesIndex].highlighted) != pass) continue;
      const float thickness = metrics.hairline * (series[seriesIndex].highlighted ? 2.8f : 1.15f);
      const ImVec4 colour = series[seriesIndex].highlighted
                                ? uiStyle::col::DataBright : series[seriesIndex].colour;
      std::vector<ImVec2> line;
      line.reserve(static_cast<std::size_t>(axisCount));
      for (int axisIndex = 0; axisIndex < axisCount; ++axisIndex) {
        line.emplace_back(axisX(axisIndex),
                          valueY(axisIndex, series[seriesIndex].values[axisIndex]));
      }
      drawList->AddPolyline(line.data(), axisCount,
                            uiStyle::u32(colour, series[seriesIndex].highlighted ? 1.0f : 0.42f),
                            ImDrawFlags_None, thickness);
      for (int axisIndex = 0; axisIndex + 1 < axisCount; ++axisIndex) {
        const float distance = pointSegmentDistanceSquared(mouse, line[axisIndex],
                                                            line[axisIndex + 1]);
        if (distance < bestDistance) {
          bestDistance = distance;
          hoveredSeries = seriesIndex;
        }
      }
    }
  }
  drawList->PopClipRect();
  const float reach = fontSize * 0.45f;
  if (!ImGui::IsItemHovered() || bestDistance > reach * reach) hoveredSeries = -1;
  if (hoveredSeries >= 0) {
    beginChartTooltip();
    if (series[hoveredSeries].label) ImGui::TextUnformatted(series[hoveredSeries].label);
    for (int axisIndex = 0; axisIndex < axisCount; ++axisIndex) {
      tooltipValue(axes[axisIndex].label ? axes[axisIndex].label : "value",
                   series[hoveredSeries].values[axisIndex], axes[axisIndex].unit);
    }
    ImGui::EndTooltip();
  }
  return hoveredSeries;
}

void waterfall(const char* id, double start, const WaterfallStep* steps, int count,
               const char* unit, ImVec2 size) {
  const DrawRect rect = reserveRect(id, size);
  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawFrame(drawList, rect);
  if (!hasArea(rect) || !steps || count <= 0) return;

  std::vector<double> before(static_cast<std::size_t>(count));
  std::vector<double> after(static_cast<std::size_t>(count));
  double running = finite(start) ? start : 0.0;
  double minimum = std::min(0.0, running);
  double maximum = std::max(0.0, running);
  for (int index = 0; index < count; ++index) {
    before[static_cast<std::size_t>(index)] = running;
    const double delta = finite(steps[index].delta) ? steps[index].delta : 0.0;
    running = steps[index].total ? delta : running + delta;
    after[static_cast<std::size_t>(index)] = running;
    minimum = std::min({minimum, before[static_cast<std::size_t>(index)], running});
    maximum = std::max({maximum, before[static_cast<std::size_t>(index)], running});
  }
  if (nearlyEqual(minimum, maximum)) return;
  const charts::Axis axis = charts::niceAxis(minimum, maximum);
  const auto& metrics = uiStyle::metrics();
  const float fontSize = ImGui::GetFontSize();
  const float labelSize = chartLabelSize(0.70f);
  const float left = rect.min.x + fontSize * 2.6f;
  const float right = rect.max.x - metrics.gap * 0.55f;
  const float top = rect.min.y + metrics.gap * 0.55f;
  const float bottom = rect.max.y - fontSize * 1.75f;
  if (right <= left || bottom <= top) return;
  auto valueY = [&](double value) {
    return bottom - static_cast<float>(axis.normalise(value)) * (bottom - top);
  };

  drawList->PushClipRect(rect.min, rect.max, true);
  for (int tick = 0; tick < axis.ticks; ++tick) {
    const double value = axis.min + axis.step * tick;
    if (value > axis.max + axis.step * 0.25) break;
    const float y = valueY(value);
    drawList->AddLine(ImVec2(left, y), ImVec2(right, y),
                      uiStyle::u32(uiStyle::col::GridLine), metrics.hairline);
    char buffer[48];
    std::snprintf(buffer, sizeof(buffer), "%.*f", axis.decimals, value);
    const ImVec2 extent = textSize(buffer, labelSize);
    addLabel(drawList, ImVec2(left - extent.x - metrics.gap * 0.28f,
                              y - extent.y * 0.5f),
             uiStyle::col::TextFaint, buffer, 0.70f);
  }

  const float slot = (right - left) / static_cast<float>(count + 1);
  const float barWidth = slot * 0.62f;
  auto barCentre = [&](int index) { return left + slot * (index + 1.0f); };
  float previousRight = left;
  float previousY = valueY(start);
  for (int index = 0; index < count; ++index) {
    const float centre = barCentre(index);
    const float x0 = centre - barWidth * 0.5f;
    const float x1 = centre + barWidth * 0.5f;
    const double lowerValue = steps[index].total ? 0.0
        : std::min(before[static_cast<std::size_t>(index)], after[static_cast<std::size_t>(index)]);
    const double upperValue = steps[index].total ? after[static_cast<std::size_t>(index)]
        : std::max(before[static_cast<std::size_t>(index)], after[static_cast<std::size_t>(index)]);
    const float y0 = valueY(upperValue);
    const float y1 = valueY(lowerValue);
    const bool gain = after[static_cast<std::size_t>(index)] >=
                      before[static_cast<std::size_t>(index)];
    const ImVec4 colour = steps[index].total ? uiStyle::col::Data
                                             : (gain ? uiStyle::col::Success
                                                     : uiStyle::col::Danger);
    if (index > 0) {
      drawList->AddLine(ImVec2(previousRight, previousY), ImVec2(x0, previousY),
                        uiStyle::u32(uiStyle::col::DataDim, 0.70f), metrics.hairline);
    }
    drawList->AddRectFilled(ImVec2(x0, std::min(y0, y1)), ImVec2(x1, std::max(y0, y1)),
                            uiStyle::u32(colour, 0.82f), metrics.radiusSm);
    drawList->AddRect(ImVec2(x0, std::min(y0, y1)), ImVec2(x1, std::max(y0, y1)),
                      uiStyle::u32(colour), metrics.radiusSm, ImDrawFlags_None,
                      metrics.hairline);

    char valueBuffer[64];
    if (unit && unit[0] != '\0') {
      std::snprintf(valueBuffer, sizeof(valueBuffer), "%+.4g %s",
                    steps[index].total ? after[static_cast<std::size_t>(index)]
                                       : steps[index].delta, unit);
    } else {
      std::snprintf(valueBuffer, sizeof(valueBuffer), "%+.4g",
                    steps[index].total ? after[static_cast<std::size_t>(index)]
                                       : steps[index].delta);
    }
    const ImVec2 valueExtent = textSize(valueBuffer, labelSize);
    if (valueExtent.x <= slot * 0.96f && valueExtent.y <= std::abs(y1 - y0) + fontSize) {
      const float labelY = std::min(y0, y1) - valueExtent.y - metrics.gap * 0.12f;
      addLabel(drawList, ImVec2(centre - valueExtent.x * 0.5f, labelY),
               uiStyle::col::TextDim, valueBuffer, 0.70f);
    }
    if (steps[index].label && steps[index].label[0] != '\0') {
      const ImVec2 extent = textSize(steps[index].label, labelSize);
      if (extent.x <= slot * 0.96f) {
        addLabel(drawList, ImVec2(centre - extent.x * 0.5f,
                                  bottom + metrics.gap * 0.28f),
                 uiStyle::col::TextDim, steps[index].label, 0.70f);
      }
    }
    previousRight = x1;
    previousY = valueY(after[static_cast<std::size_t>(index)]);
  }
  drawList->PopClipRect();
}

}  // namespace chemcad::ui::charts3d

// 2D cross-section renderer and controls for the liquid-liquid separation
// simulation (sol::Simulation): the Extraction Lab. Renders inside an
// already-open panel; never owns a top-level ImGui window.

#include "imgui.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "core/model.hpp"
#include "sol/funnel.hpp"
#include "sol/solubility.hpp"
#include "ui/app_state.hpp"
#include "ui/solubility_state.hpp"
#include "ui/theme.hpp"
#include "ui/widgets.hpp"

namespace chemcad::ui {

// Redundant forward declaration: the app shell (ui.hpp) declares the public
// signature this file implements.
void drawExtractionLab(AppState&);

namespace {

// ------------------------------------------------------------- vessel scale
// The fixed headers expose no shared "vessel height in metres" constant, so
// this picks one deterministically from the charged capacity: enough for the
// silhouette and layer bands below to stay internally consistent frame to
// frame. Clamped so tiny/huge charges never collapse the render.
float heightForVolumeMl(double vesselVolumeMl) {
  const double clampedMl = std::max(vesselVolumeMl, 1.0);
  return static_cast<float>(std::clamp(0.05 + 0.00035 * clampedMl, 0.10, 0.42));
}

constexpr int kWidthSamples = 129;  // resolution of the cached height->volume profile

// Precomputed per-vessel geometry: the wall silhouette plus a cumulative
// cross-sectional-volume profile used to turn a settled mL amount into a
// height fraction (and back) without re-sampling every layer every frame.
struct VesselGeometry {
  float heightMetres = 0.18f;
  std::vector<core::Vec2> outline;                 // vessel-local metres
  double halfWidthMetres = 0.05;                    // widest half-width, from the outline
  std::array<double, kWidthSamples> cumVolumeMl{};  // monotonic, cumVolumeMl[last] == capacityMl
};

double widthFractionAt(sol::Vessel vessel, double heightFraction) {
  return sol::vesselWidthAt(vessel, std::clamp(heightFraction, 0.0, 1.0));
}

VesselGeometry buildGeometry(sol::Vessel vessel, double vesselVolumeMl) {
  VesselGeometry geo;
  geo.heightMetres = heightForVolumeMl(vesselVolumeMl);
  geo.outline = sol::vesselOutline(vessel, static_cast<double>(geo.heightMetres));

  double maxAbsX = 1e-6;
  for (const core::Vec2& p : geo.outline) {
    maxAbsX = std::max(maxAbsX, static_cast<double>(std::fabs(p.x)));
  }
  geo.halfWidthMetres = maxAbsX;

  // Cross-sectional area at a height is treated as proportional to width^2
  // (an axisymmetric revolution of the 2D silhouette), so the cumulative
  // trapezoidal integral of width^2 gives a volume-vs-height profile.
  std::array<double, kWidthSamples> cumRaw{};
  cumRaw[0] = 0.0;
  const double dHf = 1.0 / static_cast<double>(kWidthSamples - 1);
  double prevArea = widthFractionAt(vessel, 0.0);
  prevArea *= prevArea;
  for (int i = 1; i < kWidthSamples; ++i) {
    const double hf = static_cast<double>(i) * dHf;
    const double w = widthFractionAt(vessel, hf);
    const double area = w * w;
    cumRaw[i] = cumRaw[i - 1] + 0.5 * (prevArea + area) * dHf;
    prevArea = area;
  }
  const double total = std::max(cumRaw[kWidthSamples - 1], 1e-9);
  const double capacityMl = std::max(vesselVolumeMl, 1.0);
  for (int i = 0; i < kWidthSamples; ++i) geo.cumVolumeMl[i] = capacityMl * cumRaw[i] / total;
  return geo;
}

// Rebuilds only when the vessel type or charged capacity actually changes;
// everything else (layers, droplets) reuses this per-frame without
// re-sampling the width profile.
const VesselGeometry& cachedGeometry(sol::Vessel vessel, double vesselVolumeMl) {
  static thread_local sol::Vessel lastVessel = sol::Vessel::SeparatoryFunnel;
  static thread_local double lastVolumeMl = -1.0;
  static thread_local bool initialized = false;
  static thread_local VesselGeometry cache;
  if (!initialized || lastVessel != vessel || std::fabs(lastVolumeMl - vesselVolumeMl) > 1e-6) {
    cache = buildGeometry(vessel, vesselVolumeMl);
    lastVessel = vessel;
    lastVolumeMl = vesselVolumeMl;
    initialized = true;
  }
  return cache;
}

// Inverts the cumulative-volume profile: the height fraction whose charged
// volume, measured from the bottom, equals `volumeMl`.
double heightFractionForVolume(const VesselGeometry& geo, double volumeMl) {
  const double target = std::clamp(volumeMl, 0.0, geo.cumVolumeMl[kWidthSamples - 1]);
  if (target <= 0.0) return 0.0;
  for (int i = 1; i < kWidthSamples; ++i) {
    if (geo.cumVolumeMl[i] >= target) {
      const double lo = geo.cumVolumeMl[i - 1];
      const double hi = geo.cumVolumeMl[i];
      const double t = (target - lo) / std::max(hi - lo, 1e-9);
      const double hfLo = static_cast<double>(i - 1) / (kWidthSamples - 1);
      const double hfHi = static_cast<double>(i) / (kWidthSamples - 1);
      return hfLo + t * (hfHi - hfLo);
    }
  }
  return 1.0;
}

// ---------------------------------------------------------- screen mapping
struct Transform {
  ImVec2 origin{0.0f, 0.0f};  // pixel position of vessel-local (0, 0): bottom centre
  float scale = 1.0f;         // pixels per metre
};

constexpr float kMarginLeft = 64.0f;    // room for the graduation scale
constexpr float kMarginRight = 20.0f;
constexpr float kMarginTop = 24.0f;
constexpr float kMarginBottom = 20.0f;

Transform buildTransform(const VesselGeometry& geo, ImVec2 regionMin, ImVec2 regionSize) {
  const float usableW = std::max(regionSize.x - kMarginLeft - kMarginRight, 1.0f);
  const float usableH = std::max(regionSize.y - kMarginTop - kMarginBottom, 1.0f);
  const float widthMetres = std::max(static_cast<float>(geo.halfWidthMetres) * 2.0f, 1e-4f);
  const float heightMetres = std::max(geo.heightMetres, 1e-4f);
  const float scale = std::min(usableW / widthMetres, usableH / heightMetres);

  Transform tf;
  tf.scale = scale;
  tf.origin.x = regionMin.x + kMarginLeft + usableW * 0.5f;
  tf.origin.y = regionMin.y + kMarginTop + heightMetres * scale;
  return tf;
}

ImVec2 toScreen(const Transform& tf, double x, double y) {
  return ImVec2(tf.origin.x + static_cast<float>(x) * tf.scale,
               tf.origin.y - static_cast<float>(y) * tf.scale);
}

ImVec2 toScreen(const Transform& tf, core::Vec2 p) { return toScreen(tf, p.x, p.y); }

// -------------------------------------------------------------- tick scale
double niceTickStep(double range) {
  const double roughStep = std::max(range, 1e-6) / 6.0;
  const double magnitude = std::pow(10.0, std::floor(std::log10(roughStep)));
  const double residual = roughStep / magnitude;
  double niceResidual = 10.0;
  if (residual < 1.5) niceResidual = 1.0;
  else if (residual < 3.0) niceResidual = 2.0;
  else if (residual < 7.0) niceResidual = 5.0;
  return niceResidual * magnitude;
}

void drawGraduation(ImDrawList* draw, const VesselGeometry& geo, const Transform& tf,
                    ImVec2 regionMin) {
  const double capacityMl = geo.cumVolumeMl[kWidthSamples - 1];
  const double step = niceTickStep(capacityMl);
  if (step <= 0.0) return;
  const ImU32 tickColor = style::u32(style::col::TextFaint);
  const ImU32 labelColor = style::u32(style::col::TextDim);
  const float tickX = regionMin.x + kMarginLeft - 14.0f;
  for (double v = 0.0; v <= capacityMl + 1e-6; v += step) {
    const double hf = heightFractionForVolume(geo, v);
    const ImVec2 p = toScreen(tf, 0.0, hf * geo.heightMetres);
    draw->AddLine(ImVec2(tickX, p.y), ImVec2(tickX + 10.0f, p.y), tickColor, 1.0f);
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.0f", v);
    const ImVec2 textSize = ImGui::CalcTextSize(buf);
    draw->AddText(ImVec2(tickX - 4.0f - textSize.x, p.y - textSize.y * 0.5f), labelColor, buf);
  }
}

// ------------------------------------------------------------ status text
const char* emulsionStateLabel(double fraction) {
  if (fraction > 0.35) return "EMULSIFIED";
  if (fraction < 0.02) return "SEPARATED";
  return "SETTLING";
}

ImU32 emulsionStateColor(double fraction) {
  if (fraction > 0.35) return style::u32(style::col::Danger);
  if (fraction < 0.02) return style::u32(style::col::Success);
  return style::u32(style::col::Accent);
}

void drawReadout(ImDrawList* draw, const sol::Simulation& sim, ImVec2 regionMin,
                 ImVec2 regionSize) {
  const style::Metrics& m = style::metrics();
  const float safeLeft = regionMin.x + kMarginLeft + m.gap;
  const float safeRight = regionMin.x + regionSize.x - m.gap;
  const float availableWidth = safeRight - safeLeft;
  if (availableWidth < ImGui::GetFontSize() * 5.0f || regionSize.y < m.gap * 4.0f) return;

  const double fraction = sol::emulsifiedFraction(sim);
  const float lineHeight = ImGui::GetTextLineHeightWithSpacing();
  const float padding = m.gap;
  const float boxWidth = std::min(220.0f, availableWidth);
  const size_t layerCount = std::min(sim.phases.size(), sim.settledMl.size());
  const int fixedLines = sim.shake.active ? 4 : 3;
  const float availableHeight = std::max(regionSize.y - m.gap * 2.0f, 0.0f);
  const int maxLayerLines =
      std::max(0, static_cast<int>((availableHeight - padding * 2.0f) / lineHeight) - fixedLines);
  const size_t visibleLayers =
      std::min(layerCount, static_cast<size_t>(std::max(maxLayerLines, 0)));
  const bool hasHiddenLayers = visibleLayers < layerCount;
  const int contentLines = fixedLines + static_cast<int>(visibleLayers);
  const float boxHeight =
      std::min(availableHeight, padding * 2.0f + lineHeight * static_cast<float>(contentLines));

  const ImVec2 boxMin(safeRight - boxWidth, regionMin.y + m.gap);
  const ImVec2 boxMax(safeRight, boxMin.y + boxHeight);
  draw->AddRectFilled(boxMin, boxMax, style::u32(style::col::BgSurface, 0.92f), m.radiusMd);
  draw->AddRect(boxMin, boxMax, style::u32(style::col::BorderStrong, 0.85f), m.radiusMd, 0,
                m.hairline);
  draw->PushClipRect(ImVec2(boxMin.x + m.hairline, boxMin.y + m.hairline),
                     ImVec2(boxMax.x - m.hairline, boxMax.y - m.hairline), true);

  ImVec2 cursor(boxMin.x + padding, boxMin.y + padding);
  const float wrapWidth = std::max(boxWidth - padding * 2.0f, 1.0f);
  char buf[160];

  std::snprintf(buf, sizeof(buf), "t = %.1f s", sim.elapsed);
  draw->AddText(cursor, style::u32(style::col::TextDim), buf);
  cursor.y += lineHeight;

  if (sim.shake.active) {
    std::snprintf(buf, sizeof(buf), "shaking: %.1f s left", sim.shake.remainingS);
    const float pulse =
        0.68f + 0.32f * (0.5f + 0.5f * std::sin(static_cast<float>(ImGui::GetTime()) * 5.0f));
    draw->AddText(cursor, style::u32(style::col::Accent, pulse), buf);
    cursor.y += lineHeight;
  }

  std::snprintf(buf, sizeof(buf), "emulsified: %.1f%%", fraction * 100.0);
  draw->AddText(cursor, style::u32(style::col::TextDim), buf);
  cursor.y += lineHeight;

  ImFont* bigFont = style::fonts::semibold();
  const float bigSize = ImGui::GetFontSize() * 1.10f;
  draw->AddText(bigFont, bigSize, cursor, emulsionStateColor(fraction),
                emulsionStateLabel(fraction));
  cursor.y += lineHeight;

  for (size_t i = 0; i < visibleLayers; ++i) {
    if (hasHiddenLayers && i + 1 == visibleLayers) {
      std::snprintf(buf, sizeof(buf), "+ %d more phase%s",
                    static_cast<int>(layerCount - i), layerCount - i == 1 ? "" : "s");
      draw->AddText(cursor, style::u32(style::col::TextDim), buf);
    } else {
      std::snprintf(buf, sizeof(buf), "%s: %.0f mL", sim.phases[i].label.c_str(),
                    sim.settledMl[i]);
      draw->AddText(nullptr, 0.0f, cursor, style::u32(style::col::Text), buf, nullptr,
                    wrapWidth);
    }
    cursor.y += lineHeight;
  }

  draw->PopClipRect();
}

// --------------------------------------------------------------- painting
// The vessel is rendered as curved glassware, not a flat polygon: liquid
// layers are drawn as vertical slices whose brightness rolls off toward the
// walls (cylindrical shading), each layer top carries a meniscus ellipse, and
// the glass gets a specular streak on the left wall plus rim shading on the
// right. Everything goes through `tf`, so the whole apparatus -- liquids,
// furniture, highlights -- moves as one during a shake.

// Cylindrical brightness across the vessel width: brightest left of centre
// (key light), rolling off to dark rims. xFrac in [-1, 1].
float glassShade(float xFrac) {
  const float g = std::exp(-std::pow((xFrac + 0.38f) / 0.52f, 2.0f));
  return 0.62f + 0.55f * g;
}

ImU32 phaseShade(const sol::Phase& phase, float xFrac, float alphaScale = 1.0f) {
  const float s = glassShade(xFrac);
  const ImVec4 c(std::min(phase.colour[0] * s, 1.0f), std::min(phase.colour[1] * s, 1.0f),
                 std::min(phase.colour[2] * s, 1.0f),
                 std::min(phase.colour[3] * alphaScale, 1.0f));
  return ImGui::ColorConvertFloat4ToU32(c);
}

void drawGroundShadow(ImDrawList* draw, const VesselGeometry& geo, const Transform& tf) {
  const ImVec2 c = toScreen(tf, 0.0, 0.0);
  const float rx = static_cast<float>(geo.halfWidthMetres) * tf.scale * 0.9f;
  const float ry = std::max(4.0f, rx * 0.14f);
  for (int i = 0; i < 5; ++i) {
    const float t = 1.0f - static_cast<float>(i) / 5.0f;
    draw->AddEllipseFilled(ImVec2(c.x, c.y + ry * 1.6f), ImVec2(rx * t, ry * t),
                           style::u32(style::col::BgDeep, 0.16f + 0.10f * t));
  }
}

void drawVesselGlass(ImDrawList* draw, const VesselGeometry& geo, const Transform& tf) {
  static thread_local std::vector<ImVec2> screenOutline;
  screenOutline.clear();
  screenOutline.reserve(geo.outline.size());
  for (const core::Vec2& p : geo.outline) screenOutline.push_back(toScreen(tf, p));
  if (screenOutline.size() < 3) return;
  draw->AddConcavePolyFilled(screenOutline.data(), static_cast<int>(screenOutline.size()),
                             style::u32(ImVec4(0.60f, 0.76f, 0.88f, 1.0f), 0.10f));
}

void drawVesselWall(ImDrawList* draw, const VesselGeometry& geo, const Transform& tf) {
  static thread_local std::vector<ImVec2> screenOutline;
  screenOutline.clear();
  screenOutline.reserve(geo.outline.size());
  for (const core::Vec2& p : geo.outline) screenOutline.push_back(toScreen(tf, p));
  if (screenOutline.size() < 2) return;
  // Double wall: a wide soft halo under a crisp inner line reads as glass
  // thickness far better than a single stroke.
  draw->AddPolyline(screenOutline.data(), static_cast<int>(screenOutline.size()),
                    style::u32(ImVec4(0.70f, 0.82f, 0.92f, 1.0f), 0.28f), ImDrawFlags_Closed,
                    4.0f);
  draw->AddPolyline(screenOutline.data(), static_cast<int>(screenOutline.size()),
                    style::u32(style::col::BorderStrong, 0.9f), ImDrawFlags_Closed, 1.6f);
}

// Specular streak down the left wall and a dark rim down the right: the two
// cues that make the silhouette read as a curved glass surface.
void drawGlassHighlights(ImDrawList* draw, const sol::Simulation& sim,
                         const VesselGeometry& geo, const Transform& tf) {
  static thread_local std::vector<ImVec2> streak;
  streak.clear();
  const int samples = static_cast<int>(geo.outline.size() / 2);
  for (int i = 0; i < samples; ++i) {
    const double hf = static_cast<double>(i) / std::max(1, samples - 1);
    const double halfW = widthFractionAt(sim.vessel, hf) * geo.halfWidthMetres;
    streak.push_back(toScreen(tf, -halfW * 0.80, hf * geo.heightMetres));
  }
  if (streak.size() >= 2) {
    draw->AddPolyline(streak.data(), static_cast<int>(streak.size()),
                      style::u32(ImVec4(0.92f, 0.97f, 1.0f, 1.0f), 0.10f), ImDrawFlags_None,
                      std::max(3.0f, tf.scale * 0.010f));
    draw->AddPolyline(streak.data(), static_cast<int>(streak.size()),
                      style::u32(ImVec4(0.95f, 0.99f, 1.0f, 1.0f), 0.16f), ImDrawFlags_None,
                      std::max(1.2f, tf.scale * 0.0035f));
  }

  streak.clear();
  for (int i = 0; i < samples; ++i) {
    const double hf = static_cast<double>(i) / std::max(1, samples - 1);
    const double halfW = widthFractionAt(sim.vessel, hf) * geo.halfWidthMetres;
    streak.push_back(toScreen(tf, halfW * 0.88, hf * geo.heightMetres));
  }
  if (streak.size() >= 2) {
    draw->AddPolyline(streak.data(), static_cast<int>(streak.size()),
                      style::u32(style::col::BgDeep, 0.20f), ImDrawFlags_None,
                      std::max(2.5f, tf.scale * 0.008f));
  }
}

// Stopper cap + knob at the neck, stopcock across the drain stem: the two
// pieces of furniture that make the shape unmistakably a separatory funnel.
void drawFurniture(ImDrawList* draw, const sol::Simulation& sim, const VesselGeometry& geo,
                   const Transform& tf) {
  if (sim.vessel != sol::Vessel::SeparatoryFunnel) return;

  const ImU32 glassFill = style::u32(ImVec4(0.72f, 0.84f, 0.93f, 1.0f), 0.35f);
  const ImU32 glassEdge = style::u32(style::col::BorderStrong, 0.9f);

  // Stopper: cap sitting over the neck top, with a small grip knob above it.
  const double neckW = widthFractionAt(sim.vessel, 1.0) * geo.halfWidthMetres;
  const ImVec2 topC = toScreen(tf, 0.0, geo.heightMetres);
  const float capW = static_cast<float>(neckW) * tf.scale * 2.0f + 6.0f;
  const float capH = std::max(7.0f, tf.scale * 0.012f);
  draw->AddRectFilled(ImVec2(topC.x - capW * 0.5f, topC.y - capH),
                      ImVec2(topC.x + capW * 0.5f, topC.y + 1.0f), glassFill, capH * 0.35f);
  draw->AddRect(ImVec2(topC.x - capW * 0.5f, topC.y - capH),
                ImVec2(topC.x + capW * 0.5f, topC.y + 1.0f), glassEdge, capH * 0.35f, 0, 1.2f);
  const float knobR = capH * 0.55f;
  draw->AddCircleFilled(ImVec2(topC.x, topC.y - capH - knobR * 0.8f), knobR, glassFill, 12);
  draw->AddCircle(ImVec2(topC.x, topC.y - capH - knobR * 0.8f), knobR, glassEdge, 12, 1.2f);

  // Stopcock: a barrel across the stem with a handle tab on the right.
  const double stemW = widthFractionAt(sim.vessel, 0.10) * geo.halfWidthMetres;
  const ImVec2 cockC = toScreen(tf, 0.0, 0.10 * geo.heightMetres);
  const float barW = std::max(static_cast<float>(stemW) * tf.scale * 4.2f, 26.0f);
  const float barH = std::max(7.0f, tf.scale * 0.011f);
  draw->AddRectFilled(ImVec2(cockC.x - barW * 0.5f, cockC.y - barH * 0.5f),
                      ImVec2(cockC.x + barW * 0.5f, cockC.y + barH * 0.5f), glassFill,
                      barH * 0.5f);
  draw->AddRect(ImVec2(cockC.x - barW * 0.5f, cockC.y - barH * 0.5f),
                ImVec2(cockC.x + barW * 0.5f, cockC.y + barH * 0.5f), glassEdge, barH * 0.5f, 0,
                1.2f);
  const ImVec2 tabA(cockC.x + barW * 0.5f, cockC.y);
  const ImVec2 tabB(cockC.x + barW * 0.5f + barH * 1.5f, cockC.y + barH * 0.9f);
  draw->AddLine(tabA, tabB, glassEdge, 3.0f);
  draw->AddCircleFilled(tabB, barH * 0.42f, glassFill, 10);
  draw->AddCircle(tabB, barH * 0.42f, glassEdge, 10, 1.1f);
}

// Meniscus ellipse at a liquid surface: bright rim over a soft fill, both
// squeezed to the vessel width at that height.
void drawMeniscus(ImDrawList* draw, const sol::Phase& phase, double hf,
                  const sol::Simulation& sim, const VesselGeometry& geo, const Transform& tf) {
  const double halfW = widthFractionAt(sim.vessel, hf) * geo.halfWidthMetres;
  if (halfW * tf.scale < 6.0) return;
  const ImVec2 c = toScreen(tf, 0.0, hf * geo.heightMetres);
  const float rx = static_cast<float>(halfW) * tf.scale;
  const float ry = std::max(2.5f, rx * 0.16f);
  const ImVec4 fillCol(std::min(phase.colour[0] * 1.25f, 1.0f),
                       std::min(phase.colour[1] * 1.25f, 1.0f),
                       std::min(phase.colour[2] * 1.25f, 1.0f), 0.22f);
  draw->AddEllipseFilled(c, ImVec2(rx, ry), ImGui::ColorConvertFloat4ToU32(fillCol));
  draw->AddEllipse(c, ImVec2(rx, ry), style::u32(style::col::Text, 0.45f), 0.0f, 0, 1.2f);
}

void drawLayers(ImDrawList* draw, const sol::Simulation& sim, const VesselGeometry& geo,
                const Transform& tf) {
  static thread_local std::vector<ImVec2> slice;
  const size_t n = std::min(sim.phases.size(), sim.settledMl.size());
  constexpr int kRows = 14;    // height samples per slice
  constexpr int kSlices = 26;  // vertical slices across the width

  double cursorMl = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const double volume = std::max(sim.settledMl[i], 0.0);
    if (volume <= 1e-9) continue;

    const double hfLo = heightFractionForVolume(geo, cursorMl);
    const double hfHi = heightFractionForVolume(geo, cursorMl + volume);
    cursorMl += volume;
    if (hfHi <= hfLo + 1e-6) continue;

    const sol::Phase& phase = sim.phases[i];

    // Vertical slices with cylindrical shading; each slice follows the wall
    // curvature, so the band is assembled exactly to the vessel silhouette.
    for (int j = 0; j < kSlices; ++j) {
      const float f0 = -1.0f + 2.0f * static_cast<float>(j) / kSlices;
      const float f1 = -1.0f + 2.0f * static_cast<float>(j + 1) / kSlices;
      const float fm = 0.5f * (f0 + f1);

      slice.clear();
      for (int s = 0; s <= kRows; ++s) {
        const double hf = hfLo + (hfHi - hfLo) * s / kRows;
        const double halfW = widthFractionAt(sim.vessel, hf) * geo.halfWidthMetres;
        slice.push_back(toScreen(tf, f1 * halfW, hf * geo.heightMetres));
      }
      for (int s = kRows; s >= 0; --s) {
        const double hf = hfLo + (hfHi - hfLo) * s / kRows;
        const double halfW = widthFractionAt(sim.vessel, hf) * geo.halfWidthMetres;
        slice.push_back(toScreen(tf, f0 * halfW, hf * geo.heightMetres));
      }
      if (slice.size() >= 3) {
        draw->AddConvexPolyFilled(slice.data(), static_cast<int>(slice.size()),
                                  phaseShade(phase, fm));
      }
    }

    // Meniscus at the top of every non-empty layer; the interface between two
    // settled layers reads through the two stacked ellipses.
    drawMeniscus(draw, phase, hfHi, sim, geo, tf);
  }
}

void drawDroplets(ImDrawList* draw, const sol::Simulation& sim, const Transform& tf) {
  const int phaseCount = static_cast<int>(sim.phases.size());
  for (const sol::Droplet& droplet : sim.droplets) {
    if (droplet.phase < 0 || droplet.phase >= phaseCount) continue;
    const ImVec2 center = toScreen(tf, droplet.position);
    const float radius = std::max(droplet.radius * tf.scale, 1.0f);
    const sol::Phase& phase = sim.phases[static_cast<size_t>(droplet.phase)];
    const float alpha = std::min(phase.colour[3] + 0.15f, 1.0f);
    const ImVec4 dropletColour(phase.colour[0], phase.colour[1], phase.colour[2], alpha);
    const ImU32 fill = ImGui::ColorConvertFloat4ToU32(dropletColour);
    draw->AddCircleFilled(center, radius, fill, 8);
    if (radius > 2.0f) {
      // Shaded sphere: dark rim bottom-right, specular dot top-left.
      draw->AddCircle(center, radius, style::u32(style::col::BgDeep, 0.30f), 8,
                      std::max(1.0f, radius * 0.18f));
      const ImVec2 highlight(center.x - radius * 0.35f, center.y - radius * 0.35f);
      draw->AddCircleFilled(highlight, std::max(radius * 0.25f, 0.6f),
                            IM_COL32(255, 255, 255, 90), 6);
    }
  }
}

void drawCrossSection(const sol::Simulation& sim, ImVec2 regionMin, ImVec2 regionSize) {
  ImDrawList* draw = ImGui::GetWindowDrawList();
  const ImVec2 regionMax(regionMin.x + regionSize.x, regionMin.y + regionSize.y);
  draw->PushClipRect(regionMin, regionMax, true);
  draw->AddRectFilled(regionMin, regionMax, style::u32(style::col::BgSurface));

  const VesselGeometry& geo = cachedGeometry(sim.vessel, sim.vesselVolumeMl);
  Transform tf = buildTransform(geo, regionMin, regionSize);

  // Vessel motion during a shake is real motion: oscillate horizontally at
  // the shake frequency with the stroke amplitude (soft-capped in pixels so
  // an extreme stroke never throws the vessel out of the panel). Layers,
  // droplets, glass and graduations all share `tf`, so they move as one.
  if (sim.shake.active) {
    constexpr float kTwoPi = 6.28318530718f;
    const float ampPx = std::min(static_cast<float>(sim.shake.amplitudeM) * tf.scale * 0.35f,
                                 regionSize.x * 0.10f);
    tf.origin.x += ampPx * std::sin(kTwoPi * static_cast<float>(sim.shake.frequencyHz) *
                                    static_cast<float>(sim.elapsed));
  }

  drawGroundShadow(draw, geo, tf);
  drawVesselGlass(draw, geo, tf);
  drawGraduation(draw, geo, tf, regionMin);
  drawLayers(draw, sim, geo, tf);
  drawGlassHighlights(draw, sim, geo, tf);
  drawVesselWall(draw, geo, tf);
  drawFurniture(draw, sim, geo, tf);
  drawDroplets(draw, sim, tf);
  drawReadout(draw, sim, regionMin, regionSize);

  draw->PopClipRect();
}

// ---------------------------------------------------------------- editing
void seedDefaultPhases(sol::Simulation& sim) {
  if (!sim.phases.empty()) return;

  sol::Phase aqueous;
  aqueous.label = "Aqueous";
  aqueous.volumeMl = 100.0;
  aqueous.density = 1.00;
  aqueous.viscosity = 0.89;
  aqueous.interfacialTension = 30.0;
  aqueous.emulsionStability = 0.15;
  aqueous.colour[0] = 0.55f;
  aqueous.colour[1] = 0.72f;
  aqueous.colour[2] = 0.90f;
  aqueous.colour[3] = 0.62f;

  // Dichloromethane is denser than water, so it settles as the bottom
  // layer -- the classic gotcha this suite is meant to make obvious.
  sol::Phase organic;
  organic.label = "Dichloromethane";
  organic.volumeMl = 100.0;
  organic.density = 1.33;
  organic.viscosity = 0.41;
  organic.interfacialTension = 28.0;
  organic.emulsionStability = 0.35;
  organic.colour[0] = 0.86f;
  organic.colour[1] = 0.62f;
  organic.colour[2] = 0.28f;
  organic.colour[3] = 0.68f;

  sim.phases = {aqueous, organic};
  sol::reset(sim);
}

int resizeStringInput(ImGuiInputTextCallbackData* data) {
  if (data->EventFlag != ImGuiInputTextFlags_CallbackResize) return 0;
  auto* value = static_cast<std::string*>(data->UserData);
  value->resize(static_cast<size_t>(data->BufTextLen));
  data->Buf = value->data();
  return 0;
}

bool phaseLabelInput(const char* id, std::string& value) {
  return ImGui::InputText(id, value.data(), value.capacity() + 1,
                          ImGuiInputTextFlags_CallbackResize, resizeStringInput, &value);
}

bool animatedShakeButton(bool shaking, ImVec2 size) {
  const bool clicked = widgets::primaryButton("Shake", size);
  const ImGuiID itemId = ImGui::GetItemID();
  const bool hovered = ImGui::IsItemHovered();
  const bool held = ImGui::IsItemActive();
  const float hover = widgets::hoverT(itemId, hovered || held);
  const float press = widgets::hoverT(ImGui::GetID("##shake_press_anim"), held);
  const ImVec2 min = ImGui::GetItemRectMin();
  const ImVec2 max = ImGui::GetItemRectMax();
  const style::Metrics& m = style::metrics();
  ImDrawList* draw = ImGui::GetWindowDrawList();

  const float inset = press * std::max(m.hairline, 1.0f);
  draw->AddRect(ImVec2(min.x + inset, min.y + inset), ImVec2(max.x - inset, max.y - inset),
                style::mix(style::col::AccentActive, style::col::AccentHover, hover), m.radiusMd,
                0, m.hairline * (1.0f + hover + press));
  if (shaking) {
    const float pulse =
        0.34f + 0.28f * (0.5f + 0.5f * std::sin(static_cast<float>(ImGui::GetTime()) * 5.0f));
    draw->AddRect(min, max, style::u32(style::col::AccentHover, pulse), m.radiusMd, 0,
                  m.hairline * 2.0f);
  }
  return clicked;
}

float labelledControlWidth(const char* label) {
  const float labelWidth = ImGui::CalcTextSize(label).x + ImGui::GetStyle().ItemInnerSpacing.x;
  return std::max(1.0f, std::min(138.0f, ImGui::GetContentRegionAvail().x - labelWidth));
}

void drawTransportControls(SolubilityState& s) {
  sol::Simulation& sim = s.funnel;
  static const char* kVesselNames[] = {"Separatory funnel", "Decanting flask",
                                       "Graduated cylinder"};
  const float width = ImGui::GetContentRegionAvail().x;
  const int columns = width >= 620.0f ? 3 : 1;
  constexpr ImGuiTableFlags flags =
      ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_NoSavedSettings;
  if (!ImGui::BeginTable("##transport_grid", columns, flags)) return;
  if (columns == 3) {
    ImGui::TableSetupColumn("##vessel_col", ImGuiTableColumnFlags_WidthStretch, 1.55f);
    ImGui::TableSetupColumn("##transport_col", ImGuiTableColumnFlags_WidthStretch, 1.05f);
    ImGui::TableSetupColumn("##speed_col", ImGuiTableColumnFlags_WidthStretch, 1.10f);
  }

  ImGui::TableNextColumn();
  ImGui::TextDisabled("VESSEL");
  ImGui::SetNextItemWidth(-1.0f);
  if (ImGui::Combo("##vessel", &s.funnelVessel, kVesselNames, IM_ARRAYSIZE(kVesselNames))) {
    s.funnelVessel = std::clamp(s.funnelVessel, 0, 2);
    sim.vessel = static_cast<sol::Vessel>(s.funnelVessel);
    sol::reset(sim);
    s.statusMessage = std::string("Recharged into ") + kVesselNames[s.funnelVessel];
  }

  ImGui::TableNextColumn();
  ImGui::TextDisabled("TRANSPORT");
  const float buttonsWidth = ImGui::GetContentRegionAvail().x;
  const float buttonWidth = std::max((buttonsWidth - style::metrics().gap) * 0.5f, 1.0f);
  if (widgets::ghostButton(s.funnelRunning ? "Pause" : "Run", ImVec2(buttonWidth, 0.0f)))
    s.funnelRunning = !s.funnelRunning;
  ImGui::SameLine(0.0f, style::metrics().gap);
  if (widgets::ghostButton("Reset", ImVec2(buttonWidth, 0.0f))) {
    sol::reset(sim);
    s.funnelRunning = false;
    s.statusMessage = "Vessel reset";
  }

  ImGui::TableNextColumn();
  ImGui::TextDisabled("SIMULATION SPEED");
  ImGui::SetNextItemWidth(labelledControlWidth("Speed"));
  ImGui::SliderFloat("Speed", &s.funnelSpeed, 0.1f, 10.0f, "%.1fx");
  ImGui::EndTable();
}

void drawDerivedReadout(const sol::Simulation& sim) {
  if (sim.shake.durationS <= 0.0) return;

  const double dUm = sim.shake.sauterRadiusM * 2.0e6;  // m -> um diameter
  const double weber = 1000.0 * sim.shake.peakVelocity * sim.shake.peakVelocity *
                       sim.shake.sauterRadiusM * 2.0 /
                       std::max(1e-6, sim.phases.empty() ? 0.03
                                                         : sim.phases[0].interfacialTension * 1e-3);
  char line1[160];
  char line2[160];
  std::snprintf(line1, sizeof(line1), "DERIVED   u %.2f m/s   epsilon %.1f W/kg",
                sim.shake.peakVelocity, sim.shake.specificPower);
  std::snprintf(line2, sizeof(line2), "d32 %.0f um   Weber %.1f", dUm, weber);

  const float width = ImGui::GetContentRegionAvail().x;
  if (width <= 0.0f) return;
  const style::Metrics& m = style::metrics();
  const float lineHeight = ImGui::GetTextLineHeight();
  const bool singleLine =
      ImGui::CalcTextSize(line1).x + ImGui::CalcTextSize(line2).x + m.gap * 3.0f < width;
  const ImVec2 size(width, lineHeight * (singleLine ? 1.0f : 2.0f) + m.gap * 1.5f);
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::InvisibleButton("##derived_readout", size);
  const float hover = widgets::hoverT(ImGui::GetItemID(), ImGui::IsItemHovered());
  const ImVec2 max(min.x + size.x, min.y + size.y);
  ImDrawList* draw = ImGui::GetWindowDrawList();
  draw->AddRectFilled(min, max,
                      style::mix(style::col::BgDeep, style::col::BgRaised, hover, 0.72f),
                      m.radiusSm);
  draw->AddRect(min, max, style::mix(style::col::Border, style::col::Teal, hover), m.radiusSm, 0,
                m.hairline);

  const bool mono = style::pushFont(style::fonts::mono());
  const float textY = min.y + m.gap * 0.70f;
  draw->AddText(ImVec2(min.x + m.gap, textY), style::u32(style::col::TextDim), line1);
  if (singleLine) {
    const float x = min.x + m.gap * 2.0f + ImGui::CalcTextSize(line1).x;
    draw->AddText(ImVec2(x, textY), style::u32(style::col::Teal), line2);
  } else {
    draw->AddText(ImVec2(min.x + m.gap, textY + lineHeight), style::u32(style::col::Teal),
                  line2);
  }
  style::popFont(mono);

  if (ImGui::IsItemHovered()) {
    ImGui::BeginTooltip();
    const bool tooltipMono = style::pushFont(style::fonts::mono());
    ImGui::TextUnformatted("u = 2*pi*f*A");
    ImGui::TextUnformatted("epsilon = u^2*f/2");
    ImGui::TextUnformatted("d32 = 0.725*(sigma/rho_c)^0.6*epsilon^-0.4");
    ImGui::TextUnformatted("d32 capped by We_crit = 12 (Hinze 1955)");
    ImGui::TextUnformatted("turnovers/s = u/H");
    style::popFont(tooltipMono);
    ImGui::EndTooltip();
  }
}

void drawShakeControls(SolubilityState& s) {
  sol::Simulation& sim = s.funnel;
  const float width = ImGui::GetContentRegionAvail().x;
  const int columns = width >= 700.0f ? 4 : (width >= 400.0f ? 2 : 1);
  constexpr ImGuiTableFlags flags =
      ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_NoSavedSettings;
  if (ImGui::BeginTable("##shake_grid", columns, flags)) {
    ImGui::TableNextColumn();
    ImGui::TextDisabled("ACTION");
    if (animatedShakeButton(sim.shake.active,
                            ImVec2(ImGui::GetContentRegionAvail().x, ImGui::GetFrameHeight()))) {
      sol::ShakeParams params;
      params.durationS = static_cast<double>(s.shakeDurationS);
      params.frequencyHz = static_cast<double>(s.shakeFrequencyHz);
      params.amplitudeM = static_cast<double>(s.shakeAmplitudeCm) * 0.01;
      sol::shake(sim, params);
      s.funnelRunning = true;  // a shake is motion: the clock must run
      char buf[160];
      std::snprintf(buf, sizeof(buf),
                    "Shaking %.0f s at %.1f Hz, %.0f cm: u = %.2f m/s, epsilon = %.1f W/kg",
                    params.durationS, params.frequencyHz, params.amplitudeM * 100.0,
                    sim.shake.peakVelocity, sim.shake.specificPower);
      s.statusMessage = buf;
    }

    ImGui::TableNextColumn();
    ImGui::TextDisabled("DURATION");
    ImGui::SetNextItemWidth(labelledControlWidth("Duration"));
    ImGui::DragFloat("Duration##shake", &s.shakeDurationS, 0.1f, 1.0f, 30.0f, "%.0f s");

    ImGui::TableNextColumn();
    ImGui::TextDisabled("FREQUENCY");
    ImGui::SetNextItemWidth(labelledControlWidth("Frequency"));
    ImGui::DragFloat("Frequency##shake", &s.shakeFrequencyHz, 0.05f, 0.5f, 6.0f, "%.1f Hz");

    ImGui::TableNextColumn();
    ImGui::TextDisabled("AMPLITUDE");
    ImGui::SetNextItemWidth(labelledControlWidth("Amplitude"));
    ImGui::DragFloat("Amplitude##shake", &s.shakeAmplitudeCm, 0.1f, 1.0f, 15.0f, "%.0f cm");
    ImGui::EndTable();
  }

  if (sim.shake.durationS > 0.0) {
    ImGui::Spacing();
    drawDerivedReadout(sim);
  }
}

void drawControls(SolubilityState& s) {
  if (widgets::beginCard("##transport_card", ImVec2(0.0f, 0.0f), style::col::BgSurface)) {
    widgets::sectionHeader("VESSEL & TRANSPORT", style::col::Teal);
    drawTransportControls(s);
    widgets::endCard();
  }

  ImGui::Spacing();
  if (widgets::beginCard("##shake_card", ImVec2(0.0f, 0.0f), style::col::BgSurface)) {
    widgets::sectionHeader("SHAKE", style::col::Accent);
    drawShakeControls(s);
    widgets::endCard();
  }
}

void drawPhaseTable(SolubilityState& s, bool& changed) {
  sol::Simulation& sim = s.funnel;
  constexpr ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                                    ImGuiTableFlags_SizingStretchProp |
                                    ImGuiTableFlags_NoSavedSettings | ImGuiTableFlags_PadOuterX;
  if (!ImGui::BeginTable("##phase_table", 8, flags)) return;
  ImGui::TableSetupColumn("Label", ImGuiTableColumnFlags_WidthStretch, 1.80f);
  ImGui::TableSetupColumn("mL", ImGuiTableColumnFlags_WidthStretch, 0.85f);
  ImGui::TableSetupColumn("g/mL", ImGuiTableColumnFlags_WidthStretch, 0.90f);
  ImGui::TableSetupColumn("mPa.s", ImGuiTableColumnFlags_WidthStretch, 0.90f);
  ImGui::TableSetupColumn("mN/m", ImGuiTableColumnFlags_WidthStretch, 0.90f);
  ImGui::TableSetupColumn("Stability", ImGuiTableColumnFlags_WidthStretch, 1.15f);
  ImGui::TableSetupColumn("Colour", ImGuiTableColumnFlags_WidthFixed, 46.0f);
  ImGui::TableSetupColumn("##remove", ImGuiTableColumnFlags_WidthFixed, 32.0f);
  ImGui::TableHeadersRow();

  int removeIndex = -1;
  for (size_t i = 0; i < sim.phases.size(); ++i) {
    sol::Phase& phase = sim.phases[i];
    ImGui::PushID(static_cast<int>(i));
    ImGui::TableNextRow();

    ImGui::TableNextColumn();
    ImGui::SetNextItemWidth(-1.0f);
    if (phaseLabelInput("##label", phase.label)) changed = true;

    ImGui::TableNextColumn();
    ImGui::SetNextItemWidth(-1.0f);
    {
      float v = static_cast<float>(phase.volumeMl);
      if (ImGui::DragFloat("##vol", &v, 1.0f, 0.0f, 5000.0f, "%.1f")) {
        phase.volumeMl = std::max(v, 0.0f);
        changed = true;
      }
    }

    ImGui::TableNextColumn();
    ImGui::SetNextItemWidth(-1.0f);
    {
      float v = static_cast<float>(phase.density);
      if (ImGui::DragFloat("##density", &v, 0.005f, 0.10f, 3.50f, "%.3f")) {
        phase.density = std::max(v, 0.01f);
        changed = true;
      }
    }

    ImGui::TableNextColumn();
    ImGui::SetNextItemWidth(-1.0f);
    {
      float v = static_cast<float>(phase.viscosity);
      if (ImGui::DragFloat("##visc", &v, 0.01f, 0.05f, 500.0f, "%.2f")) {
        phase.viscosity = std::max(v, 0.01f);
        changed = true;
      }
    }

    ImGui::TableNextColumn();
    ImGui::SetNextItemWidth(-1.0f);
    {
      float v = static_cast<float>(phase.interfacialTension);
      if (ImGui::DragFloat("##ift", &v, 0.5f, 0.0f, 100.0f, "%.1f")) {
        phase.interfacialTension = std::max(v, 0.0f);
        changed = true;
      }
    }

    ImGui::TableNextColumn();
    ImGui::SetNextItemWidth(-1.0f);
    {
      float v = static_cast<float>(phase.emulsionStability);
      if (ImGui::SliderFloat("##stab", &v, 0.0f, 1.0f, "%.2f")) {
        phase.emulsionStability = v;
        changed = true;
      }
    }

    ImGui::TableNextColumn();
    if (ImGui::ColorEdit4("##colour", phase.colour,
                          ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_AlphaBar)) {
      changed = true;
    }

    ImGui::TableNextColumn();
    if (widgets::ghostButton("x", ImVec2(24.0f, 0.0f))) removeIndex = static_cast<int>(i);

    ImGui::PopID();
  }
  ImGui::EndTable();

  if (removeIndex >= 0) {
    sim.phases.erase(sim.phases.begin() + removeIndex);
    changed = true;
  }
}

void drawPhaseEditor(SolubilityState& s) {
  bool changed = false;
  drawPhaseTable(s, changed);

  if (widgets::ghostButton("+ Add phase")) {
    sol::Phase phase;
    phase.label = "Phase " + std::to_string(s.funnel.phases.size() + 1);
    phase.volumeMl = 50.0;
    phase.density = 1.0;
    phase.viscosity = 1.0;
    phase.interfacialTension = 30.0;
    phase.emulsionStability = 0.3;
    phase.colour[0] = 0.40f;
    phase.colour[1] = 0.75f;
    phase.colour[2] = 0.65f;
    phase.colour[3] = 0.65f;
    s.funnel.phases.push_back(phase);
    changed = true;
  }
  const float noteWidth = ImGui::CalcTextSize("Densest phase settles to the bottom").x;
  if (ImGui::GetContentRegionAvail().x > noteWidth + style::metrics().gap) {
    ImGui::SameLine(0.0f, style::metrics().gap);
  }
  ImGui::TextDisabled("Densest phase settles to the bottom");

  if (changed) {
    sol::reset(s.funnel);
    const double total = sol::totalVolumeMl(s.funnel);
    char buf[128];
    std::snprintf(buf, sizeof(buf), "Recharged vessel: %d phase(s), %.0f mL total",
                 static_cast<int>(s.funnel.phases.size()), total);
    s.statusMessage = buf;
  }
}

void stepSimulation(SolubilityState& s) {
  if (!s.funnelRunning) return;
  // Clamp so a stalled/hitched frame cannot hand the physics a huge dt.
  const float dt = std::clamp(ImGui::GetIO().DeltaTime * s.funnelSpeed, 0.0f, 0.1f);
  if (dt <= 0.0f) return;
  sol::step(s.funnel, static_cast<double>(dt));
}

// -------------------------------------------------------- suite hand-off
// UI-only mass for the solute distribution panel, in mg. Seeded from the
// Solubility Suite's import payload; the slider owns it afterwards.
float soluteMassMgUi = 100.0f;

void colourForFamily(const std::string& family, float out[4]) {
  // Muted lab-liquid tints per solvent family; alpha stays low enough that
  // the vessel wall and interface lines still read through the layer.
  if (family == "water" || family == "alcohol") {
    out[0] = 0.42f; out[1] = 0.70f; out[2] = 0.88f;
  } else if (family == "halogenated") {
    out[0] = 0.52f; out[1] = 0.84f; out[2] = 0.58f;
  } else if (family == "alkane" || family == "aromatic") {
    out[0] = 0.88f; out[1] = 0.68f; out[2] = 0.32f;
  } else {
    out[0] = 0.70f; out[1] = 0.60f; out[2] = 0.88f;
  }
  out[3] = 0.55f;
}

sol::Phase makeImportedPhase(const sol::Solvent* solvent, double volumeMl) {
  sol::Phase phase;
  phase.label = solvent->name;
  phase.volumeMl = std::max(volumeMl, 1.0);
  phase.density = solvent->density;
  // The DB carries no viscosities; water is 0.89 mPa.s, everything else gets
  // a light-organic default close to the stock DCM phase.
  phase.viscosity = solvent->family == "water" ? 0.89 : 0.60;
  phase.interfacialTension = 30.0;
  phase.emulsionStability = 0.35;
  colourForFamily(solvent->family, phase.colour);
  return phase;
}

// The suite's "Send to Extraction Lab" button stages solvent ids + volumes;
// applied here, exactly once, on the frame after it arrives.
void consumeExtractionImport(SolubilityState& s) {
  ExtractionImport& imp = s.extractionImport;
  if (!imp.pending) return;
  imp.pending = false;  // consume-once, even when the payload is unusable

  const sol::Solvent* a = sol::findSolvent(imp.solventIdA);
  const sol::Solvent* b = sol::findSolvent(imp.solventIdB);
  if (!a || !b) {
    s.statusMessage = "Extraction import failed: unknown solvent id";
    return;
  }

  s.funnel.phases = {makeImportedPhase(a, imp.volumeMlA), makeImportedPhase(b, imp.volumeMlB)};
  sol::reset(s.funnel);
  s.funnelRunning = false;
  soluteMassMgUi = static_cast<float>(imp.soluteMassMg);
  s.statusMessage = "Imported " + a->name + " + " + b->name + " from the Solubility Suite";
}

// ------------------------------------------------------ solute distribution
// Predicts how the suite's solute splits between the aqueous and organic
// layers from logP, so the extraction is quantitative rather than just
// visual.
void drawSoluteDistribution(const SolubilityState& s) {
  if (!s.soluteValid || s.funnel.phases.size() < 2) {
    ImGui::TextDisabled("A valid solute and at least two charged phases are required.");
    return;
  }
  const sol::Simulation& sim = s.funnel;
  const size_t count = sim.phases.size();

  ImGui::TextDisabled("%s  ·  logP %.2f", s.solute.name.c_str(), s.solute.logP);

  // Aqueous phase picker: default the phase that IS water when one is
  // present; otherwise the least dense phase, which is where water sits in
  // a normal water/organic pair (halogenated solvents sink).
  static int aqueousPick = -1;  // -1 = auto
  int autoIndex = 0;
  for (size_t i = 0; i < count; ++i) {
    if (sim.phases[i].label.find("ater") != std::string::npos) {  // "Water"/"water"
      autoIndex = static_cast<int>(i);
      break;
    }
    if (sim.phases[i].density < sim.phases[static_cast<size_t>(autoIndex)].density)
      autoIndex = static_cast<int>(i);
  }
  const int aq = (aqueousPick >= 0 && static_cast<size_t>(aqueousPick) < count)
                     ? aqueousPick
                     : autoIndex;
  const int controlColumns = ImGui::GetContentRegionAvail().x >= 520.0f ? 2 : 1;
  constexpr ImGuiTableFlags controlFlags =
      ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_NoSavedSettings;
  if (ImGui::BeginTable("##distribution_controls", controlColumns, controlFlags)) {
    ImGui::TableNextColumn();
    const float aqueousLabelWidth =
        ImGui::CalcTextSize("Aqueous phase").x + ImGui::GetStyle().ItemInnerSpacing.x;
    ImGui::SetNextItemWidth(std::max(
        1.0f, std::min(220.0f, ImGui::GetContentRegionAvail().x - aqueousLabelWidth)));
    if (ImGui::BeginCombo("Aqueous phase",
                          sim.phases[static_cast<size_t>(aq)].label.c_str())) {
      for (size_t i = 0; i < count; ++i) {
        const bool selected = static_cast<int>(i) == aq;
        if (ImGui::Selectable(sim.phases[i].label.c_str(), selected))
          aqueousPick = static_cast<int>(i);
        if (selected) ImGui::SetItemDefaultFocus();
      }
      ImGui::EndCombo();
    }

    ImGui::TableNextColumn();
    const float massLabelWidth =
        ImGui::CalcTextSize("Solute mass").x + ImGui::GetStyle().ItemInnerSpacing.x;
    ImGui::SetNextItemWidth(
        std::max(1.0f, std::min(220.0f, ImGui::GetContentRegionAvail().x - massLabelWidth)));
    ImGui::SliderFloat("Solute mass", &soluteMassMgUi, 1.0f, 1000.0f, "%.0f mg");
    ImGui::EndTable();
  }

  // Everything that is not the aqueous phase counts as the organic side; in
  // the common two-phase case this is exactly the partner layer.
  double volAq = sim.phases[static_cast<size_t>(aq)].volumeMl;
  double volOrg = 0.0;
  std::string organicLabel;
  for (size_t i = 0; i < count; ++i) {
    if (static_cast<int>(i) == aq) continue;
    volOrg += sim.phases[i].volumeMl;
    if (organicLabel.empty()) organicLabel = sim.phases[i].label;
  }

  const sol::Partition p =
      sol::partition(static_cast<double>(soluteMassMgUi), s.solute.logP, volAq, volOrg);

  const float avail = ImGui::GetContentRegionAvail().x;
  if (avail <= 0.0f) return;
  const float barH = ImGui::GetFontSize() * 1.5f;
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::Dummy(ImVec2(avail, barH));
  const ImVec2 max(min.x + avail, min.y + barH);
  const float fracAq = static_cast<float>(p.mgAqueous / std::max(p.mgAqueous + p.mgOrganic, 1e-12));
  const float splitX = min.x + avail * fracAq;
  const style::Metrics& m = style::metrics();
  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddRectFilled(min, ImVec2(splitX, max.y), style::u32(style::col::Teal, 0.85f), m.radiusSm,
                    ImDrawFlags_RoundCornersLeft);
  dl->AddRectFilled(ImVec2(splitX, min.y), max, style::u32(style::col::Accent, 0.85f), m.radiusSm,
                    ImDrawFlags_RoundCornersRight);
  dl->AddRect(min, max, style::u32(style::col::BorderStrong), m.radiusSm, 0, m.hairline);

  const ImU32 onBar = style::u32(style::col::OnAccent);
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.3g mg", p.mgAqueous);
  ImVec2 textSize = ImGui::CalcTextSize(buf);
  if (fracAq * avail > textSize.x + 8.0f) {
    dl->AddText(ImVec2(min.x + 6.0f, min.y + (barH - textSize.y) * 0.5f), onBar, buf);
  }
  std::snprintf(buf, sizeof(buf), "%.3g mg", p.mgOrganic);
  textSize = ImGui::CalcTextSize(buf);
  if ((1.0f - fracAq) * avail > textSize.x + 8.0f) {
    dl->AddText(ImVec2(max.x - textSize.x - 6.0f, min.y + (barH - textSize.y) * 0.5f), onBar,
                buf);
  }

  ImGui::TextDisabled("%.1f%% extracted into %s · Neutral-species logP approximation (no pH "
                      "correction)",
                      p.fractionOrganic * 100.0, organicLabel.c_str());
}

}  // namespace

void drawExtractionLab(AppState& st) {
  SolubilityState& s = st.solubility;
  consumeExtractionImport(s);
  seedDefaultPhases(s.funnel);

  drawControls(s);
  ImGui::Spacing();

  if (widgets::beginCard("##phases_card", ImVec2(0.0f, 0.0f), style::col::BgSurface)) {
    widgets::sectionHeader("CHARGED PHASES", style::col::Violet);
    drawPhaseEditor(s);
    widgets::endCard();
  }
  ImGui::Spacing();

  if (widgets::beginCard("##distribution_card", ImVec2(0.0f, 0.0f),
                         style::col::BgSurface)) {
    widgets::sectionHeader("SOLUTE DISTRIBUTION", style::col::Teal);
    drawSoluteDistribution(s);
    widgets::endCard();
  }
  ImGui::Spacing();

  stepSimulation(s);

  const ImVec2 remaining = ImGui::GetContentRegionAvail();
  if (remaining.x <= 0.0f || remaining.y <= 0.0f) return;
  if (widgets::beginCard("##cross_section_card", remaining, style::col::BgSurface,
                         ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    widgets::sectionHeader("CROSS-SECTION", style::col::Accent);
    const ImVec2 canvasSize = ImGui::GetContentRegionAvail();
    if (canvasSize.x > 0.0f && canvasSize.y > 0.0f) {
      ImGui::InvisibleButton("##funnel_canvas", canvasSize, ImGuiButtonFlags_None);
      const ImVec2 rectMin = ImGui::GetItemRectMin();
      const ImVec2 rectMax = ImGui::GetItemRectMax();
      drawCrossSection(s.funnel, rectMin,
                       ImVec2(rectMax.x - rectMin.x, rectMax.y - rectMin.y));
    }
    widgets::endCard();
  }
}

}  // namespace chemcad::ui

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
constexpr double kPi = 3.14159265358979323846;
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

VesselGeometry buildGeometry(const sol::Simulation& sim) {
  VesselGeometry geo;
  geo.heightMetres = static_cast<float>(sol::columnHeightM(sim));
  geo.outline = sol::vesselOutline(sim.vessel, static_cast<double>(geo.heightMetres));
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
  double prevArea = widthFractionAt(sim.vessel, 0.0);
  prevArea *= prevArea;
  for (int i = 1; i < kWidthSamples; ++i) {
    const double hf = static_cast<double>(i) * dHf;
    const double w = widthFractionAt(sim.vessel, hf);
    const double area = w * w;
    cumRaw[i] = cumRaw[i - 1] + 0.5 * (prevArea + area) * dHf;
    prevArea = area;
  }
  const double total = std::max(cumRaw[kWidthSamples - 1], 1e-9);
  const double capacityMl = std::max(sim.vesselVolumeMl, 1.0);
  for (int i = 0; i < kWidthSamples; ++i) geo.cumVolumeMl[i] = capacityMl * cumRaw[i] / total;
  return geo;
}

// Rebuilds only when the vessel type or charged capacity actually changes;
// everything else (layers, droplets) reuses this per-frame without
// re-sampling the width profile.
const VesselGeometry& cachedGeometry(const sol::Simulation& sim) {
  static thread_local sol::Vessel lastVessel = sol::Vessel::SeparatoryFunnel;
  static thread_local double lastVolumeMl = -1.0;
  static thread_local bool initialized = false;
  static thread_local VesselGeometry cache;
  if (!initialized || lastVessel != sim.vessel ||
      std::fabs(lastVolumeMl - sim.vesselVolumeMl) > 1e-6) {
    cache = buildGeometry(sim);
    lastVessel = sim.vessel;
    lastVolumeMl = sim.vesselVolumeMl;
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

double volumeForHeightFraction(const VesselGeometry& geo, double heightFraction) {
  const double scaled =
      std::clamp(heightFraction, 0.0, 1.0) * static_cast<double>(kWidthSamples - 1);
  const int lo = std::clamp(static_cast<int>(std::floor(scaled)), 0, kWidthSamples - 1);
  const int hi = std::min(lo + 1, kWidthSamples - 1);
  const double t = scaled - static_cast<double>(lo);
  return geo.cumVolumeMl[lo] + (geo.cumVolumeMl[hi] - geo.cumVolumeMl[lo]) * t;
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

std::string ellipsizeText(const std::string& text, float maxWidth) {
  if (maxWidth <= 0.0f) return {};
  if (ImGui::CalcTextSize(text.c_str()).x <= maxWidth) return text;
  constexpr const char* kEllipsis = "...";
  if (ImGui::CalcTextSize(kEllipsis).x > maxWidth) return {};

  std::string fitted = text;
  while (!fitted.empty()) {
    size_t cut = fitted.size() - 1;
    while (cut > 0 &&
           (static_cast<unsigned char>(fitted[cut]) & 0xC0u) == 0x80u) {
      --cut;
    }
    fitted.resize(cut);
    const std::string candidate = fitted + kEllipsis;
    if (ImGui::CalcTextSize(candidate.c_str()).x <= maxWidth) return candidate;
  }
  return kEllipsis;
}

ImVec2 textButtonSize(const char* label) {
  const ImVec2 text = ImGui::CalcTextSize(label);
  const ImVec2 padding = ImGui::GetStyle().FramePadding;
  return ImVec2(text.x + padding.x * 2.0f, ImGui::GetFrameHeight());
}

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
  double step = niceTickStep(capacityMl);
  if (step <= 0.0) return;
  const ImU32 tickColor = style::u32(style::col::TextFaint);
  const ImU32 labelColor = style::u32(style::col::TextDim);
  const float tickX = regionMin.x + kMarginLeft - 14.0f;
  const float labelH = ImGui::GetFontSize();
  // Widen the step until no adjacent pair of labels can collide vertically.
  // Every interval is checked: the vessel neck compresses the top ticks even
  // when the stem intervals are generous.
  for (int guard = 0; guard < 8; ++guard) {
    float minPx = 1e9f;
    float prevY = 0.0f;
    for (double v = 0.0; v <= capacityMl + 1e-6; v += step) {
      const double hf = heightFractionForVolume(geo, v);
      const float y = toScreen(tf, 0.0, hf * geo.heightMetres).y;
      if (v > 0.0) minPx = std::min(minPx, std::fabs(y - prevY));
      prevY = y;
    }
    if (minPx >= labelH * 1.25f) break;
    step *= 2.0;
  }
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
  const float minimumBoxWidth =
      std::max(ImGui::CalcTextSize("emulsified: 100.0%").x,
               ImGui::CalcTextSize("EMULSIFIED").x) +
      m.gap * 2.0f;
  if (availableWidth < minimumBoxWidth || regionSize.y < m.gap * 4.0f) return;

  const double fraction = sol::emulsifiedFraction(sim);
  const float lineHeight = ImGui::GetTextLineHeightWithSpacing();
  const float padding = m.gap;
  const float boxWidth =
      std::min(std::max(ImGui::GetFontSize() * 11.0f, minimumBoxWidth), availableWidth);
  const size_t layerCount = std::min(sim.phases.size(), sim.settledMl.size());
  const int fixedLines = sim.shake.active ? 4 : 3;
  const float availableHeight = std::max(regionSize.y - m.gap * 2.0f, 0.0f);
  if (availableHeight < padding * 2.0f + lineHeight * fixedLines) return;
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
  if (bigFont) {
    draw->AddText(bigFont, bigSize, cursor, emulsionStateColor(fraction),
                  emulsionStateLabel(fraction));
  } else {
    draw->AddText(cursor, emulsionStateColor(fraction), emulsionStateLabel(fraction));
  }
  cursor.y += lineHeight;

  for (size_t i = 0; i < visibleLayers; ++i) {
    if (hasHiddenLayers && i + 1 == visibleLayers) {
      std::snprintf(buf, sizeof(buf), "+ %d more phase%s",
                    static_cast<int>(layerCount - i), layerCount - i == 1 ? "" : "s");
      draw->AddText(cursor, style::u32(style::col::TextDim), buf);
    } else {
      char suffix[40];
      std::snprintf(suffix, sizeof(suffix), ": %.0f mL", sim.settledMl[i]);
      const float labelWidth = std::max(wrapWidth - ImGui::CalcTextSize(suffix).x, 0.0f);
      const std::string line = ellipsizeText(sim.phases[i].label, labelWidth) + suffix;
      draw->AddText(cursor, style::u32(style::col::Text), line.c_str());
    }
    cursor.y += lineHeight;
  }

  draw->PopClipRect();
}

// --------------------------------------------------------------- painting
// The vessel uses the same local-metre geometry as physics. Flat mode keeps
// textbook phase colours; shaded mode adds cylindrical light rolloff. Both
// share the sloshing boundaries, parcel haze, glass, furniture and highlights.

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

void buildScreenOutline(const VesselGeometry& geo, const Transform& tf,
                        std::vector<ImVec2>& screenOutline) {
  screenOutline.clear();
  screenOutline.reserve(geo.outline.size());
  for (const core::Vec2& p : geo.outline) {
    const ImVec2 screen = toScreen(tf, p);
    if (!screenOutline.empty()) {
      const ImVec2& prev = screenOutline.back();
      const float dx = screen.x - prev.x;
      const float dy = screen.y - prev.y;
      if (dx * dx + dy * dy < 0.04f) continue;
    }
    screenOutline.push_back(screen);
  }
  if (screenOutline.size() > 2) {
    const float dx = screenOutline.front().x - screenOutline.back().x;
    const float dy = screenOutline.front().y - screenOutline.back().y;
    if (dx * dx + dy * dy < 0.04f) screenOutline.pop_back();
  }
}

void drawVesselGlass(ImDrawList* draw, const VesselGeometry& geo, const Transform& tf) {
  static thread_local std::vector<ImVec2> screenOutline;
  buildScreenOutline(geo, tf, screenOutline);
  if (screenOutline.size() < 3) return;
  draw->AddConcavePolyFilled(screenOutline.data(), static_cast<int>(screenOutline.size()),
                             style::u32(ImVec4(0.60f, 0.76f, 0.88f, 1.0f), 0.10f));
}

void drawVesselWall(ImDrawList* draw, const VesselGeometry& geo, const Transform& tf) {
  static thread_local std::vector<ImVec2> screenOutline;
  buildScreenOutline(geo, tf, screenOutline);
  if (screenOutline.size() < 3) return;
  // A soft outer stroke under a crisp closed silhouette reads as glass
  // thickness without exposing the dense analytic sampling.
  draw->AddPolyline(screenOutline.data(), static_cast<int>(screenOutline.size()),
                    style::u32(ImVec4(0.70f, 0.82f, 0.92f, 1.0f), 0.28f),
                    ImDrawFlags_Closed, 4.0f);
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

// Stopper cap + knob at the neck, stopcock across the straight drain stem:
// both are kept clear of the funnel cone so the silhouette stays readable.
void drawFurniture(ImDrawList* draw, const sol::Simulation& sim, const VesselGeometry& geo,
                   const Transform& tf) {
  if (sim.vessel != sol::Vessel::SeparatoryFunnel) return;

  const ImU32 glassFill = style::u32(ImVec4(0.72f, 0.84f, 0.93f, 1.0f), 0.35f);
  const ImU32 glassEdge = style::u32(style::col::BorderStrong, 0.9f);

  // The analytic neck is 0.20 of the maximum half-width at t = 1.
  const double neckW = widthFractionAt(sim.vessel, 1.0) * geo.halfWidthMetres;
  const ImVec2 topC = toScreen(tf, 0.0, geo.heightMetres);
  const float capW = std::max(static_cast<float>(neckW) * tf.scale * 2.0f + 4.0f, 10.0f);
  const float capH = std::min(std::max(6.0f, tf.scale * 0.009f), kMarginTop * 0.48f);
  draw->AddRectFilled(ImVec2(topC.x - capW * 0.5f, topC.y - capH),
                      ImVec2(topC.x + capW * 0.5f, topC.y + 1.0f), glassFill, capH * 0.35f);
  draw->AddRect(ImVec2(topC.x - capW * 0.5f, topC.y - capH),
                ImVec2(topC.x + capW * 0.5f, topC.y + 1.0f), glassEdge, capH * 0.35f, 0, 1.2f);
  const float knobR = capH * 0.36f;
  const ImVec2 knobC(topC.x, topC.y - capH - knobR * 0.75f);
  draw->AddCircleFilled(knobC, knobR, glassFill, 24);
  draw->AddCircle(knobC, knobR, glassEdge, 24, 1.2f);

  // t = 0.24 is within the straight stem and below the body cone.
  constexpr double kStopcockT = 0.24;
  const double stemW = widthFractionAt(sim.vessel, kStopcockT) * geo.halfWidthMetres;
  const ImVec2 cockC = toScreen(tf, 0.0, kStopcockT * geo.heightMetres);
  const float stemWidthPx = static_cast<float>(stemW) * tf.scale;
  const float maxStemHeight = std::max(5.0f, tf.scale * geo.heightMetres * 0.030f);
  const float barH = std::min(std::max(6.0f, tf.scale * 0.008f), maxStemHeight);
  const float barW = std::max(stemWidthPx * 2.0f + 10.0f, 20.0f);
  draw->AddRectFilled(ImVec2(cockC.x - barW * 0.5f, cockC.y - barH * 0.5f),
                      ImVec2(cockC.x + barW * 0.5f, cockC.y + barH * 0.5f), glassFill,
                      barH * 0.5f);
  draw->AddRect(ImVec2(cockC.x - barW * 0.5f, cockC.y - barH * 0.5f),
                ImVec2(cockC.x + barW * 0.5f, cockC.y + barH * 0.5f), glassEdge,
                barH * 0.5f, 0, 1.2f);
  const ImVec2 tabA(cockC.x + barW * 0.5f, cockC.y);
  const ImVec2 tabB(tabA.x + barH * 1.35f, tabA.y + barH * 0.75f);
  draw->AddLine(tabA, tabB, glassEdge, 3.0f);
  draw->AddCircleFilled(tabB, barH * 0.38f, glassFill, 24);
  draw->AddCircle(tabB, barH * 0.38f, glassEdge, 24, 1.1f);
}

double halfWidthAtHeight(const sol::Simulation& sim, const VesselGeometry& geo,
                         double heightM) {
  const double hf = heightM / std::max(static_cast<double>(geo.heightMetres), 1e-9);
  return widthFractionAt(sim.vessel, hf) * geo.halfWidthMetres;
}

double surfaceHeightAt(const sol::Simulation& sim, const VesselGeometry& geo,
                       double baseHeightM, double clearanceM, double xM) {
  const double energy = std::clamp(sim.shakeEnergy, 0.0, 1.0);
  if (energy <= 1e-5 || clearanceM <= 1e-9) return baseHeightM;

  const double frequency = std::max(sim.shake.frequencyHz, 0.0);
  const double omega = 2.0 * kPi * frequency;
  const double halfW = std::max(halfWidthAtHeight(sim, geo, baseHeightM), 1e-5);
  const double tilt = 0.28 * energy * std::sin(omega * sim.elapsed);
  const double amplitude = 0.16 * energy * std::max(sim.shake.amplitudeM, 0.0);
  const double maxOffset = std::max(clearanceM, 0.0);
  const double rawExtent = std::fabs(tilt) * halfW + amplitude;
  const double scale = rawExtent > maxOffset ? maxOffset / rawExtent : 1.0;
  const double waveNumber = kPi / halfW;  // one wavelength across the local width
  return baseHeightM +
         scale * (tilt * xM + amplitude * std::sin(waveNumber * xM - omega * sim.elapsed));
}

double boundaryClearance(const std::vector<double>& boundaries, size_t boundary,
                         double vesselHeightM) {
  if (boundary == 0 || boundary >= boundaries.size()) return 0.0;
  double nearestBand = 1e9;
  if (boundary > 0 && boundaries[boundary] > boundaries[boundary - 1])
    nearestBand = std::min(nearestBand, boundaries[boundary] - boundaries[boundary - 1]);
  if (boundary + 1 < boundaries.size() && boundaries[boundary + 1] > boundaries[boundary])
    nearestBand = std::min(nearestBand, boundaries[boundary + 1] - boundaries[boundary]);
  if (nearestBand == 1e9) return 0.0;
  const double edgeRoom =
      std::min(boundaries[boundary], vesselHeightM - boundaries[boundary]);
  return std::max(0.0, std::min(nearestBand * 0.24, edgeRoom * 0.80));
}

ImVec2 surfacePoint(const sol::Simulation& sim, const VesselGeometry& geo, const Transform& tf,
                    double baseHeightM, double clearanceM, double xFraction) {
  const double baseHalfW = halfWidthAtHeight(sim, geo, baseHeightM);
  const double xForWave = xFraction * baseHalfW;
  const double y = std::clamp(surfaceHeightAt(sim, geo, baseHeightM, clearanceM, xForWave),
                              0.0, static_cast<double>(geo.heightMetres));
  const double x = xFraction * halfWidthAtHeight(sim, geo, y);
  return toScreen(tf, x, y);
}

void appendSurface(std::vector<ImVec2>& points, const sol::Simulation& sim,
                   const VesselGeometry& geo, const Transform& tf, double baseHeightM,
                   double clearanceM, int samples, bool reverse) {
  for (int i = 0; i <= samples; ++i) {
    const int sample = reverse ? samples - i : i;
    const double f = -1.0 + 2.0 * static_cast<double>(sample) / samples;
    points.push_back(surfacePoint(sim, geo, tf, baseHeightM, clearanceM, f));
  }
}

void drawLayers(ImDrawList* draw, const sol::Simulation& sim, const VesselGeometry& geo,
                const Transform& tf, bool shaded) {
  static thread_local std::vector<ImVec2> polygon;
  static thread_local std::vector<ImVec2> surface;
  static thread_local std::vector<double> boundaries;
  const size_t n = std::min(sim.phases.size(), sim.settledMl.size());
  boundaries.assign(n + 1, 0.0);

  double cursorMl = 0.0;
  size_t topCarrier = n;
  for (size_t i = 0; i < n; ++i) {
    const double settled = std::max(sim.settledMl[i], 0.0);
    cursorMl += settled;
    boundaries[i + 1] =
        heightFractionForVolume(geo, cursorMl) * static_cast<double>(geo.heightMetres);
    if (settled > 1e-9) topCarrier = i;
  }
  // Dispersion changes composition, not the free-surface volume. Let the
  // upper continuous band occupy the parcel volume; haze/clouds tint that
  // volume as emulsion while the top remains at the charged-volume height.
  if (n > 0) {
    if (topCarrier == n) topCarrier = n - 1;
    const double liquidTop =
        heightFractionForVolume(geo, sol::totalVolumeMl(sim)) * geo.heightMetres;
    boundaries[topCarrier + 1] = std::max(boundaries[topCarrier + 1], liquidTop);
    for (size_t i = topCarrier + 2; i < boundaries.size(); ++i)
      boundaries[i] = boundaries[topCarrier + 1];
  }

  constexpr int kSurfaceSamples = 48;
  constexpr int kShadeSlices = 36;
  for (size_t i = 0; i < n; ++i) {
    if (boundaries[i + 1] <= boundaries[i] + 1e-9) continue;
    const sol::Phase& phase = sim.phases[i];
    const double loClear = boundaryClearance(boundaries, i, geo.heightMetres);
    const double hiClear = boundaryClearance(boundaries, i + 1, geo.heightMetres);

    if (!shaded) {
      polygon.clear();
      appendSurface(polygon, sim, geo, tf, boundaries[i + 1], hiClear, kSurfaceSamples, false);
      appendSurface(polygon, sim, geo, tf, boundaries[i], loClear, kSurfaceSamples, true);
      draw->AddConcavePolyFilled(
          polygon.data(), static_cast<int>(polygon.size()),
          ImGui::ColorConvertFloat4ToU32(
              ImVec4(phase.colour[0], phase.colour[1], phase.colour[2], phase.colour[3])));
    } else {
      for (int slice = 0; slice < kShadeSlices; ++slice) {
        const double f0 = -1.0 + 2.0 * static_cast<double>(slice) / kShadeSlices;
        const double f1 = -1.0 + 2.0 * static_cast<double>(slice + 1) / kShadeSlices;
        ImVec2 quad[4] = {
            surfacePoint(sim, geo, tf, boundaries[i], loClear, f0),
            surfacePoint(sim, geo, tf, boundaries[i], loClear, f1),
            surfacePoint(sim, geo, tf, boundaries[i + 1], hiClear, f1),
            surfacePoint(sim, geo, tf, boundaries[i + 1], hiClear, f0),
        };
        draw->AddConvexPolyFilled(quad, 4, phaseShade(phase, static_cast<float>((f0 + f1) * 0.5)));
      }
    }

    surface.clear();
    appendSurface(surface, sim, geo, tf, boundaries[i + 1], hiClear, kSurfaceSamples, false);
    draw->AddPolyline(surface.data(), static_cast<int>(surface.size()),
                      style::u32(style::col::Text, shaded ? 0.36f : 0.62f),
                      ImDrawFlags_None, shaded ? 1.1f : 1.5f);
  }
}

void drawEmulsionHaze(ImDrawList* draw, const sol::Simulation& sim,
                      const VesselGeometry& geo, const Transform& tf) {
  constexpr int kBins = 48;
  static thread_local std::array<double, kBins> parcelMl;
  static thread_local std::array<double, kBins> red;
  static thread_local std::array<double, kBins> green;
  static thread_local std::array<double, kBins> blue;
  parcelMl.fill(0.0);
  red.fill(0.0);
  green.fill(0.0);
  blue.fill(0.0);

  const double liquidHeight =
      heightFractionForVolume(geo, sol::totalVolumeMl(sim)) * geo.heightMetres;
  if (liquidHeight <= 1e-9) return;
  for (const sol::Droplet& droplet : sim.droplets) {
    if (droplet.phase < 0 || static_cast<size_t>(droplet.phase) >= sim.phases.size() ||
        droplet.parcelMl <= 0.0f) {
      continue;
    }
    const int bin = std::clamp(
        static_cast<int>(droplet.position.y / liquidHeight * kBins), 0, kBins - 1);
    const double volume = droplet.parcelMl;
    const sol::Phase& phase = sim.phases[static_cast<size_t>(droplet.phase)];
    parcelMl[bin] += volume;
    red[bin] += volume * phase.colour[0];
    green[bin] += volume * phase.colour[1];
    blue[bin] += volume * phase.colour[2];
  }

  for (int bin = 0; bin < kBins; ++bin) {
    if (parcelMl[bin] <= 1e-9) continue;
    const double yLo = liquidHeight * static_cast<double>(bin) / kBins;
    const double yHi = liquidHeight * static_cast<double>(bin + 1) / kBins;
    const double hfLo = yLo / geo.heightMetres;
    const double hfHi = yHi / geo.heightMetres;
    const double capacity =
        std::max(volumeForHeightFraction(geo, hfHi) - volumeForHeightFraction(geo, hfLo), 1e-6);
    const float fraction = static_cast<float>(std::clamp(parcelMl[bin] / capacity, 0.0, 1.0));
    const float alpha = 0.035f + 0.16f * fraction;
    const double invVolume = 1.0 / parcelMl[bin];
    const ImU32 colour = ImGui::ColorConvertFloat4ToU32(
        ImVec4(static_cast<float>(red[bin] * invVolume),
               static_cast<float>(green[bin] * invVolume),
               static_cast<float>(blue[bin] * invVolume), alpha));
    const double wLo = halfWidthAtHeight(sim, geo, yLo);
    const double wHi = halfWidthAtHeight(sim, geo, yHi);
    ImVec2 band[4] = {toScreen(tf, -wLo, yLo), toScreen(tf, wLo, yLo),
                      toScreen(tf, wHi, yHi), toScreen(tf, -wHi, yHi)};
    draw->AddConvexPolyFilled(band, 4, colour);
  }
}

void drawParcelClouds(ImDrawList* draw, const sol::Simulation& sim,
                      const VesselGeometry& geo, const Transform& tf, bool shaded) {
  const int phaseCount = static_cast<int>(sim.phases.size());
  for (const sol::Droplet& droplet : sim.droplets) {
    if (droplet.phase < 0 || droplet.phase >= phaseCount || droplet.parcelMl <= 0.0f) continue;
    const double x = droplet.position.x;
    const double y = droplet.position.y;
    if (y <= 0.0 || y >= geo.heightMetres ||
        std::fabs(x) >= halfWidthAtHeight(sim, geo, y)) {
      continue;
    }

    double radiusM = std::cbrt(3.0 * static_cast<double>(droplet.parcelMl) * 1e-6 /
                               (4.0 * kPi));
    radiusM = std::min(radiusM, 0.90 * std::min(y, static_cast<double>(geo.heightMetres) - y));
    for (int pass = 0; pass < 3 && radiusM > 0.0; ++pass) {
      const double lowerWidth = halfWidthAtHeight(sim, geo, y - radiusM);
      const double upperWidth = halfWidthAtHeight(sim, geo, y + radiusM);
      const double sideRoom = std::min(lowerWidth, upperWidth) - std::fabs(x);
      radiusM = std::min(radiusM, std::max(0.0, sideRoom * 0.88));
    }
    const float radiusPx = static_cast<float>(radiusM) * tf.scale;
    if (radiusPx < 0.6f) continue;

    const ImVec2 center = toScreen(tf, x, y);
    const sol::Phase& phase = sim.phases[static_cast<size_t>(droplet.phase)];
    const auto cloudColour = [&phase](float alphaScale) {
      return ImGui::ColorConvertFloat4ToU32(
          ImVec4(phase.colour[0], phase.colour[1], phase.colour[2],
                 std::clamp(phase.colour[3] * alphaScale, 0.0f, 1.0f)));
    };
    draw->AddCircleFilled(center, radiusPx, cloudColour(0.10f), 0);
    draw->AddCircleFilled(center, radiusPx * 0.72f, cloudColour(0.18f), 0);
    draw->AddCircleFilled(center, radiusPx * 0.43f, cloudColour(0.30f), 0);

    // Physical droplet radius controls only the fine texture within a parcel
    // cloud; the parcel's bulk volume controls the visible cloud envelope.
    const float textureRadius =
        std::clamp(droplet.radius * tf.scale, 0.55f, std::max(radiusPx * 0.14f, 0.55f));
    draw->AddCircleFilled(ImVec2(center.x + radiusPx * 0.16f, center.y + radiusPx * 0.08f),
                          textureRadius, cloudColour(0.34f), 24);
    if (shaded && radiusPx >= 6.0f) {
      const ImVec2 fleck(center.x - radiusPx * 0.31f, center.y - radiusPx * 0.29f);
      draw->AddCircleFilled(fleck, std::max(textureRadius * 0.55f, 0.6f),
                            IM_COL32(255, 255, 255, 76), 24);
    }
  }
}

void drawCrossSection(const sol::Simulation& sim, ImVec2 regionMin, ImVec2 regionSize,
                      bool shaded, float dragOffsetPx = 0.0f) {
  ImDrawList* draw = ImGui::GetWindowDrawList();
  const ImVec2 regionMax(regionMin.x + regionSize.x, regionMin.y + regionSize.y);
  draw->PushClipRect(regionMin, regionMax, true);
  draw->AddRectFilled(regionMin, regionMax, style::u32(style::col::BgSurface));

  const VesselGeometry& geo = cachedGeometry(sim);
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
  // Grab-and-shake: while the user drags the vessel it follows the mouse.
  tf.origin.x += dragOffsetPx;

  if (shaded) drawGroundShadow(draw, geo, tf);
  drawVesselGlass(draw, geo, tf);
  drawGraduation(draw, geo, tf, regionMin);
  drawLayers(draw, sim, geo, tf, shaded);
  drawEmulsionHaze(draw, sim, geo, tf);
  drawParcelClouds(draw, sim, geo, tf, shaded);
  if (shaded) drawGlassHighlights(draw, sim, geo, tf);
  drawVesselWall(draw, geo, tf);
  drawFurniture(draw, sim, geo, tf);
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
  const char* runLabel = s.funnelRunning ? "Pause" : "Run";
  if (widgets::ghostButton(runLabel, textButtonSize(runLabel)))
    s.funnelRunning = !s.funnelRunning;
  ImGui::SameLine(0.0f, ImGui::GetStyle().ItemSpacing.x);
  if (widgets::ghostButton("Reset", textButtonSize("Reset"))) {
    sol::reset(sim);
    s.funnelRunning = false;
    s.statusMessage = "Vessel reset";
  }

  ImGui::TableNextColumn();
  ImGui::TextDisabled("SIMULATION SPEED");
  ImGui::SetNextItemWidth(-1.0f);
  ImGui::SliderFloat("##speed", &s.funnelSpeed, 0.1f, 10.0f, "%.1fx");
  ImGui::EndTable();
}

void drawDerivedReadout(const sol::Simulation& sim) {
  if (sim.shake.durationS <= 0.0) return;

  const double dUm = sim.shake.sauterRadiusM * 2.0e6;  // m -> um diameter
  const double weber = 1000.0 * sim.shake.peakVelocity * sim.shake.peakVelocity *
                       sim.shake.sauterRadiusM * 2.0 /
                       std::max(1e-6, sim.phases.empty() ? 0.03
                                                         : sim.phases[0].interfacialTension * 1e-3);
  char line1[112];
  char line2[96];
  std::snprintf(line1, sizeof(line1), "DERIVED   u %.2f m/s   eps %.1f W/kg",
                sim.shake.peakVelocity, sim.shake.specificPower);
  std::snprintf(line2, sizeof(line2), "d32 %.0f um   We %.1f", dUm, weber);

  const float width = ImGui::GetContentRegionAvail().x;
  if (width <= 0.0f) return;
  const style::Metrics& m = style::metrics();
  const float lineHeight = ImGui::GetTextLineHeight();
  const ImVec2 size(width, lineHeight * 2.0f + m.gap * 1.5f);
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
  const float textX = min.x + m.gap;
  const float textY = min.y + m.gap * 0.70f;
  const float textWidth = std::max(max.x - m.gap - textX, 0.0f);
  const std::string fittedLine1 = ellipsizeText(line1, textWidth);
  const std::string fittedLine2 = ellipsizeText(line2, textWidth);
  draw->PushClipRect(ImVec2(textX, min.y), ImVec2(max.x - m.gap, max.y), true);
  draw->AddText(ImVec2(textX, textY), style::u32(style::col::TextDim), fittedLine1.c_str());
  draw->AddText(ImVec2(textX, textY + lineHeight), style::u32(style::col::Teal),
                fittedLine2.c_str());
  draw->PopClipRect();
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
    if (animatedShakeButton(sim.shake.active, textButtonSize("Shake"))) {
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
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::DragFloat("##duration", &s.shakeDurationS, 0.1f, 1.0f, 30.0f, "%.0f s");

    ImGui::TableNextColumn();
    ImGui::TextDisabled("FREQUENCY");
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::DragFloat("##frequency", &s.shakeFrequencyHz, 0.05f, 0.5f, 6.0f, "%.1f Hz");

    ImGui::TableNextColumn();
    ImGui::TextDisabled("AMPLITUDE");
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::DragFloat("##amplitude", &s.shakeAmplitudeCm, 0.1f, 1.0f, 15.0f, "%.0f cm");
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
  const float tableWidth = ImGui::GetContentRegionAvail().x;
  const bool compactHeaders = tableWidth < ImGui::GetFontSize() * 34.0f;
  if (!ImGui::BeginTable("##phase_table", 8, flags)) return;
  ImGui::TableSetupColumn("Label", ImGuiTableColumnFlags_WidthStretch, 1.80f);
  ImGui::TableSetupColumn("mL", ImGuiTableColumnFlags_WidthStretch, 0.85f);
  ImGui::TableSetupColumn(compactHeaders ? "rho" : "g/mL",
                          ImGuiTableColumnFlags_WidthStretch, 0.90f);
  ImGui::TableSetupColumn(compactHeaders ? "eta" : "mPa.s",
                          ImGuiTableColumnFlags_WidthStretch, 0.90f);
  ImGui::TableSetupColumn(compactHeaders ? "IFT" : "mN/m",
                          ImGuiTableColumnFlags_WidthStretch, 0.90f);
  ImGui::TableSetupColumn(compactHeaders ? "Stab" : "Stability",
                          ImGuiTableColumnFlags_WidthStretch, 1.15f);
  ImGui::TableSetupColumn(compactHeaders ? "Col" : "Colour",
                          ImGuiTableColumnFlags_WidthStretch, 0.72f);
  ImGui::TableSetupColumn("##remove", ImGuiTableColumnFlags_WidthFixed,
                          ImGui::GetFrameHeight() + ImGui::GetStyle().CellPadding.x * 2.0f);
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
    const float glyphButtonH = ImGui::GetFrameHeight();
    if (widgets::ghostButton("x", ImVec2(glyphButtonH, glyphButtonH)))
      removeIndex = static_cast<int>(i);

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

  if (widgets::ghostButton("+ Add phase", textButtonSize("+ Add phase"))) {
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
    ImGui::SameLine(0.0f, ImGui::GetStyle().ItemSpacing.x);
  }
  ImGui::TextWrapped("Densest phase settles to the bottom");

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

// Which phase counts as aqueous for partition maths. -1 = auto: the phase
// that IS water when one is present, else the least dense phase (where water
// sits in a normal water/organic pair -- halogenated solvents sink).
int aqueousPick = -1;

struct PartitionContext {
  bool valid = false;
  int aqueousIndex = 0;
  double volAq = 0.0;
  double volOrg = 0.0;
  double K = 1.0;  // 10^logP
  std::string organicLabel;
};

PartitionContext partitionContext(const SolubilityState& s) {
  PartitionContext ctx;
  if (!s.soluteValid || s.funnel.phases.size() < 2) return ctx;
  const sol::Simulation& sim = s.funnel;
  const size_t count = sim.phases.size();

  int autoIndex = 0;
  for (size_t i = 0; i < count; ++i) {
    if (sim.phases[i].label.find("ater") != std::string::npos) {  // "Water"/"water"
      autoIndex = static_cast<int>(i);
      break;
    }
    if (sim.phases[i].density < sim.phases[static_cast<size_t>(autoIndex)].density)
      autoIndex = static_cast<int>(i);
  }
  ctx.aqueousIndex = (aqueousPick >= 0 && static_cast<size_t>(aqueousPick) < count)
                         ? aqueousPick
                         : autoIndex;
  ctx.volAq = sim.phases[static_cast<size_t>(ctx.aqueousIndex)].volumeMl;
  for (size_t i = 0; i < count; ++i) {
    if (static_cast<int>(i) == ctx.aqueousIndex) continue;
    ctx.volOrg += sim.phases[i].volumeMl;
    if (ctx.organicLabel.empty()) ctx.organicLabel = sim.phases[i].label;
  }
  ctx.K = std::pow(10.0, s.solute.logP);
  ctx.valid = ctx.volAq > 0.0 && ctx.volOrg > 0.0;
  return ctx;
}

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
  const PartitionContext ctx = partitionContext(s);
  if (!ctx.valid) {
    ImGui::TextWrapped("A valid solute and at least two charged phases are required.");
    return;
  }
  const sol::Simulation& sim = s.funnel;
  const size_t count = sim.phases.size();

  ImGui::TextWrapped("%s  ·  logP %.2f", s.solute.name.c_str(), s.solute.logP);

  const int controlColumns = ImGui::GetContentRegionAvail().x >= 520.0f ? 2 : 1;
  constexpr ImGuiTableFlags controlFlags =
      ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_NoSavedSettings;
  if (ImGui::BeginTable("##distribution_controls", controlColumns, controlFlags)) {
    ImGui::TableNextColumn();
    const float aqueousLabelWidth =
        ImGui::CalcTextSize("Aqueous phase").x + ImGui::GetStyle().ItemInnerSpacing.x;
    const float comboWidth =
        std::max(1.0f, ImGui::GetContentRegionAvail().x - aqueousLabelWidth);
    const float previewWidth =
        std::max(comboWidth - ImGui::GetFrameHeight() -
                     ImGui::GetStyle().FramePadding.x * 2.0f,
                 0.0f);
    const std::string aqueousPreview = ellipsizeText(
        sim.phases[static_cast<size_t>(ctx.aqueousIndex)].label, previewWidth);
    ImGui::SetNextItemWidth(comboWidth);
    if (ImGui::BeginCombo("Aqueous phase", aqueousPreview.c_str())) {
      for (size_t i = 0; i < count; ++i) {
        const bool selected = static_cast<int>(i) == ctx.aqueousIndex;
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
        std::max(1.0f, ImGui::GetContentRegionAvail().x - massLabelWidth));
    ImGui::SliderFloat("Solute mass", &soluteMassMgUi, 1.0f, 1000.0f, "%.0f mg");
    ImGui::EndTable();
  }

  const sol::Partition p =
      sol::partition(static_cast<double>(soluteMassMgUi), s.solute.logP, ctx.volAq, ctx.volOrg);

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

  ImGui::TextWrapped("%.1f%% extracted into %s · Neutral-species logP approximation (no pH "
                     "correction)",
                     p.fractionOrganic * 100.0, ctx.organicLabel.c_str());
}

// ------------------------------------------------ multi-stage extraction
// The classic counter-current question: how many washes to strip the solute?
// Each wash with fresh organic removes the same fraction, so the aqueous
// remainder after n washes is q^n with q = V_aq / (K V_org + V_aq).
void drawWashReadiness(const SolubilityState& s, bool fillHeight) {
  const style::Metrics& m = style::metrics();
  const float width = ImGui::GetContentRegionAvail().x;
  if (width <= 0.0f) return;

  char phaseValue[48];
  char volumeValue[48];
  std::snprintf(phaseValue, sizeof(phaseValue), "%d charged",
                static_cast<int>(s.funnel.phases.size()));
  std::snprintf(volumeValue, sizeof(volumeValue), "%.0f mL total",
                sol::totalVolumeMl(s.funnel));
  const std::string soluteValue = s.soluteValid ? s.solute.name : "Not loaded";
  const char* labels[3] = {"SOLUTE", "PHASES", "CHARGE"};
  const std::string values[3] = {soluteValue, phaseValue, volumeValue};
  const float availableHeight =
      fillHeight ? std::max(ImGui::GetContentRegionAvail().y, ImGui::GetFrameHeight() * 3.0f)
                 : ImGui::GetFrameHeight() * 3.0f;
  const float rowHeight = availableHeight / 3.0f;
  ImDrawList* draw = ImGui::GetWindowDrawList();

  for (int row = 0; row < 3; ++row) {
    const ImVec2 min = ImGui::GetCursorScreenPos();
    ImGui::Dummy(ImVec2(width, rowHeight));
    const ImVec2 max(min.x + width, min.y + rowHeight);
    if (row > 0)
      draw->AddLine(min, ImVec2(max.x, min.y), style::u32(style::col::Border, 0.75f),
                    m.hairline);
    const ImVec2 labelSize = ImGui::CalcTextSize(labels[row]);
    const float valueWidth = std::max(width - labelSize.x - m.gap * 2.0f, 0.0f);
    const std::string value = ellipsizeText(values[row], valueWidth);
    const ImVec2 valueSize = ImGui::CalcTextSize(value.c_str());
    const float textY = min.y + (rowHeight - labelSize.y) * 0.5f;
    draw->AddText(ImVec2(min.x, textY), style::u32(style::col::TextFaint), labels[row]);
    draw->AddText(ImVec2(max.x - valueSize.x, textY),
                  style::u32(row == 0 && !s.soluteValid ? style::col::Accent
                                                       : style::col::Text),
                  value.c_str());
  }
}

void drawMultiStageExtraction(const SolubilityState& s, bool fillHeight) {
  const PartitionContext ctx = partitionContext(s);
  if (!ctx.valid) {
    ImGui::TextWrapped("Wash planning needs an imported solute and two charged phases.");
    ImGui::Spacing();
    drawWashReadiness(s, fillHeight);
    return;
  }

  const double q = ctx.volAq / (ctx.K * ctx.volOrg + ctx.volAq);  // stays in aqueous
  const double perWash = 1.0 - q;

  const char* recoveryLabel = "Per-wash recovery E = K·Vorg / (K·Vorg + Vaq)";
  const float recoveryValueWidth = ImGui::CalcTextSize("100.0%").x;
  const bool recoveryOnOneLine =
      ImGui::CalcTextSize(recoveryLabel).x + recoveryValueWidth + style::metrics().gap <=
      ImGui::GetContentRegionAvail().x;
  ImGui::TextWrapped("%s", recoveryLabel);
  if (recoveryOnOneLine) ImGui::SameLine(0.0f, style::metrics().gap);
  const bool mono = style::pushFont(style::fonts::mono());
  ImGui::TextColored(style::col::Accent, "%.1f%%", perWash * 100.0);
  style::popFont(mono);

  // Cumulative recovery bars for 1..6 washes, eased toward their targets so
  // the chart animates when volumes or logP change.
  constexpr int kMaxWashes = 6;
  static double heights[kMaxWashes] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const float dt = std::clamp(ImGui::GetIO().DeltaTime, 0.0f, 0.1f);
  const float blend = 1.0f - std::exp(-dt * 9.0f);

  const float avail = ImGui::GetContentRegionAvail().x;
  float chartH = ImGui::GetFontSize() * 5.6f;
  if (fillHeight) {
    const float resultRoom = ImGui::GetTextLineHeightWithSpacing() * 2.0f + style::metrics().gap;
    chartH = std::max(chartH, ImGui::GetContentRegionAvail().y - resultRoom);
  }
  const ImVec2 chartMin = ImGui::GetCursorScreenPos();
  ImGui::Dummy(ImVec2(avail, chartH));
  const ImVec2 chartMax(chartMin.x + avail, chartMin.y + chartH);
  ImDrawList* dl = ImGui::GetWindowDrawList();
  const style::Metrics& m = style::metrics();
  dl->AddRectFilled(chartMin, chartMax, style::u32(style::col::BgDeep, 0.55f), m.radiusSm);

  const float labelH = ImGui::GetFontSize() * 1.1f;
  const float plotH = chartH - labelH - m.gap * 1.4f;
  const float barSlot = std::max((avail - m.gap * 2.0f) / kMaxWashes, 1.0f);
  const float barW = std::max(std::min(barSlot * 0.62f, 54.0f), 1.0f);
  const float plotBase = chartMin.y + m.gap * 0.7f + plotH;

  int recommended = 0;
  for (int n = 1; n <= kMaxWashes; ++n) {
    const double recovered = 1.0 - std::pow(q, n);
    if (recommended == 0 && recovered >= 0.99) recommended = n;
    heights[n - 1] += (recovered - heights[n - 1]) * blend;

    const float x0 = chartMin.x + m.gap + (n - 1) * barSlot + (barSlot - barW) * 0.5f;
    const float h = static_cast<float>(heights[n - 1]) * plotH;
    const bool isRecommended = n == recommended;
    const ImU32 fill = isRecommended
                           ? style::u32(style::col::Accent, 0.92f)
                           : style::u32(style::col::Teal, 0.55f + 0.35f * heights[n - 1]);
    dl->AddRectFilled(ImVec2(x0, plotBase - h), ImVec2(x0 + barW, plotBase), fill,
                      m.radiusSm, ImDrawFlags_RoundCornersTop);
    if (isRecommended) {
      dl->AddRect(ImVec2(x0, plotBase - h), ImVec2(x0 + barW, plotBase),
                  style::u32(style::col::Accent), m.radiusSm, ImDrawFlags_RoundCornersTop,
                  m.hairline * 1.5f);
    }

    char buf[24];
    std::snprintf(buf, sizeof(buf), "%.0f%%", recovered * 100.0);
    const ImVec2 pctSize = ImGui::CalcTextSize(buf);
    const float pctY =
        std::max(chartMin.y + 2.0f, plotBase - h - pctSize.y - 2.0f);
    if (pctSize.x + 2.0f <= barSlot) {
      dl->AddText(ImVec2(x0 + (barW - pctSize.x) * 0.5f, pctY),
                  style::u32(style::col::TextDim), buf);
    }
    std::snprintf(buf, sizeof(buf), "%dx", n);
    const ImVec2 nSize = ImGui::CalcTextSize(buf);
    if (nSize.x + 2.0f <= barSlot) {
      dl->AddText(ImVec2(x0 + (barW - nSize.x) * 0.5f, plotBase + 3.0f),
                  style::u32(style::col::TextFaint), buf);
    }
  }

  if (recommended > 0) {
    ImGui::TextWrapped("%d wash%s with fresh %s recovers >= 99%% of the solute.",
                       recommended, recommended == 1 ? "" : "es",
                       ctx.organicLabel.c_str());
  } else {
    ImGui::TextWrapped("Even 6 washes leave > 1%% behind (q = %.3f) -- raise the organic "
                       "volume or pick a better solvent.",
                       q);
  }
}

}  // namespace

void drawExtractionLab(AppState& st) {
  SolubilityState& s = st.solubility;
  consumeExtractionImport(s);
  seedDefaultPhases(s.funnel);

  // Two-zone workspace: the cross-section gets a full-height stage on the
  // right instead of whatever scroll space the control cards left over.
  // Narrow windows fall back to the classic stacked layout.
  const float totalW = ImGui::GetContentRegionAvail().x;
  const bool twoColumns = totalW >= 860.0f;

  if (twoColumns) {
    ImGui::BeginChild("##ext_controls", ImVec2(totalW * 0.58f, 0.0f), ImGuiChildFlags_None);
  }

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

  const float multiStageHeight =
      twoColumns
          ? std::max(ImGui::GetContentRegionAvail().y - ImGui::GetStyle().ItemSpacing.y, 1.0f)
          : 0.0f;
  if (widgets::beginCard("##multi_stage_card", ImVec2(0.0f, multiStageHeight),
                         style::col::BgSurface)) {
    widgets::sectionHeader("MULTI-STAGE EXTRACTION", style::col::Violet);
    drawMultiStageExtraction(s, twoColumns);
    widgets::endCard();
  }
  if (!twoColumns) ImGui::Spacing();

  stepSimulation(s);

  if (twoColumns) {
    ImGui::EndChild();
    ImGui::SameLine(0.0f, style::metrics().gap);
  }

  const ImVec2 remaining = ImGui::GetContentRegionAvail();
  if (remaining.x <= 0.0f || remaining.y <= 0.0f) return;
  if (widgets::beginCard("##cross_section_card", remaining, style::col::BgSurface,
                         ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    widgets::sectionHeader("CROSS-SECTION", style::col::Accent);
    // Render style: flat textbook cross-section by default, shaded glassware
    // on request. Small segmented control pinned to the header line.
    {
      ImGui::PushID("##render_style");
      const style::Metrics& m = style::metrics();
      ImDrawList* dl = ImGui::GetWindowDrawList();
      const char* labels[2] = {"2D", "Shaded"};
      const ImVec2 buttonSizes[2] = {textButtonSize(labels[0]), textButtonSize(labels[1])};
      const float spacing = ImGui::GetStyle().ItemSpacing.x;
      const float totalButtonWidth = buttonSizes[0].x + spacing + buttonSizes[1].x;
      const float buttonStart =
          ImGui::GetCursorPosX() + ImGui::GetContentRegionAvail().x - totalButtonWidth;
      if (buttonStart > ImGui::GetCursorPosX()) ImGui::SameLine(buttonStart);
      for (int i = 0; i < 2; ++i) {
        if (i > 0) ImGui::SameLine(0.0f, spacing);
        ImGui::PushID(i);
        ImGui::InvisibleButton("##style", buttonSizes[i]);
        const bool clicked = ImGui::IsItemClicked();
        const bool active = (i == 1) == s.funnelRender3D;
        const bool hovered = ImGui::IsItemHovered();
        const ImVec2 bMin = ImGui::GetItemRectMin();
        const ImVec2 bMax = ImGui::GetItemRectMax();
        dl->AddRectFilled(bMin, bMax,
                          active ? style::u32(style::col::Accent, 0.85f)
                                 : style::u32(style::col::BgRaised, hovered ? 1.0f : 0.6f),
                          m.radiusSm);
        const ImVec2 tSize = ImGui::CalcTextSize(labels[i]);
        dl->AddText(ImVec2(bMin.x + (buttonSizes[i].x - tSize.x) * 0.5f,
                           bMin.y + (buttonSizes[i].y - tSize.y) * 0.5f),
                    style::u32(active ? style::col::OnAccent : style::col::TextDim), labels[i]);
        if (clicked) s.funnelRender3D = (i == 1);
        ImGui::PopID();
      }
      ImGui::PopID();
    }
    const ImVec2 canvasSize = ImGui::GetContentRegionAvail();
    if (canvasSize.x > 0.0f && canvasSize.y > 0.0f) {
      ImGui::InvisibleButton("##funnel_canvas", canvasSize, ImGuiButtonFlags_MouseButtonLeft);
      const ImVec2 rectMin = ImGui::GetItemRectMin();
      const ImVec2 rectMax = ImGui::GetItemRectMax();

      // ---- grab-and-shake ----------------------------------------------
      // Dragging the vessel with the mouse IS the shake: the smoothed drag
      // velocity becomes the slosh velocity the physics disperses with, and
      // the vessel follows the pointer while held.
      sol::Simulation& sim = s.funnel;
      const VesselGeometry& geo = cachedGeometry(sim);
      const Transform tf0 = buildTransform(geo, rectMin, canvasSize);
      const float dt = std::clamp(ImGui::GetIO().DeltaTime, 0.0f, 0.1f);

      if (ImGui::IsItemActivated()) {
        s.funnelGrabbed = true;
        s.funnelGrabAnchorX = ImGui::GetMousePos().x;
        s.funnelDragOffsetPx = 0.0f;
      }
      if (s.funnelGrabbed && !ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
        s.funnelGrabbed = false;  // released: the shake state decays naturally
      }
      if (s.funnelGrabbed) {
        const float mouseX = ImGui::GetMousePos().x;
        const float rawOffset = mouseX - s.funnelGrabAnchorX;
        const float maxOffset = canvasSize.x * 0.12f;
        const float offset = std::clamp(rawOffset, -maxOffset, maxOffset);

        if (dt > 0.0f) {
          const float instantVel = std::fabs(offset - s.funnelDragOffsetPx) / dt / tf0.scale;
          const float blend = 1.0f - std::exp(-dt * 14.0f);
          s.funnelMouseVel += (instantVel - s.funnelMouseVel) * blend;
        }
        s.funnelDragOffsetPx = offset;

        // Drive the physics with the real mouse motion: peak slosh velocity
        // = smoothed drag velocity; a nominal 5 cm stroke converts it to an
        // equivalent frequency for the power readout.
        const double u = static_cast<double>(s.funnelMouseVel);
        sim.shake.active = true;
        sim.shake.remainingS = 0.12;  // refreshed per frame; lapses on release
        sim.shake.durationS = std::max(sim.shake.durationS, 0.12);
        sim.shake.peakVelocity = u;
        sim.shake.frequencyHz = u / (2.0 * kPi * 0.05);
        sim.shake.amplitudeM = 0.05;
        sim.shake.specificPower = 0.5 * u * u * sim.shake.frequencyHz;
        sol::step(sim, static_cast<double>(dt * s.funnelSpeed));
      } else {
        // Ease the vessel back to centre after release.
        s.funnelDragOffsetPx *= std::exp(-dt * 8.0f);
        if (std::fabs(s.funnelDragOffsetPx) < 0.5f) s.funnelDragOffsetPx = 0.0f;
        s.funnelMouseVel *= std::exp(-dt * 10.0f);
      }

      drawCrossSection(sim, rectMin, ImVec2(rectMax.x - rectMin.x, rectMax.y - rectMin.y),
                       s.funnelRender3D, s.funnelDragOffsetPx);

      // Hint, bottom-centre of the stage, ellipsised inside the canvas.
      const std::string hint =
          ellipsizeText("grab the funnel and shake it", std::max(canvasSize.x - 12.0f, 0.0f));
      if (!hint.empty()) {
        const ImVec2 hintSize = ImGui::CalcTextSize(hint.c_str());
        ImGui::GetWindowDrawList()->AddText(
            ImVec2(rectMin.x + (canvasSize.x - hintSize.x) * 0.5f,
                   rectMax.y - hintSize.y - 6.0f),
            style::u32(
                style::col::TextFaint,
                s.funnelGrabbed ? 0.9f
                                : 0.45f + 0.25f * std::sin(ImGui::GetTime() * 2.0f)),
            hint.c_str());
      }
    }
    widgets::endCard();
  }
}

}  // namespace chemcad::ui

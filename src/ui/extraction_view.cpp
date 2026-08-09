// Extraction Calculator controls and fluid stage. The default view composites
// the real 3D particle simulation; the schematic is an x-z cut through that
// same immutable snapshot, never a second liquid model.

#include "imgui.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "core/model.hpp"
#include "fluid/frame.hpp"
#include "fluid/simulation.hpp"
#include "gfx/fluid_stage.hpp"
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
// the particle cut and measured bulk bands reuse this profile every frame.
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


// --------------------------------------------------------------- painting
// The schematic uses the same vessel geometry as the solver. Its phase bands
// and particle cut are measurements of the immutable fluid snapshot; glass and
// graduation furniture remain explanatory drawing aids.

// Cylindrical brightness across the vessel width: brightest left of centre
// (key light), rolling off to dark rims. xFrac in [-1, 1].
float glassShade(float xFrac) {
  const float g = std::exp(-std::pow((xFrac + 0.38f) / 0.52f, 2.0f));
  return 0.62f + 0.55f * g;
}

ImU32 phaseShade(const fluid::PhaseMaterial& phase, float xFrac,
                 float alphaScale = 1.0f) {
  const float s = glassShade(xFrac);
  const ImVec4 c(std::min(phase.colour[0] * s, 1.0f),
                 std::min(phase.colour[1] * s, 1.0f),
                 std::min(phase.colour[2] * s, 1.0f),
                 std::min(phase.colour[3] * alphaScale, 1.0f));
  return ImGui::ColorConvertFloat4ToU32(c);
}
ImU32 phaseShade(const sol::Phase& phase, float xFrac,
                 float alphaScale = 1.0f) {
  const float shade = glassShade(xFrac);
  const ImVec4 colour(std::min(phase.colour[0] * shade, 1.0f),
                      std::min(phase.colour[1] * shade, 1.0f),
                      std::min(phase.colour[2] * shade, 1.0f),
                      std::min(phase.colour[3] * alphaScale, 1.0f));
  return ImGui::ColorConvertFloat4ToU32(colour);
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
                             style::u32(style::col::Text, 0.10f));
}

void drawVesselWall(ImDrawList* draw, const VesselGeometry& geo, const Transform& tf) {
  static thread_local std::vector<ImVec2> screenOutline;
  buildScreenOutline(geo, tf, screenOutline);
  if (screenOutline.size() < 3) return;
  // A soft outer stroke under a crisp closed silhouette reads as glass
  // thickness without exposing the dense analytic sampling.
  draw->AddPolyline(screenOutline.data(), static_cast<int>(screenOutline.size()),
                    style::u32(style::col::Text, 0.28f), ImDrawFlags_Closed, 4.0f);
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
                      style::u32(style::col::Text, 0.10f), ImDrawFlags_None,
                      std::max(3.0f, tf.scale * 0.010f));
    draw->AddPolyline(streak.data(), static_cast<int>(streak.size()),
                      style::u32(style::col::Text, 0.16f), ImDrawFlags_None,
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

  const ImU32 glassFill = style::u32(style::col::TextDim, 0.35f);
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

void drawBulkBands(ImDrawList* draw, const fluid::Snapshot& snapshot,
                   const sol::Simulation& charge, const VesselGeometry& geo,
                   const Transform& tf) {
  const fluid::Diagnostics& diagnostics = snapshot.diagnostics;
  const size_t count = std::min(snapshot.phases.size(), diagnostics.phases.size());
  if (!diagnostics.valid || count == 0) return;

  // Diagnostics reports each phase's measured bulk top. Sorting these measured
  // heights reconstructs the vertical stack without assuming that edit-table
  // order is density order.
  static thread_local std::vector<size_t> order;
  static thread_local std::vector<ImVec2> polygon;
  order.resize(count);
  for (size_t i = 0; i < count; ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&diagnostics](size_t a, size_t b) {
    return diagnostics.phases[a].layerTopM < diagnostics.phases[b].layerTopM;
  });
  polygon.reserve(66);

  constexpr int kWallSamples = 32;
  double lowerM = 0.0;
  for (size_t phaseIndex : order) {
    const fluid::PhaseDiagnostics& phaseDiagnostics = diagnostics.phases[phaseIndex];
    if (!phaseDiagnostics.bulkResolved) continue;
    const double upperM =
        std::clamp(phaseDiagnostics.layerTopM, lowerM, static_cast<double>(geo.heightMetres));
    if (upperM <= lowerM + 1e-9) continue;

    polygon.clear();
    for (int sample = 0; sample <= kWallSamples; ++sample) {
      const double t = static_cast<double>(sample) / kWallSamples;
      const double z = lowerM + (upperM - lowerM) * t;
      polygon.push_back(toScreen(tf, -halfWidthAtHeight(charge, geo, z), z));
    }
    for (int sample = kWallSamples; sample >= 0; --sample) {
      const double t = static_cast<double>(sample) / kWallSamples;
      const double z = lowerM + (upperM - lowerM) * t;
      polygon.push_back(toScreen(tf, halfWidthAtHeight(charge, geo, z), z));
    }
    draw->AddConcavePolyFilled(polygon.data(), static_cast<int>(polygon.size()),
                               phaseShade(snapshot.phases[phaseIndex], -0.15f, 0.30f));

    const double halfWidth = halfWidthAtHeight(charge, geo, upperM);
    draw->AddLine(toScreen(tf, -halfWidth, upperM), toScreen(tf, halfWidth, upperM),
                  phaseShade(snapshot.phases[phaseIndex], -0.15f, 0.72f), 1.2f);
    lowerM = upperM;
  }
}

void drawParticleCut(ImDrawList* draw, const fluid::Snapshot& snapshot,
                     const sol::Simulation& charge, const VesselGeometry& geo,
                     const Transform& tf) {
  const size_t count =
      std::min({snapshot.px.size(), snapshot.py.size(), snapshot.pz.size(),
                snapshot.phase.size()});
  const double slabHalfWidthM = snapshot.particleRadiusM * 1.5;
  const float radiusPx =
      std::max(static_cast<float>(snapshot.particleRadiusM) * tf.scale, 1.0f);
  const int segments = radiusPx < 3.0f ? 8 : 0;

  // This is a true x-z section, not a silhouette: only particles whose centres
  // lie in the near-axis slab are painted.
  for (size_t i = 0; i < count; ++i) {
    if (std::fabs(static_cast<double>(snapshot.py[i])) >= slabHalfWidthM) continue;
    const size_t phaseIndex = snapshot.phase[i];
    if (phaseIndex >= snapshot.phases.size()) continue;

    const double x = snapshot.px[i];
    const double z = snapshot.pz[i];
    if (z < 0.0 || z > geo.heightMetres) continue;
    const double halfWidth = halfWidthAtHeight(charge, geo, z);
    if (std::fabs(x) > halfWidth) continue;
    const float xFraction =
        static_cast<float>(x / std::max(halfWidth, static_cast<double>(1e-9)));
    const ImVec2 centre = toScreen(tf, x, z);
    const fluid::PhaseMaterial& phase = snapshot.phases[phaseIndex];
    draw->AddCircleFilled(centre, radiusPx, phaseShade(phase, xFraction, 0.96f), segments);
    draw->AddCircle(centre, radiusPx, phaseShade(phase, xFraction, 1.0f), segments,
                    std::max(style::metrics().hairline, radiusPx * 0.12f));
    if (radiusPx >= 2.4f) {
      draw->AddCircleFilled(ImVec2(centre.x - radiusPx * 0.28f,
                                   centre.y - radiusPx * 0.28f),
                            std::max(radiusPx * 0.20f, 0.6f),
                            style::u32(style::col::Text, 0.30f), 6);
    }
  }
}

void drawCrossSection(const fluid::Snapshot& snapshot, const sol::Simulation& charge,
                      ImVec2 regionMin, ImVec2 regionSize) {
  ImDrawList* draw = ImGui::GetWindowDrawList();
  const ImVec2 regionMax(regionMin.x + regionSize.x, regionMin.y + regionSize.y);
  draw->PushClipRect(regionMin, regionMax, true);
  draw->AddRectFilled(regionMin, regionMax, style::u32(style::col::BgSurface));

  const VesselGeometry& geo = cachedGeometry(charge);
  const Transform tf = buildTransform(geo, regionMin, regionSize);
  drawGroundShadow(draw, geo, tf);
  drawVesselGlass(draw, geo, tf);
  drawGraduation(draw, geo, tf, regionMin);
  drawBulkBands(draw, snapshot, charge, geo, tf);
  drawParticleCut(draw, snapshot, charge, geo, tf);
  drawGlassHighlights(draw, charge, geo, tf);
  drawVesselWall(draw, geo, tf);
  drawFurniture(draw, charge, geo, tf);
  draw->PopClipRect();
}
void drawAnalyticFallback(const sol::Simulation& charge, ImVec2 regionMin,
                          ImVec2 regionSize) {
  ImDrawList* draw = ImGui::GetWindowDrawList();
  const ImVec2 regionMax(regionMin.x + regionSize.x, regionMin.y + regionSize.y);
  draw->PushClipRect(regionMin, regionMax, true);
  draw->AddRectFilled(regionMin, regionMax, style::u32(style::col::BgSurface));

  const VesselGeometry& geometry = cachedGeometry(charge);
  const Transform transform = buildTransform(geometry, regionMin, regionSize);
  drawGroundShadow(draw, geometry, transform);
  drawVesselGlass(draw, geometry, transform);
  drawGraduation(draw, geometry, transform, regionMin);

  static thread_local std::vector<ImVec2> band;
  constexpr int kWallSamples = 32;
  const size_t phaseCount = std::min(charge.phases.size(), charge.settledMl.size());
  double lowerMl = 0.0;
  for (size_t phaseIndex = 0; phaseIndex < phaseCount; ++phaseIndex) {
    const double upperMl =
        lowerMl + std::max(static_cast<double>(charge.settledMl[phaseIndex]), 0.0);
    const double lowerM =
        heightFractionForVolume(geometry, lowerMl) * geometry.heightMetres;
    const double upperM =
        heightFractionForVolume(geometry, upperMl) * geometry.heightMetres;
    lowerMl = upperMl;
    if (upperM <= lowerM) continue;

    band.clear();
    for (int sample = 0; sample <= kWallSamples; ++sample) {
      const double t = static_cast<double>(sample) / kWallSamples;
      const double height = lowerM + (upperM - lowerM) * t;
      band.push_back(
          toScreen(transform, -halfWidthAtHeight(charge, geometry, height), height));
    }
    for (int sample = kWallSamples; sample >= 0; --sample) {
      const double t = static_cast<double>(sample) / kWallSamples;
      const double height = lowerM + (upperM - lowerM) * t;
      band.push_back(
          toScreen(transform, halfWidthAtHeight(charge, geometry, height), height));
    }
    draw->AddConcavePolyFilled(
        band.data(), static_cast<int>(band.size()),
        phaseShade(charge.phases[phaseIndex], -0.15f, 0.62f));
  }

  for (const sol::Droplet& droplet : charge.droplets) {
    if (droplet.phase < 0 ||
        static_cast<size_t>(droplet.phase) >= charge.phases.size())
      continue;
    const double x = droplet.position.x;
    const double height = droplet.position.y;
    if (height <= 0.0 || height >= geometry.heightMetres ||
        std::fabs(x) >= halfWidthAtHeight(charge, geometry, height))
      continue;
    const double parcelRadiusM =
        std::cbrt(3.0 * std::max(static_cast<double>(droplet.parcelMl), 0.0) *
                  1.0e-6 / (4.0 * kPi));
    const float radius =
        std::max(static_cast<float>(parcelRadiusM) * transform.scale * 0.42f,
                 style::metrics().hairline);
    const ImVec2 centre = toScreen(transform, x, height);
    const sol::Phase& phase = charge.phases[static_cast<size_t>(droplet.phase)];
    draw->AddCircleFilled(centre, radius, phaseShade(phase, 0.0f, 0.76f),
                          radius < 3.0f ? 8 : 0);
    draw->AddCircle(centre, radius, phaseShade(phase, 0.0f), radius < 3.0f ? 8 : 0,
                    style::metrics().hairline);
  }

  drawGlassHighlights(draw, charge, geometry, transform);
  drawVesselWall(draw, geometry, transform);
  drawFurniture(draw, charge, geometry, transform);
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

constexpr std::array<double, 3> kFluidSpacingM{6.0e-3, 4.0e-3, 3.0e-3};
constexpr std::array<const char*, 3> kFluidResolutionNames{"Coarse", "Normal", "Fine"};

double selectedFluidSpacing(const SolubilityState& s) {
  const size_t index = static_cast<size_t>(s.fluidResolution);
  return kFluidSpacingM[std::min(index, kFluidSpacingM.size() - 1)];
}
template <typename T>
void hashFluidValue(size_t& seed, const T& value) {
  seed ^= std::hash<T>{}(value) + 0x9e3779b9u + (seed << 6u) + (seed >> 2u);
}

size_t fluidConfigurationSignature(const SolubilityState& s) {
  size_t signature = 0;
  hashFluidValue(signature, static_cast<int>(s.funnel.vessel));
  hashFluidValue(signature, s.funnel.vesselVolumeMl);
  hashFluidValue(signature, static_cast<int>(s.fluidResolution));
  hashFluidValue(signature, s.funnel.phases.size());
  for (const sol::Phase& phase : s.funnel.phases) {
    hashFluidValue(signature, phase.label);
    hashFluidValue(signature, phase.volumeMl);
    hashFluidValue(signature, phase.density);
    hashFluidValue(signature, phase.viscosity);
    hashFluidValue(signature, phase.interfacialTension);
    for (float channel : phase.colour) hashFluidValue(signature, channel);
  }
  return signature;
}


std::vector<fluid::PhaseMaterial> fluidMaterials(const sol::Simulation& charge) {
  std::vector<fluid::PhaseMaterial> materials;
  materials.reserve(charge.phases.size());
  for (const sol::Phase& phase : charge.phases) {
    fluid::PhaseMaterial material;
    material.label = phase.label;
    material.restDensity = phase.density * 1000.0;          // g/mL -> kg/m^3
    material.dynamicViscosity = phase.viscosity * 1.0e-3;  // mPa.s -> Pa.s
    material.volumeMl = phase.volumeMl;
    for (size_t channel = 0; channel < 4; ++channel)
      material.colour[channel] = phase.colour[channel];
    materials.push_back(std::move(material));
  }
  return materials;
}

std::vector<double> interfaceTensions(const sol::Simulation& charge) {
  const size_t count = charge.phases.size();
  std::vector<double> sigma(count * count, 0.0);
  for (size_t i = 0; i < count; ++i) {
    for (size_t j = i + 1; j < count; ++j) {
      // The legacy editor stores one representative IFT per phase. Their
      // symmetric mean is the only order-independent pair value available.
      const double pair =
          0.5 * (charge.phases[i].interfacialTension +
                 charge.phases[j].interfacialTension) *
          1.0e-3;  // mN/m -> N/m
      sigma[i * count + j] = pair;
      sigma[j * count + i] = pair;
    }
  }
  return sigma;
}

struct FluidBoundaryState {
  const SolubilityState* owner = nullptr;
  fluid::Simulation* observedSimulation = nullptr;
  size_t configuration = 0;
  bool unavailable = false;
  std::string reason;
};

FluidBoundaryState& fluidBoundaryState(SolubilityState& s) {
  static FluidBoundaryState state;
  const size_t configuration = fluidConfigurationSignature(s);
  const bool replacedState =
      state.owner != &s ||
      (state.observedSimulation != nullptr && s.fluid.get() != state.observedSimulation);
  if (replacedState) {
    state = {};
    state.owner = &s;
    state.configuration = configuration;
  } else if (state.configuration != configuration) {
    // A changed charge is a new attempt. The failed instance stays quarantined
    // until this exact configuration boundary changes or Retry is pressed.
    state.configuration = configuration;
    state.unavailable = false;
    state.reason.clear();
  }
  state.observedSimulation = s.fluid.get();
  return state;
}

void recordFluidFailure(SolubilityState& s, const std::exception& error) {
  FluidBoundaryState& state = fluidBoundaryState(s);
  state.unavailable = true;
  state.reason = error.what();
  if (state.reason.empty()) state.reason = "The fluid solver reported an unknown error.";
  state.observedSimulation = s.fluid.get();
  s.funnelRunning = false;
  s.fluidGrabActive = false;
  s.fluidManualAcceleration = {0.0, 0.0, 0.0};
  s.extractionRenderMode = ExtractionRenderMode::Schematic2D;
  s.statusMessage = "Physics unavailable: " + state.reason;
}

template <typename Interaction>
bool runFluidInteraction(SolubilityState& s, Interaction&& interaction,
                         bool forceAttempt = false) {
  FluidBoundaryState& state = fluidBoundaryState(s);
  if (state.unavailable && !forceAttempt) return false;
  try {
    interaction();
    state.observedSimulation = s.fluid.get();
    return true;
  } catch (const std::exception& error) {
    recordFluidFailure(s, error);
    return false;
  }
}

bool rechargeFluid(SolubilityState& s, bool forceAttempt = false) {
  const bool charged = runFluidInteraction(
      s,
      [&s] {
        // Configure a fresh candidate so a failed calibration cannot leave the
        // panel holding a partially configured Simulation.
        auto candidate = std::make_unique<fluid::Simulation>();
        candidate->setVessel(s.funnel.vessel, s.funnel.vesselVolumeMl);
        candidate->setResolution(selectedFluidSpacing(s));
        candidate->setPhases(fluidMaterials(s.funnel), interfaceTensions(s.funnel));
        candidate->charge();
        s.fluidManualAcceleration = {0.0, 0.0, 0.0};
        candidate->setManualAcceleration(s.fluidManualAcceleration);
        s.fluid = std::move(candidate);
      },
      forceAttempt);
  if (charged) {
    FluidBoundaryState& state = fluidBoundaryState(s);
    state.unavailable = false;
    state.reason.clear();
  }
  return charged;
}

fluid::Simulation* availableFluid(SolubilityState& s) {
  FluidBoundaryState& state = fluidBoundaryState(s);
  if (state.unavailable) return nullptr;
  if (!s.fluid && !rechargeFluid(s)) return nullptr;
  return s.fluid.get();
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


constexpr std::array<const char*, 3> kShakeAxisNames{"Vertical", "Horizontal",
                                                     "Diagonal"};

std::array<double, 3> selectedShakeAxis(FluidShakeAxis axis) {
  constexpr double kInvSqrtTwo = 0.70710678118654752440;
  switch (axis) {
    case FluidShakeAxis::Vertical:
      return {0.0, 0.0, 1.0};
    case FluidShakeAxis::Horizontal:
      return {1.0, 0.0, 0.0};
    case FluidShakeAxis::Diagonal:
      return {kInvSqrtTwo, 0.0, kInvSqrtTwo};
  }
  return {0.0, 0.0, 1.0};
}

fluid::VesselMotion shakeReportMotion(const SolubilityState& s) {
  fluid::VesselMotion motion;
  motion.shaking = true;
  motion.shakeAxis = selectedShakeAxis(s.shakeAxis);
  motion.shakeFrequencyHz = s.shakeFrequencyHz;
  motion.shakeAmplitudeM = static_cast<double>(s.shakeAmplitudeCm) * 0.01;
  return motion;
}

void drawFluidDiagnostics(SolubilityState& s) {
  const size_t resolutionIndex =
      std::min(static_cast<size_t>(s.fluidResolution),
               kFluidResolutionNames.size() - 1);

  ImGui::Spacing();
  ImGui::Separator();
  ImGui::TextDisabled("REAL FLUID READOUT");

  fluid::Simulation* simulation = availableFluid(s);
  std::shared_ptr<const fluid::Snapshot> snapshot;
  fluid::Solver::Stats stats;
  double realTimeFactor = std::numeric_limits<double>::quiet_NaN();
  bool stepping = false;
  std::string status;
  const bool readable =
      simulation &&
      runFluidInteraction(s, [&] {
        snapshot = simulation->snapshot();
        stats = simulation->solverStats();
        realTimeFactor = simulation->realTimeFactor();
        stepping = simulation->stepping();
        status = simulation->statusLine();
      });
  if (!readable || !snapshot) {
    const FluidBoundaryState& state = fluidBoundaryState(s);
    ImGui::TextColored(style::col::Danger, "Physics unavailable");
    ImGui::TextWrapped("%s", state.reason.empty() ? "Fluid setup did not complete."
                                                  : state.reason.c_str());
    return;
  }

  const char* runState = s.funnelRunning ? "RUNNING" : "PAUSED";
  ImGui::TextColored(s.funnelRunning ? style::col::Success : style::col::TextDim,
                     "%s%s", runState, stepping ? " - physics step active" : "");
  ImGui::TextWrapped("%s", status.c_str());
  ImGui::TextWrapped("Particle-resolved %s-resolution estimates (dx %.0f mm)",
                     kFluidResolutionNames[resolutionIndex],
                     selectedFluidSpacing(s) * 1000.0);

  const bool completedStep = snapshot->elapsedS > 0.0;
  const fluid::Diagnostics& diagnostics = snapshot->diagnostics;
  const int columns = ImGui::GetContentRegionAvail().x >= 500.0f ? 3 : 1;
  constexpr ImGuiTableFlags flags =
      ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_NoSavedSettings;
  if (ImGui::BeginTable("##fluid_diagnostics", columns, flags)) {
    ImGui::TableNextColumn();
    ImGui::TextDisabled("SUBSTEPS");
    if (completedStep && stats.substeps > 0)
      ImGui::Text("%d", stats.substeps);
    else
      ImGui::TextUnformatted("Waiting for completed step");

    ImGui::TableNextColumn();
    ImGui::TextDisabled("WORST DENSITY ERROR");
    if (completedStep && stats.substeps > 0)
      ImGui::Text("%.2f%%", stats.maxDensityError * 100.0);
    else
      ImGui::TextUnformatted("Unavailable");

    ImGui::TableNextColumn();
    ImGui::TextDisabled("REAL-TIME FACTOR");
    if (completedStep && std::isfinite(realTimeFactor))
      ImGui::Text("%.2fx", realTimeFactor);
    else
      ImGui::TextUnformatted("Waiting for completed step");

    ImGui::TableNextColumn();
    ImGui::TextDisabled("DISPERSED");
    if (completedStep && diagnostics.valid)
      ImGui::Text("%.1f%%", diagnostics.dispersedFraction * 100.0);
    else
      ImGui::TextUnformatted("Unavailable");

    ImGui::TableNextColumn();
    ImGui::TextDisabled("SAUTER d32");
    if (!completedStep || !diagnostics.valid)
      ImGui::TextUnformatted("Unavailable");
    else if (diagnostics.sauterDiameterM > 0.0)
      ImGui::Text("%.0f um", diagnostics.sauterDiameterM * 1.0e6);
    else
      ImGui::TextUnformatted("No resolved drops");

    ImGui::TableNextColumn();
    ImGui::TextDisabled("INTERFACIAL AREA");
    if (completedStep && diagnostics.valid)
      ImGui::Text("%.2f cm^2", diagnostics.interfacialAreaM2 * 1.0e4);
    else
      ImGui::TextUnformatted("Unavailable");

    ImGui::TableNextColumn();
    ImGui::TextDisabled("FREE SURFACE");
    if (completedStep && diagnostics.valid)
      ImGui::Text("%.1f mm", diagnostics.freeSurfaceM * 1000.0);
    else
      ImGui::TextUnformatted("Unavailable");

    if (completedStep && diagnostics.valid) {
      const size_t phaseCount =
          std::min(diagnostics.phases.size(), s.funnel.phases.size());
      for (size_t i = 0; i < phaseCount; ++i) {
        ImGui::TableNextColumn();
        ImGui::TextDisabled("%s BULK", s.funnel.phases[i].label.c_str());
        if (diagnostics.phases[i].bulkResolved)
          ImGui::Text("%.1f mL", diagnostics.phases[i].bulkMl);
        else
          ImGui::TextUnformatted("Bulk unresolved");
      }
    }
    ImGui::EndTable();
  }
}

void drawTransportControls(SolubilityState& s) {
  sol::Simulation& charge = s.funnel;
  static const char* kVesselNames[] = {"Separatory funnel", "Decanting flask",
                                       "Graduated cylinder"};
  const float width = ImGui::GetContentRegionAvail().x;
  const int columns = width >= 720.0f ? 4 : (width >= 420.0f ? 2 : 1);
  constexpr ImGuiTableFlags flags =
      ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_NoSavedSettings;
  if (!ImGui::BeginTable("##transport_grid", columns, flags)) return;

  ImGui::TableNextColumn();
  ImGui::TextDisabled("VESSEL");
  ImGui::SetNextItemWidth(-1.0f);
  if (ImGui::Combo("##vessel", &s.funnelVessel, kVesselNames, IM_ARRAYSIZE(kVesselNames))) {
    s.funnelVessel = std::clamp(s.funnelVessel, 0, 2);
    charge.vessel = static_cast<sol::Vessel>(s.funnelVessel);
    sol::reset(charge);
    if (rechargeFluid(s))
      s.statusMessage = std::string("Recharged into ") + kVesselNames[s.funnelVessel];
  }

  ImGui::TableNextColumn();
  ImGui::TextDisabled("TRANSPORT");
  const char* runLabel = s.funnelRunning ? "Pause" : "Run";
  const ImVec2 runSize = textButtonSize(runLabel);
  const ImVec2 resetSize = textButtonSize("Reset");
  const float actionSpacing = ImGui::GetStyle().ItemSpacing.x;
  const bool actionsFit =
      runSize.x + actionSpacing + resetSize.x <= ImGui::GetContentRegionAvail().x;
  if (widgets::ghostButton(runLabel, runSize))
    s.funnelRunning = !s.funnelRunning;
  if (actionsFit) ImGui::SameLine(0.0f, actionSpacing);
  if (widgets::ghostButton("Reset", resetSize)) {
    sol::reset(charge);
    const bool recharged = rechargeFluid(s);
    s.funnelRunning = false;
    s.fluidTiltTargetDeg = 0.0f;
    s.fluidTiltCurrentDeg = 0.0f;
    s.fluidTiltAngularVelocityRadS = 0.0f;
    if (recharged) s.statusMessage = "Particle vessel recharged";
  }
  ImGui::TextColored(s.funnelRunning ? style::col::Success : style::col::TextDim,
                     "%s", s.funnelRunning ? "Running" : "Paused");

  ImGui::TableNextColumn();
  ImGui::TextDisabled("SIMULATION SPEED");
  ImGui::SetNextItemWidth(-1.0f);
  ImGui::SliderFloat("##speed", &s.funnelSpeed, 0.1f, 10.0f, "%.1fx");

  ImGui::TableNextColumn();
  ImGui::TextDisabled("RESOLUTION");
  int resolution = static_cast<int>(s.fluidResolution);
  ImGui::SetNextItemWidth(-1.0f);
  if (ImGui::Combo("##resolution", &resolution, kFluidResolutionNames.data(),
                   static_cast<int>(kFluidResolutionNames.size()))) {
    resolution = std::clamp(resolution, 0,
                            static_cast<int>(kFluidResolutionNames.size()) - 1);
    s.fluidResolution = static_cast<FluidResolution>(resolution);
    if (rechargeFluid(s)) {
      s.statusMessage = std::string("Recharged at ") + kFluidResolutionNames[resolution] +
                        " particle resolution";
    }
  }
  ImGui::EndTable();
}

void drawDerivedReadout(const SolubilityState& s) {
  const fluid::VesselMotion motion = shakeReportMotion(s);
  const double peakVelocity = fluid::shakePeakVelocity(motion);
  const double specificPower = fluid::shakeSpecificPower(motion);
  char line[144];
  std::snprintf(line, sizeof(line), "PHYSICS USED   peak velocity %.2f m/s   specific power %.1f W/kg",
                peakVelocity, specificPower);

  const float width = ImGui::GetContentRegionAvail().x;
  if (width <= 0.0f) return;
  const style::Metrics& metrics = style::metrics();
  const ImVec2 size(width, ImGui::GetTextLineHeight() + metrics.gap * 1.5f);
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::InvisibleButton("##derived_readout", size);
  const float hover = widgets::hoverT(ImGui::GetItemID(), ImGui::IsItemHovered());
  const ImVec2 max(min.x + size.x, min.y + size.y);
  ImDrawList* draw = ImGui::GetWindowDrawList();
  draw->AddRectFilled(min, max,
                      style::mix(style::col::BgDeep, style::col::BgRaised, hover, 0.72f),
                      metrics.radiusSm);
  draw->AddRect(min, max, style::mix(style::col::Border, style::col::Teal, hover),
                metrics.radiusSm, 0, metrics.hairline);

  const bool mono = style::pushFont(style::fonts::mono());
  const std::string fitted =
      ellipsizeText(line, std::max(width - metrics.gap * 2.0f, 0.0f));
  draw->AddText(ImVec2(min.x + metrics.gap, min.y + metrics.gap * 0.70f),
                style::u32(style::col::Teal), fitted.c_str());
  style::popFont(mono);
  if (ImGui::IsItemHovered()) {
    ImGui::BeginTooltip();
    ImGui::TextUnformatted("u_peak = 2*pi*f*A");
    ImGui::TextUnformatted("specific power = u_peak^2*f/2");
    ImGui::TextUnformatted("Computed by fluid::shakePeakVelocity and fluid::shakeSpecificPower.");
    ImGui::EndTooltip();
  }
}

void drawShakeControls(SolubilityState& s) {
  fluid::Simulation* simulation = availableFluid(s);
  bool shaking = false;
  if (simulation)
    runFluidInteraction(s, [&] { shaking = simulation->shaking(); });
  const float width = ImGui::GetContentRegionAvail().x;
  const int columns = width >= 800.0f ? 5 : (width >= 420.0f ? 2 : 1);
  constexpr ImGuiTableFlags flags =
      ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_NoSavedSettings;
  if (ImGui::BeginTable("##shake_grid", columns, flags)) {
    ImGui::TableNextColumn();
    ImGui::TextDisabled("ACTION");
    ImGui::BeginDisabled(simulation == nullptr);
    const bool shakeRequested =
        animatedShakeButton(shaking, textButtonSize("Shake"));
    ImGui::EndDisabled();
    if (shakeRequested && simulation) {
      const std::array<double, 3> axis = selectedShakeAxis(s.shakeAxis);
      const double amplitudeM = static_cast<double>(s.shakeAmplitudeCm) * 0.01;
      const bool started = runFluidInteraction(s, [&] {
        simulation->shake(axis, s.shakeDurationS, s.shakeFrequencyHz, amplitudeM);
      });
      if (started) {
        s.funnelRunning = true;
        const fluid::VesselMotion report = shakeReportMotion(s);
        char message[176];
        std::snprintf(
            message, sizeof(message),
            "%s shake: %.0f s at %.1f Hz, %.0f cm; u %.2f m/s, power %.1f W/kg",
            kShakeAxisNames[static_cast<size_t>(s.shakeAxis)], s.shakeDurationS,
            s.shakeFrequencyHz, s.shakeAmplitudeCm, fluid::shakePeakVelocity(report),
            fluid::shakeSpecificPower(report));
        s.statusMessage = message;
      }
    }

    ImGui::TableNextColumn();
    ImGui::TextDisabled("AXIS");
    int axis = static_cast<int>(s.shakeAxis);
    ImGui::SetNextItemWidth(-1.0f);
    if (ImGui::Combo("##shake_axis", &axis, kShakeAxisNames.data(),
                     static_cast<int>(kShakeAxisNames.size()))) {
      s.shakeAxis = static_cast<FluidShakeAxis>(
          std::clamp(axis, 0, static_cast<int>(kShakeAxisNames.size()) - 1));
    }

    ImGui::TableNextColumn();
    ImGui::TextDisabled("DURATION");
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::DragFloat("##duration", &s.shakeDurationS, 0.1f, 1.0f, 30.0f, "%.0f s");

    ImGui::TableNextColumn();
    ImGui::TextDisabled("FREQUENCY");
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::DragFloat("##frequency", &s.shakeFrequencyHz, 0.05f, 0.5f, 6.0f,
                     "%.1f Hz");

    ImGui::TableNextColumn();
    ImGui::TextDisabled("AMPLITUDE");
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::DragFloat("##amplitude", &s.shakeAmplitudeCm, 0.1f, 1.0f, 15.0f,
                     "%.0f cm");
    ImGui::EndTable();
  }

  ImGui::Spacing();
  drawDerivedReadout(s);
}

void setTiltTarget(SolubilityState& s, float degrees) {
  s.fluidTiltTargetDeg = std::clamp(degrees, 0.0f, 180.0f);
  s.funnelRunning = true;
}

void drawTiltControls(SolubilityState& s) {
  ImGui::Separator();
  ImGui::TextDisabled("TILT / INVERT");
  const int columns = ImGui::GetContentRegionAvail().x >= 460.0f ? 2 : 1;
  constexpr ImGuiTableFlags flags =
      ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_NoSavedSettings;
  if (!ImGui::BeginTable("##tilt_grid", columns, flags)) return;
  ImGui::TableNextColumn();
  ImGui::SetNextItemWidth(-1.0f);
  if (ImGui::SliderFloat("##tilt", &s.fluidTiltTargetDeg, 0.0f, 180.0f, "%.0f deg"))
    s.funnelRunning = true;

  ImGui::TableNextColumn();
  const float actionWidth = ImGui::GetContentRegionAvail().x;
  const ImVec2 uprightSize = textButtonSize("Upright");
  const ImVec2 ventSize = textButtonSize("Vent");
  const ImVec2 invertSize = textButtonSize("Invert");
  const float spacing = ImGui::GetStyle().ItemSpacing.x;
  const bool oneLine =
      uprightSize.x + ventSize.x + invertSize.x + spacing * 2.0f <= actionWidth;
  if (widgets::ghostButton("Upright", uprightSize)) setTiltTarget(s, 0.0f);
  if (oneLine) ImGui::SameLine(0.0f, spacing);
  if (widgets::ghostButton("Vent", ventSize)) setTiltTarget(s, 135.0f);
  if (oneLine) ImGui::SameLine(0.0f, spacing);
  if (widgets::ghostButton("Invert", invertSize)) setTiltTarget(s, 180.0f);
  ImGui::EndTable();
}

void drawControls(SolubilityState& s) {
  if (widgets::beginCard("##transport_card", ImVec2(0.0f, 0.0f),
                         style::col::BgSurface)) {
    widgets::sectionHeader("VESSEL & TRANSPORT", style::col::Teal);
    drawTransportControls(s);
    drawFluidDiagnostics(s);
    widgets::endCard();
  }

  ImGui::Spacing();
  if (widgets::beginCard("##shake_card", ImVec2(0.0f, 0.0f),
                         style::col::BgSurface)) {
    widgets::sectionHeader("SHAKE", style::col::Accent);
    drawShakeControls(s);
    ImGui::Spacing();
    drawTiltControls(s);
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
    const bool recharged = rechargeFluid(s);
    if (recharged) {
      const double total = sol::totalVolumeMl(s.funnel);
      char buf[128];
      std::snprintf(buf, sizeof(buf), "Recharged vessel: %d phase(s), %.0f mL total",
                    static_cast<int>(s.funnel.phases.size()), total);
      s.statusMessage = buf;
    }
  }
}

void updateFluidPose(SolubilityState& s, fluid::Simulation& simulation, double dt) {
  if (dt <= 0.0) return;
  double angleRad = s.fluidTiltCurrentDeg * kPi / 180.0;
  const double targetRad = s.fluidTiltTargetDeg * kPi / 180.0;
  const double previousRate = s.fluidTiltAngularVelocityRadS;

  // A rate-limited, exponentially driven hand motion avoids an impulsive Euler
  // term when the target changes. The rate and its derivative come from the
  // exact trajectory drawn, so frameAcceleration receives consistent
  // Coriolis/Euler terms (frame.hpp).
  constexpr double kMaximumRate = 2.0 * kPi / 3.0;  // 120 degrees/s
  const double desiredRate =
      std::clamp((targetRad - angleRad) * 4.5, -kMaximumRate, kMaximumRate);
  const double rateBlend = 1.0 - std::exp(-dt * 10.0);
  double angularRate = previousRate + (desiredRate - previousRate) * rateBlend;
  const double previousError = targetRad - angleRad;
  angleRad += angularRate * dt;
  const bool crossedTarget = previousError != 0.0 &&
                             previousError * (targetRad - angleRad) <= 0.0;
  if (crossedTarget ||
      (std::fabs(targetRad - angleRad) < 1.0e-4 &&
       std::fabs(angularRate) < 1.0e-3)) {
    angleRad = targetRad;
    angularRate = 0.0;
  }
  const double angularAcceleration = (angularRate - previousRate) / dt;
  s.fluidTiltCurrentDeg = static_cast<float>(angleRad * 180.0 / kPi);
  s.fluidTiltAngularVelocityRadS = static_cast<float>(angularRate);

  fluid::Pose pose;
  pose.orientation = {std::cos(angleRad * 0.5), 0.0, std::sin(angleRad * 0.5), 0.0};
  simulation.setPose(pose, {0.0, angularRate, 0.0},
                     {0.0, angularAcceleration, 0.0});
}

void stepSimulation(SolubilityState& s) {
  fluid::Simulation* simulation = availableFluid(s);
  if (!simulation) return;
  const double frameDelta = static_cast<double>(ImGui::GetIO().DeltaTime);
  if (!(frameDelta > 0.0) || !std::isfinite(frameDelta)) return;
  const double requestedSeconds = frameDelta * s.funnelSpeed;
  const double poseSeconds = std::min(requestedSeconds, 0.1);
  runFluidInteraction(s, [&] {
    updateFluidPose(s, *simulation, poseSeconds);
    if (s.funnelRunning)
      simulation->requestAdvance(frameDelta * s.funnelSpeed);
  });
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
  const bool recharged = rechargeFluid(s);
  s.funnelRunning = false;
  soluteMassMgUi = static_cast<float>(imp.soluteMassMg);
  if (recharged)
    s.statusMessage =
        "Imported " + a->name + " + " + b->name + " from the Solubility Suite";
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

  ImGui::TextWrapped("%s  |  logP %.2f", s.solute.name.c_str(), s.solute.logP);

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

  ImGui::TextWrapped("%.1f%% extracted into %s | Neutral-species logP approximation (no pH "
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

  const char* recoveryLabel = "Per-wash recovery E = K*Vorg / (K*Vorg + Vaq)";
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


bool stageChoiceButton(const char* id, const char* label, bool selected,
                       ImVec2 requestedSize = ImVec2(0.0f, 0.0f)) {
  const ImVec2 size =
      requestedSize.x > 0.0f && requestedSize.y > 0.0f
          ? requestedSize
          : ImVec2(std::max(ImGui::GetContentRegionAvail().x, 1.0f),
                   ImGui::GetFrameHeight());
  ImGui::InvisibleButton(id, size);
  const bool clicked = ImGui::IsItemClicked();
  const bool hovered = ImGui::IsItemHovered();
  const ImVec2 min = ImGui::GetItemRectMin();
  const ImVec2 max = ImGui::GetItemRectMax();
  const style::Metrics& metrics = style::metrics();
  ImDrawList* draw = ImGui::GetWindowDrawList();
  draw->AddRectFilled(min, max,
                      selected ? style::u32(style::col::Accent, 0.88f)
                               : style::u32(style::col::BgRaised,
                                            hovered ? 1.0f : 0.62f),
                      metrics.radiusSm);
  draw->AddRect(min, max,
                selected ? style::u32(style::col::AccentHover)
                         : style::u32(style::col::Border),
                metrics.radiusSm, 0, metrics.hairline);
  const ImVec2 text = ImGui::CalcTextSize(label);
  draw->AddText(ImVec2(min.x + (size.x - text.x) * 0.5f,
                       min.y + (size.y - text.y) * 0.5f),
                style::u32(selected ? style::col::OnAccent : style::col::TextDim),
                label);
  return clicked;
}

void drawStageToolbar(SolubilityState& s) {
  const int columns = ImGui::GetContentRegionAvail().x >= 450.0f ? 3 : 1;
  constexpr ImGuiTableFlags flags =
      ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_NoSavedSettings;
  if (ImGui::BeginTable("##stage_modes", columns, flags)) {
    ImGui::TableNextColumn();
    if (stageChoiceButton("##fluid_3d", "3D fluid",
                          s.extractionRenderMode == ExtractionRenderMode::Fluid3D)) {
      s.extractionRenderMode = ExtractionRenderMode::Fluid3D;
    }
    ImGui::TableNextColumn();
    if (stageChoiceButton(
            "##schematic_2d", "2D schematic",
            s.extractionRenderMode == ExtractionRenderMode::Schematic2D)) {
      s.extractionRenderMode = ExtractionRenderMode::Schematic2D;
      s.fluidGrabMode = false;
      s.fluidGrabActive = false;
    }
    ImGui::TableNextColumn();
    ImGui::BeginDisabled(s.extractionRenderMode != ExtractionRenderMode::Fluid3D);
    const float grabButtonHeight = ImGui::GetFrameHeight();
    if (stageChoiceButton("##grab_shake", "G", s.fluidGrabMode,
                          ImVec2(grabButtonHeight, grabButtonHeight))) {
      s.fluidGrabMode = !s.fluidGrabMode;
      if (!s.fluidGrabMode) s.fluidGrabActive = false;
    }
    ImGui::SameLine(0.0f, ImGui::GetStyle().ItemSpacing.x);
    ImGui::TextUnformatted(s.fluidGrabMode ? "GRAB active" : "GRAB off");
    ImGui::EndDisabled();
    ImGui::EndTable();
  }

  // A dedicated square toggle makes the gesture discoverable and prevents
  // left-drag orbit from competing with hand forcing.
  if (s.extractionRenderMode == ExtractionRenderMode::Schematic2D) {
    ImGui::TextWrapped("Particle cut: |y| < 1.5 particle radii.");
  } else if (s.fluidGrabMode) {
    ImGui::TextWrapped(
        "Drag in any direction; vertical drag drives world +z / -z.");
  } else {
    ImGui::TextWrapped("Left-drag orbit | wheel zoom | double-click re-frame.");
  }
}

void updateStageInput(SolubilityState& s, fluid::Simulation& simulation,
                      gfx::FluidStage& stage, const fluid::Snapshot& snapshot,
                      ImVec2 canvasSize) {
  const ImGuiIO& io = ImGui::GetIO();
  const double dt = std::clamp(static_cast<double>(io.DeltaTime), 0.0, 0.1);
  const bool hovered = ImGui::IsItemHovered();
  const bool active = ImGui::IsItemActive();
  const bool reframe =
      hovered && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left);

  if (reframe) {
    stage.camera.frame(snapshot.vesselHeightM, snapshot.maxRadiusM,
                       canvasSize.x / std::max(canvasSize.y, 1.0f));
    s.fluidGrabActive = false;
  } else if (hovered && io.MouseWheel != 0.0f) {
    stage.camera.zoom(io.MouseWheel);
  }

  if (s.fluidGrabMode && !reframe && ImGui::IsItemActivated()) {
    s.fluidGrabActive = true;
    s.fluidGrabAnchorPx = {io.MousePos.x, io.MousePos.y};
    s.funnelRunning = true;
  }
  if (!s.fluidGrabMode || !ImGui::IsMouseDown(ImGuiMouseButton_Left))
    s.fluidGrabActive = false;

  std::array<double, 3> targetAcceleration{0.0, 0.0, 0.0};
  if (s.fluidGrabMode && s.fluidGrabActive && active && dt > 0.0) {
    const std::array<float, 3> eye = stage.camera.eye();
    const double horizontalLength =
        std::max(std::hypot(static_cast<double>(eye[0]), static_cast<double>(eye[1])),
                 1.0e-9);
    const std::array<double, 3> cameraRight{
        -static_cast<double>(eye[1]) / horizontalLength,
        static_cast<double>(eye[0]) / horizontalLength, 0.0};

    // A full-height one-frame gesture maps to 4 g and is magnitude-clamped
    // there. This makes screen scale irrelevant while bounding a noisy input
    // device before it reaches the pressure solve.
    constexpr double kGravity = 9.80665;
    const double accelerationPerPixel =
        4.0 * kGravity / std::max(static_cast<double>(canvasSize.y), 1.0);
    const double pointerDeltaX =
        static_cast<double>(io.MousePos.x - s.fluidGrabAnchorPx[0]);
    const double pointerDeltaY =
        static_cast<double>(io.MousePos.y - s.fluidGrabAnchorPx[1]);
    s.fluidGrabAnchorPx = {io.MousePos.x, io.MousePos.y};
    const double horizontal = pointerDeltaX * accelerationPerPixel;
    const double vertical = -pointerDeltaY * accelerationPerPixel;
    targetAcceleration = {cameraRight[0] * horizontal,
                          cameraRight[1] * horizontal, vertical};
    const double magnitude =
        std::sqrt(targetAcceleration[0] * targetAcceleration[0] +
                  targetAcceleration[1] * targetAcceleration[1] +
                  targetAcceleration[2] * targetAcceleration[2]);
    if (magnitude > 4.0 * kGravity) {
      const double scale = 4.0 * kGravity / magnitude;
      for (double& component : targetAcceleration) component *= scale;
    }
  } else if (!s.fluidGrabMode && active && !reframe) {
    stage.camera.orbit(io.MouseDelta.x, io.MouseDelta.y);
  }

  if (dt > 0.0) {
    // Exponential filtering is independent of frame rate. Release targets zero
    // with a slightly slower decay so the fictitious acceleration does not
    // acquire a discontinuity at the end of a hand stroke.
    const double timeConstant = s.fluidGrabActive ? 0.055 : 0.090;
    const double blend = 1.0 - std::exp(-dt / timeConstant);
    for (size_t component = 0; component < 3; ++component) {
      s.fluidManualAcceleration[component] +=
          (targetAcceleration[component] - s.fluidManualAcceleration[component]) * blend;
      if (!s.fluidGrabActive &&
          std::fabs(s.fluidManualAcceleration[component]) < 1.0e-4) {
        s.fluidManualAcceleration[component] = 0.0;
      }
    }
  }
  runFluidInteraction(
      s, [&] { simulation.setManualAcceleration(s.fluidManualAcceleration); });
}

bool drawPhysicsUnavailableState(SolubilityState& s, gfx::FluidStage& stage) {
  FluidBoundaryState& boundary = fluidBoundaryState(s);
  stage.requested = false;
  s.extractionRenderMode = ExtractionRenderMode::Schematic2D;
  s.fluidGrabMode = false;
  s.fluidGrabActive = false;
  stage.snapshot.reset();

  ImGui::TextColored(style::col::Danger, "Physics unavailable");
  ImGui::TextWrapped("The 3D stage is missing because: %s",
                     boundary.reason.empty() ? "fluid setup did not complete."
                                             : boundary.reason.c_str());
  ImGui::TextDisabled(
      "The calculator remains available; showing the 2D charge schematic instead.");

  if (widgets::ghostButton("Retry physics", textButtonSize("Retry physics"))) {
    if (rechargeFluid(s, true)) {
      s.statusMessage = "Fluid physics restored; showing the 2D schematic";
      return false;
    }
  }

  const ImVec2 canvasSize = ImGui::GetContentRegionAvail();
  if (canvasSize.x <= 0.0f || canvasSize.y <= 0.0f) return true;
  ImGui::InvisibleButton("##physics_unavailable_stage", canvasSize);
  drawAnalyticFallback(s.funnel, ImGui::GetItemRectMin(), canvasSize);
  return true;
}

void drawFluidStage(AppState& st, SolubilityState& s) {
  gfx::FluidStage& stage = st.fluidStage;
  drawStageToolbar(s);

  fluid::Simulation* simulation = availableFluid(s);
  if (!simulation) {
    drawPhysicsUnavailableState(s, stage);
    return;
  }

  const bool wants3D =
      s.extractionRenderMode == ExtractionRenderMode::Fluid3D;
  const bool rendererReady = stage.available && stage.texture != 0;
  if (wants3D && !rendererReady) {
    if (!stage.status.empty())
      ImGui::TextWrapped("3D unavailable (%s); showing the particle schematic.",
                         stage.status.c_str());
    else if (!stage.available)
      ImGui::TextWrapped(
          "3D unavailable (headless session or no OpenGL); showing the particle schematic.");
    else
      ImGui::TextWrapped(
          "3D renderer has no completed frame yet; showing the particle schematic.");
  }

  const ImVec2 canvasSize = ImGui::GetContentRegionAvail();
  if (canvasSize.x <= 0.0f || canvasSize.y <= 0.0f) {
    stage.requested = false;
    return;
  }

  std::shared_ptr<const fluid::Snapshot> snapshot;
  if (!runFluidInteraction(s, [&] { snapshot = simulation->snapshot(); })) {
    drawPhysicsUnavailableState(s, stage);
    return;
  }
  if (!snapshot) {
    stage.requested = false;
    stage.snapshot.reset();
    ImGui::TextDisabled("Fluid snapshot unavailable; recharge the vessel.");
    return;
  }

  if (wants3D) {
    const ImVec2 framebufferScale = ImGui::GetIO().DisplayFramebufferScale;
    stage.requested = true;
    stage.width = std::max(
        1, static_cast<int>(std::lround(canvasSize.x * framebufferScale.x)));
    stage.height = std::max(
        1, static_cast<int>(std::lround(canvasSize.y * framebufferScale.y)));
    stage.snapshot = snapshot;
  } else {
    stage.requested = false;
    stage.snapshot.reset();
  }

  const ImVec2 cursor = ImGui::GetCursorPos();
  ImGui::InvisibleButton("##fluid_stage_input", canvasSize,
                         ImGuiButtonFlags_MouseButtonLeft);
  const ImVec2 rectMin = ImGui::GetItemRectMin();
  if (wants3D) {
    updateStageInput(s, *simulation, stage, *snapshot, canvasSize);
  } else {
    s.fluidGrabActive = false;
    const double dt =
        std::clamp(static_cast<double>(ImGui::GetIO().DeltaTime), 0.0, 0.1);
    const double decay = dt > 0.0 ? std::exp(-dt / 0.090) : 1.0;
    for (double& component : s.fluidManualAcceleration) component *= decay;
    runFluidInteraction(
        s, [&] { simulation->setManualAcceleration(s.fluidManualAcceleration); });
  }

  if (fluidBoundaryState(s).unavailable) {
    stage.requested = false;
    stage.snapshot.reset();
    drawAnalyticFallback(s.funnel, rectMin, canvasSize);
  } else if (wants3D && rendererReady) {
    // The behaviour item is submitted first and owns all input; Image is then
    // placed over exactly the same rectangle using the previous FBO frame.
    ImGui::SetCursorPos(cursor);
    ImGui::Image(ImTextureRef(static_cast<ImTextureID>(stage.texture)), canvasSize,
                 ImVec2(0.0f, 1.0f), ImVec2(1.0f, 0.0f));
  } else {
    drawCrossSection(*snapshot, s.funnel, rectMin, canvasSize);
  }
}
}  // namespace

void drawExtractionLab(AppState& st) {
  SolubilityState& s = st.solubility;
  consumeExtractionImport(s);
  seedDefaultPhases(s.funnel);
  availableFluid(s);
  st.fluidStage.requested = false;

  if (!s.statusMessage.empty()) {
    const bool physicsUnavailable = fluidBoundaryState(s).unavailable;
    ImGui::PushStyleColor(
        ImGuiCol_Text,
        physicsUnavailable ? style::col::Danger : style::col::TextDim);
    ImGui::TextWrapped("%s", s.statusMessage.c_str());
    ImGui::PopStyleColor();
    ImGui::Spacing();
  }

  // Wide panels reserve a persistent stage column; narrow panels stack a
  // minimum-height stage into the parent scroll region rather than clipping it
  // below the calculator cards.
  const float totalWidth = ImGui::GetContentRegionAvail().x;
  const bool twoColumns = totalWidth >= 860.0f;
  if (twoColumns) {
    ImGui::BeginChild("##ext_controls", ImVec2(totalWidth * 0.58f, 0.0f),
                      ImGuiChildFlags_None);
  }

  drawControls(s);
  ImGui::Spacing();

  if (widgets::beginCard("##phases_card", ImVec2(0.0f, 0.0f),
                         style::col::BgSurface)) {
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
          ? std::max(ImGui::GetContentRegionAvail().y -
                         ImGui::GetStyle().ItemSpacing.y,
                     1.0f)
          : 0.0f;
  if (widgets::beginCard("##multi_stage_card", ImVec2(0.0f, multiStageHeight),
                         style::col::BgSurface)) {
    widgets::sectionHeader("MULTI-STAGE EXTRACTION", style::col::Violet);
    drawMultiStageExtraction(s, twoColumns);
    widgets::endCard();
  }
  if (!twoColumns) ImGui::Spacing();

  // Exactly one advance is submitted per visible panel frame. Stage input from
  // the prior frame is already filtered into Simulation::setManualAcceleration.
  stepSimulation(s);

  if (twoColumns) {
    ImGui::EndChild();
    ImGui::SameLine(0.0f, style::metrics().gap);
  }

  const ImVec2 remaining = ImGui::GetContentRegionAvail();
  if (remaining.x <= 0.0f) return;
  const float minimumStageHeight = ImGui::GetFrameHeight() * 18.0f;
  const float stageHeight =
      twoColumns ? std::max(remaining.y, ImGui::GetFrameHeight() * 6.0f)
                 : std::max(remaining.y, minimumStageHeight);
  if (widgets::beginCard(
          "##fluid_stage_card", ImVec2(remaining.x, stageHeight),
          style::col::BgSurface,
          ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    widgets::sectionHeader("VESSEL STAGE", style::col::Accent);
    drawFluidStage(st, s);
    widgets::endCard();
  }
}

}  // namespace chemcad::ui

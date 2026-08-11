// Extraction Calculator controls and fluid stage. The default view composites
// the real 3D particle simulation; the schematic is an x-z cut through that
// same immutable snapshot, never a second liquid model.

#include "imgui.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "core/model.hpp"
#include "fluid/frame.hpp"
#include "fluid/simulation.hpp"
#include "gfx/fluid_stage.hpp"
#include "sol/funnel.hpp"
#include "sol/solubility.hpp"
#include "ui/app_state.hpp"
#include "ui/charts.hpp"
#include "ui/icons.hpp"
#include "ui/layout.hpp"
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

double halfWidthAtHeight(const sol::Simulation& sim, const VesselGeometry& geo,
                         double heightM);

// The vessel body is a surface of revolution, so its section is the region
// |x| <= w(z): a stack of trapezoids. Handing the concave outline to
// AddConcavePolyFilled instead used to fill something close to its convex hull,
// which is a kite far larger than the funnel and read as a glow around it.
void drawVesselGlass(ImDrawList* draw, const sol::Simulation& sim,
                     const VesselGeometry& geo, const Transform& tf) {
  constexpr int kStrips = 96;
  const ImU32 body = style::u32(style::col::Text, 0.10f);
  const ImDrawListFlags savedFlags = draw->Flags;
  // Shared strip edges would otherwise show as seams through their AA fringes.
  draw->Flags &= ~ImDrawListFlags_AntiAliasedFill;
  double lowerZ = 0.0;
  double lowerHalf = halfWidthAtHeight(sim, geo, 0.0);
  for (int strip = 1; strip <= kStrips; ++strip) {
    const double upperZ = geo.heightMetres * static_cast<double>(strip) / kStrips;
    const double upperHalf = halfWidthAtHeight(sim, geo, upperZ);
    if (lowerHalf > 1.0e-9 || upperHalf > 1.0e-9) {
      const ImVec2 quad[4] = {
          toScreen(tf, -lowerHalf, lowerZ), toScreen(tf, lowerHalf, lowerZ),
          toScreen(tf, upperHalf, upperZ), toScreen(tf, -upperHalf, upperZ)};
      draw->AddConvexPolyFilled(quad, 4, body);
    }
    lowerZ = upperZ;
    lowerHalf = upperHalf;
  }
  draw->Flags = savedFlags;
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

// The particles are point samples of a continuum, not beads of liquid. Drawing
// one translucent circle each is what made the section read as a bag of orbs:
// at dx = 8 mm there are only about ten particles across the funnel, so every
// one of them was individually legible. This instead reconstructs the field the
// particles sample and fills its iso-contour, which is the 2D counterpart of
// what the screen-space renderer does in 3D, and the two views therefore agree.
//
// Marching TRIANGLES rather than squares: splitting each cell on a diagonal
// removes the saddle ambiguity entirely and every emitted polygon is convex,
// which is exactly what ImDrawList::AddConvexPolyFilled requires.
void drawLiquidSection(ImDrawList* draw, const fluid::Snapshot& snapshot,
                       const sol::Simulation& charge, const VesselGeometry& geo,
                       const Transform& tf) {
  const size_t count =
      std::min({snapshot.px.size(), snapshot.py.size(), snapshot.pz.size(),
                snapshot.phase.size()});
  if (count == 0 || snapshot.phases.empty()) return;

  const double slabHalfWidthM = snapshot.particleSpacingM * 0.75;
  const float radiusPx =
      std::max(static_cast<float>(snapshot.particleRadiusM) * tf.scale, 1.0f);
  // A lone particle must resolve to a disc of exactly radiusPx, so the influence
  // radius and the iso level are chosen together: with w(d) = (1 - (d/R)^2)^2
  // and R = 1.6 r, w(r) is the level below.
  const float influencePx = radiusPx * 1.6f;
  constexpr float kIsoLevel = 0.372f;  // (1 - (1/1.6)^2)^2

  const ImVec2 topLeft = toScreen(tf, -geo.halfWidthMetres, geo.heightMetres);
  const ImVec2 bottomRight = toScreen(tf, geo.halfWidthMetres, 0.0);
  const float padding = influencePx + 2.0f;
  const float minX = topLeft.x - padding, maxX = bottomRight.x + padding;
  const float minY = topLeft.y - padding, maxY = bottomRight.y + padding;
  if (!(maxX > minX) || !(maxY > minY)) return;

  // Cell size follows the particle radius so the contour stays smooth at any
  // zoom, and the node budget bounds the cost when the stage is large.
  constexpr int kMaxNodes = 24000;
  float cell = std::clamp(radiusPx * 0.5f, 3.0f, 16.0f);
  int nx = static_cast<int>((maxX - minX) / cell) + 2;
  int ny = static_cast<int>((maxY - minY) / cell) + 2;
  if (nx * ny > kMaxNodes) {
    cell *= std::sqrt(static_cast<float>(nx * ny) / static_cast<float>(kMaxNodes));
    nx = static_cast<int>((maxX - minX) / cell) + 2;
    ny = static_cast<int>((maxY - minY) / cell) + 2;
  }
  if (nx < 2 || ny < 2) return;

  // Reused across frames: the section is redrawn every frame and this must not
  // allocate in the draw path.
  static std::vector<float> field;
  static std::vector<float> phaseWeight;
  field.assign(static_cast<size_t>(nx) * ny, 0.0f);
  phaseWeight.assign(static_cast<size_t>(nx) * ny, 0.0f);

  const float inverseCell = 1.0f / cell;
  const float influence2 = influencePx * influencePx;
  for (size_t i = 0; i < count; ++i) {
    // A true x-z section: only particles near the cutting plane contribute.
    if (std::fabs(static_cast<double>(snapshot.py[i])) >= slabHalfWidthM) continue;
    const size_t phaseIndex = snapshot.phase[i];
    if (phaseIndex >= snapshot.phases.size()) continue;
    const ImVec2 centre = toScreen(tf, snapshot.px[i], snapshot.pz[i]);
    const float phaseB = phaseIndex == 0 ? 0.0f : 1.0f;

    const int i0 = std::max(0, static_cast<int>((centre.x - influencePx - minX) * inverseCell));
    const int i1 = std::min(nx - 1, static_cast<int>((centre.x + influencePx - minX) * inverseCell) + 1);
    const int j0 = std::max(0, static_cast<int>((centre.y - influencePx - minY) * inverseCell));
    const int j1 = std::min(ny - 1, static_cast<int>((centre.y + influencePx - minY) * inverseCell) + 1);
    for (int j = j0; j <= j1; ++j) {
      const float ny_ = minY + static_cast<float>(j) * cell - centre.y;
      for (int k = i0; k <= i1; ++k) {
        const float nx_ = minX + static_cast<float>(k) * cell - centre.x;
        const float d2 = nx_ * nx_ + ny_ * ny_;
        if (d2 >= influence2) continue;
        const float t = 1.0f - d2 / influence2;
        const float w = t * t;
        const size_t node = static_cast<size_t>(j) * nx + k;
        field[node] += w;
        phaseWeight[node] += w * phaseB;
      }
    }
  }

  // Confining the liquid to the glass by zeroing whole nodes quantises the wall
  // to the cell size and leaves a visible staircase. Each emitted vertex is
  // projected onto the wall instead, so the liquid meets the glass exactly.
  const auto clampToVessel = [&](ImVec2 point) {
    const double z = static_cast<double>(tf.origin.y - point.y) / tf.scale;
    const double halfWidth =
        (z < 0.0 || z > geo.heightMetres) ? 0.0 : halfWidthAtHeight(charge, geo, z);
    const float limit = static_cast<float>(halfWidth) * tf.scale;
    point.x = std::clamp(point.x, tf.origin.x - limit, tf.origin.x + limit);
    return point;
  };

  const auto nodePosition = [&](int k, int j) {
    return ImVec2(minX + static_cast<float>(k) * cell, minY + static_cast<float>(j) * cell);
  };
  const auto crossing = [&](int ka, int ja, int kb, int jb) {
    const float fa = field[static_cast<size_t>(ja) * nx + ka];
    const float fb = field[static_cast<size_t>(jb) * nx + kb];
    const float denominator = fb - fa;
    const float t = std::fabs(denominator) < 1.0e-6f
                        ? 0.5f
                        : std::clamp((kIsoLevel - fa) / denominator, 0.0f, 1.0f);
    const ImVec2 a = nodePosition(ka, ja);
    const ImVec2 b = nodePosition(kb, jb);
    return ImVec2(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t);
  };

  // Anti-aliased fills put a blended fringe around every polygon. Adjacent
  // triangles share their edges exactly, so those fringes would print the
  // triangulation onto the liquid as a visible mesh. The outer silhouette is
  // covered by the drawn glass wall, so plain fills are the right trade.
  const ImDrawListFlags savedFlags = draw->Flags;
  draw->Flags &= ~ImDrawListFlags_AntiAliasedFill;

  const fluid::PhaseMaterial& phaseA = snapshot.phases[0];
  const fluid::PhaseMaterial& phaseB =
      snapshot.phases[snapshot.phases.size() > 1 ? 1 : 0];
  ImVec2 polygon[4];
  const int corners[2][3][2] = {{{0, 0}, {1, 0}, {1, 1}}, {{0, 0}, {1, 1}, {0, 1}}};
  for (int j = 0; j + 1 < ny; ++j) {
    for (int k = 0; k + 1 < nx; ++k) {
      for (const auto& triangle : corners) {
        int used = 0;
        float phaseSum = 0.0f;
        float weightSum = 0.0f;
        for (int v = 0; v < 3; ++v) {
          const int ka = k + triangle[v][0], ja = j + triangle[v][1];
          const int kb = k + triangle[(v + 1) % 3][0], jb = j + triangle[(v + 1) % 3][1];
          const size_t nodeA = static_cast<size_t>(ja) * nx + ka;
          const size_t nodeB = static_cast<size_t>(jb) * nx + kb;
          const bool insideA = field[nodeA] >= kIsoLevel;
          const bool insideB = field[nodeB] >= kIsoLevel;
          if (insideA && used < 4) {
            polygon[used++] = clampToVessel(nodePosition(ka, ja));
            phaseSum += phaseWeight[nodeA];
            weightSum += field[nodeA];
          }
          if (insideA != insideB && used < 4)
            polygon[used++] = clampToVessel(crossing(ka, ja, kb, jb));
        }
        if (used < 3) continue;

        // Interpolate the two liquids by their local mixture. Picking a single
        // dominant phase instead would flip between colours wherever the mix
        // sits near a half, which prints the triangulation as a checkerboard.
        const float mix =
            weightSum > 0.0f ? std::clamp(phaseSum / weightSum, 0.0f, 1.0f) : 0.0f;
        fluid::PhaseMaterial shaded = phaseA;
        for (int channel = 0; channel < 4; ++channel) {
          shaded.colour[channel] = phaseA.colour[channel] +
                                   (phaseB.colour[channel] - phaseA.colour[channel]) * mix;
        }
        const float xFraction = std::clamp(
            (polygon[0].x - tf.origin.x) /
                std::max(1.0f, static_cast<float>(geo.halfWidthMetres) * tf.scale),
            -1.0f, 1.0f);
        draw->AddConvexPolyFilled(polygon, used, phaseShade(shaded, xFraction, 0.96f));
      }
    }
  }
  draw->Flags = savedFlags;
}

void drawCrossSection(const fluid::Snapshot& snapshot, const sol::Simulation& charge,
                      const std::array<double, 3>& vesselPositionM, ImVec2 regionMin,
                      ImVec2 regionSize) {
  ImDrawList* draw = ImGui::GetWindowDrawList();
  const ImVec2 regionMax(regionMin.x + regionSize.x, regionMin.y + regionSize.y);
  draw->PushClipRect(regionMin, regionMax, true);
  draw->AddRectFilled(regionMin, regionMax, style::u32(style::col::BgSurface));

  const VesselGeometry& geo = cachedGeometry(charge);
  const Transform stationary = buildTransform(geo, regionMin, regionSize);
  // The solver works in vessel coordinates, so a shaken vessel is stationary in
  // its own frame and only the contents appear to move. Offsetting the
  // transform by the vessel's world displacement puts the glassware back in the
  // hand: the vessel travels and the liquid lags inside it, which is what a
  // shaken separatory funnel actually looks like. The caller supplies the
  // position rather than the snapshot doing so, because the hand-driven part is
  // live every frame while the snapshot is only as fresh as the last step.
  Transform tf = stationary;
  tf.origin.x += static_cast<float>(vesselPositionM[0]) * tf.scale;
  tf.origin.y -= static_cast<float>(vesselPositionM[2]) * tf.scale;

  drawGroundShadow(draw, geo, stationary);
  drawVesselGlass(draw, charge, geo, tf);
  // The graduation is a fixed scale read against the bench, not paint on the
  // glass, so it stays put while the vessel moves past it.
  drawGraduation(draw, geo, stationary, regionMin);
  drawBulkBands(draw, snapshot, charge, geo, tf);
  drawLiquidSection(draw, snapshot, charge, geo, tf);
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
  drawVesselGlass(draw, charge, geometry, transform);
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

struct FluidResolutionPreset {
  const char* name;
  const char* choiceLabel;
  fluid::QualityProfile quality;
};

constexpr std::array<FluidResolutionPreset, 3> kFluidPresets{{
    {"Interactive", "Interactive preview", fluid::QualityProfile::interactive()},
    {"Balanced", "Balanced", fluid::QualityProfile::balanced()},
    {"Quality", "Quality - slower than real time", fluid::QualityProfile::quality()},
}};

size_t selectedFluidPreset(const SolubilityState& s) {
  return std::min(static_cast<size_t>(s.fluidResolution),
                  kFluidPresets.size() - 1);
}

double selectedFluidSpacing(const SolubilityState& s) {
  return kFluidPresets[selectedFluidPreset(s)].quality.spacing;
}

uint64_t estimatedParticleCount(const SolubilityState& s, size_t presetIndex) {
  const double volumeM3 = std::max(sol::totalVolumeMl(s.funnel), 0.0) * 1.0e-6;
  const double spacing =
      kFluidPresets[std::min(presetIndex, kFluidPresets.size() - 1)].quality.spacing;
  return static_cast<uint64_t>(
      std::llround(volumeM3 / (spacing * spacing * spacing)));
}

void rememberFluidRate(SolubilityState& s, double realTimeFactor) {
  if (!std::isfinite(realTimeFactor) || realTimeFactor < 0.0) return;
  const size_t preset = selectedFluidPreset(s);
  s.fluidPresetRealTimeFactor[preset] = realTimeFactor;
  s.fluidPresetRealTimeFactorValid[preset] = true;
  s.fluidPresetMeasuredParticles[preset] =
      estimatedParticleCount(s, preset);
}

std::string fluidPresetTradeText(const SolubilityState& s,
                                 size_t presetIndex) {
  const size_t preset = std::min(presetIndex, kFluidPresets.size() - 1);
  const uint64_t estimatedParticles = estimatedParticleCount(s, preset);
  char text[256];
  const int written = std::snprintf(
      text, sizeof(text),
      "%s | dx %.0f mm | compression limit %.1f%% | ~%llu particles",
      kFluidPresets[preset].choiceLabel,
      kFluidPresets[preset].quality.spacing * 1000.0,
      kFluidPresets[preset].quality.densityTolerance * 100.0,
      static_cast<unsigned long long>(estimatedParticles));
  // rememberFluidRate has been recording this per preset all along so the
  // picker could show an observed cost next to the extrapolated particle count
  // (solubility_state.hpp). An observed rate beats an estimate; the count it
  // was measured at travels with it, because the charge may have changed since.
  if (written > 0 && written < static_cast<int>(sizeof(text)) &&
      s.fluidPresetRealTimeFactorValid[preset]) {
    std::snprintf(text + written, sizeof(text) - static_cast<size_t>(written),
                  " | measured %.2fx at ~%llu",
                  s.fluidPresetRealTimeFactor[preset],
                  static_cast<unsigned long long>(
                      s.fluidPresetMeasuredParticles[preset]));
  }
  return text;
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
  s.fluidHand.reset();
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

// Everything the worker needs, copied off the UI thread so the build touches no
// panel state while it runs.
struct FluidBuildRequest {
  sol::Vessel vessel = sol::Vessel::SeparatoryFunnel;
  double vesselVolumeMl = 250.0;
  fluid::QualityProfile quality;
  std::vector<fluid::PhaseMaterial> materials;
  std::vector<double> tensions;
};

struct FluidBuildResult {
  std::shared_ptr<fluid::Simulation> simulation;
  std::string error;
};

// Starts a background build. The panel keeps drawing the analytic schematic
// until it lands, which is what stops the vessel SDF build and the first charge
// from stalling the workspace the first time it is opened.
void startFluidBuild(SolubilityState& s) {
  if (s.fluidBuildPending || s.fluidTasks == nullptr) return;

  FluidBuildRequest request;
  request.vessel = s.funnel.vessel;
  request.vesselVolumeMl = s.funnel.vesselVolumeMl;
  request.quality = kFluidPresets[selectedFluidPreset(s)].quality;
  request.materials = fluidMaterials(s.funnel);
  request.tensions = interfaceTensions(s.funnel);

  const std::size_t signature = fluidConfigurationSignature(s);
  s.fluidBuildPending = true;
  s.fluidBuildSignature = signature;
  s.fluidHand.reset();
  s.fluidShakeProgressValid = false;
  s.fluidShakeStartElapsedS = 0.0;
  s.fluidShakeEndElapsedS = 0.0;
  s.fluidTasks->run<FluidBuildResult>(
      [request] {
        FluidBuildResult result;
        try {
          // Configured fully before it is published, so a failed setup can
          // never hand the panel a half-built Simulation.
          auto candidate = std::make_shared<fluid::Simulation>();
          candidate->setVessel(request.vessel, request.vesselVolumeMl);
          candidate->setQuality(request.quality);
          candidate->setPhases(request.materials, request.tensions);
          candidate->setManualMotion({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
          result.simulation = std::move(candidate);
        } catch (const std::exception& error) {
          result.error = error.what();
        } catch (...) {
          result.error = "fluid setup failed";
        }
        return result;
      },
      [&s, signature](FluidBuildResult result) {
        s.fluidBuildPending = false;
        // The user may have changed vessel, volumes, solvents or resolution
        // while this was in flight; that build is for a configuration that no
        // longer exists, so it is dropped and the next frame asks again.
        if (signature != fluidConfigurationSignature(s)) return;
        if (!result.simulation) {
          recordFluidFailure(s, std::runtime_error(result.error.empty()
                                                       ? "fluid setup failed"
                                                       : result.error));
          return;
        }
        s.fluid = std::move(result.simulation);
        // Every rebuild gets it: a Simulation is replaced whenever the vessel,
        // the solvents or the resolution change, and a device left behind on
        // the old one would quietly stop being used.
        s.fluid->setAccelerator(s.fluidAccelerator);
        FluidBoundaryState& state = fluidBoundaryState(s);
        state.unavailable = false;
        state.reason.clear();
        state.observedSimulation = s.fluid.get();
      });
}

// Drops the current simulation and requests a replacement. Returns whether a
// build was started; the caller's own charge edit has already taken effect on
// the analytic funnel, so the panel stays usable while the particles catch up.
bool rechargeFluid(SolubilityState& s, bool forceAttempt = false) {
  FluidBoundaryState& state = fluidBoundaryState(s);
  if (state.unavailable && !forceAttempt) return false;
  if (forceAttempt) {
    state.unavailable = false;
    state.reason.clear();
  }
  s.fluid.reset();
  state.observedSimulation = nullptr;
  startFluidBuild(s);
  return s.fluidBuildPending;
}

fluid::Simulation* availableFluid(SolubilityState& s) {
  FluidBoundaryState& state = fluidBoundaryState(s);
  if (state.unavailable) return nullptr;
  // Deferred until the workspace is active; see fluidConstructionAllowed.
  if (!s.fluid && !s.fluidBuildPending && s.fluidConstructionAllowed) startFluidBuild(s);
  return s.fluid.get();
}

// Why the particle solver cannot be driven right now, in the user's terms.
// Deferred construction is NOT a failure: the panel only builds the simulation
// once the workspace is demonstrably on top, and calling that "unavailable"
// reports a fault where there is none.
const char* fluidAbsenceReason(SolubilityState& s) {
  const FluidBoundaryState& state = fluidBoundaryState(s);
  if (state.unavailable)
    return state.reason.empty() ? "The fluid solver reported an unknown error."
                                : state.reason.c_str();
  if (s.fluid) return "The particle solver has not published a frame yet.";
  if (s.fluidBuildPending) return "The particle solver is still being built.";
  return "The particle solver starts once the Extraction workspace is on top.";
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

template <typename Stats>
double maxDensityCompression(const Stats& stats) {
  if constexpr (requires(const Stats& value) { value.maxDensityCompression; }) {
    return stats.maxDensityCompression;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

template <typename Stats>
double maxDensityDeficit(const Stats& stats) {
  if constexpr (requires(const Stats& value) { value.maxDensityDeficit; }) {
    return stats.maxDensityDeficit;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

struct FluidDiagnosticTraces {
  charts::Trace compression;
  charts::Trace deficit;
  charts::Trace realTimeFactor;
  charts::Trace dispersedFraction;
  charts::Trace interfacialArea;
  charts::Trace sauterDiameter;
  fluid::Simulation* source = nullptr;
  uint64_t revision = 0;
  bool hasRevision = false;
  std::vector<charts::StackSegment> phaseSegments;

  void clear() {
    compression.clear();
    deficit.clear();
    realTimeFactor.clear();
    dispersedFraction.clear();
    interfacialArea.clear();
    sauterDiameter.clear();
  }
};

FluidDiagnosticTraces& fluidDiagnosticTraces() {
  // This panel has one active fluid simulation, so retaining the ring buffers
  // here avoids coupling frame-rate-only presentation state to SolubilityState.
  static FluidDiagnosticTraces traces;
  return traces;
}

void sampleFluidDiagnostics(FluidDiagnosticTraces& traces,
                            fluid::Simulation* simulation,
                            const fluid::Snapshot& snapshot,
                            const fluid::Solver::Stats& stats,
                            double compression, double deficit,
                            double realTimeFactor) {
  if (traces.source != simulation) {
    traces.clear();
    traces.source = simulation;
    traces.hasRevision = false;
  }

  const bool completedStep = snapshot.elapsedS > 0.0;
  if (!completedStep) {
    if (!traces.hasRevision || traces.revision != snapshot.revision) {
      traces.clear();
      traces.revision = snapshot.revision;
      traces.hasRevision = true;
    }
    return;
  }
  // UI frames commonly repeat one immutable publication; only revisions are
  // solver-time samples, so frame rate cannot distort the trace history.
  if (traces.hasRevision && traces.revision == snapshot.revision) return;
  if (traces.hasRevision && snapshot.revision < traces.revision) traces.clear();

  traces.revision = snapshot.revision;
  traces.hasRevision = true;
  if (stats.substeps > 0) {
    traces.compression.push(compression);
    traces.deficit.push(deficit);
  }
  traces.realTimeFactor.push(realTimeFactor);
  if (snapshot.diagnostics.valid) {
    traces.dispersedFraction.push(snapshot.diagnostics.dispersedFraction);
    traces.interfacialArea.push(snapshot.diagnostics.interfacialAreaM2 * 1.0e4);
    traces.sauterDiameter.push(snapshot.diagnostics.sauterDiameterM * 1.0e6);
  }
}


void drawFluidDiagnostics(SolubilityState& s, float height) {
  const size_t resolutionIndex = selectedFluidPreset(s);
  fluid::Simulation* simulation = availableFluid(s);
  std::shared_ptr<const fluid::Snapshot> snapshot;
  fluid::Solver::Stats stats;
  double realTimeFactor = std::numeric_limits<double>::quiet_NaN();
  const bool readable =
      simulation &&
      runFluidInteraction(s, [&] {
        snapshot = simulation->snapshot();
        stats = simulation->solverStats();
        realTimeFactor = simulation->realTimeFactor();
      });
  if (!readable || !snapshot) {
    // A solver that has not been built yet is not a broken one. Claiming a
    // fault here sent users looking for a problem that did not exist.
    const bool failed = fluidBoundaryState(s).unavailable;
    widgets::emptyState(failed ? icons::Icon::Warning : icons::Icon::Timer,
                        failed ? "Physics unavailable" : "Physics starting",
                        fluidAbsenceReason(s));
    return;
  }

  const bool completedStep = snapshot->elapsedS > 0.0;
  if (completedStep) rememberFluidRate(s, realTimeFactor);
  const double compression = maxDensityCompression(stats);
  const double deficit = maxDensityDeficit(stats);
  const fluid::Diagnostics& diagnostics = snapshot->diagnostics;
  FluidDiagnosticTraces& traces = fluidDiagnosticTraces();
  sampleFluidDiagnostics(traces, simulation, *snapshot, stats, compression,
                         deficit, realTimeFactor);

  const bool compressionAvailable =
      completedStep && stats.substeps > 0 && std::isfinite(compression);
  const bool deficitAvailable =
      completedStep && stats.substeps > 0 && std::isfinite(deficit);
  const bool rateAvailable = completedStep && std::isfinite(realTimeFactor);
  const bool dispersedAvailable =
      completedStep && diagnostics.valid &&
      std::isfinite(diagnostics.dispersedFraction);
  const bool sauterAvailable =
      completedStep && diagnostics.valid &&
      std::isfinite(diagnostics.sauterDiameterM);
  const bool areaAvailable =
      completedStep && diagnostics.valid &&
      std::isfinite(diagnostics.interfacialAreaM2);
  const bool surfaceAvailable =
      completedStep && diagnostics.valid &&
      std::isfinite(diagnostics.freeSurfaceM);

  char compressionValue[32] = "--";
  char deficitValue[32] = "--";
  char rateValue[32] = "--";
  char dispersedValue[32] = "--";
  char sauterValue[32] = "--";
  char areaValue[32] = "--";
  if (compressionAvailable)
    std::snprintf(compressionValue, sizeof(compressionValue), "%.2f",
                  compression * 100.0);
  if (deficitAvailable)
    std::snprintf(deficitValue, sizeof(deficitValue), "%.2f", deficit * 100.0);
  if (rateAvailable)
    std::snprintf(rateValue, sizeof(rateValue), "%.2f", realTimeFactor);
  if (dispersedAvailable)
    std::snprintf(dispersedValue, sizeof(dispersedValue), "%.1f",
                  diagnostics.dispersedFraction * 100.0);
  if (sauterAvailable && diagnostics.sauterDiameterM > 0.0)
    std::snprintf(sauterValue, sizeof(sauterValue), "%.0f",
                  diagnostics.sauterDiameterM * 1.0e6);
  else if (sauterAvailable)
    std::snprintf(sauterValue, sizeof(sauterValue), "No drops");
  if (areaAvailable)
    std::snprintf(areaValue, sizeof(areaValue), "%.2f",
                  diagnostics.interfacialAreaM2 * 1.0e4);

  charts::SparklineStyle compressionStyle;
  compressionStyle.accent = style::col::Data;
  compressionStyle.ceilingValue =
      kFluidPresets[resolutionIndex].quality.densityTolerance;
  charts::SparklineStyle deficitStyle;
  deficitStyle.accent = style::col::DataDim;
  deficitStyle.ceilingValue = 1.0;
  charts::SparklineStyle rateStyle;
  rateStyle.accent = style::col::DataBright;
  charts::SparklineStyle dispersedStyle;
  dispersedStyle.accent = style::col::DataBright;
  dispersedStyle.ceilingValue = 1.0;
  charts::SparklineStyle sauterStyle;
  sauterStyle.accent = style::col::DataDim;
  charts::SparklineStyle areaStyle;
  areaStyle.accent = style::col::Data;

  const float width = ImGui::GetContentRegionAvail().x;
  const layout::Frame frame = layout::measure(ImVec2(width, height));
  constexpr int kInstrumentCount = 6;
  const int instrumentColumns =
      std::clamp(layout::columnsThatFit(frame, 10.0f), 1, kInstrumentCount);
  const int instrumentRows =
      (kInstrumentCount + instrumentColumns - 1) / instrumentColumns;
  const float supportHeight = frame.row * 4.35f + frame.gap * 3.0f;
  const float gridHeight =
      std::max(height - supportHeight, frame.row);
  const float meterHeight = frame.em * 0.35f;
  const float tileHeight = std::max(
      (gridHeight -
       static_cast<float>(instrumentRows) * (meterHeight + frame.gap)) /
          static_cast<float>(instrumentRows),
      frame.em);
  constexpr ImGuiTableFlags instrumentFlags =
      ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_NoSavedSettings;
  if (ImGui::BeginTable("##fluid_instruments", instrumentColumns,
                        instrumentFlags)) {
    ImGui::TableNextColumn();
    charts::instrument(
        "##compression", "Worst compression", compressionValue,
        compressionAvailable ? "%" : "Unavailable", traces.compression,
        ImVec2(ImGui::GetContentRegionAvail().x, tileHeight), compressionStyle);
    charts::MeterStyle compressionMeter;
    compressionMeter.accent = style::col::Data;
    const double compressionLimit =
        kFluidPresets[resolutionIndex].quality.densityTolerance;
    compressionMeter.warnAt = 0.75;
    compressionMeter.dangerAt = 1.0;
    charts::meter(
        "##compression_headroom",
        compressionAvailable && compressionLimit > 0.0
            ? compression / compressionLimit
            : 0.0,
        ImVec2(ImGui::GetContentRegionAvail().x, meterHeight),
        compressionMeter);

    ImGui::TableNextColumn();
    charts::instrument(
        "##deficit", "Free-surface deficit", deficitValue,
        deficitAvailable ? "%" : "Unavailable", traces.deficit,
        ImVec2(ImGui::GetContentRegionAvail().x, tileHeight), deficitStyle);
    charts::MeterStyle deficitMeter;
    deficitMeter.accent = style::col::DataDim;
    deficitMeter.warnAt = 0.75;
    deficitMeter.dangerAt = 1.0;
    charts::meter("##deficit_headroom", deficitAvailable ? deficit : 0.0,
                  ImVec2(ImGui::GetContentRegionAvail().x, meterHeight),
                  deficitMeter);

    ImGui::TableNextColumn();
    charts::instrument(
        "##real_time_factor", "Real-time factor", rateValue,
        rateAvailable ? "x" : "Unavailable", traces.realTimeFactor,
        ImVec2(ImGui::GetContentRegionAvail().x, tileHeight), rateStyle);

    ImGui::TableNextColumn();
    charts::instrument(
        "##dispersed", "Dispersed phase", dispersedValue,
        dispersedAvailable ? "%" : "Unavailable", traces.dispersedFraction,
        ImVec2(ImGui::GetContentRegionAvail().x, tileHeight), dispersedStyle);

    ImGui::TableNextColumn();
    charts::instrument(
        "##sauter", "Sauter d32", sauterValue,
        !sauterAvailable
            ? "Unavailable"
            : (diagnostics.sauterDiameterM > 0.0 ? "um" : nullptr),
        traces.sauterDiameter,
        ImVec2(ImGui::GetContentRegionAvail().x, tileHeight), sauterStyle);

    ImGui::TableNextColumn();
    charts::instrument(
        "##interfacial_area", "Interfacial area", areaValue,
        areaAvailable ? "cm^2" : "Unavailable", traces.interfacialArea,
        ImVec2(ImGui::GetContentRegionAvail().x, tileHeight), areaStyle);
    ImGui::EndTable();
  }

  char substepsValue[32] = "Unavailable";
  if (completedStep && stats.substeps > 0)
    std::snprintf(substepsValue, sizeof(substepsValue), "%d", stats.substeps);
  char surfaceValue[32] = "Unavailable";
  if (surfaceAvailable)
    std::snprintf(surfaceValue, sizeof(surfaceValue), "%.1f mm",
                  diagnostics.freeSurfaceM * 1000.0);
  constexpr widgets::Column supportColumns[] = {
      {"Substeps", true, false, nullptr, 5.0f},
      {"Free surface", true, false, nullptr, 7.0f}};
  if (widgets::beginDataTable("##fluid_support", supportColumns, 2,
                              ImVec2(0.0f, frame.row * 2.0f))) {
    widgets::dataRow(style::col::DataDim);
    widgets::dataCell(substepsValue);
    widgets::dataCell(surfaceValue);
    widgets::endDataTable();
  }

  traces.phaseSegments.clear();
  const size_t phaseCount =
      completedStep && diagnostics.valid
          ? std::min(diagnostics.phases.size(), s.funnel.phases.size())
          : 0;
  traces.phaseSegments.reserve(phaseCount);
  for (size_t i = 0; i < phaseCount; ++i) {
    ImVec4 phaseColour = style::col::Data;
    if (i < snapshot->phases.size()) {
      const fluid::PhaseMaterial& phase = snapshot->phases[i];
      phaseColour =
          ImVec4(phase.colour[0], phase.colour[1], phase.colour[2],
                 phase.colour[3]);
    }
    const double bulkMl =
        diagnostics.phases[i].bulkResolved &&
                std::isfinite(diagnostics.phases[i].bulkMl)
            ? std::max(0.0, diagnostics.phases[i].bulkMl)
            : 0.0;
    traces.phaseSegments.push_back(
        {s.funnel.phases[i].label.c_str(), bulkMl, phaseColour});
  }
  charts::stackedBar(
      "##phase_bulk_chart", traces.phaseSegments.data(),
      static_cast<int>(traces.phaseSegments.size()),
      ImVec2(ImGui::GetContentRegionAvail().x, frame.row * 1.35f));

  const fluid::VesselMotion motion = shakeReportMotion(s);
  char shakeValue[96];
  std::snprintf(shakeValue, sizeof(shakeValue), "%.2f m/s peak  |  %.1f W/kg",
                fluid::shakePeakVelocity(motion),
                fluid::shakeSpecificPower(motion));
  widgets::keyValue("Driven shake", shakeValue, style::col::DataBright);
}

void drawTransportControls(SolubilityState& s) {
  sol::Simulation& charge = s.funnel;
  static const char* kVesselNames[] = {"Separatory funnel", "Decanting flask",
                                       "Graduated cylinder"};
  static const icons::Icon kVesselGlyphs[] = {
      icons::Icon::SepFunnel, icons::Icon::Flask, icons::Icon::Beaker};
  const layout::Frame frame = layout::measure();
  const int columns =
      std::clamp(layout::columnsThatFit(frame, 13.0f), 1, 4);
  constexpr ImGuiTableFlags flags =
      ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_NoSavedSettings;
  if (!ImGui::BeginTable("##transport_grid", columns, flags)) return;

  ImGui::TableNextColumn();
  ImGui::TextDisabled("Vessel");
  int vessel = std::clamp(s.funnelVessel, 0, 2);
  if (widgets::segmentedIcons("##vessel", kVesselGlyphs, kVesselNames, 3, vessel)) {
    s.funnelVessel = vessel;
    charge.vessel = static_cast<sol::Vessel>(vessel);
    sol::reset(charge);
    if (rechargeFluid(s))
      s.statusMessage = std::string("Recharged into ") + kVesselNames[vessel];
  }

  ImGui::TableNextColumn();
  ImGui::TextDisabled("Transport");
  const float actionSpacing = ImGui::GetStyle().ItemSpacing.x;
  const float actionWidth =
      std::max((ImGui::GetContentRegionAvail().x - actionSpacing) * 0.5f,
               ImGui::GetFrameHeight() * 2.0f);
  if (widgets::actionButton("##run_solver",
                            s.funnelRunning ? icons::Icon::Pause : icons::Icon::Play,
                            s.funnelRunning ? "Pause" : "Run",
                            ImVec2(actionWidth, 0.0f), !s.funnelRunning,
                            "Advance or hold the particle solver")) {
    s.funnelRunning = !s.funnelRunning;
  }
  ImGui::SameLine(0.0f, actionSpacing);
  if (widgets::actionButton("##reset_solver", icons::Icon::Rewind, "Reset",
                            ImVec2(actionWidth, 0.0f), false,
                            "Recharge the vessel and stand it upright")) {
    sol::reset(charge);
    const bool recharged = rechargeFluid(s);
    s.funnelRunning = false;
    s.fluidTiltTargetDeg = 0.0f;
    s.fluidTiltCurrentDeg = 0.0f;
    s.fluidTiltAngularVelocityRadS = 0.0f;
    if (recharged) s.statusMessage = "Particle vessel recharged";
  }
  widgets::statusDot(s.funnelRunning ? "Solver running" : "Solver paused",
                     s.funnelRunning, style::col::Data);
  if (!s.statusMessage.empty() && ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", s.statusMessage.c_str());

  ImGui::TableNextColumn();
  ImGui::TextDisabled("Simulation rate");
  widgets::glyphSlider("##speed", icons::Icon::Timer, "rate", s.funnelSpeed, 0.1f,
                       10.0f, "%.1fx",
                       "Simulated seconds per wall-clock second");

  ImGui::TableNextColumn();
  ImGui::TextDisabled("Solver quality");
  if (fluid::Simulation* simulation = s.fluid.get()) {
    double measured = std::numeric_limits<double>::quiet_NaN();
    if (runFluidInteraction(s, [&] { measured = simulation->realTimeFactor(); }))
      rememberFluidRate(s, measured);
  }

  size_t resolution = selectedFluidPreset(s);
  ImGui::SetNextItemWidth(-1.0f);
  if (ImGui::BeginCombo("##resolution", kFluidPresets[resolution].name)) {
    for (size_t candidate = 0; candidate < kFluidPresets.size(); ++candidate) {
      const std::string choice = fluidPresetTradeText(s, candidate);
      const bool selected = candidate == resolution;
      if (ImGui::Selectable(choice.c_str(), selected) && !selected) {
        s.fluidResolution = static_cast<FluidResolution>(candidate);
        resolution = candidate;
        if (rechargeFluid(s))
          s.statusMessage =
              std::string(kFluidPresets[candidate].name) +
              " solver quality charged";
      }
      if (selected) ImGui::SetItemDefaultFocus();
    }
    ImGui::EndCombo();
  }
  ImGui::EndTable();
}


void drawShakeControls(SolubilityState& s) {
  fluid::Simulation* simulation = availableFluid(s);
  bool shaking = false;
  double elapsedS = 0.0;
  double realTimeFactor = std::numeric_limits<double>::quiet_NaN();
  if (simulation) {
    runFluidInteraction(s, [&] {
      shaking = simulation->shaking();
      elapsedS = simulation->elapsedS();
      realTimeFactor = simulation->realTimeFactor();
    });
    rememberFluidRate(s, realTimeFactor);
    if (shaking && !s.fluidShakeProgressValid) {
      s.fluidShakeStartElapsedS = elapsedS;
      s.fluidShakeEndElapsedS =
          elapsedS + static_cast<double>(s.shakeDurationS);
      s.fluidShakeProgressValid = true;
    }
  }
  const layout::Frame frame = layout::measure();
  const int columns =
      std::clamp(layout::columnsThatFit(frame, 10.0f), 1, 5);
  constexpr ImGuiTableFlags flags =
      ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_NoSavedSettings;
  if (ImGui::BeginTable("##shake_grid", columns, flags)) {
    ImGui::TableNextColumn();
    ImGui::TextDisabled("Action");
    ImGui::BeginDisabled(simulation == nullptr);
    const bool shakeRequested =
        animatedShakeButton(shaking, textButtonSize("Shake"));
    ImGui::EndDisabled();
    // ImGui suppresses hover on a disabled item, so the button cannot carry
    // its own explanation here; without one, a greyed Shake is indistinguishable
    // from a broken one.
    if (simulation == nullptr &&
        ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
      ImGui::SetTooltip("%s", fluidAbsenceReason(s));
    if (shakeRequested && simulation) {
      const std::array<double, 3> axis = selectedShakeAxis(s.shakeAxis);
      const double amplitudeM = static_cast<double>(s.shakeAmplitudeCm) * 0.01;
      const bool started = runFluidInteraction(s, [&] {
        simulation->shake(axis, s.shakeDurationS, s.shakeFrequencyHz, amplitudeM);
      });
      if (started) {
        s.funnelRunning = true;
        shaking = true;
        s.fluidShakeStartElapsedS = elapsedS;
        s.fluidShakeEndElapsedS =
            elapsedS + static_cast<double>(s.shakeDurationS);
        s.fluidShakeProgressValid = true;
        s.statusMessage =
            std::string(kShakeAxisNames[static_cast<size_t>(s.shakeAxis)]) +
            " driven shake started";
      }
    }

    ImGui::TableNextColumn();
    ImGui::TextDisabled("Axis");
    int axis = static_cast<int>(s.shakeAxis);
    if (widgets::segmented("##shake_axis", kShakeAxisNames.data(),
                           static_cast<int>(kShakeAxisNames.size()), axis)) {
      s.shakeAxis = static_cast<FluidShakeAxis>(
          std::clamp(axis, 0, static_cast<int>(kShakeAxisNames.size()) - 1));
    }

    ImGui::TableNextColumn();
    ImGui::TextDisabled("Duration");
    widgets::glyphSlider("##duration", icons::Icon::Timer, "t", s.shakeDurationS, 1.0f,
                         30.0f, "%.0f s", "How long the driven shake runs");

    ImGui::TableNextColumn();
    ImGui::TextDisabled("Frequency");
    widgets::glyphSlider("##frequency", icons::Icon::Shake, "f", s.shakeFrequencyHz, 0.5f,
                         6.0f, "%.1f Hz", "Strokes per second; peak speed is 2*pi*f*A");

    ImGui::TableNextColumn();
    ImGui::TextDisabled("Amplitude");
    widgets::glyphSlider("##amplitude", icons::Icon::Ruler, "A", s.shakeAmplitudeCm, 1.0f,
                         15.0f, "%.0f cm", "Half-stroke of the driven oscillation");
    ImGui::EndTable();
  }
  if (simulation == nullptr) {
    // The rate message is about a solver that exists; when there is none, the
    // whole shake row is inert and that is what has to be said.
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextDisabled("Shake is unavailable: %s", fluidAbsenceReason(s));
    ImGui::PopTextWrapPos();
  } else if (!std::isfinite(realTimeFactor)) {
    ImGui::TextDisabled("Physics rate appears after the first completed step.");
  }

  if (widgets::onlyWhen(
          shaking && s.fluidShakeProgressValid,
          "Shake progress appears while a driven shake is running")) {
    const double duration =
        std::max(s.fluidShakeEndElapsedS - s.fluidShakeStartElapsedS, 0.0);
    const double remaining =
        std::clamp(s.fluidShakeEndElapsedS - elapsedS, 0.0, duration);
    const float progress =
        duration > 0.0
            ? static_cast<float>(std::clamp(1.0 - remaining / duration,
                                            0.0, 1.0))
            : 1.0f;
    char progressLabel[64];
    std::snprintf(progressLabel, sizeof(progressLabel),
                  "%.1f simulated s remaining", remaining);
    ImGui::ProgressBar(progress, ImVec2(-1.0f, 0.0f), progressLabel);
  }
}

void setTiltTarget(SolubilityState& s, float degrees) {
  s.fluidTiltTargetDeg = std::clamp(degrees, 0.0f, 180.0f);
  s.funnelRunning = true;
}

void drawTiltControls(SolubilityState& s) {
  ImGui::Separator();
  ImGui::TextDisabled("Tilt and invert");
  const layout::Frame frame = layout::measure();
  const int columns = std::min(layout::columnsThatFit(frame, 18.0f), 2);
  constexpr ImGuiTableFlags flags =
      ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_NoSavedSettings;
  if (!ImGui::BeginTable("##tilt_grid", columns, flags)) return;
  ImGui::TableNextColumn();
  if (widgets::glyphSlider("##tilt", icons::Icon::Gauge, "tilt",
                           s.fluidTiltTargetDeg, 0.0f, 180.0f, "%.0f deg",
                           "Vessel attitude: 0 upright, 180 fully inverted"))
    s.funnelRunning = true;

  ImGui::TableNextColumn();
  const float spacing = ImGui::GetStyle().ItemSpacing.x;
  const float presetWidth =
      std::max((ImGui::GetContentRegionAvail().x - spacing * 2.0f) / 3.0f,
               ImGui::GetFrameHeight() * 1.8f);
  const ImVec2 presetSize(presetWidth, 0.0f);
  if (widgets::actionButton("##tilt_upright", icons::Icon::ChevronUp, "Upright",
                            presetSize, false, "Stand the vessel up (0 deg)"))
    setTiltTarget(s, 0.0f);
  ImGui::SameLine(0.0f, spacing);
  if (widgets::actionButton("##tilt_vent", icons::Icon::ChevronRight, "Vent",
                            presetSize, false,
                            "Tip to the venting attitude (135 deg)"))
    setTiltTarget(s, 135.0f);
  ImGui::SameLine(0.0f, spacing);
  if (widgets::actionButton("##tilt_invert", icons::Icon::ChevronDown, "Invert",
                            presetSize, false,
                            "Fully invert to drain from the neck (180 deg)"))
    setTiltTarget(s, 180.0f);
  ImGui::EndTable();
}

void drawRunTab(SolubilityState& s, float height) {
  const layout::Frame frame =
      layout::measure(ImVec2(ImGui::GetContentRegionAvail().x, height));
  const float weights[] = {0.42f, 0.58f};
  const float minimums[] = {frame.control * 2.0f, frame.control * 3.0f};
  float rows[2]{};
  layout::distribute(height, weights, minimums, 2, frame.gap, rows);
  const float startY = ImGui::GetCursorPosY();

  if (widgets::beginCard(
          "##vessel_card", ImVec2(0.0f, rows[0]), style::col::BgSurface,
          ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    widgets::sectionHeader("Vessel", style::col::Data);
    drawTransportControls(s);
    widgets::endCard();
  }
  layout::nextRow(startY + rows[0] + frame.gap);
  if (widgets::beginCard(
          "##motion_card", ImVec2(0.0f, rows[1]), style::col::BgSurface,
          ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    widgets::sectionHeader("Shake and tilt", style::col::Data);
    drawShakeControls(s);
    drawTiltControls(s);
    widgets::endCard();
  }
}

// Defined with the partition maths further down: removing a phase renumbers
// every phase after it, so a manual aqueous pick made against the old numbering
// has to be dropped rather than silently re-pointed at a different liquid.
void forgetAqueousPick();

void drawPhaseTable(SolubilityState& s, bool& changed, float height) {
  sol::Simulation& sim = s.funnel;
  constexpr widgets::Column columns[] = {
      {"Phase", false, true, nullptr, 8.0f},
      {"Volume", true, false, "mL", 5.0f},
      {"Density", true, false, "g/mL", 5.0f},
      {"Viscosity", true, false, "mPa.s", 5.0f},
      {"Interfacial tension", true, false, "mN/m", 7.0f},
      {"Stability", true, false, nullptr, 5.0f},
      {"", false, false, nullptr, 3.0f}};
  if (!widgets::beginDataTable("##phase_table", columns, 7,
                               ImVec2(0.0f, height)))
    return;

  // Turns the cell just drawn into a click target. The cursor is deliberately
  // NOT restored afterwards: widgets::dataCell opens the next cell with
  // TableNextColumn, which repositions it anyway, and a bare SetCursorScreenPos
  // with no item behind it is what ImGui reports as "code uses SetCursorPos()
  // to extend window/parent boundaries" -- an assert in a debug build.
  auto cellActivator = [](const char* id) {
    const ImVec2 min = ImGui::GetItemRectMin();
    const ImVec2 max = ImGui::GetItemRectMax();
    ImGui::SetCursorScreenPos(min);
    ImGui::InvisibleButton(id, ImVec2(max.x - min.x, max.y - min.y));
    const bool activated = ImGui::IsItemClicked();
    if (ImGui::IsItemHovered())
      ImGui::SetMouseCursor(ImGuiMouseCursor_Hand);
    return activated;
  };

  int removeIndex = -1;
  for (size_t i = 0; i < sim.phases.size(); ++i) {
    sol::Phase& phase = sim.phases[i];
    ImGui::PushID(static_cast<int>(i));
    widgets::dataRow(
        ImVec4(phase.colour[0], phase.colour[1], phase.colour[2], phase.colour[3]));

    widgets::dataCell(phase.label.c_str());
    if (cellActivator("##edit_phase_name")) ImGui::OpenPopup("##phase_editor");
    widgets::dataCellf("%.1f", phase.volumeMl);
    if (cellActivator("##edit_phase_volume")) ImGui::OpenPopup("##phase_editor");
    widgets::dataCellf("%.3f", phase.density);
    if (cellActivator("##edit_phase_density")) ImGui::OpenPopup("##phase_editor");
    widgets::dataCellf("%.2f", phase.viscosity);
    if (cellActivator("##edit_phase_viscosity")) ImGui::OpenPopup("##phase_editor");
    widgets::dataCellf("%.1f", phase.interfacialTension);
    if (cellActivator("##edit_phase_tension")) ImGui::OpenPopup("##phase_editor");
    widgets::dataCellf("%.2f", phase.emulsionStability);
    if (cellActivator("##edit_phase_stability")) ImGui::OpenPopup("##phase_editor");

    widgets::dataCell("");
    const ImVec2 cellMin = ImGui::GetItemRectMin();
    const float controlSize = ImGui::GetFontSize();
    // Same reason as cellActivator: the row is closed by dataRow/endDataTable,
    // both of which reposition the cursor, so restoring it here would only
    // trip ImGui's parent-boundary check.
    ImGui::SetCursorScreenPos(cellMin);
    ImGui::InvisibleButton("##phase_swatch", ImVec2(controlSize, controlSize));
    const ImVec2 swatchMin = ImGui::GetItemRectMin();
    const ImVec2 swatchMax = ImGui::GetItemRectMax();
    ImGui::GetWindowDrawList()->AddRectFilled(
        swatchMin, swatchMax,
        style::u32(ImVec4(phase.colour[0], phase.colour[1], phase.colour[2],
                          phase.colour[3])),
        style::metrics().radiusSm);
    ImGui::GetWindowDrawList()->AddRect(
        swatchMin, swatchMax, style::u32(style::col::GridLine),
        style::metrics().radiusSm, 0, style::metrics().hairline);
    if (ImGui::IsItemClicked()) ImGui::OpenPopup("##phase_editor");
    ImGui::SameLine(0.0f, style::metrics().gap * 0.25f);
    // Emptying the table is not a thing the panel can honour: drawExtractionLab
    // re-seeds its default charge whenever no phase is left, so removing the
    // last one would look like the button swapped the liquid for a different
    // pair rather than removing it.
    const bool lastPhase = sim.phases.size() <= 1;
    ImGui::BeginDisabled(lastPhase);
    if (widgets::iconButton("##remove_phase", icons::Icon::Trash,
                            ImVec2(controlSize, controlSize), false,
                            "Remove phase"))
      removeIndex = static_cast<int>(i);
    ImGui::EndDisabled();
    if (lastPhase && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
      ImGui::SetTooltip("The vessel needs at least one charged phase.");

    if (ImGui::BeginPopup("##phase_editor")) {
      ImGui::TextUnformatted("Edit charged phase");
      ImGui::Separator();
      ImGui::SetNextItemWidth(ImGui::GetFontSize() * 18.0f);
      if (phaseLabelInput("Phase", phase.label)) changed = true;
      float value = static_cast<float>(phase.volumeMl);
      if (ImGui::DragFloat("Volume (mL)", &value, 1.0f, 0.0f, 5000.0f,
                           "%.1f")) {
        phase.volumeMl = std::max(value, 0.0f);
        changed = true;
      }
      value = static_cast<float>(phase.density);
      if (ImGui::DragFloat("Density (g/mL)", &value, 0.005f, 0.10f, 3.50f,
                           "%.3f")) {
        phase.density = std::max(value, 0.01f);
        changed = true;
      }
      value = static_cast<float>(phase.viscosity);
      if (ImGui::DragFloat("Viscosity (mPa.s)", &value, 0.01f, 0.05f, 500.0f,
                           "%.2f")) {
        phase.viscosity = std::max(value, 0.01f);
        changed = true;
      }
      value = static_cast<float>(phase.interfacialTension);
      if (ImGui::DragFloat("Interfacial tension (mN/m)", &value, 0.5f, 0.0f,
                           100.0f, "%.1f")) {
        phase.interfacialTension = std::max(value, 0.0f);
        changed = true;
      }
      value = static_cast<float>(phase.emulsionStability);
      if (ImGui::SliderFloat("Stability", &value, 0.0f, 1.0f, "%.2f")) {
        phase.emulsionStability = value;
        changed = true;
      }
      if (ImGui::ColorEdit4(
              "Colour", phase.colour,
              ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_AlphaBar))
        changed = true;
      ImGui::EndPopup();
    }
    ImGui::PopID();
  }
  widgets::endDataTable();

  if (removeIndex >= 0) {
    sim.phases.erase(sim.phases.begin() + removeIndex);
    forgetAqueousPick();
    changed = true;
  }
}

void drawPhaseEditor(SolubilityState& s, float height) {
  bool changed = false;
  const layout::Frame frame =
      layout::measure(ImVec2(ImGui::GetContentRegionAvail().x, height));
  std::vector<charts::StackSegment> segments;
  segments.reserve(s.funnel.phases.size());
  for (const sol::Phase& phase : s.funnel.phases) {
    segments.push_back(
        {phase.label.c_str(), phase.volumeMl,
         ImVec4(phase.colour[0], phase.colour[1], phase.colour[2],
                phase.colour[3])});
  }
  const float barHeight = frame.row * 1.35f;
  charts::stackedBar("##charge_volume_split", segments.data(),
                     static_cast<int>(segments.size()),
                     ImVec2(ImGui::GetContentRegionAvail().x, barHeight));

  const float tableHeight = std::max(
      height - barHeight - frame.control - frame.gap * 2.0f, frame.row * 2.0f);
  drawPhaseTable(s, changed, tableHeight);

  if (widgets::actionButton("##add_phase", icons::Icon::Plus, "Add phase",
                            ImVec2(0.0f, 0.0f), false)) {
    sol::Phase phase;
    phase.label = "Phase " + std::to_string(s.funnel.phases.size() + 1);
    phase.volumeMl = 50.0;
    phase.density = 1.0;
    phase.viscosity = 1.0;
    phase.interfacialTension = 30.0;
    phase.emulsionStability = 0.3;
    phase.colour[0] = style::col::Teal.x;
    phase.colour[1] = style::col::Teal.y;
    phase.colour[2] = style::col::Teal.z;
    phase.colour[3] = style::col::Teal.w;
    s.funnel.phases.push_back(phase);
    changed = true;
  }
  ImGui::SameLine(0.0f, frame.gap);
  ImGui::TextDisabled("Densest phase settles to the bottom");

  if (changed) {
    sol::reset(s.funnel);
    const bool recharged = rechargeFluid(s);
    if (recharged) {
      s.statusMessage = "Charged phases updated";
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

void forgetAqueousPick() { aqueousPick = -1; }

struct PartitionContext {
  bool valid = false;
  int aqueousIndex = 0;
  double volAq = 0.0;
  double volOrg = 0.0;
  // Fraction of the solute that one equilibration against `volOrg` of fresh
  // organic leaves in the aqueous layer, V_aq / (D V_org + V_aq). Read out of
  // sol::partition rather than re-derived from 10^logP here: partition() bounds
  // D to six decades, so deriving it twice made the wash chart and the split
  // bar report different distributions for a solute with |logP| > 6.
  double aqueousRemainder = 1.0;
  std::string organicLabel;
};

PartitionContext partitionContext(const SolubilityState& s) {
  PartitionContext ctx;
  if (!s.soluteValid || s.funnel.phases.size() < 2) return ctx;
  // A solute whose logP failed to resolve has no distribution to report; the
  // callers' onlyWhen() lines then explain the absence instead of the panel
  // printing "nan%" across the wash chart.
  if (!std::isfinite(s.solute.logP)) return ctx;
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
  int organicCount = 0;
  for (size_t i = 0; i < count; ++i) {
    if (static_cast<int>(i) == ctx.aqueousIndex) continue;
    ctx.volOrg += sim.phases[i].volumeMl;
    if (++organicCount == 1) ctx.organicLabel = sim.phases[i].label;
  }
  // Every non-aqueous layer is pooled into one organic volume, so naming that
  // pool after the first of them would put the whole extracted mass in a layer
  // that only holds part of it.
  if (organicCount > 1)
    ctx.organicLabel = "Organic (" + std::to_string(organicCount) + " phases)";
  ctx.valid = ctx.volAq > 0.0 && ctx.volOrg > 0.0;
  if (ctx.valid) {
    const sol::Partition unitMass =
        sol::partition(1.0, s.solute.logP, ctx.volAq, ctx.volOrg);
    ctx.aqueousRemainder = std::clamp(unitMass.mgAqueous, 0.0, 1.0);
  }
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
  // The charge is a different pair of liquids now, so a pick made against the
  // old one would keep pointing at a phase index that means something else.
  forgetAqueousPick();
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
  if (!ctx.valid) return;
  const sol::Simulation& sim = s.funnel;
  const size_t count = sim.phases.size();
  const layout::Frame frame = layout::measure();

  const int controlColumns =
      std::min(layout::columnsThatFit(frame, 18.0f), 2);
  constexpr ImGuiTableFlags controlFlags =
      ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_NoSavedSettings;
  if (ImGui::BeginTable("##distribution_controls", controlColumns,
                        controlFlags)) {
    ImGui::TableNextColumn();
    const float aqueousLabelWidth =
        ImGui::CalcTextSize("Aqueous phase").x +
        ImGui::GetStyle().ItemInnerSpacing.x;
    const float comboWidth =
        std::max(frame.em, ImGui::GetContentRegionAvail().x - aqueousLabelWidth);
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
    widgets::glyphSlider("##solute_mass", icons::Icon::Molecule, "solute mass",
                         soluteMassMgUi, 1.0f, 1000.0f, "%.0f mg");
    ImGui::EndTable();
  }

  const sol::Partition p =
      sol::partition(static_cast<double>(soluteMassMgUi), s.solute.logP,
                     ctx.volAq, ctx.volOrg);
  charts::StackSegment split[] = {
      {sim.phases[static_cast<size_t>(ctx.aqueousIndex)].label.c_str(),
       p.mgAqueous, style::col::Teal},
      {ctx.organicLabel.c_str(), p.mgOrganic, style::col::Data}};
  charts::stackedBar("##solute_split", split, 2,
                     ImVec2(ImGui::GetContentRegionAvail().x, frame.row * 1.6f));

  char partitionValue[128];
  std::snprintf(partitionValue, sizeof(partitionValue),
                "%.1f%% into %s  |  logP %.2f", p.fractionOrganic * 100.0,
                ctx.organicLabel.c_str(), s.solute.logP);
  widgets::keyValue(s.solute.name.c_str(), partitionValue, style::col::DataBright);
  ImGui::TextDisabled("Neutral-species logP approximation; no pH correction");
}

// ------------------------------------------------ multi-stage extraction
// The classic counter-current question: how many washes to strip the solute?
// Each wash with fresh organic removes the same fraction, so the aqueous
// remainder after n washes is q^n with q = V_aq / (D V_org + V_aq).

void drawMultiStageExtraction(const SolubilityState& s, float height) {
  const PartitionContext ctx = partitionContext(s);
  if (!ctx.valid) return;
  const double q = ctx.aqueousRemainder;
  const double perWash = 1.0 - q;
  const layout::Frame frame =
      layout::measure(ImVec2(ImGui::GetContentRegionAvail().x, height));

  char recoveryValue[32];
  std::snprintf(recoveryValue, sizeof(recoveryValue), "%.1f%%",
                perWash * 100.0);
  widgets::keyValue("Per-wash recovery", recoveryValue,
                    style::col::DataBright);

  constexpr int kMaxWashes = 6;
  charts::BarRow rows[kMaxWashes]{};
  char labels[kMaxWashes][16]{};
  char annotations[kMaxWashes][24]{};
  int recommended = 0;
  for (int n = 1; n <= kMaxWashes; ++n) {
    const double recovered = 1.0 - std::pow(q, n);
    if (recommended == 0 && recovered >= 0.99) recommended = n;
    std::snprintf(labels[n - 1], sizeof(labels[n - 1]), "%d wash%s", n,
                  n == 1 ? "" : "es");
    std::snprintf(annotations[n - 1], sizeof(annotations[n - 1]), "%.1f%%",
                  recovered * 100.0);
    rows[n - 1] = {labels[n - 1], recovered, annotations[n - 1],
                   n == recommended ? style::col::DataBright
                                    : style::col::Data,
                   false};
  }

  const float resultHeight = frame.row * 2.0f + frame.gap;
  const float chartHeight =
      std::max(height - resultHeight - frame.row, frame.row * 3.0f);
  charts::rankedBars("##wash_recovery", rows, kMaxWashes,
                     ImVec2(ImGui::GetContentRegionAvail().x, chartHeight));

  if (recommended > 0) {
    ImGui::TextWrapped("%d wash%s with fresh %s recovers >= 99%% of the solute.",
                       recommended, recommended == 1 ? "" : "es",
                       ctx.organicLabel.c_str());
  } else {
    ImGui::TextWrapped(
        "Even 6 washes leave > 1%% behind (q = %.3f); raise the organic "
        "volume or choose a stronger partitioning solvent.",
        q);
  }
}

void drawChargeTab(SolubilityState& s, float height) {
  const layout::Frame frame =
      layout::measure(ImVec2(ImGui::GetContentRegionAvail().x, height));
  const float weights[] = {0.62f, 0.38f};
  const float minimums[] = {frame.control * 4.0f, frame.control * 3.0f};
  float rows[2]{};
  layout::distribute(height, weights, minimums, 2, frame.gap, rows);
  const float startY = ImGui::GetCursorPosY();

  if (widgets::beginCard(
          "##charge_card", ImVec2(0.0f, rows[0]), style::col::BgSurface,
          ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    widgets::sectionHeader("Charge", style::col::Data);
    drawPhaseEditor(s, layout::pageHeight());
    widgets::endCard();
  }
  layout::nextRow(startY + rows[0] + frame.gap);
  if (widgets::beginCard(
          "##partition_card", ImVec2(0.0f, rows[1]), style::col::BgSurface,
          ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    widgets::sectionHeader("Partitioning", style::col::Teal);
    const PartitionContext ctx = partitionContext(s);
    if (widgets::onlyWhen(
            ctx.valid,
            "Partitioning needs a valid solute and two charged phases"))
      drawSoluteDistribution(s);
    widgets::endCard();
  }
}

void drawAnalyseTab(SolubilityState& s, float height) {
  const layout::Frame frame =
      layout::measure(ImVec2(ImGui::GetContentRegionAvail().x, height));
  const float weights[] = {0.66f, 0.34f};
  const float minimums[] = {frame.control * 6.0f, frame.control * 4.0f};
  float rows[2]{};
  layout::distribute(height, weights, minimums, 2, frame.gap, rows);
  const float startY = ImGui::GetCursorPosY();

  if (widgets::beginCard(
          "##physics_card", ImVec2(0.0f, rows[0]), style::col::BgSurface,
          ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    widgets::sectionHeader("Live physics", style::col::Data);
    drawFluidDiagnostics(s, layout::pageHeight());
    widgets::endCard();
  }
  layout::nextRow(startY + rows[0] + frame.gap);
  if (widgets::beginCard(
          "##wash_card", ImVec2(0.0f, rows[1]), style::col::BgSurface,
          ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    widgets::sectionHeader("Wash planning", style::col::Data);
    const PartitionContext ctx = partitionContext(s);
    if (widgets::onlyWhen(
            ctx.valid,
            "Wash planning needs a valid solute and two charged phases"))
      drawMultiStageExtraction(s, layout::pageHeight());
    widgets::endCard();
  }
}


void drawStageToolbar(SolubilityState& s) {
  // The two stages answer different questions -- a solid is what you show
  // someone, a section is what you measure -- so they are one exclusive switch
  // rather than two buttons that happen to be mutually exclusive.
  static const icons::Icon kStageGlyphs[2] = {icons::Icon::Cube, icons::Icon::Layers};
  static const char* kStageTips[2] = {
      "3D fluid: screen-space surface reconstruction of the particle field",
      "2D schematic: exact x-z section through the vessel axis"};
  int mode = s.extractionRenderMode == ExtractionRenderMode::Schematic2D ? 1 : 0;

  widgets::beginToolbar("##stage_modes");
  widgets::segmentedIcons("##stage_mode", kStageGlyphs, kStageTips, 2, mode,
                          ImGui::GetFontSize() * 7.0f);
  s.extractionRenderMode =
      mode == 1 ? ExtractionRenderMode::Schematic2D : ExtractionRenderMode::Fluid3D;
  ImGui::SameLine();
  widgets::toolbarSeparator();
  ImGui::SameLine();
  const bool reframable = s.extractionRenderMode == ExtractionRenderMode::Fluid3D;
  ImGui::BeginDisabled(!reframable);
  if (widgets::actionButton("##reframe_stage", icons::Icon::Crosshair, "Re-frame",
                            ImVec2(0.0f, 0.0f), false,
                            "Fit the vessel to the stage; double-clicking the stage "
                            "does the same thing")) {
    s.fluidReframeRequested = true;
  }
  ImGui::EndDisabled();
  // A disabled item is never "hovered", so actionButton's own tooltip cannot
  // reach the user here -- and a greyed control whose reason is unreadable is
  // indistinguishable from a broken one.
  if (!reframable && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
    ImGui::SetTooltip(
        "Re-frame moves the 3D camera. The schematic is a fixed section with no "
        "camera, and always fits the vessel to the stage.");
  }
  widgets::endToolbar();

  // Grab is no longer a mode: the left button always shakes, in both views, and
  // the vessel is free to leave the stage entirely while it is held.
  const char* hint =
      s.extractionRenderMode == ExtractionRenderMode::Schematic2D
          ? "Drag to shake -- fling it off the stage and let go to watch it "
            "swing back | particle cut: |y| < 0.75 particle spacings."
          : "Drag to shake -- fling it off the stage and let go to watch it "
            "swing back | right-drag orbit | wheel zoom | double-click re-frame.";
  const std::string fitted =
      ellipsizeText(hint, ImGui::GetContentRegionAvail().x);
  ImGui::TextDisabled("%s", fitted.c_str());
  if (ImGui::IsItemHovered() && fitted != hint) ImGui::SetTooltip("%s", hint);
}

// Pixels a metre of world motion covers on the stage, from the perspective
// projection at the orbit distance. Dragging then moves the vessel with the
// pointer instead of by an arbitrary screen-size-dependent amount.
double stagePixelsPerMetre(const gfx::Camera3D& camera, float canvasHeightPx) {
  constexpr double kPi = 3.14159265358979323846;
  const double halfHeightM = std::max(
      1.0e-6, static_cast<double>(camera.distanceM) *
                  std::tan(0.5 * static_cast<double>(camera.fovDeg) * kPi / 180.0));
  return 0.5 * static_cast<double>(canvasHeightPx) / halfHeightM;
}

// Advances the hand-follower one frame and publishes the resulting vessel
// motion. `handDelta` is this frame's commanded hand movement in world metres;
// it is zero on the frames the user is not dragging. The mechanics live in
// fluid::HandFollower, which is physics and is tested as such.
void advanceVesselShake(SolubilityState& s, fluid::Simulation& simulation,
                        const std::array<double, 3>& handDelta, double dt) {
  if (!(dt > 0.0)) return;
  s.fluidHand.advance(handDelta, s.fluidGrabActive, dt);
  runFluidInteraction(s, [&] {
    simulation.setManualMotion(s.fluidHand.position, s.fluidHand.acceleration);
  });
}

// Tracks the pointer while the left button is held and reports this frame's
// hand movement in screen pixels. Shared by both stage views.
std::array<double, 2> pollGrabDeltaPixels(SolubilityState& s, bool allowStart) {
  const ImGuiIO& io = ImGui::GetIO();
  if (allowStart && ImGui::IsItemActivated()) {
    s.fluidGrabActive = true;
    s.fluidGrabAnchorPx = {io.MousePos.x, io.MousePos.y};
    s.funnelRunning = true;  // a shake is meaningless against a paused solver
  }
  if (!ImGui::IsMouseDown(ImGuiMouseButton_Left)) s.fluidGrabActive = false;
  if (!s.fluidGrabActive) return {0.0, 0.0};
  const std::array<double, 2> delta{
      static_cast<double>(io.MousePos.x - s.fluidGrabAnchorPx[0]),
      static_cast<double>(io.MousePos.y - s.fluidGrabAnchorPx[1])};
  s.fluidGrabAnchorPx = {io.MousePos.x, io.MousePos.y};
  return delta;
}

void updateStageInput(SolubilityState& s, fluid::Simulation& simulation,
                      gfx::FluidStage& stage, const fluid::Snapshot& snapshot,
                      ImVec2 canvasSize) {
  const ImGuiIO& io = ImGui::GetIO();
  const double dt = std::clamp(static_cast<double>(io.DeltaTime), 0.0, 0.1);
  const bool hovered = ImGui::IsItemHovered();
  const bool active = ImGui::IsItemActive();
  const bool reframe =
      s.fluidReframeRequested ||
      (hovered && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left));
  s.fluidReframeRequested = false;

  if (reframe) {
    stage.camera.frame(snapshot.vesselHeightM, snapshot.maxRadiusM,
                       canvasSize.x / std::max(canvasSize.y, 1.0f));
    s.fluidGrabActive = false;
  } else if (hovered && io.MouseWheel != 0.0f) {
    stage.camera.zoom(io.MouseWheel);
  }

  // Orbit is on the right button so the left button can do the thing the panel
  // is actually about: picking the vessel up and shaking it.
  if ((hovered || active) && ImGui::IsMouseDragging(ImGuiMouseButton_Right)) {
    stage.camera.orbit(io.MouseDelta.x, io.MouseDelta.y);
  }

  const std::array<double, 2> deltaPx = pollGrabDeltaPixels(s, !reframe);
  std::array<double, 3> handDelta{0.0, 0.0, 0.0};
  if (s.fluidGrabActive && active) {
    const double metresPerPixel = 1.0 / stagePixelsPerMetre(stage.camera, canvasSize.y);
    // Screen right is the camera's horizontal axis in the world; screen up is
    // world +z, so a vertical drag shakes along the funnel axis.
    const std::array<float, 3> eye = stage.camera.eye();
    const double horizontalLength = std::max(
        std::hypot(static_cast<double>(eye[0]), static_cast<double>(eye[1])), 1.0e-9);
    const std::array<double, 3> cameraRight{-static_cast<double>(eye[1]) / horizontalLength,
                                            static_cast<double>(eye[0]) / horizontalLength,
                                            0.0};
    const double horizontal = deltaPx[0] * metresPerPixel;
    const double vertical = -deltaPx[1] * metresPerPixel;
    handDelta = {cameraRight[0] * horizontal, cameraRight[1] * horizontal, vertical};
  }
  advanceVesselShake(s, simulation, handDelta, dt);
}

bool drawPhysicsUnavailableState(SolubilityState& s, gfx::FluidStage& stage) {
  FluidBoundaryState& boundary = fluidBoundaryState(s);
  stage.requested = false;
  s.extractionRenderMode = ExtractionRenderMode::Schematic2D;
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
    if (fluidBoundaryState(s).unavailable) {
      drawPhysicsUnavailableState(s, stage);
    } else {
      // Construction is merely deferred, so nothing is wrong and the user's
      // render-mode preference must survive: drawPhysicsUnavailableState would
      // force the 2D schematic on and silently lose it. Draw the analytic
      // vessel meanwhile rather than a grey line of text -- the charge, the
      // layers and the emulsion are all known without the particle solver, so
      // there is no reason to show the user an empty box while it builds.
      stage.requested = false;
      stage.snapshot.reset();
      const ImVec2 canvasSize = ImGui::GetContentRegionAvail();
      if (canvasSize.x > 0.0f && canvasSize.y > 0.0f) {
        const ImVec2 rectMin = ImGui::GetCursorScreenPos();
        ImGui::Dummy(canvasSize);
        drawAnalyticFallback(s.funnel, rectMin, canvasSize);
      }
    }
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

  // The stage rect is where the vessel LIVES; the overlay rect is where it is
  // allowed to travel. While the vessel is in hand -- and until it has settled
  // back afterwards -- the render target grows to the whole application window
  // and its texture is drawn on the foreground draw list, so the funnel passes
  // over the neighbouring panels instead of being clipped at its dock edge.
  const ImVec2 stageMinPx = ImGui::GetCursorScreenPos();
  const ImGuiViewport* viewport = ImGui::GetMainViewport();
  const bool overlay = wants3D && (s.fluidGrabActive || !s.fluidHand.atRest());
  const ImVec2 targetMin = overlay ? viewport->Pos : stageMinPx;
  const ImVec2 targetSize = overlay ? viewport->Size : canvasSize;
  s.fluidOverlayActive = overlay;

  // A resized stage must re-fit, not crop: shrinking the window should zoom the
  // apparatus out until it fits again, which is what the eye expects from a
  // viewport and what keeps the vessel usable on a small display. Only a real
  // change counts, so a one-pixel dock jitter does not fight the user's zoom.
  if (wants3D && !s.fluidGrabActive) {
    const float widthChange = std::fabs(canvasSize.x - s.fluidStageSizePx[0]);
    const float heightChange = std::fabs(canvasSize.y - s.fluidStageSizePx[1]);
    const float threshold = std::max(2.0f, ImGui::GetFontSize() * 0.25f);
    if (widthChange > threshold || heightChange > threshold) {
      stage.camera.frame(snapshot->vesselHeightM, snapshot->maxRadiusM,
                         canvasSize.x / std::max(canvasSize.y, 1.0f));
      s.fluidStageSizePx = {canvasSize.x, canvasSize.y};
    }
  }

  if (wants3D) {
    const ImVec2 framebufferScale = ImGui::GetIO().DisplayFramebufferScale;
    stage.requested = true;
    stage.width = std::max(
        1, static_cast<int>(std::lround(targetSize.x * framebufferScale.x)));
    stage.height = std::max(
        1, static_cast<int>(std::lround(targetSize.y * framebufferScale.y)));
    stage.snapshot = snapshot;
    // Docked, the stage is a window onto a bench and wants its backdrop. In
    // hand it is stretched to the whole application, so an opaque backdrop
    // would black out every panel the funnel passes over -- which is the entire
    // screen, since the target IS the viewport. Dropping the backdrop's alpha
    // leaves the vessel and its contents composited over the workspace.
    stage.settings.backgroundAlpha = overlay ? 0.0f : 1.0f;

    // Growing the target must not resize the apparatus, so the vertical field
    // of view narrows by exactly the height ratio -- that holds pixels-per-metre
    // constant -- and the principal point shifts so the vessel's rest position
    // stays over the dock rect rather than jumping to the middle of the screen.
    if (overlay && canvasSize.y > 1.0f) {
      constexpr double kPi = 3.14159265358979323846;
      const double half = 0.5 * static_cast<double>(s.fluidStageFovDeg) * kPi / 180.0;
      const double scaled =
          std::tan(half) * static_cast<double>(canvasSize.y / targetSize.y);
      stage.camera.fovDeg =
          static_cast<float>(2.0 * std::atan(scaled) * 180.0 / kPi);
      const ImVec2 stageCentre(stageMinPx.x + canvasSize.x * 0.5f,
                               stageMinPx.y + canvasSize.y * 0.5f);
      stage.camera.shiftX =
          2.0f * (stageCentre.x - targetMin.x) / std::max(targetSize.x, 1.0f) - 1.0f;
      stage.camera.shiftY =
          1.0f - 2.0f * (stageCentre.y - targetMin.y) / std::max(targetSize.y, 1.0f);
    } else {
      stage.camera.fovDeg = s.fluidStageFovDeg;
      stage.camera.shiftX = 0.0f;
      stage.camera.shiftY = 0.0f;
    }
  } else {
    stage.requested = false;
    stage.snapshot.reset();
  }

  const ImVec2 cursor = ImGui::GetCursorPos();
  ImGui::InvisibleButton("##fluid_stage_input", canvasSize,
                         ImGuiButtonFlags_MouseButtonLeft |
                             ImGuiButtonFlags_MouseButtonRight);
  const ImVec2 rectMin = ImGui::GetItemRectMin();
  if (wants3D) {
    updateStageInput(s, *simulation, stage, *snapshot, canvasSize);
  } else {
    // The schematic is an exact x-z section, so the pointer maps straight onto
    // the vessel axes and shaking works here too -- there is no camera to
    // orbit, which makes this the most direct way to shake the funnel.
    const double dt =
        std::clamp(static_cast<double>(ImGui::GetIO().DeltaTime), 0.0, 0.1);
    const VesselGeometry& geometry = cachedGeometry(s.funnel);
    const Transform transform = buildTransform(geometry, rectMin, canvasSize);
    const double metresPerPixel = 1.0 / std::max(transform.scale, 1.0e-6f);
    const std::array<double, 2> deltaPx = pollGrabDeltaPixels(s, true);
    std::array<double, 3> handDelta{0.0, 0.0, 0.0};
    if (s.fluidGrabActive && ImGui::IsItemActive()) {
      handDelta = {deltaPx[0] * metresPerPixel, 0.0, -deltaPx[1] * metresPerPixel};
    }
    advanceVesselShake(s, *simulation, handDelta, dt);
  }

  // Input has just moved the hand, so the vessel's rigid translation is known
  // now, this frame. The driven shake stays physics-timed and comes from the
  // snapshot; only the hand is substituted. Both stages use the same value, so
  // the 3D and schematic views agree on where the glassware is.
  fluid::Pose livePose = snapshot->pose;
  for (std::size_t axis = 0; axis < 3; ++axis) {
    livePose.position[axis] = snapshot->shakeOffset[axis] + s.fluidHand.position[axis];
  }
  stage.pose = livePose;

  if (fluidBoundaryState(s).unavailable) {
    stage.requested = false;
    stage.snapshot.reset();
    s.fluidPresentedValid = false;
    drawAnalyticFallback(s.funnel, rectMin, canvasSize);
  } else if (wants3D && rendererReady) {
    // The texture on hand was rendered for LAST frame's target, so it is drawn
    // into last frame's rect. Blitting it into this frame's rect would stretch
    // the image for one frame every time the overlay engages or releases, which
    // is exactly when the eye is on it.
    if (s.fluidPresentedValid) {
      const ImVec2 presentedMin(s.fluidPresentedRectPx[0], s.fluidPresentedRectPx[1]);
      const ImVec2 presentedMax(s.fluidPresentedRectPx[2], s.fluidPresentedRectPx[3]);
      const bool presentedOverlay =
          presentedMin.x < stageMinPx.x - 0.5f || presentedMin.y < stageMinPx.y - 0.5f;
      if (presentedOverlay) {
        // Over the whole application: the funnel has left its dock, so it is
        // drawn on the foreground list, above every other panel. There is no
        // scrim. One used to be painted here at 55% BgDeep across the entire
        // presented rect -- and the presented rect while in hand IS the
        // viewport -- so grabbing the vessel blacked out the workspace behind
        // it. The stage's own backdrop is transparent in this mode, so what
        // travels over the bench is the apparatus and nothing else.
        ImDrawList* front = ImGui::GetForegroundDrawList();
        front->AddImage(ImTextureRef(static_cast<ImTextureID>(stage.texture)),
                        presentedMin, presentedMax, ImVec2(0.0f, 1.0f),
                        ImVec2(1.0f, 0.0f));
      } else {
        ImGui::SetCursorPos(cursor);
        ImGui::Image(ImTextureRef(static_cast<ImTextureID>(stage.texture)), canvasSize,
                     ImVec2(0.0f, 1.0f), ImVec2(1.0f, 0.0f));
      }
    }
    s.fluidPresentedRectPx = {targetMin.x, targetMin.y, targetMin.x + targetSize.x,
                              targetMin.y + targetSize.y};
    s.fluidPresentedValid = true;
  } else {
    s.fluidPresentedValid = false;
    drawCrossSection(*snapshot, s.funnel, livePose.position, rectMin, canvasSize);
  }
}
}  // namespace

void drawExtractionLab(AppState& st) {
  SolubilityState& s = st.solubility;
  s.fluidTasks = &st.tasks;
  // ImGui submits every docked panel on the frame the dock layout is built,
  // before it knows which tab is on top. Two consecutive draws are the first
  // trustworthy evidence that this expensive workspace is actually visible.
  const int frameCount = ImGui::GetFrameCount();
  const bool drawnLastFrame =
      s.extractionLastDrawnFrame == frameCount - 1;
  s.extractionLastDrawnFrame = frameCount;
  s.fluidConstructionAllowed =
      drawnLastFrame && st.tab == MainTab::Extraction;
  consumeExtractionImport(s);
  seedDefaultPhases(s.funnel);
  availableFluid(s);
  st.fluidStage.requested = false;

  const layout::Frame frame = layout::measure();
  int& selectedTab = s.extractionTab;
  static const char* tabLabels[] = {"Run", "Charge", "Analyse"};
  static const icons::Icon tabGlyphs[] = {
      icons::Icon::SepFunnel, icons::Icon::Droplet, icons::Icon::ChartBars};

  auto drawConsole = [&](ImVec2 size) {
    const bool visible = ImGui::BeginChild(
        "##extraction_console", size, ImGuiChildFlags_None,
        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
    if (visible) {
      widgets::subTabs("##extraction_tabs", tabLabels, tabGlyphs, 3,
                       selectedTab);
      const float contentHeight = layout::pageHeight();
      switch (std::clamp(selectedTab, 0, 2)) {
        case 0:
          drawRunTab(s, contentHeight);
          break;
        case 1:
          drawChargeTab(s, contentHeight);
          break;
        default:
          drawAnalyseTab(s, contentHeight);
          break;
      }
    }
    ImGui::EndChild();
  };

  auto drawStage = [&](ImVec2 size) {
    if (widgets::beginCard(
            "##fluid_stage_card", size, style::col::BgSurface,
            ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
      widgets::sectionHeader("Vessel stage", style::col::Data);
      drawFluidStage(st, s);
      widgets::endCard();
    }
  };

  const bool stageBeside = frame.wide && !frame.tall;
  if (stageBeside) {
    const float weights[] = {0.44f, 0.56f};
    const float minimums[] = {frame.em * 32.0f, frame.em * 30.0f};
    float columns[2]{};
    layout::distribute(frame.size.x, weights, minimums, 2, frame.gap,
                       columns);
    drawConsole(ImVec2(columns[0], frame.size.y));
    ImGui::SameLine(0.0f, frame.gap);
    stepSimulation(s);
    drawStage(ImVec2(columns[1], frame.size.y));
  } else {
    const float weights[] = {0.48f, 0.52f};
    const float minimums[] = {frame.control * 8.0f, frame.control * 9.0f};
    float rows[2]{};
    layout::distribute(frame.size.y, weights, minimums, 2, frame.gap, rows);
    const float startY = ImGui::GetCursorPosY();
    layout::nextRow(startY + rows[0] + frame.gap);
    drawConsole(ImVec2(frame.size.x, rows[1]));
    stepSimulation(s);
    ImGui::SetCursorPosY(startY);
    drawStage(ImVec2(frame.size.x, rows[0]));
  }

}

// Called once at startup so the first visit to the workspace finds the
// simulation already built. The build is the same one the panel would start;
// doing it early pays for the vessel SDF, the boundary tables and the initial
// charge on a worker thread while the user is still on the sketch canvas.
void warmExtractionPhysics(AppState& st) {
  SolubilityState& s = st.solubility;
  s.fluidTasks = &st.tasks;
  seedDefaultPhases(s.funnel);
  if (!s.fluid && !s.fluidBuildPending) startFluidBuild(s);
}

// The physics worker may hold a second OpenGL context, which it releases when
// its thread exits -- that is, when the last reference to the Simulation goes.
// AppState is a stack local in main() and outlives glfwTerminate(), so waiting
// for its destructor would destroy every context while another thread still
// had one current. Drain the in-flight build first: a queued build holds its
// own reference and would keep the simulation alive past this point.
void shutdownExtractionPhysics(AppState& st) {
  SolubilityState& s = st.solubility;
  while (st.tasks.busy()) {
    st.tasks.pump();
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  }
  st.tasks.pump();
  if (s.fluid) s.fluid->waitForIdle();
  s.fluid.reset();
  s.fluidTasks = nullptr;
  s.fluidBuildPending = false;
}

}  // namespace chemcad::ui

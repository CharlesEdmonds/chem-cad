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
  const double fraction = sol::emulsifiedFraction(sim);
  const float boxWidth = 190.0f;
  ImVec2 cursor(regionMin.x + regionSize.x - boxWidth, regionMin.y + 8.0f);
  const float lineHeight = ImGui::GetTextLineHeightWithSpacing();
  char buf[160];

  std::snprintf(buf, sizeof(buf), "t = %.1f s", sim.elapsed);
  draw->AddText(cursor, style::u32(style::col::TextDim), buf);
  cursor.y += lineHeight;

  std::snprintf(buf, sizeof(buf), "emulsified: %.1f%%", fraction * 100.0);
  draw->AddText(cursor, style::u32(style::col::TextDim), buf);
  cursor.y += lineHeight;

  ImFont* bigFont = style::fonts::semibold();
  const float bigSize = ImGui::GetFontSize() * 1.15f;
  draw->AddText(bigFont, bigSize, cursor, emulsionStateColor(fraction),
               emulsionStateLabel(fraction));
  cursor.y += lineHeight * 1.4f;

  const size_t layerCount = std::min(sim.phases.size(), sim.settledMl.size());
  for (size_t i = 0; i < layerCount; ++i) {
    std::snprintf(buf, sizeof(buf), "%s: %.0f mL", sim.phases[i].label.c_str(), sim.settledMl[i]);
    draw->AddText(cursor, style::u32(style::col::Text), buf);
    cursor.y += lineHeight;
  }
}

// --------------------------------------------------------------- painting
void drawVesselGlass(ImDrawList* draw, const VesselGeometry& geo, const Transform& tf) {
  static thread_local std::vector<ImVec2> screenOutline;
  screenOutline.clear();
  screenOutline.reserve(geo.outline.size());
  for (const core::Vec2& p : geo.outline) screenOutline.push_back(toScreen(tf, p));
  if (screenOutline.size() < 3) return;
  draw->AddConcavePolyFilled(screenOutline.data(), static_cast<int>(screenOutline.size()),
                             style::u32(ImVec4(0.55f, 0.72f, 0.85f, 1.0f), 0.10f));
}

void drawVesselWall(ImDrawList* draw, const VesselGeometry& geo, const Transform& tf) {
  static thread_local std::vector<ImVec2> screenOutline;
  screenOutline.clear();
  screenOutline.reserve(geo.outline.size());
  for (const core::Vec2& p : geo.outline) screenOutline.push_back(toScreen(tf, p));
  if (screenOutline.size() < 2) return;
  draw->AddPolyline(screenOutline.data(), static_cast<int>(screenOutline.size()),
                    style::u32(style::col::BorderStrong), ImDrawFlags_Closed, 2.0f);
}

void drawLayers(ImDrawList* draw, const sol::Simulation& sim, const VesselGeometry& geo,
                const Transform& tf) {
  static thread_local std::vector<ImVec2> band;
  const size_t n = std::min(sim.phases.size(), sim.settledMl.size());
  constexpr int kBandSamples = 20;

  double cursorMl = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const double volume = std::max(sim.settledMl[i], 0.0);
    if (volume <= 1e-9) continue;

    const double hfLo = heightFractionForVolume(geo, cursorMl);
    const double hfHi = heightFractionForVolume(geo, cursorMl + volume);
    cursorMl += volume;
    if (hfHi <= hfLo + 1e-6) continue;

    band.clear();
    for (int s = 0; s <= kBandSamples; ++s) {
      const double hf = hfLo + (hfHi - hfLo) * s / kBandSamples;
      const double halfW = widthFractionAt(sim.vessel, hf) * geo.halfWidthMetres;
      band.push_back(toScreen(tf, halfW, hf * geo.heightMetres));
    }
    for (int s = kBandSamples; s >= 0; --s) {
      const double hf = hfLo + (hfHi - hfLo) * s / kBandSamples;
      const double halfW = widthFractionAt(sim.vessel, hf) * geo.halfWidthMetres;
      band.push_back(toScreen(tf, -halfW, hf * geo.heightMetres));
    }

    const sol::Phase& phase = sim.phases[i];
    const ImU32 fill = ImGui::ColorConvertFloat4ToU32(
        ImVec4(phase.colour[0], phase.colour[1], phase.colour[2], phase.colour[3]));
    if (band.size() >= 3) {
      draw->AddConcavePolyFilled(band.data(), static_cast<int>(band.size()), fill);
    }

    // Bright interface line only where a non-empty layer actually sits above.
    if (i + 1 < n && sim.settledMl[i + 1] > 1e-9) {
      const double halfWTop = widthFractionAt(sim.vessel, hfHi) * geo.halfWidthMetres;
      const ImVec2 a = toScreen(tf, -halfWTop, hfHi * geo.heightMetres);
      const ImVec2 b = toScreen(tf, halfWTop, hfHi * geo.heightMetres);
      draw->AddLine(a, b, style::u32(style::col::Text, 0.7f), 1.5f);
    }
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
      const ImVec2 highlight(center.x - radius * 0.35f, center.y - radius * 0.35f);
      draw->AddCircleFilled(highlight, std::max(radius * 0.25f, 0.6f),
                            IM_COL32(255, 255, 255, 70), 6);
    }
  }
}

void drawCrossSection(const sol::Simulation& sim, ImVec2 regionMin, ImVec2 regionSize) {
  ImDrawList* draw = ImGui::GetWindowDrawList();
  const ImVec2 regionMax(regionMin.x + regionSize.x, regionMin.y + regionSize.y);
  draw->PushClipRect(regionMin, regionMax, true);
  draw->AddRectFilled(regionMin, regionMax, style::u32(style::col::BgSurface));

  const VesselGeometry& geo = cachedGeometry(sim.vessel, sim.vesselVolumeMl);
  const Transform tf = buildTransform(geo, regionMin, regionSize);

  drawVesselGlass(draw, geo, tf);
  drawGraduation(draw, geo, tf, regionMin);
  drawLayers(draw, sim, geo, tf);
  drawVesselWall(draw, geo, tf);
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

void drawControls(SolubilityState& s) {
  sol::Simulation& sim = s.funnel;
  static const char* kVesselNames[] = {"Separatory funnel", "Decanting flask",
                                       "Graduated cylinder"};

  ImGui::SetNextItemWidth(200.0f);
  if (ImGui::Combo("##vessel", &s.funnelVessel, kVesselNames, IM_ARRAYSIZE(kVesselNames))) {
    s.funnelVessel = std::clamp(s.funnelVessel, 0, 2);
    sim.vessel = static_cast<sol::Vessel>(s.funnelVessel);
    sol::reset(sim);
    s.statusMessage = std::string("Recharged into ") + kVesselNames[s.funnelVessel];
  }
  ImGui::SameLine();
  if (widgets::primaryButton("Shake")) {
    sol::shake(sim, static_cast<double>(s.shakeVigour));
    char buf[96];
    std::snprintf(buf, sizeof(buf), "Shaken at %.0f%% vigour", s.shakeVigour * 100.0f);
    s.statusMessage = buf;
  }
  ImGui::SameLine();
  ImGui::SetNextItemWidth(140.0f);
  ImGui::SliderFloat("Vigour", &s.shakeVigour, 0.0f, 1.0f, "%.2f");

  if (widgets::ghostButton(s.funnelRunning ? "Pause" : "Run")) s.funnelRunning = !s.funnelRunning;
  ImGui::SameLine();
  if (widgets::ghostButton("Reset")) {
    sol::reset(sim);
    s.funnelRunning = false;
    s.statusMessage = "Vessel reset";
  }
  ImGui::SameLine();
  ImGui::SetNextItemWidth(160.0f);
  ImGui::SliderFloat("Speed", &s.funnelSpeed, 0.1f, 10.0f, "%.1fx");
}

void drawPhaseTable(SolubilityState& s, bool& changed) {
  sol::Simulation& sim = s.funnel;
  constexpr ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                                    ImGuiTableFlags_SizingStretchProp |
                                    ImGuiTableFlags_NoSavedSettings;
  if (!ImGui::BeginTable("##phase_table", 8, flags)) return;
  ImGui::TableSetupColumn("Label", ImGuiTableColumnFlags_WidthStretch, 0.22f);
  ImGui::TableSetupColumn("mL", ImGuiTableColumnFlags_WidthStretch, 0.11f);
  ImGui::TableSetupColumn("g/mL", ImGuiTableColumnFlags_WidthStretch, 0.11f);
  ImGui::TableSetupColumn("mPa.s", ImGuiTableColumnFlags_WidthStretch, 0.11f);
  ImGui::TableSetupColumn("mN/m", ImGuiTableColumnFlags_WidthStretch, 0.11f);
  ImGui::TableSetupColumn("Stability", ImGuiTableColumnFlags_WidthStretch, 0.13f);
  ImGui::TableSetupColumn("Colour", ImGuiTableColumnFlags_WidthFixed, 40.0f);
  ImGui::TableSetupColumn("##remove", ImGuiTableColumnFlags_WidthFixed, 28.0f);
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
  widgets::sectionHeader("Charged phases");
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
  ImGui::SameLine();
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
  if (!s.soluteValid || s.funnel.phases.size() < 2) return;
  const sol::Simulation& sim = s.funnel;
  const size_t count = sim.phases.size();

  widgets::sectionHeader("Solute distribution");
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

  ImGui::SetNextItemWidth(200.0f);
  if (ImGui::BeginCombo("Aqueous phase", sim.phases[static_cast<size_t>(aq)].label.c_str())) {
    for (size_t i = 0; i < count; ++i) {
      const bool selected = static_cast<int>(i) == aq;
      if (ImGui::Selectable(sim.phases[i].label.c_str(), selected))
        aqueousPick = static_cast<int>(i);
      if (selected) ImGui::SetItemDefaultFocus();
    }
    ImGui::EndCombo();
  }
  ImGui::SameLine();
  ImGui::SetNextItemWidth(220.0f);
  ImGui::SliderFloat("Solute mass", &soluteMassMgUi, 1.0f, 1000.0f, "%.0f mg");

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

  // Stacked bar: aqueous share teal, organic share amber.
  const float avail = std::max(ImGui::GetContentRegionAvail().x, 60.0f);
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
  drawPhaseEditor(s);
  ImGui::Spacing();
  drawSoluteDistribution(s);
  ImGui::Spacing();

  stepSimulation(s);

  widgets::sectionHeader("Cross-section");
  ImVec2 size = ImGui::GetContentRegionAvail();
  size.x = std::max(size.x, 1.0f);
  size.y = std::max(size.y, 320.0f);
  ImGui::InvisibleButton("##funnel_canvas", size, ImGuiButtonFlags_None);
  const ImVec2 rectMin = ImGui::GetItemRectMin();
  const ImVec2 rectMax = ImGui::GetItemRectMax();

  drawCrossSection(s.funnel, rectMin, ImVec2(rectMax.x - rectMin.x, rectMax.y - rectMin.y));
}

}  // namespace chemcad::ui

#include "imgui.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

#include "chem/bridge.hpp"
#include "core/model.hpp"
#include "sol/solubility.hpp"
#include "ui/solubility_state.hpp"
#include "ui/theme.hpp"
#include "ui/ui.hpp"
#include "ui/widgets.hpp"

namespace chemcad::ui {
namespace {

// Sized from the current font rather than fixed pixels -- the app runs at a
// 1.25 UI scale by default, and hardcoded widths clip their own labels.
inline float uiScale() { return ImGui::GetFontSize() / 13.0f; }

int resizeStringInput(ImGuiInputTextCallbackData* data) {
  if (data->EventFlag != ImGuiInputTextFlags_CallbackResize) return 0;
  auto* value = static_cast<std::string*>(data->UserData);
  value->resize(static_cast<size_t>(data->BufTextLen));
  data->Buf = value->data();
  return 0;
}

bool stringInputWithHint(const char* label, const char* hint, std::string& value,
                         ImGuiInputTextFlags flags, bool mono = false) {
  flags |= ImGuiInputTextFlags_CallbackResize;
  const bool pushed = mono ? style::pushFont(style::fonts::mono()) : false;
  const bool result = ImGui::InputTextWithHint(label, hint, value.data(), value.capacity() + 1,
                                               flags, resizeStringInput, &value);
  style::popFont(pushed);
  return result;
}

bool containsCaseInsensitive(const std::string& haystack, const std::string& needle) {
  if (needle.empty()) return true;
  const auto lower = [](char c) { return std::tolower(static_cast<unsigned char>(c)); };
  const auto it = std::search(haystack.begin(), haystack.end(), needle.begin(), needle.end(),
                              [&](char a, char b) { return lower(a) == lower(b); });
  return it != haystack.end();
}

// Formats a g/mL value across the many orders of magnitude a solubility
// prediction can span: scientific notation below 1e-3, fixed otherwise.
std::string formatSolubility(double value) {
  char buf[64];
  const double magnitude = std::fabs(value);
  if (value == 0.0) {
    return "0";
  } else if (magnitude < 1e-3 || magnitude >= 1e6) {
    std::snprintf(buf, sizeof(buf), "%.3e", value);
  } else if (magnitude < 1.0) {
    std::snprintf(buf, sizeof(buf), "%.5f", value);
  } else {
    std::snprintf(buf, sizeof(buf), "%.4g", value);
  }
  return buf;
}

// A small perceptually-ordered ramp (violet -> blue -> teal -> green ->
// yellow), used for both the ternary fill and its legend.
ImVec4 colorRampF(float t) {
  static constexpr ImVec4 kStops[] = {
      {0.231f, 0.114f, 0.494f, 1.0f},
      {0.129f, 0.404f, 0.674f, 1.0f},
      {0.204f, 0.635f, 0.573f, 1.0f},
      {0.612f, 0.792f, 0.220f, 1.0f},
      {0.992f, 0.906f, 0.145f, 1.0f},
  };
  constexpr int n = static_cast<int>(sizeof(kStops) / sizeof(kStops[0]));
  t = std::clamp(t, 0.0f, 1.0f);
  const float scaled = t * static_cast<float>(n - 1);
  const int idx = std::min(n - 2, static_cast<int>(scaled));
  const float frac = scaled - static_cast<float>(idx);
  const ImVec4& c0 = kStops[idx];
  const ImVec4& c1 = kStops[idx + 1];
  return ImVec4(c0.x + (c1.x - c0.x) * frac, c0.y + (c1.y - c0.y) * frac,
               c0.z + (c1.z - c0.z) * frac, 1.0f);
}
ImU32 colorRamp(float t) { return style::u32(colorRampF(t)); }

// -------------------------------------------------------------- solute cache

void applySoluteOverrides(SolubilityState& sb) {
  if (!sb.soluteValid) return;
  sb.solute.meltingPoint = static_cast<double>(sb.meltingPointC);
  sb.solute.interactionRadius = static_cast<double>(sb.interactionRadius);
}

void recomputeSoluteFromSketch(AppState& st) {
  SolubilityState& sb = st.solubility;
  sb.soluteError.clear();
  sb.soluteValid = false;
  try {
    const core::Molecule* found = nullptr;
    for (const core::Molecule& mol : st.doc.molecules) {
      if (!mol.empty()) {
        found = &mol;
        break;
      }
    }
    if (!found) throw sol::SolError("Sketch is empty -- draw a structure first.");
    sb.solute = sol::describeSolute(*found);
    sb.soluteValid = true;
  } catch (const sol::SolError& err) {
    sb.soluteError = err.what();
  } catch (const chem::ChemError& err) {
    sb.soluteError = err.what();
  } catch (const std::exception& err) {
    sb.soluteError = err.what();
  }
  sb.sourceRevision = st.docRevision;
  if (sb.overrideSolute) applySoluteOverrides(sb);
  ++sb.soluteVersion;
}

void recomputeSoluteFromSmiles(AppState& st) {
  SolubilityState& sb = st.solubility;
  sb.soluteError.clear();
  sb.soluteValid = false;
  if (sb.manualSmiles.empty()) {
    sb.soluteError = "Enter a SMILES string.";
  } else {
    try {
      const core::Molecule mol = chem::fromSmiles(sb.manualSmiles);
      sb.solute = sol::describeSolute(mol);
      sb.soluteValid = true;
    } catch (const sol::SolError& err) {
      sb.soluteError = err.what();
    } catch (const chem::ChemError& err) {
      sb.soluteError = err.what();
    } catch (const std::exception& err) {
      sb.soluteError = err.what();
    }
  }
  if (sb.overrideSolute) applySoluteOverrides(sb);
  ++sb.soluteVersion;
}

// -------------------------------------------------------------- sweep cache

void recomputeSweepIfNeeded(SolubilityState& sb, const std::vector<const sol::Solvent*>& chosen) {
  std::string signature;
  signature.reserve(64);
  for (const sol::Solvent* solvent : chosen) {
    signature += solvent->id;
    signature += '|';
  }
  signature += std::to_string(sb.sweepSteps);
  signature += '|';
  signature += std::to_string(static_cast<int>(sb.temperatureC * 100.0f));
  signature += '|';
  signature += std::to_string(sb.soluteVersion);
  if (signature == sb.sweepSignature) return;
  sb.sweepSignature = signature;
  sb.sweepPeakIndex = -1;
  try {
    sb.sweep = sol::sweep(sb.solute, chosen, sb.sweepSteps, static_cast<double>(sb.temperatureC));
    // Cache the maximum-solubility sample once per sweep -- this is the
    // co-solvency peak the ratio plot marks, and finding it is work we do
    // not want to redo on every draw call of an unchanged frame.
    double best = -1.0;
    for (size_t i = 0; i < sb.sweep.size(); ++i) {
      const double v = sb.sweep[i].prediction.gramsPerMillilitre;
      if (v > best) {
        best = v;
        sb.sweepPeakIndex = static_cast<int>(i);
      }
    }
  } catch (const sol::SolError& err) {
    sb.sweep.clear();
    sb.statusMessage = std::string("Ratio sweep failed: ") + err.what();
  }
}

// ------------------------------------------------------------------ widgets

// Minimum height for widgets::statCard's stacked "value above caption"
// layout (mono value line, then the dim label line beneath it), mirroring
// the padding statCard itself uses internally so a card fed this height can
// never let the label baseline land under the value's at any
// CHEMCAD_UI_SCALE -- every term below scales in lockstep with the font, so
// the inequality that keeps them apart is scale-invariant.
float statCardHeight() {
  const style::Metrics& m = style::metrics();
  const float pad = m.gap * 0.9f;
  const float valueLineH = ImGui::GetFontSize();
  const float labelLineH = ImGui::GetFontSize() * 0.82f;
  return pad * 2.0f + valueLineH + labelLineH + m.gap * 0.5f;  // + a visible gap between lines
}

// Minimum height for drawHeadlineReadout's stacked "big value / sub-value* /
// caption" layout. Same reasoning as statCardHeight above.
float headlineCardHeight(int subLines) {
  const style::Metrics& m = style::metrics();
  const float pad = m.gap * 1.1f;
  const float valueFontSize = ImGui::GetFontSize() * 2.0f;
  const float subFontSize = ImGui::GetFontSize() * 0.95f;
  const float labelSize = ImGui::GetFontSize() * 0.82f;
  const float lineGap = m.gap * 0.4f;
  float h = pad * 2.0f + valueFontSize + lineGap + labelSize;
  h += static_cast<float>(subLines) * (subFontSize + lineGap);
  return h;
}

void drawHeadlineReadout(const char* label, const std::string& value,
                         const std::string& subValue, ImVec4 accent, ImVec2 size) {
  const style::Metrics& m = style::metrics();
  const ImVec2 min = ImGui::GetCursorScreenPos();
  ImGui::Dummy(size);
  const ImVec2 max(min.x + size.x, min.y + size.y);

  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddRectFilled(min, max, style::u32(style::col::BgRaised), m.radiusLg);
  dl->AddRect(min, max, style::u32(accent, 0.6f), m.radiusLg, 0, m.hairline);

  const float pad = m.gap * 1.1f;
  const float valueFontSize = ImGui::GetFontSize() * 2.0f;
  const float subFontSize = ImGui::GetFontSize() * 0.95f;
  const float labelSize = ImGui::GetFontSize() * 0.82f;
  const float lineGap = m.gap * 0.4f;

  float y = min.y + pad;
  dl->AddText(style::fonts::mono(), valueFontSize, ImVec2(min.x + pad, y),
             style::u32(style::col::Text), value.c_str());
  y += valueFontSize + lineGap;
  if (!subValue.empty()) {
    dl->AddText(style::fonts::mono(), subFontSize, ImVec2(min.x + pad, y), style::u32(accent),
               subValue.c_str());
  }
  dl->AddText(nullptr, labelSize, ImVec2(min.x + pad, max.y - pad - labelSize),
             style::u32(style::col::TextFaint), label);
}

void drawResultCard(SolubilityState& sb) {
  const sol::Prediction& p = sb.prediction;
  const float avail = ImGui::GetContentRegionAvail().x;
  const float spacing = ImGui::GetStyle().ItemSpacing.x;
  const float headlineW = std::min(300.0f * uiScale(), avail * 0.42f);

  // Size both the headline card and the stat grid from font metrics -- not
  // a hand-picked multiple of font size -- so the two rows of stat cards
  // always have room for a value line, a gap, and a caption line.
  const float statH = statCardHeight();
  const float rowGap = ImGui::GetStyle().ItemSpacing.y;
  const float statGridH = statH * 2.0f + rowGap;
  const float headlineH = std::max(headlineCardHeight(1), statGridH);

  const std::string headline = formatSolubility(p.gramsPerMillilitre) + " g/mL";
  const std::string molar = formatSolubility(p.molesPerLitre) + " mol/L";
  drawHeadlineReadout("SOLUBILITY", headline, molar, style::col::Accent,
                      ImVec2(headlineW, headlineH));

  ImGui::SameLine(0.0f, spacing);
  ImGui::BeginGroup();
  const float statW = (avail - headlineW - spacing * 3.0f) / 3.0f;
  char moleFrac[32];
  std::snprintf(moleFrac, sizeof(moleFrac), "%.4g", p.moleFraction);
  char idealFrac[32];
  std::snprintf(idealFrac, sizeof(idealFrac), "%.4g", p.idealMoleFraction);
  char ra[32];
  std::snprintf(ra, sizeof(ra), "%.2f", p.ra);
  char red[32];
  std::snprintf(red, sizeof(red), "%.2f", p.relativeEnergyDifference);
  char gamma[32];
  std::snprintf(gamma, sizeof(gamma), "%.3g", p.activityCoefficient);
  char chi[32];
  std::snprintf(chi, sizeof(chi), "%.3f", p.chi);
  widgets::statCard("MOLE FRACTION", moleFrac, ImVec2(statW, statH));
  ImGui::SameLine(0.0f, spacing);
  widgets::statCard("IDEAL MOLE FRAC", idealFrac, ImVec2(statW, statH));
  ImGui::SameLine(0.0f, spacing);
  widgets::statCard("Ra  MPa^0.5", ra, ImVec2(statW, statH));
  widgets::statCard("RED (Ra / R0)", red, ImVec2(statW, statH));
  ImGui::SameLine(0.0f, spacing);
  widgets::statCard("ACTIVITY COEF", gamma, ImVec2(statW, statH));
  ImGui::SameLine(0.0f, spacing);
  widgets::statCard("CHI (F-H)", chi, ImVec2(statW, statH));
  ImGui::EndGroup();

  if (p.outsideSphere) {
    ImGui::TextColored(style::col::Danger,
                       "Outside the Hansen sphere (RED %.2f > 1) -- result is extrapolated.",
                       p.relativeEnergyDifference);
  }
  if (!p.converged) {
    ImGui::TextColored(style::col::Danger,
                       "Saturation composition did not converge -- treat this result as "
                       "unreliable.");
  }
}

void drawSolventSlot(SolubilityState& sb, const std::vector<sol::Solvent>& all, int index,
                     std::string& filter) {
  static const char* kSlotLabel[3] = {"Solvent A", "Solvent B", "Solvent C"};
  std::string& id = sb.solventIds[static_cast<size_t>(index)];
  const sol::Solvent* current = sol::findSolvent(id);

  ImGui::TextUnformatted(kSlotLabel[index]);
  ImGui::SameLine(110.0f * uiScale());
  ImGui::SetNextItemWidth(220.0f * uiScale());
  if (ImGui::BeginCombo("##solvent_pick", current ? current->name.c_str() : "Select solvent...")) {
    ImGui::SetNextItemWidth(-1.0f);
    stringInputWithHint("##filter", "Search...", filter, 0);
    ImGui::Separator();
    for (const sol::Solvent& sv : all) {
      if (!containsCaseInsensitive(sv.name, filter)) continue;
      const bool selected = sv.id == id;
      if (ImGui::Selectable(sv.name.c_str(), selected)) id = sv.id;
      if (selected) ImGui::SetItemDefaultFocus();
    }
    ImGui::EndCombo();
  }
  ImGui::SameLine();
  ImGui::SetNextItemWidth(150.0f * uiScale());
  ImGui::SliderFloat("##ratio", &sb.ratios[static_cast<size_t>(index)], 0.0f, 10.0f,
                     "%.2f parts");
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", "Relative volume parts of this solvent in the blend");
}

// ------------------------------------------------------- binary ratio plot

// The ratio graph is the headline feature: give it real screen room at any
// CHEMCAD_UI_SCALE, and reserve the same height for the ternary plot below
// so switching solvent count doesn't jump the layout around.
constexpr float kPlotMinHeight = 260.0f;

void drawBinarySweepPlot(SolubilityState& sb, const sol::Solvent& a, const sol::Solvent& b) {
  if (sb.sweep.empty()) {
    ImGui::TextDisabled("Not enough data for a ratio plot.");
    return;
  }

  std::vector<const sol::SweepPoint*> sorted;
  sorted.reserve(sb.sweep.size());
  for (const sol::SweepPoint& sp : sb.sweep) sorted.push_back(&sp);
  std::sort(sorted.begin(), sorted.end(), [](const sol::SweepPoint* x, const sol::SweepPoint* y) {
    return x->fractions[0] < y->fractions[0];
  });

  double lo = sorted.front()->prediction.gramsPerMillilitre;
  double hi = lo;
  for (const sol::SweepPoint* sp : sorted) {
    lo = std::min(lo, sp->prediction.gramsPerMillilitre);
    hi = std::max(hi, sp->prediction.gramsPerMillilitre);
  }
  // A flat curve (or an all-zero one) would divide by zero below; fall back
  // to a small sensible span around the value so the frame still renders.
  if (hi - lo < 1e-12) hi = lo + std::max(1e-9, std::fabs(lo) * 0.1 + 1e-9);

  const sol::SweepPoint* peak =
      (sb.sweepPeakIndex >= 0 && static_cast<size_t>(sb.sweepPeakIndex) < sb.sweep.size())
          ? &sb.sweep[static_cast<size_t>(sb.sweepPeakIndex)]
          : nullptr;

  const ImVec2 canvasSize(ImGui::GetContentRegionAvail().x, kPlotMinHeight * uiScale());
  const ImVec2 origin = ImGui::GetCursorScreenPos();
  ImGui::InvisibleButton("##binary_plot", canvasSize);
  const bool hovered = ImGui::IsItemHovered();
  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddRectFilled(origin, ImVec2(origin.x + canvasSize.x, origin.y + canvasSize.y),
                    style::u32(style::col::BgSurface), style::metrics().radiusMd);

  const float padL = 62.0f * uiScale();
  const float padR = 20.0f * uiScale();
  const float padT = 16.0f * uiScale();
  const float padB = 46.0f * uiScale();
  const ImVec2 plotMin(origin.x + padL, origin.y + padT);
  const ImVec2 plotMax(origin.x + canvasSize.x - padR, origin.y + canvasSize.y - padB);
  dl->AddRectFilled(plotMin, plotMax, style::u32(style::col::BgPanel));
  dl->AddRect(plotMin, plotMax, style::u32(style::col::BorderStrong));

  const auto toScreen = [&](double fracA, double value) {
    const float x = plotMin.x + static_cast<float>(fracA) * (plotMax.x - plotMin.x);
    const float t = static_cast<float>((value - lo) / (hi - lo));
    const float y = plotMax.y - t * (plotMax.y - plotMin.y);
    return ImVec2(x, y);
  };

  // x-axis ticks: volume % of solvent A, 0/25/50/75/100.
  for (int i = 0; i <= 4; ++i) {
    const float t = static_cast<float>(i) / 4.0f;
    const float x = plotMin.x + t * (plotMax.x - plotMin.x);
    dl->AddLine(ImVec2(x, plotMin.y), ImVec2(x, plotMax.y), style::u32(style::col::Border, 0.35f));
    char buf[8];
    std::snprintf(buf, sizeof(buf), "%.0f%%", t * 100.0f);
    const ImVec2 tickSize = ImGui::CalcTextSize(buf);
    dl->AddText(ImVec2(x - tickSize.x * 0.5f, plotMax.y + 4.0f * uiScale()),
               style::u32(style::col::TextDim), buf);
  }
  // y-axis ticks: g/mL across [lo, hi], scientific notation below 1e-3 (via
  // formatSolubility).
  for (int i = 0; i <= 4; ++i) {
    const float t = static_cast<float>(i) / 4.0f;
    const double value = lo + (hi - lo) * (1.0 - t);
    const float y = plotMin.y + t * (plotMax.y - plotMin.y);
    dl->AddLine(ImVec2(plotMin.x, y), ImVec2(plotMax.x, y), style::u32(style::col::Border, 0.2f));
    const std::string label = formatSolubility(value);
    const ImVec2 labelSize = ImGui::CalcTextSize(label.c_str());
    dl->AddText(ImVec2(plotMin.x - labelSize.x - 6.0f * uiScale(), y - labelSize.y * 0.5f),
               style::u32(style::col::TextDim), label.c_str());
  }

  // Axis titles: the pure endpoints of the blend (fracA=0 is 100% b, the
  // opposite of the x-axis's "% of solvent a" ticks).
  const std::string leftTitle = "100% " + b.name;
  const std::string rightTitle = "100% " + a.name;
  const ImVec2 rightTitleSize = ImGui::CalcTextSize(rightTitle.c_str());
  dl->AddText(ImVec2(plotMin.x, plotMax.y + 18.0f * uiScale()), style::u32(style::col::TextFaint),
             leftTitle.c_str());
  dl->AddText(ImVec2(plotMax.x - rightTitleSize.x, plotMax.y + 18.0f * uiScale()),
             style::u32(style::col::TextFaint), rightTitle.c_str());

  // Shaded area under the curve -- one quad per segment rather than a single
  // convex fill, so the rise-then-fall co-solvency shape fills correctly
  // even though the region under a peaked curve is not convex.
  for (size_t i = 0; i + 1 < sorted.size(); ++i) {
    const ImVec2 p0 = toScreen(sorted[i]->fractions[0], sorted[i]->prediction.gramsPerMillilitre);
    const ImVec2 p1 =
        toScreen(sorted[i + 1]->fractions[0], sorted[i + 1]->prediction.gramsPerMillilitre);
    const ImVec2 b0(p0.x, plotMax.y);
    const ImVec2 b1(p1.x, plotMax.y);
    dl->AddQuadFilled(p0, p1, b1, b0, style::u32(style::col::Teal, 0.16f));
  }

  std::vector<ImVec2> points;
  points.reserve(sorted.size());
  for (const sol::SweepPoint* sp : sorted)
    points.push_back(toScreen(sp->fractions[0], sp->prediction.gramsPerMillilitre));
  if (points.size() >= 2)
    dl->AddPolyline(points.data(), static_cast<int>(points.size()), style::u32(style::col::Teal),
                    0, 2.4f);
  else if (points.size() == 1)
    dl->AddCircleFilled(points.front(), 3.5f, style::u32(style::col::Teal));

  // Mark the co-solvency maximum: drop-line, filled marker, and a
  // composition + value label.
  if (peak) {
    const ImVec2 peakPt = toScreen(peak->fractions[0], peak->prediction.gramsPerMillilitre);
    dl->AddLine(ImVec2(peakPt.x, plotMax.y), peakPt, style::u32(style::col::Accent, 0.7f), 1.5f);
    dl->AddCircleFilled(peakPt, 5.5f, style::u32(style::col::Accent));
    dl->AddCircle(peakPt, 8.5f, style::u32(style::col::Accent), 0, 1.5f);

    char peakLabel[96];
    std::snprintf(peakLabel, sizeof(peakLabel), "peak %s g/mL at %.0f%% %s",
                 formatSolubility(peak->prediction.gramsPerMillilitre).c_str(),
                 peak->fractions[0] * 100.0, a.name.c_str());
    const ImVec2 peakLabelSize = ImGui::CalcTextSize(peakLabel);
    const float lxRaw = peakPt.x - peakLabelSize.x * 0.5f;
    float lx = std::clamp(lxRaw, plotMin.x, plotMax.x - peakLabelSize.x);
    float ly = peakPt.y - peakLabelSize.y - 10.0f * uiScale();
    if (ly < plotMin.y) ly = peakPt.y + 10.0f * uiScale();
    dl->AddRectFilled(ImVec2(lx - 4.0f, ly - 2.0f),
                      ImVec2(lx + peakLabelSize.x + 4.0f, ly + peakLabelSize.y + 2.0f),
                      style::u32(style::col::BgRaised, 0.9f), 3.0f);
    dl->AddText(ImVec2(lx, ly), style::u32(style::col::Accent), peakLabel);
  }

  // Mark the current working ratio -- a distinct colour from the peak so
  // the two are never confused at a glance.
  const double totalRatio = static_cast<double>(sb.ratios[0]) + static_cast<double>(sb.ratios[1]);
  if (totalRatio > 1e-9) {
    const double fracA = static_cast<double>(sb.ratios[0]) / totalRatio;
    const ImVec2 marker = toScreen(fracA, sb.prediction.gramsPerMillilitre);
    dl->AddLine(ImVec2(marker.x, plotMin.y), ImVec2(marker.x, plotMax.y),
               style::u32(style::col::Violet, 0.55f), 1.5f);
    dl->AddCircleFilled(marker, 4.5f, style::u32(style::col::Violet));
    dl->AddText(ImVec2(marker.x + 6.0f, plotMin.y + 2.0f), style::u32(style::col::Violet),
               "current");
  }

  if (hovered) {
    const ImVec2 mouse = ImGui::GetMousePos();
    float t = (mouse.x - plotMin.x) / std::max(1.0f, plotMax.x - plotMin.x);
    t = std::clamp(t, 0.0f, 1.0f);
    size_t nearest = 0;
    double best = 1e18;
    for (size_t i = 0; i < sorted.size(); ++i) {
      const double d = std::fabs(sorted[i]->fractions[0] - t);
      if (d < best) {
        best = d;
        nearest = i;
      }
    }
    const sol::SweepPoint& sp = *sorted[nearest];
    dl->AddLine(ImVec2(mouse.x, plotMin.y), ImVec2(mouse.x, plotMax.y),
               style::u32(style::col::TextFaint, 0.5f));
    char tip[128];
    std::snprintf(tip, sizeof(tip), "%s %.0f%%  /  %s %.0f%%   ->   %s g/mL", a.name.c_str(),
                  sp.fractions[0] * 100.0, b.name.c_str(), sp.fractions[1] * 100.0,
                  formatSolubility(sp.prediction.gramsPerMillilitre).c_str());
    ImGui::SetTooltip("%s", tip);
  }

  if (peak) {
    ImGui::Text("Co-solvency peak: %s g/mL at %.0f%% %s (%.0f%% %s).",
               formatSolubility(peak->prediction.gramsPerMillilitre).c_str(),
               peak->fractions[0] * 100.0, a.name.c_str(), peak->fractions[1] * 100.0,
               b.name.c_str());
    ImGui::SameLine();
    if (widgets::ghostButton("Apply peak ratio")) {
      sb.ratios[0] = static_cast<float>(peak->fractions[0]);
      sb.ratios[1] = static_cast<float>(peak->fractions[1]);
    }
  }
}

// ------------------------------------------------------------ ternary plot

struct TernaryGrid {
  int steps = 0;
  std::unordered_map<int64_t, const sol::SweepPoint*> byIndex;
};

TernaryGrid buildTernaryGrid(const std::vector<sol::SweepPoint>& sweep, int steps) {
  TernaryGrid grid;
  grid.steps = steps;
  grid.byIndex.reserve(sweep.size() * 2);
  for (const sol::SweepPoint& sp : sweep) {
    const int i = static_cast<int>(std::lround(sp.fractions[0] * steps));
    const int j = static_cast<int>(std::lround(sp.fractions[1] * steps));
    grid.byIndex[static_cast<int64_t>(i) * 4096 + j] = &sp;
  }
  return grid;
}

const sol::SweepPoint* gridAt(const TernaryGrid& grid, int i, int j) {
  if (i < 0 || j < 0 || i + j > grid.steps) return nullptr;
  const auto it = grid.byIndex.find(static_cast<int64_t>(i) * 4096 + j);
  return it == grid.byIndex.end() ? nullptr : it->second;
}

void drawTernarySweepPlot(SolubilityState& sb, const sol::Solvent& a, const sol::Solvent& b,
                          const sol::Solvent& c) {
  if (sb.sweep.size() < 3) {
    ImGui::TextDisabled("Not enough data for a ternary plot.");
    return;
  }

  double lo = sb.sweep.front().prediction.gramsPerMillilitre;
  double hi = lo;
  for (const sol::SweepPoint& sp : sb.sweep) {
    lo = std::min(lo, sp.prediction.gramsPerMillilitre);
    hi = std::max(hi, sp.prediction.gramsPerMillilitre);
  }
  const double range = hi - lo;
  const bool degenerate = range < 1e-12;

  const float side =
      std::clamp(ImGui::GetContentRegionAvail().x - 170.0f * uiScale(), 140.0f * uiScale(),
                340.0f * uiScale());
  const float height = side * 0.8660254f;  // sqrt(3)/2
  // Reserve at least as much vertical room as the binary plot so switching
  // between 2- and 3-solvent blends doesn't jump the layout around.
  const float canvasH = std::max(kPlotMinHeight * uiScale(), height + 60.0f * uiScale());
  const ImVec2 canvasSize(side + 190.0f * uiScale(), canvasH);
  const ImVec2 origin = ImGui::GetCursorScreenPos();
  ImGui::InvisibleButton("##ternary_plot", canvasSize);
  const bool hovered = ImGui::IsItemHovered();
  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddRectFilled(origin, ImVec2(origin.x + canvasSize.x, origin.y + canvasSize.y),
                    style::u32(style::col::BgSurface), style::metrics().radiusMd);

  const float cx = origin.x + side * 0.5f + 24.0f * uiScale();
  const float top = origin.y + std::max(24.0f * uiScale(), (canvasH - height) * 0.5f);
  const ImVec2 vA(cx, top);                            // 100% solvent A (top)
  const ImVec2 vB(cx - side * 0.5f, top + height);      // 100% solvent B (bottom-left)
  const ImVec2 vC(cx + side * 0.5f, top + height);      // 100% solvent C (bottom-right)

  const auto toScreen = [&](double wa, double wb, double wc) {
    return ImVec2(static_cast<float>(wa) * vA.x + static_cast<float>(wb) * vB.x +
                      static_cast<float>(wc) * vC.x,
                  static_cast<float>(wa) * vA.y + static_cast<float>(wb) * vB.y +
                      static_cast<float>(wc) * vC.y);
  };
  const auto valueColor = [&](double v) {
    const float t = degenerate ? 0.5f : static_cast<float>((v - lo) / range);
    return colorRamp(t);
  };

  const int steps = sb.sweepSteps;
  const TernaryGrid grid = buildTernaryGrid(sb.sweep, steps);
  for (int i = 0; i < steps; ++i) {
    for (int j = 0; i + j < steps; ++j) {
      const sol::SweepPoint* p00 = gridAt(grid, i, j);
      const sol::SweepPoint* p10 = gridAt(grid, i + 1, j);
      const sol::SweepPoint* p01 = gridAt(grid, i, j + 1);
      if (p00 && p10 && p01) {
        const double avg = (p00->prediction.gramsPerMillilitre +
                            p10->prediction.gramsPerMillilitre +
                            p01->prediction.gramsPerMillilitre) /
                           3.0;
        dl->AddTriangleFilled(toScreen(p00->fractions[0], p00->fractions[1], p00->fractions[2]),
                              toScreen(p10->fractions[0], p10->fractions[1], p10->fractions[2]),
                              toScreen(p01->fractions[0], p01->fractions[1], p01->fractions[2]),
                              valueColor(avg));
      }
      const sol::SweepPoint* p11 = gridAt(grid, i + 1, j + 1);
      if (p10 && p01 && p11) {
        const double avg = (p10->prediction.gramsPerMillilitre +
                            p01->prediction.gramsPerMillilitre +
                            p11->prediction.gramsPerMillilitre) /
                           3.0;
        dl->AddTriangleFilled(toScreen(p10->fractions[0], p10->fractions[1], p10->fractions[2]),
                              toScreen(p01->fractions[0], p01->fractions[1], p01->fractions[2]),
                              toScreen(p11->fractions[0], p11->fractions[1], p11->fractions[2]),
                              valueColor(avg));
      }
    }
  }

  dl->AddTriangle(vA, vB, vC, style::u32(style::col::BorderStrong), 1.6f);

  const ImVec2 aLabel = ImGui::CalcTextSize(a.name.c_str());
  dl->AddText(ImVec2(vA.x - aLabel.x * 0.5f, vA.y - 18.0f * uiScale()),
             style::u32(style::col::Text), a.name.c_str());
  const ImVec2 bLabel = ImGui::CalcTextSize(b.name.c_str());
  dl->AddText(ImVec2(vB.x - bLabel.x - 4.0f, vB.y + 4.0f), style::u32(style::col::Text),
             b.name.c_str());
  dl->AddText(ImVec2(vC.x + 4.0f, vC.y + 4.0f), style::u32(style::col::Text), c.name.c_str());

  const double totalRatio = static_cast<double>(sb.ratios[0]) + static_cast<double>(sb.ratios[1]) +
                            static_cast<double>(sb.ratios[2]);
  if (totalRatio > 1e-9) {
    const ImVec2 marker = toScreen(sb.ratios[0] / totalRatio, sb.ratios[1] / totalRatio,
                                   sb.ratios[2] / totalRatio);
    dl->AddCircleFilled(marker, 5.0f, style::u32(style::col::Accent));
    dl->AddCircle(marker, 8.0f, style::u32(style::col::Accent), 0, 1.5f);
  }

  // Legend: vertical gradient bar with min/max labels.
  const float legendX = origin.x + side + 76.0f * uiScale();
  const float legendW = 14.0f * uiScale();
  constexpr int kStripes = 24;
  for (int s = 0; s < kStripes; ++s) {
    const float t0 = static_cast<float>(s) / kStripes;
    const float t1 = static_cast<float>(s + 1) / kStripes;
    const float y0 = top + height * (1.0f - t1);
    const float y1 = top + height * (1.0f - t0);
    dl->AddRectFilled(ImVec2(legendX, y0), ImVec2(legendX + legendW, y1),
                      colorRamp((t0 + t1) * 0.5f));
  }
  dl->AddRect(ImVec2(legendX, top), ImVec2(legendX + legendW, top + height),
             style::u32(style::col::Border));
  dl->AddText(ImVec2(legendX + legendW + 6.0f, top - 2.0f), style::u32(style::col::TextDim),
             formatSolubility(hi).c_str());
  dl->AddText(ImVec2(legendX + legendW + 6.0f, top + height - 12.0f),
             style::u32(style::col::TextDim), formatSolubility(lo).c_str());
  dl->AddText(ImVec2(legendX, top + height + 6.0f), style::u32(style::col::TextFaint), "g/mL");

  if (hovered) {
    const ImVec2 mouse = ImGui::GetMousePos();
    const float denom = (vB.y - vC.y) * (vA.x - vC.x) + (vC.x - vB.x) * (vA.y - vC.y);
    if (std::fabs(denom) > 1e-6f) {
      const float wa =
          ((vB.y - vC.y) * (mouse.x - vC.x) + (vC.x - vB.x) * (mouse.y - vC.y)) / denom;
      const float wb =
          ((vC.y - vA.y) * (mouse.x - vC.x) + (vA.x - vC.x) * (mouse.y - vC.y)) / denom;
      const float wc = 1.0f - wa - wb;
      if (wa >= -0.02f && wb >= -0.02f && wc >= -0.02f) {
        const sol::SweepPoint* nearest = nullptr;
        double best = 1e18;
        for (const sol::SweepPoint& sp : sb.sweep) {
          const double d = std::fabs(sp.fractions[0] - wa) + std::fabs(sp.fractions[1] - wb) +
                           std::fabs(sp.fractions[2] - wc);
          if (d < best) {
            best = d;
            nearest = &sp;
          }
        }
        if (nearest) {
          char tip[192];
          std::snprintf(tip, sizeof(tip), "%s %.0f%%  /  %s %.0f%%  /  %s %.0f%%   ->   %s g/mL",
                        a.name.c_str(), nearest->fractions[0] * 100.0, b.name.c_str(),
                        nearest->fractions[1] * 100.0, c.name.c_str(),
                        nearest->fractions[2] * 100.0,
                        formatSolubility(nearest->prediction.gramsPerMillilitre).c_str());
          ImGui::SetTooltip("%s", tip);
        }
      }
    }
  }
}

}  // namespace

void drawSolubilitySuite(AppState& st) {
  SolubilityState& sb = st.solubility;
  static bool prevUseSketch = sb.useSketch;
  static bool prevOverride = sb.overrideSolute;
  static std::array<std::string, 3> filterBuf;

  // ------------------------------------------------------------- solute ---
  widgets::sectionHeader("Solute");
  if (ImGui::RadioButton("From sketch", sb.useSketch)) sb.useSketch = true;
  ImGui::SameLine();
  if (ImGui::RadioButton("From SMILES", !sb.useSketch)) sb.useSketch = false;

  const bool modeToggled = sb.useSketch != prevUseSketch;
  prevUseSketch = sb.useSketch;

  if (sb.useSketch) {
    if (modeToggled || st.docRevision != sb.sourceRevision) recomputeSoluteFromSketch(st);
  } else {
    if (modeToggled) {
      sb.soluteValid = false;
      sb.soluteError.clear();
    }
    ImGui::SetNextItemWidth(-1.0f);
    const bool enter = stringInputWithHint("##solute_smiles", "SMILES, e.g. CC(=O)Oc1ccccc1C(=O)O",
                                           sb.manualSmiles, ImGuiInputTextFlags_EnterReturnsTrue,
                                           true);
    if (enter || ImGui::IsItemDeactivatedAfterEdit()) recomputeSoluteFromSmiles(st);
  }

  if (!sb.soluteError.empty()) {
    ImGui::TextColored(style::col::Danger, "%s", sb.soluteError.c_str());
  } else if (!sb.soluteValid) {
    ImGui::TextDisabled("%s", sb.useSketch ? "Draw a structure to describe the solute."
                                           : "Enter a SMILES string and press Enter.");
  } else {
    const float avail = ImGui::GetContentRegionAvail().x;
    const float spacing = ImGui::GetStyle().ItemSpacing.x;
    const float cardW = (avail - spacing * 4.0f) / 5.0f;
    const float cardH = statCardHeight();
    char mw[32];
    std::snprintf(mw, sizeof(mw), "%.1f", sb.solute.molarMass);
    char mv[32];
    std::snprintf(mv, sizeof(mv), "%.1f", sb.solute.molarVolume);
    char logp[32];
    std::snprintf(logp, sizeof(logp), "%.2f", sb.solute.logP);
    char hansen[64];
    std::snprintf(hansen, sizeof(hansen), "%.1f / %.1f / %.1f", sb.solute.hansen.dispersion,
                 sb.solute.hansen.polar, sb.solute.hansen.hydrogenBond);
    char meltPt[32];
    if (sb.solute.meltingPointEstimated) {
      std::snprintf(meltPt, sizeof(meltPt), "%.0f (est)", sb.solute.meltingPoint);
    } else {
      std::snprintf(meltPt, sizeof(meltPt), "%.0f", sb.solute.meltingPoint);
    }
    widgets::statCard("MW g/mol", mw, ImVec2(cardW, cardH));
    ImGui::SameLine(0.0f, spacing);
    widgets::statCard("MOLAR VOL cm3/mol", mv, ImVec2(cardW, cardH));
    ImGui::SameLine(0.0f, spacing);
    widgets::statCard("LOGP", logp, ImVec2(cardW, cardH));
    ImGui::SameLine(0.0f, spacing);
    widgets::statCard("dD / dP / dH", hansen, ImVec2(cardW, cardH));
    ImGui::SameLine(0.0f, spacing);
    widgets::statCard("MELTING PT C", meltPt, ImVec2(cardW, cardH));
    if (sb.solute.meltingPointEstimated) {
      ImGui::SameLine();
      widgets::badge("Tm ESTIMATED", style::col::Violet);
    }
  }

  if (sb.soluteValid && sb.solute.meltingPointEstimated && !sb.overrideSolute) {
    // An estimated Tm is the single biggest source of error in the
    // prediction (a Joback group-contribution guess vs. a literature value
    // can be tens of degrees off) -- make the fix impossible to miss.
    ImGui::TextColored(style::col::Violet,
                       "Melting point is a Joback group-contribution estimate (%.0f C), not a "
                       "measured value -- override it below for an accurate prediction.",
                       sb.solute.meltingPoint);
  }

  const bool overrideClicked = ImGui::Checkbox("Override melting point / Hansen radius R0",
                                               &sb.overrideSolute);
  if (overrideClicked && sb.overrideSolute && !prevOverride && sb.soluteValid) {
    sb.meltingPointC = static_cast<float>(sb.solute.meltingPoint);
    sb.interactionRadius = static_cast<float>(sb.solute.interactionRadius);
  }
  prevOverride = sb.overrideSolute;
  if (sb.overrideSolute) {
    bool changed = false;
    ImGui::SetNextItemWidth(220.0f * uiScale());
    changed |= ImGui::SliderFloat("Melting point (C)", &sb.meltingPointC, -50.0f, 300.0f, "%.1f");
    ImGui::SetNextItemWidth(220.0f * uiScale());
    changed |= ImGui::SliderFloat("Hansen radius R0", &sb.interactionRadius, 1.0f, 30.0f, "%.1f");
    if (sb.soluteValid) {
      applySoluteOverrides(sb);
      if (changed) ++sb.soluteVersion;
    }
  }

  ImGui::Spacing();

  // ------------------------------------------------------------ solvents --
  widgets::sectionHeader("Solvents");
  static const char* kCountItems[] = {"1  (pure solvent)", "2  (binary blend)",
                                      "3  (ternary blend)"};
  int countIdx = sb.solventCount - 1;
  ImGui::SetNextItemWidth(220.0f * uiScale());
  if (ImGui::Combo("Number of solvents", &countIdx, kCountItems, 3)) sb.solventCount = countIdx + 1;

  bool solventsOk = false;
  std::vector<const sol::Solvent*> chosenAll(static_cast<size_t>(sb.solventCount), nullptr);
  std::vector<sol::Component> components;
  try {
    const std::vector<sol::Solvent>& all = sol::solvents();
    for (int i = 0; i < sb.solventCount; ++i) {
      ImGui::PushID(i);
      drawSolventSlot(sb, all, i, filterBuf[static_cast<size_t>(i)]);
      ImGui::PopID();
      chosenAll[static_cast<size_t>(i)] = sol::findSolvent(sb.solventIds[static_cast<size_t>(i)]);
    }
    solventsOk = true;
    for (int i = 0; i < sb.solventCount; ++i) {
      const sol::Solvent* found = chosenAll[static_cast<size_t>(i)];
      if (found && sb.ratios[static_cast<size_t>(i)] > 0.0f)
        components.push_back({found, static_cast<double>(sb.ratios[static_cast<size_t>(i)])});
    }
  } catch (const sol::SolError& err) {
    ImGui::TextColored(style::col::Danger, "%s", err.what());
    sb.statusMessage = std::string("Solvent database error: ") + err.what();
  }

  sol::Mixture mixture;
  if (solventsOk && !components.empty()) {
    mixture = sol::blend(components);
    ImGui::TextColored(style::col::TextDim,
                       "Blend Hansen: dD %.1f  dP %.1f  dH %.1f MPa^0.5   density %.3f g/mL",
                       mixture.hansen.dispersion, mixture.hansen.polar, mixture.hansen.hydrogenBond,
                       mixture.density);
  } else if (solventsOk) {
    ImGui::TextDisabled("Select at least one solvent.");
  }

  ImGui::Spacing();
  ImGui::SetNextItemWidth(260.0f * uiScale());
  ImGui::SliderFloat("Temperature (C)", &sb.temperatureC, -20.0f, 150.0f, "%.1f");

  ImGui::Spacing();

  // ------------------------------------------------------------- result ---
  widgets::sectionHeader("Predicted solubility");
  const bool canPredict = sb.soluteValid && solventsOk && !components.empty();
  if (canPredict) {
    sb.prediction = sol::predict(sb.solute, components, static_cast<double>(sb.temperatureC));
    drawResultCard(sb);
    sb.statusMessage = "Predicted " + formatSolubility(sb.prediction.gramsPerMillilitre) + " g/mL.";
  } else {
    sb.prediction = sol::Prediction{};
    ImGui::TextDisabled("Provide a valid solute and at least one solvent to see a prediction.");
  }

  // -------------------------------------------------------- ratio graph ---
  const bool haveAllSlots =
      solventsOk &&
      std::all_of(chosenAll.begin(), chosenAll.end(),
                 [](const sol::Solvent* s) { return s != nullptr; });
  if (sb.solventCount >= 2) {
    widgets::sectionHeader("Solubility vs. ratio");
    ImGui::SetNextItemWidth(220.0f * uiScale());
    ImGui::SliderInt("Grid resolution", &sb.sweepSteps, 2, 64);

    if (sb.soluteValid && haveAllSlots) {
      recomputeSweepIfNeeded(sb, chosenAll);
      if (sb.solventCount == 2) {
        drawBinarySweepPlot(sb, *chosenAll[0], *chosenAll[1]);
      } else {
        drawTernarySweepPlot(sb, *chosenAll[0], *chosenAll[1], *chosenAll[2]);
      }
    } else {
      sb.sweep.clear();
      sb.sweepSignature.clear();
      ImGui::TextDisabled("Choose a solute and every solvent slot to see the ratio plot.");
    }
    ImGui::Spacing();
  } else {
    sb.sweep.clear();
    sb.sweepSignature.clear();
  }

  // ------------------------------------------------------------- funnel ---
  widgets::sectionHeader("Separatory funnel");
  drawFunnelView(st);
}

}  // namespace chemcad::ui

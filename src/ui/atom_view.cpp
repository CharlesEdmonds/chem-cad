// Animated atomic-structure renderer for the periodic-table hover card.
//
// Two representations of the same electron configuration:
//   * Orbital clouds (primary): a rotating pseudo-3D point cloud where every
//     occupied subshell (s/p/d/f) contributes dots sampled from the shape of
//     its orbitals, clustered around a proton+neutron nucleus.
//   * Bohr shell model (backup): tilted shell rings with electrons orbiting.
//
// Everything is drawn through ImDrawList -- no GL, headless-safe. The point
// cloud is generated once per element (seeded by Z) and cached; only the
// rotation is animated, so frame cost is one sort + dots while hovering.

#include "ui/atom_view.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "ui/theme.hpp"

namespace chemcad::ui {
namespace {

constexpr double kPi = 3.14159265358979323846;

// ------------------------------------------------------------ configuration
struct Subshell {
  int n = 1;
  int l = 0;  // 0 s, 1 p, 2 d, 3 f
  int electrons = 0;
};

// Standard Aufbau filling order.
constexpr Subshell kAufbau[] = {
    {1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {3, 0, 0}, {3, 1, 0}, {4, 0, 0}, {3, 2, 0},
    {4, 1, 0}, {5, 0, 0}, {4, 2, 0}, {5, 1, 0}, {6, 0, 0}, {4, 3, 0}, {5, 2, 0},
    {6, 1, 0}, {7, 0, 0}, {5, 3, 0}, {6, 2, 0}, {7, 1, 0},
};

int subshellCapacity(int l) { return 2 * (2 * l + 1); }

// Known anomalous ground-state configurations (s->d / s->f promotions).
// {z, fromN, fromL, toN, toL, electronsMoved}
struct ConfigPatch {
  uint8_t z;
  int fromN, fromL, toN, toL, moved;
};
constexpr ConfigPatch kExceptions[] = {
    {24, 4, 0, 3, 2, 1},  // Cr  4s1 3d5
    {29, 4, 0, 3, 2, 1},  // Cu  4s1 3d10
    {41, 5, 0, 4, 2, 1},  // Nb  5s1 4d4
    {42, 5, 0, 4, 2, 1},  // Mo  5s1 4d5
    {44, 5, 0, 4, 2, 1},  // Ru  5s1 4d7
    {45, 5, 0, 4, 2, 1},  // Rh  5s1 4d8
    {46, 5, 0, 4, 2, 2},  // Pd  4d10
    {47, 5, 0, 4, 2, 1},  // Ag  5s1 4d10
    {64, 6, 0, 5, 2, 1},  // Gd  4f7 5d1 6s2
    {78, 6, 0, 5, 2, 1},  // Pt  6s1 5d9
    {79, 6, 0, 5, 2, 1},  // Au  6s1 5d10
    {89, 7, 0, 6, 2, 2},  // Th  6d2 7s2
    {96, 7, 0, 6, 2, 1},  // Cm  5f7 6d1 7s2
    {103, 7, 0, 6, 2, 1}, // Lr  5f14 6d1 7s2 (vs 7p1)
};

std::vector<Subshell> electronConfig(uint8_t z) {
  std::vector<Subshell> config(std::begin(kAufbau), std::end(kAufbau));
  int remaining = z;
  for (Subshell& sub : config) {
    const int take = std::min(remaining, subshellCapacity(sub.l));
    sub.electrons = take;
    remaining -= take;
    if (remaining <= 0) break;
  }
  for (const ConfigPatch& patch : kExceptions) {
    if (patch.z != z) continue;
    for (Subshell& from : config) {
      if (from.n == patch.fromN && from.l == patch.fromL) from.electrons -= patch.moved;
    }
    bool found = false;
    for (Subshell& to : config) {
      if (to.n == patch.toN && to.l == patch.toL) {
        to.electrons += patch.moved;
        found = true;
      }
    }
    if (!found) config.push_back({patch.toN, patch.toL, patch.moved});
  }
  config.erase(std::remove_if(config.begin(), config.end(),
                              [](const Subshell& s) { return s.electrons <= 0; }),
               config.end());
  return config;
}

// Noble-gas-core shorthand, e.g. [Ar] 3d10 4s2 4p4.
std::string configString(const std::vector<Subshell>& config, uint8_t z) {
  static const char* kCoreSymbols[] = {"He", "Ne", "Ar", "Kr", "Xe", "Rn", "Og"};
  static const int kCoreZ[] = {2, 10, 18, 36, 54, 86, 118};
  static const char* kL[] = {"s", "p", "d", "f"};

  int coreIndex = -1;
  for (int i = 0; i < 7; ++i) {
    if (kCoreZ[i] < z) coreIndex = i;
  }
  const int coreZ = coreIndex >= 0 ? kCoreZ[coreIndex] : 0;

  std::string out;
  if (coreIndex >= 0) out = std::string("[") + kCoreSymbols[coreIndex] + "] ";
  int cum = 0;
  for (const Subshell& sub : config) {
    const int start = cum;
    cum += sub.electrons;
    if (cum <= coreZ) continue;  // entirely inside the noble-gas core
    const int shown = cum - std::max(start, coreZ);
    char buf[24];
    std::snprintf(buf, sizeof(buf), "%d%s%d  ", sub.n, kL[sub.l], shown);
    out += buf;
  }
  while (!out.empty() && out.back() == ' ') out.pop_back();
  return out;
}

// ------------------------------------------------------------- point cloud
struct CloudPoint {
  float x, y, z;
  uint8_t l;      // subshell type, for colour
  bool core = false;  // nucleus sphere
  bool proton = false;
};

struct AtomModel {
  std::vector<CloudPoint> points;   // orbital dots
  std::vector<CloudPoint> nucleus;  // nucleon spheres (unit space)
  int shells[7] = {0, 0, 0, 0, 0, 0, 0};  // Bohr shell totals
};

float gauss(std::mt19937& rng, float sigma) {
  std::normal_distribution<float> d(0.0f, sigma);
  return d(rng);
}

// Unit vector on the sphere, uniform direction.
void sphereDir(std::mt19937& rng, float& x, float& y, float& z) {
  std::uniform_real_distribution<float> u(-1.0f, 1.0f);
  for (;;) {
    x = u(rng);
    y = u(rng);
    z = u(rng);
    const float n = x * x + y * y + z * z;
    if (n > 1e-4f && n <= 1.0f) {
      const float inv = 1.0f / std::sqrt(n);
      x *= inv;
      y *= inv;
      z *= inv;
      return;
    }
  }
}

void sampleS(std::mt19937& rng, float rn, CloudPoint& p) {
  std::uniform_real_distribution<float> u(0.0f, 1.0f);
  float dx, dy, dz;
  sphereDir(rng, dx, dy, dz);
  const float r = rn * (0.55f + 0.45f * std::cbrt(u(rng)));
  p.x = dx * r;
  p.y = dy * r;
  p.z = dz * r;
}

void sampleLobe(std::mt19937& rng, float rn, float ax, float ay, float az, float spread,
                CloudPoint& p) {
  std::uniform_real_distribution<float> u(0.0f, 1.0f);
  const float t = rn * (0.28f + 0.72f * u(rng));
  p.x = ax * t + gauss(rng, spread * rn);
  p.y = ay * t + gauss(rng, spread * rn);
  p.z = az * t + gauss(rng, spread * rn);
}

const AtomModel& modelFor(const ElementData& element) {
  static const AtomModel* cache = nullptr;
  static uint8_t cacheZ = 0;
  if (cache && cacheZ == element.z) return *cache;

  static AtomModel model;
  model = AtomModel{};
  cacheZ = element.z;

  const std::vector<Subshell> config = electronConfig(element.z);
  std::mt19937 rng(static_cast<unsigned>(element.z) * 2654435761u + 7u);
  std::uniform_real_distribution<float> u(0.0f, 1.0f);
  std::uniform_int_distribution<int> sign(0, 1);

  // Dots per electron scale down for heavy atoms so the stage never saturates.
  const float dotsPerElectron =
      std::clamp(1150.0f / static_cast<float>(std::max(1, static_cast<int>(element.z))),
                 7.0f, 52.0f);

  for (const Subshell& sub : config) {
    model.shells[sub.n - 1] += sub.electrons;
    const float rn = 0.34f + 0.17f * static_cast<float>(sub.n);  // unit-space radius
    const int dots = static_cast<int>(dotsPerElectron * sub.electrons + 0.5f);
    for (int i = 0; i < dots; ++i) {
      CloudPoint p{};
      p.l = static_cast<uint8_t>(sub.l);
      // Which orbital within the subshell this dot belongs to (pairs fill in
      // order, mirroring how the electrons are counted).
      const int pair = sub.electrons > 0 ? (i * (subshellCapacity(sub.l) / 2)) / dots : 0;
      const float s = sign(rng) ? 1.0f : -1.0f;
      switch (sub.l) {
        case 0:
          sampleS(rng, rn, p);
          break;
        case 1: {
          const int orbital = std::min(pair, 2);
          const float axis[3] = {orbital == 0 ? s : 0.0f, orbital == 1 ? s : 0.0f,
                                 orbital == 2 ? s : 0.0f};
          sampleLobe(rng, rn, axis[0], axis[1], axis[2], 0.16f, p);
          break;
        }
        case 2: {
          const int orbital = std::min(pair, 4);
          const float k = 0.70710678f;
          if (orbital == 0) {  // dz2: two z lobes + a small xy collar
            if (u(rng) < 0.18f) {
              const float a = u(rng) * 2.0f * static_cast<float>(kPi);
              const float r = rn * 0.40f;
              p.x = std::cos(a) * r + gauss(rng, 0.05f * rn);
              p.y = std::sin(a) * r + gauss(rng, 0.05f * rn);
              p.z = gauss(rng, 0.05f * rn);
            } else {
              sampleLobe(rng, rn, 0.0f, 0.0f, s, 0.13f, p);
            }
          } else if (orbital == 1) {  // dxz
            sampleLobe(rng, rn, s * k, 0.0f, (sign(rng) ? 1.0f : -1.0f) * k, 0.14f, p);
          } else if (orbital == 2) {  // dyz
            sampleLobe(rng, rn, 0.0f, s * k, (sign(rng) ? 1.0f : -1.0f) * k, 0.14f, p);
          } else if (orbital == 3) {  // dxy
            sampleLobe(rng, rn, s * k, (sign(rng) ? 1.0f : -1.0f) * k, 0.0f, 0.14f, p);
          } else {  // dx2-y2
            if (sign(rng))
              sampleLobe(rng, rn, s, 0.0f, 0.0f, 0.14f, p);
            else
              sampleLobe(rng, rn, 0.0f, s, 0.0f, 0.14f, p);
          }
          break;
        }
        default: {  // f: stylised eight-lobe cubic set
          const float k = 0.57735027f;
          sampleLobe(rng, rn, s * k, (sign(rng) ? 1.0f : -1.0f) * k,
                     (sign(rng) ? 1.0f : -1.0f) * k, 0.13f, p);
          break;
        }
      }
      model.points.push_back(p);
    }
  }

  // Nucleus: protons + neutrons on a Fibonacci ball. Capped for heavy atoms;
  // the exact counts are printed next to the stage either way. A real nucleus
  // is ~10^-4 of the atom's radius; we keep it small enough to read as a
  // compact core inside the cloud, not a rival to it.
  const int neutrons = std::max(0, static_cast<int>(std::lround(element.mass)) - element.z);
  const int nucleons = element.z + neutrons;
  const int shown = std::min(nucleons, element.z <= 20 ? nucleons : 64);
  const float nucR = 0.045f + 0.012f * std::cbrt(static_cast<float>(std::max(shown, 1)));
  const float golden = 2.39996323f;
  for (int i = 0; i < shown; ++i) {
    const float y = 1.0f - 2.0f * (i + 0.5f) / shown;
    const float rr = std::sqrt(std::max(0.0f, 1.0f - y * y));
    const float th = golden * i;
    CloudPoint p{};
    p.core = true;
    // Interleave protons/neutrons with a golden-ratio shuffle for a mixed core.
    const int idx = (i * 89) % shown;
    p.proton = idx < element.z * shown / std::max(nucleons, 1);
    const float jitter = 0.72f + 0.28f * ((i * 37 % 13) / 12.0f);
    p.x = std::cos(th) * rr * nucR * jitter;
    p.y = y * nucR * jitter;
    p.z = std::sin(th) * rr * nucR * jitter;
    model.nucleus.push_back(p);
  }

  cache = &model;
  return *cache;
}

ImVec4 cloudColor(int l) {
  switch (l) {
    case 0: return {0.45f, 0.83f, 0.78f, 1.0f};  // s: teal
    case 1: return {0.96f, 0.72f, 0.32f, 1.0f};  // p: amber
    case 2: return {0.74f, 0.57f, 0.96f, 1.0f};  // d: violet
    default: return {0.58f, 0.87f, 0.55f, 1.0f}; // f: green
  }
}

struct Projected {
  float x, y, depth, persp;
  const CloudPoint* src;
};

void rotateProject(const std::vector<CloudPoint>& in, std::vector<Projected>& out, float yaw,
                   float pitch, ImVec2 centre, float scale) {
  const float cy = std::cos(yaw), sy = std::sin(yaw);
  const float cp = std::cos(pitch), sp = std::sin(pitch);
  out.resize(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    const CloudPoint& p = in[i];
    const float x1 = p.x * cy + p.z * sy;
    const float z1 = -p.x * sy + p.z * cy;
    const float y2 = p.y * cp - z1 * sp;
    const float z2 = p.y * sp + z1 * cp;
    const float persp = 1.0f / (1.0f + z2 * 0.22f);
    out[i] = Projected{centre.x + x1 * scale * persp, centre.y - y2 * scale * persp, z2,
                       persp, &p};
  }
}

void shadedSphere(ImDrawList* dl, ImVec2 c, float r, ImVec4 base, float fog) {
  const ImVec4 bg(0.043f, 0.055f, 0.075f, 1.0f);
  const ImVec4 col(base.x + (bg.x - base.x) * fog, base.y + (bg.y - base.y) * fog,
                   base.z + (bg.z - base.z) * fog, 1.0f);
  dl->AddCircleFilled(c, r, ImGui::ColorConvertFloat4ToU32(col), 12);
  const ImVec4 hi(col.x + (1.0f - col.x) * 0.5f, col.y + (1.0f - col.y) * 0.5f,
                  col.z + (1.0f - col.z) * 0.5f, 0.8f);
  dl->AddCircleFilled(ImVec2(c.x - r * 0.32f, c.y - r * 0.34f), r * 0.34f,
                      ImGui::ColorConvertFloat4ToU32(hi), 8);
}

void drawNucleus(ImDrawList* dl, const AtomModel& model, float yaw, float pitch,
                 ImVec2 centre, float scale, float unitR) {
  std::vector<Projected> pts;
  rotateProject(model.nucleus, pts, yaw, pitch, centre, scale);
  std::sort(pts.begin(), pts.end(),
            [](const Projected& a, const Projected& b) { return a.depth > b.depth; });
  const float r = std::max(1.1f, unitR * scale * 0.055f);
  for (const Projected& p : pts) {
    const float fog = std::clamp((p.depth / 1.2f + 1.0f) * 0.5f * 0.45f, 0.0f, 0.5f);
    shadedSphere(dl, ImVec2(p.x, p.y), r * p.persp,
                 p.src->proton ? ImVec4(0.87f, 0.34f, 0.31f, 1.0f)
                               : ImVec4(0.76f, 0.76f, 0.81f, 1.0f),
                 fog);
  }
}

void drawOrbitalCloud(const ElementData& element, ImVec2 min, ImVec2 max) {
  const AtomModel& model = modelFor(element);
  ImDrawList* dl = ImGui::GetWindowDrawList();
  const double t = ImGui::GetTime();
  const float yaw = static_cast<float>(t * 0.55);
  const float pitch = 0.44f + 0.10f * std::sin(static_cast<float>(t) * 0.30f);
  const ImVec2 centre((min.x + max.x) * 0.5f, (min.y + max.y) * 0.5f);
  const float scale = std::min(max.x - min.x, max.y - min.y) * 0.5f / 1.85f;

  std::vector<Projected> pts;
  rotateProject(model.points, pts, yaw, pitch, centre, scale);
  std::sort(pts.begin(), pts.end(),
            [](const Projected& a, const Projected& b) { return a.depth > b.depth; });

  for (const Projected& p : pts) {
    const ImVec4 base = cloudColor(p.src->l);
    const float depthT = std::clamp((p.depth / 1.6f + 1.0f) * 0.5f, 0.0f, 1.0f);
    const float alpha = 0.16f + 0.72f * (1.0f - depthT);
    const float r = (1.05f + 1.05f * (1.0f - depthT)) * p.persp;
    dl->AddCircleFilled(ImVec2(p.x, p.y), r,
                        ImGui::ColorConvertFloat4ToU32(ImVec4(base.x, base.y, base.z, alpha)),
                        6);
  }

  drawNucleus(dl, model, yaw, pitch, centre, scale, 1.0f);
}

void drawBohr(const ElementData& element, ImVec2 min, ImVec2 max) {
  const AtomModel& model = modelFor(element);
  ImDrawList* dl = ImGui::GetWindowDrawList();
  const double t = ImGui::GetTime();
  const ImVec2 centre((min.x + max.x) * 0.5f, (min.y + max.y) * 0.5f);
  const float maxR = std::min(max.x - min.x, max.y - min.y) * 0.5f / 1.15f;
  const float tilt = 0.38f;
  const float yaw = static_cast<float>(t * 0.35);

  int outer = 7;
  while (outer > 1 && model.shells[outer - 1] == 0) --outer;

  for (int n = outer; n >= 1; --n) {
    const int count = model.shells[n - 1];
    if (count <= 0) continue;
    const float r = maxR * (0.24f + 0.76f * static_cast<float>(n) / outer);
    dl->AddEllipse(centre, ImVec2(r, r * tilt), style::u32(style::col::BorderStrong, 0.35f),
                   yaw * (n % 2 == 0 ? 0.12f : -0.09f), 0, 1.0f);
    const float speed = 1.9f / (0.6f + 0.4f * n);
    for (int e = 0; e < count; ++e) {
      const float a = static_cast<float>(t * speed * (n % 2 == 0 ? 1.0 : -1.0)) +
                      2.0f * static_cast<float>(kPi) * e / count;
      const float ca = std::cos(a + yaw * 0.4f), sa = std::sin(a + yaw * 0.4f);
      const float depth = 0.5f + 0.5f * sa;  // front when sa > 0
      const float er = 2.6f * (0.75f + 0.45f * depth);
      const ImVec4 col(0.55f + 0.35f * depth, 0.82f, 0.95f, 0.35f + 0.6f * depth);
      dl->AddCircleFilled(ImVec2(centre.x + ca * r, centre.y + sa * r * tilt), er,
                          ImGui::ColorConvertFloat4ToU32(col), 8);
    }
  }

  drawNucleus(dl, model, yaw, 0.35f, centre, maxR * 0.62f, 1.0f);
}

}  // namespace

std::string elementConfigString(const ElementData& element) {
  return configString(electronConfig(element.z), element.z);
}

void drawAtomModel(const ElementData& element, ImVec2 min, ImVec2 max, bool bohrMode) {
  ImDrawList* dl = ImGui::GetWindowDrawList();
  const style::Metrics& m = style::metrics();
  dl->AddRectFilled(min, max, style::u32(style::col::BgDeep), m.radiusMd);
  dl->AddRect(min, max, style::u32(style::col::Border), m.radiusMd, 0, m.hairline);
  dl->PushClipRect(min, max, true);
  if (bohrMode) {
    drawBohr(element, min, max);
  } else {
    drawOrbitalCloud(element, min, max);
  }
  dl->PopClipRect();
}

}  // namespace chemcad::ui

#include "sol/funnel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <unordered_map>
#include <utility>

// Design note: `Simulation` (the fixed public contract) carries no leftover-
// time accumulator, so `step` cannot maintain a persistent fixed-timestep
// carry across calls. Instead every call divides its own `dt` into the
// integer number of ~1/240 s slices closest to `dt`, so the *slice length* is
// pinned near the target substep regardless of frame rate and the total
// integrated time always equals `dt` exactly (no drift, no carry state).

namespace chemcad::sol {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kGravity = 9.80665;          // m/s^2
constexpr double kSubstepDt = 1.0 / 240.0;    // s, target physics substep
constexpr double kMlToM3 = 1e-6;

// Internal physics length scale: the column height a fully-charged vessel
// maps to. Purely an implementation constant for the Stokes/Weber maths
// below; it is independent of whatever height a renderer asks
// `vesselOutline` to draw at, so the two never need to agree.
constexpr double kReferenceHeightM = 0.22;
constexpr int kProfileSamples = 128;  // resolution of the height<->volume map

// Weber-style Sauter-radius scaling constant for `shake`. Calibrated so a
// typical organic/water interface (~30 mN/m) at a moderate shake (vigour ~
// 0.5) lands near 0.2 mm droplets, with size falling as shaking intensifies.
constexpr double kSauterK = 0.002;
constexpr double kMinDropletRadius = 20e-6;  // m
constexpr double kMaxDropletRadius = 3e-3;   // m
constexpr int kMaxDroplets = 4000;           // hard cap so a session never explodes

constexpr double kDropletRelaxTau = 0.12;   // s, velocity-toward-terminal time constant
constexpr double kShakeDecayRate = 0.15;    // 1/s, shakeEnergy exponential decay
constexpr double kStirAmplitude = 0.05;     // m/s, churn velocity scale at shakeEnergy = 1
constexpr double kCoalesceRate = 3.0;       // droplet-droplet merge attempts, see `step`
// Bulk re-absorption attempt rate at stability 0, see `integrateSubstep`. A
// droplet overlapping its own settled layer is offered this many absorption
// rolls per second, scaled by (1 - emulsionStability), instead of being
// deleted the instant it geometrically overlaps -- so a stable emulsion can
// actually retain droplets that `shake` happened to scatter inside their own
// layer band, and an unstable one still clears out well within tens of
// seconds.
constexpr double kBulkCoalesceRate = 4.0;   // 1/s

// -- Vessel cross-section profiles ------------------------------------------
// Each profile is the single source of truth for both `vesselWidthAt` (used
// by physics to clamp droplets) and `vesselOutline` (used to draw the same
// vessel), so the two are consistent by construction.

struct WidthPoint {
  double t;      // height fraction, 0 bottom .. 1 top
  double width;  // half-width fraction of the widest half-width, (0, 1]
};

double interpolateWidth(const WidthPoint* pts, size_t n, double t) {
  if (n == 0) return 1.0;
  if (t <= pts[0].t) return pts[0].width;
  if (t >= pts[n - 1].t) return pts[n - 1].width;
  for (size_t i = 0; i + 1 < n; ++i) {
    if (t < pts[i].t || t > pts[i + 1].t) continue;
    double span = pts[i + 1].t - pts[i].t;
    double localT = span > 0.0 ? (t - pts[i].t) / span : 0.0;
    double eased = localT * localT * (3.0 - 2.0 * localT);  // smoothstep
    return pts[i].width + (pts[i + 1].width - pts[i].width) * eased;
  }
  return pts[n - 1].width;
}

// Narrow stopcock stem, a pear-shaped bulb widening to about 60% of the
// height, tapering into a neck, and a slight ground-glass stopper flare.
constexpr std::array<WidthPoint, 12> kFunnelProfile = {{
    {0.00, 0.05}, {0.05, 0.09}, {0.10, 0.10}, {0.22, 0.16}, {0.35, 0.55}, {0.50, 0.90},
    {0.60, 1.00}, {0.70, 0.85}, {0.80, 0.42}, {0.85, 0.17}, {0.93, 0.17}, {1.00, 0.27},
}};

// Flat, wide conical base (Erlenmeyer) tapering to a straight narrow neck.
constexpr std::array<WidthPoint, 5> kFlaskProfile = {{
    {0.00, 1.00}, {0.06, 0.96}, {0.55, 0.34}, {0.66, 0.20}, {1.00, 0.20},
}};

// Flared foot, straight-walled tube, and a small pour-spout flare at the top.
constexpr std::array<WidthPoint, 7> kCylinderProfile = {{
    {0.00, 0.62}, {0.03, 1.00}, {0.09, 0.55}, {0.10, 0.55},
    {0.92, 0.55}, {0.96, 0.62}, {1.00, 0.66},
}};

// Aspect ratio (widest half-width / height) used only to give `vesselOutline`
// concrete metres; physics never reads this, it solves its own radius from
// charged volume.
double referenceMaxHalfWidth(Vessel vessel) {
  switch (vessel) {
    case Vessel::SeparatoryFunnel: return 0.30;
    case Vessel::DecantingFlask: return 0.42;
    case Vessel::GraduatedCylinder: return 0.22;
  }
  return 0.30;
}

double dropletVolumeM3(double radiusMetres) {
  return 4.0 / 3.0 * kPi * radiusMetres * radiusMetres * radiusMetres;
}
double dropletVolumeMl(double radiusMetres) { return dropletVolumeM3(radiusMetres) * 1e6; }

}  // namespace

double vesselWidthAt(Vessel vessel, double heightFraction) {
  double t = std::clamp(heightFraction, 0.0, 1.0);
  switch (vessel) {
    case Vessel::SeparatoryFunnel:
      return interpolateWidth(kFunnelProfile.data(), kFunnelProfile.size(), t);
    case Vessel::DecantingFlask:
      return interpolateWidth(kFlaskProfile.data(), kFlaskProfile.size(), t);
    case Vessel::GraduatedCylinder:
      return interpolateWidth(kCylinderProfile.data(), kCylinderProfile.size(), t);
  }
  return 1.0;
}

std::vector<core::Vec2> vesselOutline(Vessel vessel, double heightMetres) {
  heightMetres = std::max(heightMetres, 1e-6);
  double maxHalfWidth = referenceMaxHalfWidth(vessel) * heightMetres;
  constexpr int kSegments = 48;

  std::vector<core::Vec2> outline;
  outline.reserve(static_cast<size_t>(kSegments + 1) * 2);

  // Right wall, bottom to top, then the mirrored left wall back down: this
  // traces the boundary counter-clockwise for a y-up, x-right frame.
  for (int i = 0; i <= kSegments; ++i) {
    double t = double(i) / kSegments;
    double x = vesselWidthAt(vessel, t) * maxHalfWidth;
    outline.push_back({static_cast<float>(x), static_cast<float>(t * heightMetres)});
  }
  for (int i = kSegments; i >= 0; --i) {
    double t = double(i) / kSegments;
    double x = -vesselWidthAt(vessel, t) * maxHalfWidth;
    outline.push_back({static_cast<float>(x), static_cast<float>(t * heightMetres)});
  }
  return outline;
}

namespace {

// Cumulative volume (m^3) a unit-maxRadius vessel of this shape holds from
// the bottom up to each sampled height fraction, assuming the cross-section
// is a circle of radius `maxRadius * vesselWidthAt(...)` (axisymmetric
// revolution of the same profile the outline draws). Rescaling by the real
// maxRadius^2 turns this into real volume for any charged amount.
struct HeightProfile {
  std::array<double, kProfileSamples + 1> cumulativeVolume{};
};

const HeightProfile& heightProfile(Vessel vessel) {
  static const std::array<HeightProfile, 3> kProfiles = [] {
    std::array<HeightProfile, 3> result;
    for (int v = 0; v < 3; ++v) {
      Vessel vessel = static_cast<Vessel>(v);
      double acc = 0.0;
      result[v].cumulativeVolume[0] = 0.0;
      for (int i = 1; i <= kProfileSamples; ++i) {
        double h0 = double(i - 1) / kProfileSamples;
        double h1 = double(i) / kProfileSamples;
        double w0 = vesselWidthAt(vessel, h0);
        double w1 = vesselWidthAt(vessel, h1);
        double wm = vesselWidthAt(vessel, (h0 + h1) * 0.5);
        // Simpson's rule per slab; area is pi*r^2 for a unit maxRadius.
        double area = (w0 * w0 + 4.0 * wm * wm + w1 * w1) / 6.0 * kPi;
        acc += area * (h1 - h0) * kReferenceHeightM;
        result[v].cumulativeVolume[i] = acc;
      }
    }
    return result;
  }();
  return kProfiles[static_cast<int>(vessel)];
}

// Inverts the profile's cumulative-volume curve: what height fraction holds
// exactly `volumeOverMaxRadius2` (i.e. real volume / maxRadius^2)?
double heightFractionForVolume(const HeightProfile& profile, double volumeOverMaxRadius2) {
  const auto& c = profile.cumulativeVolume;
  if (volumeOverMaxRadius2 <= 0.0) return 0.0;
  if (volumeOverMaxRadius2 >= c.back()) return 1.0;
  int lo = 0, hi = kProfileSamples;
  while (lo + 1 < hi) {
    int mid = (lo + hi) / 2;
    if (c[mid] <= volumeOverMaxRadius2) lo = mid; else hi = mid;
  }
  double segSpan = c[hi] - c[lo];
  double segFrac = segSpan > 0.0 ? (volumeOverMaxRadius2 - c[lo]) / segSpan : 0.0;
  return (double(lo) + segFrac) / kProfileSamples;
}

struct ColumnGeometry {
  double maxRadius = 0.0;   // m, solved so the vessel's rated capacity fills kReferenceHeightM
  double fillHeight = 0.0;  // m, height of the total charged liquid column
};

ColumnGeometry columnGeometry(const Simulation& sim) {
  const HeightProfile& profile = heightProfile(sim.vessel);
  double capacityM3 = std::max(1e-9, sim.vesselVolumeMl) * kMlToM3;
  double capacityUnit = std::max(1e-18, profile.cumulativeVolume.back());

  ColumnGeometry geo;
  geo.maxRadius = std::sqrt(capacityM3 / capacityUnit);
  double totalUnit = totalVolumeMl(sim) * kMlToM3 / std::max(1e-18, geo.maxRadius * geo.maxRadius);
  geo.fillHeight = heightFractionForVolume(profile, totalUnit) * kReferenceHeightM;
  return geo;
}

// Height (m) of the top of each bulk layer, bottom to top, index 0 is always
// 0. boundaries[phases.size()] is the top of the settled stack, which is at
// or below `geo.fillHeight` -- the gap between the two is where dispersed
// droplets currently roam.
std::vector<double> layerBoundaries(const Simulation& sim, const ColumnGeometry& geo) {
  std::vector<double> boundaries(sim.phases.size() + 1, 0.0);
  const HeightProfile& profile = heightProfile(sim.vessel);
  double maxR2 = std::max(1e-18, geo.maxRadius * geo.maxRadius);
  double cumMl = 0.0;
  for (size_t i = 0; i < sim.phases.size(); ++i) {
    cumMl += sim.settledMl[i];
    double unit = cumMl * kMlToM3 / maxR2;
    boundaries[i + 1] = heightFractionForVolume(profile, unit) * kReferenceHeightM;
  }
  return boundaries;
}

// Which phase index a droplet at height `y` is currently surrounded by. Above
// the whole bulk stack (the dispersed zone) counts as the topmost, lightest
// phase, matching where that space actually drains into as it settles.
int continuousPhaseAt(const std::vector<double>& boundaries, double y) {
  int n = static_cast<int>(boundaries.size()) - 1;
  if (n <= 0) return 0;
  for (int i = 0; i < n; ++i) {
    if (y < boundaries[i + 1] || i == n - 1) return i;
  }
  return n - 1;
}

// A `float` radius cannot exactly represent every reachable combined volume:
// `cbrt` runs in double precision, but rounding its result down to the
// droplet's `float` storage silently drops (or adds) a sliver of volume that
// then has nowhere else to go. Left alone it just vanishes from the
// simulation, and thousands of merges per run compound that into a visible
// drift (this is the bug: mm-scale droplets lose enough per merge that it
// used to be invisible at micron scale but is not anymore). Parking the
// exact residual in `settledMl` -- a `double` with no such precision
// ceiling -- keeps the total exact by construction instead of merely close.
// Returns false (no mutation of either droplet) in the rare case where even
// shrinking the merged radius by one ULP would still drive settledMl
// negative; the caller must then leave both droplets unmerged rather than
// lose volume.
bool mergeDroplets(Droplet& survivor, const Droplet& absorbed, std::vector<double>& settledMl) {
  double va = dropletVolumeM3(survivor.radius);
  double vb = dropletVolumeM3(absorbed.radius);
  double vt = va + vb;
  if (vt <= 0.0) return false;

  double combinedVolumeMl = vt * 1e6;
  float mergedRadius = static_cast<float>(std::cbrt(3.0 * vt / (4.0 * kPi)));
  double residualMl = combinedVolumeMl - dropletVolumeMl(static_cast<double>(mergedRadius));
  double& settled = settledMl[static_cast<size_t>(survivor.phase)];
  if (residualMl < 0.0 && settled + residualMl < 0.0) {
    float shrunk = std::nextafter(mergedRadius, 0.0f);
    double shrunkResidual = combinedVolumeMl - dropletVolumeMl(static_cast<double>(shrunk));
    if (settled + shrunkResidual < 0.0) return false;
    mergedRadius = shrunk;
    residualMl = shrunkResidual;
  }

  double wa = va / vt, wb = vb / vt;
  survivor.position.x = static_cast<float>(survivor.position.x * wa + absorbed.position.x * wb);
  survivor.position.y = static_cast<float>(survivor.position.y * wa + absorbed.position.y * wb);
  survivor.velocity.x = static_cast<float>(survivor.velocity.x * wa + absorbed.velocity.x * wb);
  survivor.velocity.y = static_cast<float>(survivor.velocity.y * wa + absorbed.velocity.y * wb);
  survivor.radius = mergedRadius;
  settled += residualMl;
  return true;
}

// Droplet-droplet coalescence via a uniform grid keyed by cell: only droplets
// in the same or an adjacent cell are ever compared, so this stays O(n) on
// average (O(n log n) worst case via the hash map) instead of O(n^2).
void coalesceDroplets(Simulation& sim, double dt, std::mt19937& rng) {
  size_t n = sim.droplets.size();
  if (n < 2 || dt <= 0.0) return;

  float maxRadius = 0.0f;
  for (const auto& d : sim.droplets) maxRadius = std::max(maxRadius, d.radius);
  double cell = std::max(double(maxRadius) * 2.0, 1e-5);

  auto cellKey = [cell](const core::Vec2& p) -> uint64_t {
    auto cx = static_cast<int64_t>(std::floor(double(p.x) / cell));
    auto cy = static_cast<int64_t>(std::floor(double(p.y) / cell));
    return (static_cast<uint64_t>(cx) << 32) ^ (static_cast<uint64_t>(cy) & 0xffffffffull);
  };

  std::unordered_map<uint64_t, std::vector<uint32_t>> grid;
  grid.reserve(n * 2);
  for (uint32_t i = 0; i < n; ++i) grid[cellKey(sim.droplets[i].position)].push_back(i);

  std::vector<bool> removed(n, false);
  std::uniform_real_distribution<double> unit(0.0, 1.0);

  for (uint32_t i = 0; i < n; ++i) {
    if (removed[i]) continue;
    auto cx = static_cast<int64_t>(std::floor(double(sim.droplets[i].position.x) / cell));
    auto cy = static_cast<int64_t>(std::floor(double(sim.droplets[i].position.y) / cell));
    for (int64_t ny = cy - 1; ny <= cy + 1 && !removed[i]; ++ny) {
      for (int64_t nx = cx - 1; nx <= cx + 1 && !removed[i]; ++nx) {
        uint64_t key =
            (static_cast<uint64_t>(nx) << 32) ^ (static_cast<uint64_t>(ny) & 0xffffffffull);
        auto it = grid.find(key);
        if (it == grid.end()) continue;
        for (uint32_t j : it->second) {
          if (j <= i || removed[j] || removed[i]) continue;
          const Droplet& di = sim.droplets[i];
          const Droplet& dj = sim.droplets[j];
          if (dj.phase != di.phase) continue;
          double dx = double(di.position.x) - double(dj.position.x);
          double dy = double(di.position.y) - double(dj.position.y);
          double touch = double(di.radius) + double(dj.radius);
          if (dx * dx + dy * dy > touch * touch) continue;
          double stability = sim.phases[di.phase].emulsionStability;
          double p = dt * (1.0 - stability) * kCoalesceRate;
          if (unit(rng) >= p) continue;
          if (mergeDroplets(sim.droplets[i], sim.droplets[j], sim.settledMl)) {
            removed[j] = true;
          }
        }
      }
    }
  }

  if (std::none_of(removed.begin(), removed.end(), [](bool b) { return b; })) return;
  std::vector<Droplet> kept;
  kept.reserve(n);
  for (uint32_t i = 0; i < n; ++i) {
    if (!removed[i]) kept.push_back(sim.droplets[i]);
  }
  sim.droplets = std::move(kept);
}

// One fixed-length physics slice: Stokes relaxation, stirring, wall clamp and
// return-to-bulk coalescence for every droplet.
void integrateSubstep(Simulation& sim, double dt, std::mt19937& rng) {
  if (sim.phases.empty()) return;

  sim.shakeEnergy *= std::exp(-kShakeDecayRate * dt);

  ColumnGeometry geo = columnGeometry(sim);
  std::vector<double> boundaries = layerBoundaries(sim, geo);
  double blend = 1.0 - std::exp(-dt / kDropletRelaxTau);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  for (size_t di = 0; di < sim.droplets.size();) {
    Droplet& d = sim.droplets[di];
    int cont = continuousPhaseAt(boundaries, d.position.y);
    double rhoDrop = sim.phases[d.phase].density * 1000.0;
    double rhoCont = sim.phases[cont].density * 1000.0;
    double muCont = std::max(1e-4, sim.phases[cont].viscosity) * 1e-3;
    // A stable emulsion resists the very buoyancy that would break it (the
    // surfactant film blocks drainage), so stability directly damps mobility.
    double mobility = std::max(0.0, 1.0 - sim.phases[d.phase].emulsionStability);
    double r = double(d.radius);
    double vTerminal = -mobility * 2.0 * r * r * (rhoDrop - rhoCont) * kGravity / (9.0 * muCont);
    d.velocity.y = static_cast<float>(d.velocity.y + (vTerminal - double(d.velocity.y)) * blend);

    double stir = 0.0;
    if (sim.shakeEnergy > 1e-4) {
      stir = kStirAmplitude * sim.shakeEnergy *
             std::sin(2.0 * kPi * (sim.elapsed * 1.7 + double(d.position.x) * 37.0 +
                                    double(d.position.y) * 53.0));
    }
    d.velocity.x = static_cast<float>(d.velocity.x + (stir - double(d.velocity.x)) * blend);

    d.position.x = static_cast<float>(double(d.position.x) + double(d.velocity.x) * dt);
    d.position.y = static_cast<float>(double(d.position.y) + double(d.velocity.y) * dt);

    double vesselFrac = std::clamp(double(d.position.y) / kReferenceHeightM, 0.0, 1.0);
    double halfWidth = geo.maxRadius * vesselWidthAt(sim.vessel, vesselFrac);
    double xLimit = std::max(0.0, halfWidth - r);
    d.position.x = static_cast<float>(std::clamp(double(d.position.x), -xLimit, xLimit));
    double yLow = r;
    double yHigh = std::max(yLow, geo.fillHeight - r);
    d.position.y = static_cast<float>(std::clamp(double(d.position.y), yLow, yHigh));

    double ownLow = boundaries[static_cast<size_t>(d.phase)];
    double ownHigh = boundaries[static_cast<size_t>(d.phase) + 1];
    bool reachedOwnLayer =
        double(d.position.y) - r <= ownHigh && double(d.position.y) + r >= ownLow;
    if (reachedOwnLayer) {
      // Instant geometric absorption would ignore stability entirely (and
      // would even fire on droplets `shake` placed inside their own layer at
      // t=0); gate it by a per-substep probability instead, exactly like
      // droplet-droplet coalescence below.
      double stability = sim.phases[d.phase].emulsionStability;
      double rate = kBulkCoalesceRate * (1.0 - stability);
      double p = rate > 0.0 ? 1.0 - std::exp(-rate * dt) : 0.0;
      if (unit(rng) < p) {
        // Exact by construction, unlike a merge: this recomputes the volume
        // forward from the same `float` radius the droplet already carries
        // (the same value `dropletVolumeMl(sim)` will no longer count once
        // the droplet is erased below), so nothing is rounded away.
        sim.settledMl[static_cast<size_t>(d.phase)] += dropletVolumeMl(double(d.radius));
        sim.droplets[di] = sim.droplets.back();
        sim.droplets.pop_back();
        continue;  // re-examine the swapped-in element at this index
      }
    }
    ++di;
  }
}

}  // namespace

void reset(Simulation& sim) {
  std::stable_sort(sim.phases.begin(), sim.phases.end(),
                    [](const Phase& a, const Phase& b) { return a.density > b.density; });
  sim.settledMl.resize(sim.phases.size());
  // Straight double-to-double copy of the charged amount, no cbrt/float
  // round-trip involved, and droplets are wiped rather than reconciled --
  // nothing here can lose volume the way a merge can.
  for (size_t i = 0; i < sim.phases.size(); ++i) sim.settledMl[i] = sim.phases[i].volumeMl;
  sim.droplets.clear();
  sim.shakeEnergy = 0.0;
  sim.elapsed = 0.0;
}

void shake(Simulation& sim, double vigour) {
  vigour = std::clamp(vigour, 0.0, 1.0);
  sim.shakeEnergy = std::min(1.0, sim.shakeEnergy + vigour);
  if (vigour <= 0.0 || sim.phases.empty()) return;

  std::mt19937 rng(sim.seed);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  std::uniform_real_distribution<double> jitter(-0.02, 0.02);
  ColumnGeometry geo = columnGeometry(sim);

  // First pass: which phases actually disperse anything this call? The
  // droplet budget below is split evenly across only those phases, so a
  // phase early in `sim.phases` cannot claim the whole cap and starve a
  // later one.
  std::vector<size_t> dispersingPhases;
  std::vector<double> dispersedMl(sim.phases.size(), 0.0);
  for (size_t i = 0; i < sim.phases.size(); ++i) {
    double dispersed = sim.settledMl[i] * vigour * (1.0 - std::exp(-vigour * 3.0));
    dispersed = std::clamp(dispersed, 0.0, sim.settledMl[i]);
    dispersedMl[i] = dispersed;
    if (dispersed > 0.0) dispersingPhases.push_back(i);
  }
  if (dispersingPhases.empty()) return;

  int perPhaseBudget = kMaxDroplets / static_cast<int>(dispersingPhases.size());

  for (size_t i : dispersingPhases) {
    double dispersed = dispersedMl[i];

    double tension = sim.phases[i].interfacialTension;
    double denom = 1000.0 * (0.05 + vigour) * (0.05 + vigour);
    double sauterRadiusM = kSauterK * tension / denom;

    int capacityLeft =
        std::min(perPhaseBudget, kMaxDroplets - static_cast<int>(sim.droplets.size()));
    capacityLeft = std::max(0, capacityLeft);
    if (capacityLeft <= 0) continue;

    // A real emulsion shatters into billions of droplets; this simulation
    // caps the droplet count at kMaxDroplets for performance, so the
    // droplet budget wins over the physical Sauter radius whenever the
    // physical size would need more droplets than the budget allows to
    // carry the dispersed volume. Growing the radius to fit the budget
    // trades droplet count for a dispersed volume that is still visually
    // and physically representative, instead of clamping to the physical
    // radius and silently dispersing almost nothing.
    double dispersedM3 = dispersed * kMlToM3;
    double budgetRadiusM = std::cbrt(3.0 * dispersedM3 / (4.0 * kPi * capacityLeft));
    double radiusM =
        std::clamp(std::max(sauterRadiusM, budgetRadiusM), kMinDropletRadius, kMaxDropletRadius);
    // Round to the float precision the droplet will actually store so the
    // volume subtracted here always matches the volume recovered later.
    // Unlike `mergeDroplets`, this needs no residual: the rounding happens
    // once, going forward (radius chosen -> volume derived from it), so the
    // same double value gets subtracted here and re-summed from the droplet
    // later -- there is no separate "true" volume it could disagree with.
    float radiusF = static_cast<float>(radiusM);
    double dropVolMl = dropletVolumeMl(static_cast<double>(radiusF));
    if (dropVolMl <= 0.0) continue;

    // If even the budget at kMaxDropletRadius cannot carry the full
    // dispersed volume, count falls short of capacityLeft here and the
    // undispersed remainder below simply stays settled.
    int count = std::min(static_cast<int>(dispersed / dropVolMl), capacityLeft);
    if (count <= 0) continue;

    sim.settledMl[i] -= count * dropVolMl;  // the undispersed remainder simply stays settled

    for (int k = 0; k < count; ++k) {
      Droplet d;
      d.phase = static_cast<int>(i);
      d.radius = radiusF;
      double placementFrac = unit(rng);
      double dropletY = placementFrac * geo.fillHeight;
      double vesselFrac = std::clamp(dropletY / kReferenceHeightM, 0.0, 1.0);
      double halfWidth = geo.maxRadius * vesselWidthAt(sim.vessel, vesselFrac);
      double xLimit = std::max(0.0, halfWidth - radiusM);
      d.position.y = static_cast<float>(dropletY);
      d.position.x = static_cast<float>((unit(rng) * 2.0 - 1.0) * xLimit);
      d.velocity.x = static_cast<float>(jitter(rng));
      d.velocity.y = static_cast<float>(jitter(rng));
      sim.droplets.push_back(d);
    }
  }

  sim.seed = rng();
}

void step(Simulation& sim, double dt) {
  dt = std::clamp(dt, 0.0, 0.1);
  if (dt <= 0.0) return;

  std::mt19937 rng(sim.seed);

  int substeps = std::max(1, static_cast<int>(std::lround(dt / kSubstepDt)));
  double subDt = dt / substeps;
  for (int s = 0; s < substeps; ++s) integrateSubstep(sim, subDt, rng);

  coalesceDroplets(sim, dt, rng);
  sim.seed = rng();

  sim.elapsed += dt;
}

double totalVolumeMl(const Simulation& sim) {
  double total = 0.0;
  for (const auto& phase : sim.phases) total += phase.volumeMl;
  return total;
}

double emulsifiedFraction(const Simulation& sim) {
  double total = totalVolumeMl(sim);
  if (total <= 0.0) return 0.0;
  double dispersed = 0.0;
  for (const auto& d : sim.droplets) dispersed += dropletVolumeMl(static_cast<double>(d.radius));
  return std::clamp(dispersed / total, 0.0, 1.0);
}

}  // namespace chemcad::sol

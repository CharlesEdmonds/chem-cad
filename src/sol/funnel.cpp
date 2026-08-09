#include "sol/funnel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <unordered_map>
#include <utility>

// `Simulation` intentionally carries no fixed-step accumulator. Each call to
// step divides its own duration into slices near 1/240 s, preserving the exact
// requested elapsed time without hidden state.

namespace chemcad::sol {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kGravity = 9.80665;
constexpr double kSubstepDt = 1.0 / 240.0;
constexpr double kMlToM3 = 1e-6;
constexpr int kProfileSamples = 1024;

// Hinze (1955): d32 = C_H (sigma/rho_c)^0.6 epsilon^-0.4. The physical
// droplet size controls drag; parcel volume is deliberately independent.
constexpr double kHinzeC = 0.725;
constexpr double kWeberCrit = 12.0;
constexpr double kMinDropletRadius = 20e-6;
constexpr double kMaxDropletRadius = 3e-3;
constexpr int kTargetParcels = 900;
constexpr int kMaxDroplets = 4000;

constexpr double kDropletRelaxTau = 0.12;
constexpr double kShakeDecayRate = 0.15;
constexpr double kCoalesceRate = 3.0;
constexpr double kBulkCoalesceRate = 4.0;

inline double safeRoot(double value) { return std::sqrt(std::max(0.0, value)); }

// The Squibb profile is constructed from exact circle equations in the
// (height fraction t, half-width fraction w) plane. The stem fillet has centre
// (0.34, 0.1636815104), radius 0.0636815104, and meets the cone with slope
// 2.8118091744. A small crown circle, centre (0.72, 0.9787728299), makes that
// straight cone tangent to the horizontal at the true widest point. The steep
// shoulder uses two equal radius 0.0226556444 arcs joined by their common
// slope -8 tangent; its final arc is horizontal where it meets the neck. This
// arc-line-arc construction is necessary for a high, broad pear shoulder to
// lose 0.8 width in only 0.14 height without a spline or a kink.
double funnelWidth(double t) {
  constexpr double kFilletCentreW = 0.163681510413514;
  constexpr double kFilletRadius = 0.063681510413514;
  constexpr double kConeStartW = 0.14234293575290247;
  constexpr double kConeSlope = 2.811809174431545;
  constexpr double kCrownCentreW = 0.978772829862162;
  constexpr double kCrownRadius = 0.021227170137838036;
  constexpr double kShoulderRadius = 0.022655644370746332;
  constexpr double kShoulderLineT0 = 0.7424806946917841;
  constexpr double kShoulderLineW0 = 0.9801544424657267;
  constexpr double kShoulderLineT1 = 0.8375193053082158;

  if (t <= 0.34) return 0.10;
  if (t <= 0.40) {
    const double dt = t - 0.34;
    return kFilletCentreW - safeRoot(kFilletRadius * kFilletRadius - dt * dt);
  }
  if (t <= 0.70) return kConeStartW + kConeSlope * (t - 0.40);
  if (t <= 0.72) {
    const double dt = t - 0.72;
    return kCrownCentreW + safeRoot(kCrownRadius * kCrownRadius - dt * dt);
  }
  if (t <= kShoulderLineT0) {
    const double dt = t - 0.72;
    return 1.0 - kShoulderRadius +
           safeRoot(kShoulderRadius * kShoulderRadius - dt * dt);
  }
  if (t <= kShoulderLineT1) return kShoulderLineW0 - 8.0 * (t - kShoulderLineT0);
  if (t <= 0.86) {
    const double dt = t - 0.86;
    return 0.20 + kShoulderRadius -
           safeRoot(kShoulderRadius * kShoulderRadius - dt * dt);
  }
  return 0.20;
}

// Flat base cut, straight Erlenmeyer cone, then a circle tangent to both the
// cone and the vertical neck. Circle centre (0.66, 0.3248079393), radius
// 0.1248079393; its cone tangent is (0.56, 0.2501284796).
double flaskWidth(double t) {
  constexpr double kTangentW = 0.25012847955050965;
  constexpr double kCentreW = 0.324807939261731;
  constexpr double kRadius = 0.124807939261731;
  if (t <= 0.56) return 1.0 + (kTangentW - 1.0) * (t / 0.56);
  if (t <= 0.66) {
    const double dt = t - 0.66;
    return kCentreW - safeRoot(kRadius * kRadius - dt * dt);
  }
  return 0.20;
}

// The cylinder foot is a straight flare followed by a true circle tangent to
// the tube. At the top, a quarter-circle flare leaves the straight tube
// tangentially and forms the pour-spout lip without a control-point bump.
double cylinderWidth(double t) {
  constexpr double kFootTangentW = 0.5846370932138486;
  constexpr double kFootCentreW = 0.5904151729623988;
  constexpr double kFootRadius = 0.040415172962398715;
  if (t <= 0.06) return 1.0 + (kFootTangentW - 1.0) * (t / 0.06);
  if (t <= 0.10) {
    const double dt = t - 0.10;
    return kFootCentreW - safeRoot(kFootRadius * kFootRadius - dt * dt);
  }
  if (t <= 0.94) return 0.55;
  const double dt = t - 0.94;
  return 0.61 - safeRoot(0.06 * 0.06 - dt * dt);
}

// Widest half-width / full height. With the profile integral, these are the
// only vessel calibration constants. The 250 mL Squibb value gives H=0.190 m
// and a 90.4 mm maximum body diameter.
double referenceMaxHalfWidth(Vessel vessel) {
  switch (vessel) {
    case Vessel::SeparatoryFunnel: return 0.238;
    case Vessel::DecantingFlask: return 0.32;
    case Vessel::GraduatedCylinder: return 0.17;
  }
  return 0.238;
}

double parcelCloudRadius(double parcelMl) {
  return std::cbrt(3.0 * std::max(0.0, parcelMl) * kMlToM3 / (4.0 * kPi));
}

}  // namespace

double vesselWidthAt(Vessel vessel, double heightFraction) {
  const double t = std::clamp(heightFraction, 0.0, 1.0);
  switch (vessel) {
    case Vessel::SeparatoryFunnel: return funnelWidth(t);
    case Vessel::DecantingFlask: return flaskWidth(t);
    case Vessel::GraduatedCylinder: return cylinderWidth(t);
  }
  return 1.0;
}

std::vector<core::Vec2> vesselOutline(Vessel vessel, double heightMetres) {
  heightMetres = std::max(heightMetres, 1e-6);
  const double maxHalfWidth = referenceMaxHalfWidth(vessel) * heightMetres;
  constexpr int kSegmentsPerSide = 256;

  std::vector<core::Vec2> outline;
  outline.reserve(static_cast<size_t>(kSegmentsPerSide + 1) * 2 + 1);
  for (int i = 0; i <= kSegmentsPerSide; ++i) {
    const double t = double(i) / kSegmentsPerSide;
    outline.push_back({static_cast<float>(vesselWidthAt(vessel, t) * maxHalfWidth),
                       static_cast<float>(t * heightMetres)});
  }
  for (int i = kSegmentsPerSide; i >= 0; --i) {
    const double t = double(i) / kSegmentsPerSide;
    outline.push_back({static_cast<float>(-vesselWidthAt(vessel, t) * maxHalfWidth),
                       static_cast<float>(t * heightMetres)});
  }
  outline.push_back(outline.front());
  return outline;
}

namespace {

// Dimensionless cumulative integral of w(t)^2. Multiplication by
// pi*maxRadius^2*height gives the axisymmetric volume.
struct HeightProfile {
  std::array<double, kProfileSamples + 1> cumulativeIntegral{};
};

const HeightProfile& heightProfile(Vessel vessel) {
  static const std::array<HeightProfile, 3> kProfiles = [] {
    std::array<HeightProfile, 3> result;
    for (int v = 0; v < 3; ++v) {
      const Vessel shape = static_cast<Vessel>(v);
      double integral = 0.0;
      for (int i = 1; i <= kProfileSamples; ++i) {
        const double t0 = double(i - 1) / kProfileSamples;
        const double t1 = double(i) / kProfileSamples;
        const double tm = 0.5 * (t0 + t1);
        const double w0 = vesselWidthAt(shape, t0);
        const double wm = vesselWidthAt(shape, tm);
        const double w1 = vesselWidthAt(shape, t1);
        integral += (w0 * w0 + 4.0 * wm * wm + w1 * w1) *
                    (t1 - t0) / 6.0;
        result[v].cumulativeIntegral[i] = integral;
      }
    }
    return result;
  }();
  return kProfiles[static_cast<int>(vessel)];
}

double heightFractionForIntegral(const HeightProfile& profile, double target) {
  const auto& c = profile.cumulativeIntegral;
  if (target <= 0.0) return 0.0;
  if (target >= c.back()) return 1.0;
  int lo = 0;
  int hi = kProfileSamples;
  while (lo + 1 < hi) {
    const int mid = (lo + hi) / 2;
    if (c[mid] <= target) lo = mid;
    else hi = mid;
  }
  const double span = c[hi] - c[lo];
  const double local = span > 0.0 ? (target - c[lo]) / span : 0.0;
  return (double(lo) + local) / kProfileSamples;
}

}  // namespace

double columnHeightM(const Simulation& sim) {
  const double volumeM3 = std::max(1e-12, sim.vesselVolumeMl * kMlToM3);
  const double a = referenceMaxHalfWidth(sim.vessel);
  const double integral = std::max(1e-12, heightProfile(sim.vessel).cumulativeIntegral.back());
  return std::cbrt(volumeM3 / (kPi * a * a * integral));
}

namespace {

struct ColumnGeometry {
  double height = 0.0;
  double maxRadius = 0.0;
  double fillHeight = 0.0;
};

ColumnGeometry columnGeometry(const Simulation& sim) {
  ColumnGeometry geo;
  geo.height = columnHeightM(sim);
  geo.maxRadius = referenceMaxHalfWidth(sim.vessel) * geo.height;
  const double volumeScale =
      std::max(1e-18, kPi * geo.maxRadius * geo.maxRadius * geo.height);
  const double chargedIntegral = totalVolumeMl(sim) * kMlToM3 / volumeScale;
  geo.fillHeight = heightFractionForIntegral(heightProfile(sim.vessel), chargedIntegral) * geo.height;
  return geo;
}

std::vector<double> layerBoundaries(const Simulation& sim, const ColumnGeometry& geo) {
  std::vector<double> boundaries(sim.phases.size() + 1, 0.0);
  const double volumeScale =
      std::max(1e-18, kPi * geo.maxRadius * geo.maxRadius * geo.height);
  double cumulativeMl = 0.0;
  for (size_t i = 0; i < sim.phases.size(); ++i) {
    cumulativeMl += sim.settledMl[i];
    boundaries[i + 1] =
        heightFractionForIntegral(heightProfile(sim.vessel),
                                  cumulativeMl * kMlToM3 / volumeScale) *
        geo.height;
  }
  return boundaries;
}

int continuousPhaseAt(const std::vector<double>& boundaries, double y) {
  const int n = static_cast<int>(boundaries.size()) - 1;
  if (n <= 0) return 0;
  for (int i = 0; i < n; ++i) {
    if (y < boundaries[i + 1] || i == n - 1) return i;
  }
  return n - 1;
}

double halfWidthAt(const Simulation& sim, const ColumnGeometry& geo, double y) {
  return geo.maxRadius *
         vesselWidthAt(sim.vessel, y / std::max(geo.height, 1e-12));
}

bool parcelFitsColumn(const Simulation& sim, const ColumnGeometry& geo, double parcelMl) {
  const double cloud = parcelCloudRadius(parcelMl);
  if (geo.fillHeight < 2.0 * cloud) return false;
  constexpr int kFitSamples = 512;
  for (int i = 0; i <= kFitSamples; ++i) {
    const double y = cloud + (geo.fillHeight - 2.0 * cloud) * double(i) / kFitSamples;
    if (halfWidthAt(sim, geo, y) >= cloud) return true;
  }
  return false;
}

void clampParcelToVessel(const Simulation& sim, const ColumnGeometry& geo, Droplet& d) {
  const double cloud = parcelCloudRadius(d.parcelMl);
  const double low = cloud;
  const double high = std::max(low, geo.fillHeight - cloud);
  double y = std::clamp(double(d.position.y), low, high);

  if (halfWidthAt(sim, geo, y) < cloud) {
    constexpr int kFitSamples = 512;
    double nearest = y;
    double nearestDistance = std::numeric_limits<double>::infinity();
    for (int i = 0; i <= kFitSamples; ++i) {
      const double candidate = low + (high - low) * double(i) / kFitSamples;
      if (halfWidthAt(sim, geo, candidate) < cloud) continue;
      const double distance = std::abs(candidate - y);
      if (distance < nearestDistance) {
        nearest = candidate;
        nearestDistance = distance;
      }
    }
    y = nearest;
  }

  d.position.y = static_cast<float>(y);
  const double xLimit = std::max(0.0, halfWidthAt(sim, geo, y) - cloud);
  d.position.x = static_cast<float>(std::clamp(double(d.position.x), -xLimit, xLimit));
}

double continuousDensity(const Simulation& sim, size_t excludeIndex) {
  double mass = 0.0;
  double volume = 0.0;
  for (size_t i = 0; i < sim.phases.size(); ++i) {
    if (i == excludeIndex) continue;
    mass += sim.phases[i].density * sim.phases[i].volumeMl;
    volume += sim.phases[i].volumeMl;
  }
  if (volume <= 0.0) return sim.phases[excludeIndex].density * 1000.0;
  return mass / volume * 1000.0;
}

// Volume-weighted viscosity of everything except `excludeIndex`: the drag the
// continuous medium exerts, in mPa.s.
double continuousViscosity(const Simulation& sim, size_t excludeIndex) {
  double weighted = 0.0;
  double volume = 0.0;
  for (size_t i = 0; i < sim.phases.size(); ++i) {
    if (i == excludeIndex) continue;
    weighted += sim.phases[i].viscosity * sim.phases[i].volumeMl;
    volume += sim.phases[i].volumeMl;
  }
  if (volume <= 0.0) return sim.phases[excludeIndex].viscosity;
  return weighted / volume;
}

double weberRadius(const Simulation& sim, size_t i, double u) {
  const double sigma = std::max(1e-6, sim.phases[i].interfacialTension * 1e-3);
  const double rhoC = std::max(1.0, continuousDensity(sim, i));
  return 0.5 * kWeberCrit * sigma / (rhoC * std::max(u * u, 1e-9));
}

double physicalDropletRadius(const Simulation& sim, size_t i, double epsilon, double u) {
  const double sigma = std::max(1e-6, sim.phases[i].interfacialTension * 1e-3);
  const double rhoC = std::max(1.0, continuousDensity(sim, i));
  const double hinze = 0.5 * kHinzeC * std::pow(sigma / rhoC, 0.6) *
                       std::pow(std::max(epsilon, 1e-9), -0.4);
  return std::clamp(std::min(hinze, weberRadius(sim, i, u)),
                    kMinDropletRadius, kMaxDropletRadius);
}

// Coalescence grows DROPLETS, not parcels. A parcel is a statistical carrier
// for a slice of dispersed volume; merging parcels was both wrong and ugly --
// it collapsed 900 carriers into a few dozen mL-scale blobs (which the
// renderer then had to draw as beach balls) while leaving the physical droplet
// radius, the quantity that actually governs settling, almost untouched.
//
// Smoluchowski coagulation for a nearly monodisperse emulsion: the collision
// rate per droplet scales with the local dispersed volume fraction, and since
// each collision doubles a droplet's volume, the population mean radius obeys
//
//     d(r^3)/dt = k (1 - stability) phi r^3   =>   dr/dt = k/3 (1 - stability) phi r
//
// so a dispersion coarsens exponentially at first and then creams out as the
// growing radius drives the Stokes velocity up (v ~ r^2). While the shake is
// still running the Weber criterion caps the radius: turbulent shear breaks up
// anything bigger, which is why an emulsion only coarsens once you stop
// shaking. Parcel volume never changes here, so total volume is conserved by
// construction rather than by bookkeeping.
void coarsenDroplets(Simulation& sim, double dt) {
  if (sim.droplets.empty() || dt <= 0.0) return;

  const double total = totalVolumeMl(sim);
  if (total <= 0.0) return;

  // Dispersed volume fraction per phase: the collision partner density.
  std::vector<double> dispersedMl(sim.phases.size(), 0.0);
  for (const Droplet& d : sim.droplets) {
    if (d.phase < 0 || static_cast<size_t>(d.phase) >= dispersedMl.size()) continue;
    dispersedMl[static_cast<size_t>(d.phase)] += double(d.parcelMl);
  }

  for (Droplet& d : sim.droplets) {
    if (d.phase < 0 || static_cast<size_t>(d.phase) >= sim.phases.size()) continue;
    const size_t phase = static_cast<size_t>(d.phase);
    const double stability = std::clamp(sim.phases[phase].emulsionStability, 0.0, 1.0);
    const double phi = std::clamp(dispersedMl[phase] / total, 0.0, 1.0);
    const double growth = kCoalesceRate * (1.0 - stability) * phi / 3.0;
    double radius = double(d.radius) * std::exp(growth * dt);

    double radiusCap = kMaxDropletRadius;
    if (sim.shake.active) {
      radiusCap = std::clamp(weberRadius(sim, phase, sim.shake.peakVelocity),
                             kMinDropletRadius, kMaxDropletRadius);
    }
    // A droplet can never be bigger than the parcel that carries it.
    radiusCap = std::min(radiusCap, parcelCloudRadius(d.parcelMl));
    d.radius = static_cast<float>(std::clamp(radius, kMinDropletRadius,
                                             std::max(radiusCap, kMinDropletRadius)));
  }
}

void disperseDuringShake(Simulation& sim, double dt, std::mt19937& rng) {
  if (!sim.shake.active || dt <= 0.0 || sim.phases.empty()) return;

  const ColumnGeometry geo = columnGeometry(sim);
  const double u = sim.shake.peakVelocity;
  const double turnoverHz = u / std::max(geo.fillHeight, 0.02);
  const double fraction = 1.0 - std::exp(-turnoverHz * dt);
  const double targetParcelMl = std::max(totalVolumeMl(sim) / kTargetParcels, 1e-6);
  std::uniform_real_distribution<double> unit(0.0, 1.0);

  for (size_t i = 0; i < sim.phases.size(); ++i) {
    const double want = sim.settledMl[i] * fraction;
    int capacityLeft = kMaxDroplets - static_cast<int>(sim.droplets.size());
    if (want <= 0.0 || capacityLeft <= 0) continue;

    const double expected = std::min(want / targetParcelMl, double(capacityLeft));
    int count = static_cast<int>(expected);
    if (unit(rng) < expected - count) ++count;
    count = std::min(count, capacityLeft);
    if (count <= 0) continue;

    const float physicalRadius = static_cast<float>(
        physicalDropletRadius(sim, i, sim.shake.specificPower, u));
    for (int parcel = 0; parcel < count && sim.settledMl[i] > 0.0; ++parcel) {
      float parcelMl = static_cast<float>(std::min(targetParcelMl, sim.settledMl[i]));
      if (double(parcelMl) > sim.settledMl[i]) parcelMl = std::nextafter(parcelMl, 0.0f);
      if (parcelMl <= 0.0f || !parcelFitsColumn(sim, geo, parcelMl)) break;
      sim.settledMl[i] -= double(parcelMl);

      Droplet d;
      d.phase = static_cast<int>(i);
      d.radius = physicalRadius;
      d.parcelMl = parcelMl;
      // Seed inside the wall AT THAT HEIGHT, not across the widest section:
      // in a tapering funnel most of the column is far narrower than
      // maxRadius, so a uniform-x seed put the majority of parcels outside the
      // glass, where the wall clamp then stacked them all on the same two
      // vertical lines. That pile-up read as vertical stripes on screen.
      const double y = unit(rng) * geo.fillHeight;
      const double cloud = parcelCloudRadius(double(parcelMl));
      const double xLimit = std::max(0.0, halfWidthAt(sim, geo, y) - cloud);
      d.position.y = static_cast<float>(y);
      d.position.x = static_cast<float>((unit(rng) * 2.0 - 1.0) * xLimit);
      clampParcelToVessel(sim, geo, d);
      sim.droplets.push_back(d);
    }
  }
}

// A streamfunction makes the large-scale slosh divergence-free by
// construction: psi=(U L/pi) cos(pi x/2W(y)) sin(pi y/Hfill) sin(2 pi f t).
// Its curl is one vessel-filling recirculation cell. W'(y) is taken from the
// same analytic wall profile, so the flow follows a narrowing pear rather
// than creating the independent vertical stripes produced by per-drop noise.
core::Vec2 coherentVelocity(const Simulation& sim, const ColumnGeometry& geo,
                            const Droplet& d) {
  if (sim.shakeEnergy <= 1e-4 || geo.fillHeight <= 1e-9) return {};
  const double y = std::clamp(double(d.position.y), 0.0, geo.fillHeight);
  const double width = std::max(halfWidthAt(sim, geo, y), 1e-6);
  const double x = std::clamp(double(d.position.x), -width, width);
  const double t = y / std::max(geo.height, 1e-12);
  constexpr double kDerivativeStep = 1e-4;
  const double ta = std::max(0.0, t - kDerivativeStep);
  const double tb = std::min(1.0, t + kDerivativeStep);
  const double dwDt = (vesselWidthAt(sim.vessel, tb) - vesselWidthAt(sim.vessel, ta)) /
                      std::max(tb - ta, 1e-12);
  const double widthPrime = geo.maxRadius * dwDt / std::max(geo.height, 1e-12);

  const double q = kPi * x / (2.0 * width);
  const double py = kPi * y / geo.fillHeight;
  const double temporal = std::sin(2.0 * kPi * sim.shake.frequencyHz * sim.elapsed);
  const double peak = sim.shake.peakVelocity * sim.shakeEnergy;
  const double length = std::min(geo.fillHeight, 2.0 * geo.maxRadius);
  const double psiScale = peak * length / kPi * temporal;
  const double vx = psiScale *
      (std::sin(q) * q * widthPrime / width * std::sin(py) +
       std::cos(q) * std::cos(py) * kPi / geo.fillHeight);
  const double vy = psiScale * std::sin(q) * kPi / (2.0 * width) * std::sin(py);
  return {static_cast<float>(vx), static_cast<float>(vy)};
}

void integrateSubstep(Simulation& sim, double dt, std::mt19937& rng) {
  if (sim.phases.empty()) return;
  if (sim.shake.active) {
    sim.shake.remainingS -= dt;
    sim.shakeEnergy = 1.0;
    if (sim.shake.remainingS <= 0.0) {
      sim.shake.remainingS = 0.0;
      sim.shake.active = false;
    }
  } else {
    sim.shakeEnergy *= std::exp(-kShakeDecayRate * dt);
  }

  disperseDuringShake(sim, dt, rng);
  const ColumnGeometry geo = columnGeometry(sim);
  const std::vector<double> boundaries = layerBoundaries(sim, geo);
  const double blend = 1.0 - std::exp(-dt / kDropletRelaxTau);
  std::normal_distribution<double> normal(0.0, 1.0);
  std::uniform_real_distribution<double> unit(0.0, 1.0);

  for (size_t di = 0; di < sim.droplets.size();) {
    Droplet& d = sim.droplets[di];
    // The medium a parcel travels through. Inside the settled stack that is
    // the resident bulk phase; ABOVE it, in the dispersed zone, it is the
    // emulsion of the OTHER phases. Reading the dispersed zone as "the
    // parcel's own phase" made every droplet neutrally buoyant in itself, so
    // a fully emulsified column had zero driving force and simply froze.
    const int resident = continuousPhaseAt(boundaries, d.position.y);
    const bool inSettledStack = double(d.position.y) < boundaries.back();
    const bool ownMedium = resident == d.phase;
    const double rhoDrop = sim.phases[static_cast<size_t>(d.phase)].density * 1000.0;
    double rhoContinuous = 0.0;
    double muContinuous = 0.0;
    if (inSettledStack && !ownMedium) {
      rhoContinuous = sim.phases[static_cast<size_t>(resident)].density * 1000.0;
      muContinuous = std::max(1e-4, sim.phases[static_cast<size_t>(resident)].viscosity) * 1e-3;
    } else {
      rhoContinuous = continuousDensity(sim, static_cast<size_t>(d.phase));
      muContinuous = std::max(1e-4, continuousViscosity(sim, static_cast<size_t>(d.phase))) * 1e-3;
    }
    const double radius = double(d.radius);

    // Stokes terminal speed for a spherical physical droplet. Parcel/cloud
    // size never enters this equation, so renderer density cannot accelerate
    // separation by orders of magnitude.
    const double terminal =
        -2.0 * radius * radius * (rhoDrop - rhoContinuous) * kGravity /
        (9.0 * muContinuous);
    const core::Vec2 coherent = coherentVelocity(sim, geo, d);
    const double turbulent = sim.shakeEnergy > 1e-4
                                 ? 0.08 * std::sqrt(std::max(0.0,
                                       sim.shake.specificPower * radius))
                                 : 0.0;
    const double targetX = double(coherent.x) + turbulent * normal(rng);
    const double targetY = terminal + double(coherent.y) + turbulent * normal(rng);
    d.velocity.x = static_cast<float>(double(d.velocity.x) +
                                      (targetX - double(d.velocity.x)) * blend);
    d.velocity.y = static_cast<float>(double(d.velocity.y) +
                                      (targetY - double(d.velocity.y)) * blend);
    d.position.x = static_cast<float>(double(d.position.x) + double(d.velocity.x) * dt);
    d.position.y = static_cast<float>(double(d.position.y) + double(d.velocity.y) * dt);
    clampParcelToVessel(sim, geo, d);

    const double cloud = parcelCloudRadius(d.parcelMl);
    const double ownLow = boundaries[static_cast<size_t>(d.phase)];
    const double ownHigh = boundaries[static_cast<size_t>(d.phase) + 1];
    const bool inOwnBand = double(d.position.y) - cloud <= ownHigh &&
                           double(d.position.y) + cloud >= ownLow;
    // A parcel that has run out of column in its buoyancy direction has
    // reached its own bulk even when that bulk is still a zero-thickness band:
    // at full dispersion every settled volume is 0, so every layer boundary
    // sits at 0 and a band test alone can never fire for the LIGHT phase --
    // its notional layer is pinned to the bottom of an empty stack while its
    // droplets float at the surface. That is a deadlock: nothing settles, so
    // no band ever grows. Floor contact for a sinker and free-surface contact
    // for a floater break it, and they are also what physically happens --
    // coalescence with the film already collecting there.
    const bool sinks = terminal < 0.0;
    const bool atTravelLimit = sinks
                                   ? double(d.position.y) <= cloud + 1e-9
                                   : double(d.position.y) >= geo.fillHeight - cloud - 1e-9;
    if (!sim.shake.active && (inOwnBand || atTravelLimit)) {
      const double stability = sim.phases[static_cast<size_t>(d.phase)].emulsionStability;
      const double rate = kBulkCoalesceRate * (1.0 - stability);
      const double probability = rate > 0.0 ? 1.0 - std::exp(-rate * dt) : 0.0;
      if (unit(rng) < probability) {
        sim.settledMl[static_cast<size_t>(d.phase)] += double(d.parcelMl);
        sim.droplets[di] = sim.droplets.back();
        sim.droplets.pop_back();
        continue;
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
  for (size_t i = 0; i < sim.phases.size(); ++i) sim.settledMl[i] = sim.phases[i].volumeMl;
  sim.droplets.clear();
  sim.shake = ShakeState{};
  sim.shakeEnergy = 0.0;
  sim.elapsed = 0.0;
}

void shake(Simulation& sim, const ShakeParams& params) {
  if (sim.phases.empty()) return;
  ShakeState state;
  state.durationS = std::clamp(params.durationS, 0.1, 120.0);
  state.remainingS = state.durationS;
  state.frequencyHz = std::clamp(params.frequencyHz, 0.2, 12.0);
  state.amplitudeM = std::clamp(params.amplitudeM, 0.005, 0.30);
  state.active = true;
  state.peakVelocity = 2.0 * kPi * state.frequencyHz * state.amplitudeM;
  state.specificPower = 0.5 * state.peakVelocity * state.peakVelocity * state.frequencyHz;

  double weighted = 0.0;
  double total = 0.0;
  for (size_t i = 0; i < sim.phases.size(); ++i) {
    const double volume = std::max(sim.phases[i].volumeMl, 0.0);
    weighted += physicalDropletRadius(sim, i, state.specificPower, state.peakVelocity) * volume;
    total += volume;
  }
  state.sauterRadiusM = total > 0.0 ? weighted / total : 0.0;
  sim.shake = state;
  sim.shakeEnergy = 1.0;
}

void step(Simulation& sim, double dt) {
  dt = std::clamp(dt, 0.0, 0.1);
  if (dt <= 0.0) return;
  std::mt19937 rng(sim.seed);
  const int substeps = std::max(1, static_cast<int>(std::lround(dt / kSubstepDt)));
  const double subDt = dt / substeps;
  for (int i = 0; i < substeps; ++i) {
    integrateSubstep(sim, subDt, rng);
    sim.elapsed += subDt;
  }
  coarsenDroplets(sim, dt);
  sim.seed = rng();
}

double totalVolumeMl(const Simulation& sim) {
  double total = 0.0;
  for (const Phase& phase : sim.phases) total += phase.volumeMl;
  return total;
}

double emulsifiedFraction(const Simulation& sim) {
  const double total = totalVolumeMl(sim);
  if (total <= 0.0) return 0.0;
  double dispersed = 0.0;
  for (const Droplet& d : sim.droplets) dispersed += double(d.parcelMl);
  return std::clamp(dispersed / total, 0.0, 1.0);
}

}  // namespace chemcad::sol

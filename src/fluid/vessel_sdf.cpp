#include "fluid/vessel_sdf.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <vector>

#include "fluid/kernels.hpp"

namespace chemcad::fluid {
namespace {

constexpr std::size_t kProfileIntervals = 1024;
constexpr std::size_t kBoundaryIntervals = 4096;
constexpr int kGaussOrder = 64;

struct GaussLegendre64 {
  std::array<double, kGaussOrder / 2> roots{};
  std::array<double, kGaussOrder / 2> weights{};
};

GaussLegendre64 makeGaussLegendre64() {
  GaussLegendre64 rule;
  for (int i = 0; i < kGaussOrder / 2; ++i) {
    double root = std::cos(kPi * (static_cast<double>(i) + 0.75) /
                           (static_cast<double>(kGaussOrder) + 0.5));
    double derivative = 0.0;
    for (int iteration = 0; iteration < 16; ++iteration) {
      double previous = 1.0;
      double current = root;
      for (int degree = 2; degree <= kGaussOrder; ++degree) {
        const double next = ((2.0 * degree - 1.0) * root * current -
                             (degree - 1.0) * previous) /
                            degree;
        previous = current;
        current = next;
      }
      derivative = kGaussOrder * (root * current - previous) /
                   (root * root - 1.0);
      const double correction = current / derivative;
      root -= correction;
      if (std::abs(correction) <= 2.0 * std::numeric_limits<double>::epsilon()) break;
    }
    rule.roots[i] = root;
    rule.weights[i] = 2.0 / ((1.0 - root * root) * derivative * derivative);
  }
  return rule;
}

// Integrating the radially symmetric Wendland kernel over a plane a distance
// |u| from its centre reduces to 2*pi*integral(W(r)*r dr). Expanding the C2
// polynomial gives this closed form, leaving only the wall-normal integral for
// the 64-point Gauss-Legendre rule below.
double wendlandPlaneIntegral(double u, double support) {
  const double a = std::abs(u) / support;
  if (a >= 1.0) return 0.0;
  const double a2 = a * a;
  const double a4 = a2 * a2;
  const double a5 = a4 * a;
  const double a6 = a5 * a;
  const double a7 = a6 * a;
  const double primitive = 0.5 * a2 - 2.5 * a4 + 4.0 * a5 - 2.5 * a6 +
                           (4.0 / 7.0) * a7;
  return std::max(0.0, 21.0 * (1.0 / 14.0 - primitive) / support);
}

double integrateGlassHalfSpace(double lower, double support,
                               const GaussLegendre64& rule) {
  if (lower >= support) return 0.0;
  lower = std::max(lower, -support);
  const double midpoint = 0.5 * (lower + support);
  const double halfWidth = 0.5 * (support - lower);
  double sum = 0.0;
  for (int i = 0; i < kGaussOrder / 2; ++i) {
    const double offset = halfWidth * rule.roots[i];
    sum += rule.weights[i] *
           (wendlandPlaneIntegral(midpoint - offset, support) +
            wendlandPlaneIntegral(midpoint + offset, support));
  }
  return halfWidth * sum;
}

}  // namespace

void VesselBoundary::build(sol::Vessel vessel, double heightM, double support,
                           double spacing) {
  vessel_ = vessel;
  heightM_ = std::max(0.0, heightM);
  support_ = std::max(0.0, support);
  spacing_ = std::max(0.0, spacing);

  profileZ_.resize(kProfileIntervals + 1);
  profileR_.resize(kProfileIntervals + 1);
  cumulativeVol_.assign(kProfileIntervals + 1, 0.0);

  // vesselOutline is the public, physically-scaled representation of the same
  // analytic profile. Its first point recovers the vessel's calibrated
  // widest-radius/height aspect ratio without duplicating private constants.
  const std::vector<core::Vec2> outline = sol::vesselOutline(vessel_, heightM_);
  const double widthAtBottom = sol::vesselWidthAt(vessel_, 0.0);
  const double radialScale = !outline.empty() && widthAtBottom > 0.0
                                 ? std::abs(static_cast<double>(outline.front().x)) /
                                       widthAtBottom
                                 : 0.0;

  maxRadiusM_ = 0.0;
  for (std::size_t i = 0; i <= kProfileIntervals; ++i) {
    const double t = static_cast<double>(i) / kProfileIntervals;
    profileZ_[i] = t * heightM_;
    profileR_[i] = radialScale * sol::vesselWidthAt(vessel_, t);
    maxRadiusM_ = std::max(maxRadiusM_, profileR_[i]);
  }

  // Per-panel Simpson integration uses the analytic midpoint radius, not a
  // midpoint reconstructed from the endpoints. Thus the volume table and the
  // rendered analytic vessel share the same profile rather than a trapezoid.
  for (std::size_t i = 1; i <= kProfileIntervals; ++i) {
    const double tMid = (static_cast<double>(i) - 0.5) / kProfileIntervals;
    const double rMid = radialScale * sol::vesselWidthAt(vessel_, tMid);
    const double dz = profileZ_[i] - profileZ_[i - 1];
    const double panel = kPi * dz / 6.0 *
                         (profileR_[i - 1] * profileR_[i - 1] +
                          4.0 * rMid * rMid + profileR_[i] * profileR_[i]);
    cumulativeVol_[i] = cumulativeVol_[i - 1] + panel;
  }

  boundaryDelta_.assign(kBoundaryIntervals + 1, 0.0);
  boundarySlope_.assign(kBoundaryIntervals + 1, 0.0);
  if (!(support_ > 0.0) || !(spacing_ > 0.0)) return;

  const GaussLegendre64 quadrature = makeGaussLegendre64();
  // The continuum half-space integral is a fraction in [0,1]. Scale it by
  // the solver's discrete rest-lattice delta0, rather than by 1/dx^3, so a
  // planar half-filled kernel is corrected to the exact PCISPH density target.
  const double restDelta = restNumberDensity(spacing_, support_);
  for (std::size_t i = 0; i <= kBoundaryIntervals; ++i) {
    const double phi = -support_ + 2.0 * support_ *
                                       static_cast<double>(i) / kBoundaryIntervals;
    // Plane-cut boundary handling replaces the locally curved wall by its
    // tangent half-space. This is the standard kernel-deficiency correction;
    // its neglected terms are O(H/R_curvature), so it is valid where the local
    // glass curvature radius is much larger than the kernel support H.
    boundaryDelta_[i] =
        integrateGlassHalfSpace(-phi, support_, quadrature) * restDelta;
    // Leibniz' rule gives d/dphi integral[-phi,H] A(u)du = A(-phi).
    // Storing that analytic derivative avoids differencing a steep table near
    // contact while remaining exactly the derivative of the same plane cut.
    boundarySlope_[i] = (i == 0 || i == kBoundaryIntervals)
                            ? 0.0
                            : wendlandPlaneIntegral(-phi, support_) * restDelta;
  }
  // Compact support makes both limiting cases exact, independent of the
  // quadrature round-off in the final table entry.
  boundaryDelta_.front() = 0.0;
  boundaryDelta_.back() = restDelta;
}

double VesselBoundary::radiusAt(double z) const {
  if (profileR_.empty()) return 0.0;
  if (z <= 0.0 || !(heightM_ > 0.0)) return profileR_.front();
  if (z >= heightM_) return profileR_.back();
  const double scaled = z * kProfileIntervals / heightM_;
  const std::size_t lower = static_cast<std::size_t>(scaled);
  const double fraction = scaled - static_cast<double>(lower);
  return profileR_[lower] + fraction * (profileR_[lower + 1] - profileR_[lower]);
}

SurfaceQuery VesselBoundary::query(double x, double y, double z) const {
  SurfaceQuery result;
  if (profileZ_.size() < 2 || profileR_.size() != profileZ_.size()) return result;

  const double radial = std::sqrt(x * x + y * y);
  double bestSquared = std::numeric_limits<double>::infinity();
  double bestNormalS = 0.0;
  double bestNormalZ = 1.0;

  const auto considerSegment = [&](double s0, double z0, double s1, double z1,
                                   double normalS, double normalZ) {
    const double ds = s1 - s0;
    const double dz = z1 - z0;
    const double lengthSquared = ds * ds + dz * dz;
    double parameter = 0.0;
    if (lengthSquared > 0.0) {
      parameter = std::clamp(((radial - s0) * ds + (z - z0) * dz) /
                                 lengthSquared,
                             0.0, 1.0);
    }
    const double nearestS = s0 + parameter * ds;
    const double nearestZ = z0 + parameter * dz;
    const double errorS = radial - nearestS;
    const double errorZ = z - nearestZ;
    const double distanceSquared = errorS * errorS + errorZ * errorZ;
    // Segments are presented in a fixed geometric-index order. Strict less
    // therefore implements the required lowest-index tie break exactly.
    if (distanceSquared < bestSquared) {
      bestSquared = distanceSquared;
      const double normalLength = std::sqrt(normalS * normalS + normalZ * normalZ);
      bestNormalS = normalLength > 0.0 ? normalS / normalLength : 0.0;
      bestNormalZ = normalLength > 0.0 ? normalZ / normalLength : 1.0;
    }
  };

  // Index 0 is the bottom cap, followed by ascending wall segments and finally
  // the top cap. Normals point from glass into the fluid as SurfaceQuery
  // specifies: +z at the bottom, inward radially at the wall, -z at the top.
  considerSegment(0.0, 0.0, profileR_.front(), 0.0, 0.0, 1.0);
  for (std::size_t i = 0; i + 1 < profileZ_.size(); ++i) {
    const double ds = profileR_[i + 1] - profileR_[i];
    const double dz = profileZ_[i + 1] - profileZ_[i];
    considerSegment(profileR_[i], profileZ_[i], profileR_[i + 1],
                    profileZ_[i + 1], -dz, ds);
  }
  considerSegment(0.0, heightM_, profileR_.back(), heightM_, 0.0, -1.0);

  const bool inside = z >= 0.0 && z <= heightM_ && radial <= radiusAt(z);
  const double unsignedDistance = std::sqrt(std::max(0.0, bestSquared));
  result.distance = inside ? -unsignedDistance : unsignedDistance;
  const double radialX = radial >= 1.0e-12 ? x / radial : 1.0;
  const double radialY = radial >= 1.0e-12 ? y / radial : 0.0;
  result.nx = bestNormalS * radialX;
  result.ny = bestNormalS * radialY;
  result.nz = bestNormalZ;
  return result;
}

double VesselBoundary::boundaryDensity(double phi) const {
  if (boundaryDelta_.empty() || !(support_ > 0.0)) return 0.0;
  if (phi <= -support_) return 0.0;
  if (phi >= support_) return boundaryDelta_.back();
  const double scaled = (phi + support_) * kBoundaryIntervals / (2.0 * support_);
  const std::size_t lower = static_cast<std::size_t>(scaled);
  const double fraction = scaled - static_cast<double>(lower);
  return boundaryDelta_[lower] +
         fraction * (boundaryDelta_[lower + 1] - boundaryDelta_[lower]);
}

double VesselBoundary::boundaryDensitySlope(double phi) const {
  if (boundarySlope_.empty() || !(support_ > 0.0) || phi <= -support_ ||
      phi >= support_) {
    return 0.0;
  }
  const double scaled = (phi + support_) * kBoundaryIntervals / (2.0 * support_);
  const std::size_t lower = static_cast<std::size_t>(scaled);
  const double fraction = scaled - static_cast<double>(lower);
  return boundarySlope_[lower] +
         fraction * (boundarySlope_[lower + 1] - boundarySlope_[lower]);
}

double VesselBoundary::volumeBelow(double z) const {
  if (cumulativeVol_.empty() || z <= 0.0 || !(heightM_ > 0.0)) return 0.0;
  if (z >= heightM_) return cumulativeVol_.back();
  const double scaled = z * kProfileIntervals / heightM_;
  const std::size_t lower = static_cast<std::size_t>(scaled);
  const double fraction = scaled - static_cast<double>(lower);
  return cumulativeVol_[lower] +
         fraction * (cumulativeVol_[lower + 1] - cumulativeVol_[lower]);
}

double VesselBoundary::heightForVolume(double volume) const {
  if (cumulativeVol_.empty() || volume <= 0.0) return 0.0;
  if (volume >= cumulativeVol_.back()) return heightM_;
  const auto upper = std::upper_bound(cumulativeVol_.begin(), cumulativeVol_.end(), volume);
  const std::size_t hi = static_cast<std::size_t>(upper - cumulativeVol_.begin());
  const std::size_t lo = hi - 1;
  const double span = cumulativeVol_[hi] - cumulativeVol_[lo];
  const double fraction = span > 0.0 ? (volume - cumulativeVol_[lo]) / span : 0.0;
  return heightM_ * (static_cast<double>(lo) + fraction) / kProfileIntervals;
}

std::size_t VesselBoundary::chargeLattice(const std::vector<PhaseMaterial>& phases,
                                          double spacing,
                                          Particles& particles) const {
  particles.clear();
  if (!(spacing > 0.0) || phases.empty() || !(heightM_ > 0.0)) return 0;

  struct Charge {
    std::size_t phase = 0;
    std::size_t requested = 0;
  };
  std::vector<Charge> charges;
  charges.reserve(phases.size());
  const double siteVolume = spacing * spacing * spacing;
  std::size_t requested = 0;
  for (std::size_t phase = 0; phase < phases.size(); ++phase) {
    const double volume = std::max(0.0, phases[phase].volumeMl) * 1.0e-6;
    const double exactSites = volume / siteVolume;
    // A sub-half-site positive charge still represents a real phase. Giving it
    // one site prevents volume rounding from silently deleting that phase.
    const std::size_t sites =
        volume > 0.0
            ? std::max<std::size_t>(
                  1, static_cast<std::size_t>(std::llround(exactSites)))
            : 0;
    charges.push_back(Charge{phase, sites});
    requested += sites;
  }
  std::stable_sort(charges.begin(), charges.end(), [&](const Charge& a, const Charge& b) {
    if (phases[a.phase].restDensity != phases[b.phase].restDensity) {
      return phases[a.phase].restDensity > phases[b.phase].restDensity;
    }
    return a.phase < b.phase;
  });
  particles.reserve(requested);

  struct Site {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
  };
  std::vector<Site> sites;
  sites.reserve(requested);

  // The solver projects particles to phi=-0.35*dx (SolverConfig's contact
  // radius), so charging at that same margin is safe without discarding the
  // only columns in a narrow stem. Keep a simple cubic lattice, with number
  // density 1/dx^3: BCC or HCP at nearest-neighbour distance dx would increase
  // that density and invalidate the solver's cubic rest-number-density target.
  constexpr double kContactRadiusFactor = 0.35;
  const double clearance = kContactRadiusFactor * spacing;
  const int radialExtent = static_cast<int>(std::ceil(maxRadiusM_ / spacing));
  const int verticalCount = static_cast<int>(std::floor(heightM_ / spacing));
  for (int iz = 0; iz < verticalCount && sites.size() < requested; ++iz) {
    const double z = (static_cast<double>(iz) + 0.5) * spacing;
    for (int iy = -radialExtent; iy <= radialExtent && sites.size() < requested;
         ++iy) {
      const double y = static_cast<double>(iy) * spacing;
      for (int ix = -radialExtent; ix <= radialExtent && sites.size() < requested;
           ++ix) {
        const double x = static_cast<double>(ix) * spacing;
        if (query(x, y, z).distance > -clearance) continue;
        sites.push_back(Site{x, y, z});
      }
    }
  }

  std::size_t positiveCharges = static_cast<std::size_t>(std::count_if(
      charges.begin(), charges.end(),
      [](const Charge& charge) { return charge.requested > 0; }));
  std::size_t nextSite = 0;
  for (const Charge& charge : charges) {
    if (charge.requested == 0) continue;
    const std::size_t available = sites.size() - nextSite;
    if (available == 0) break;

    // Reserve one site for each later positive charge whenever capacity permits.
    // This preserves bottom-up density ordering without allowing an early phase
    // shortfall to consume every site and silently erase all lighter phases.
    const std::size_t allocated =
        available >= positiveCharges
            ? std::min(charge.requested, available - (positiveCharges - 1))
            : 1;
    for (std::size_t i = 0; i < allocated; ++i) {
      const Site& site = sites[nextSite++];
      particles.add(static_cast<float>(site.x), static_cast<float>(site.y),
                    static_cast<float>(site.z),
                    static_cast<uint8_t>(charge.phase));
    }
    --positiveCharges;
  }

  // The return is the actual allocation. Callers can compare it with the sum of
  // rounded volume/dx^3 requests above to observe a vessel-capacity shortfall.
  return particles.size();
}

}  // namespace chemcad::fluid

#include "fluid/solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "fluid/kernels.hpp"

namespace chemcad::fluid {
namespace {

constexpr double kTiny = 1.0e-12;

struct Vec3d {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

Vec3d operator+(const Vec3d& a, const Vec3d& b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}
Vec3d operator-(const Vec3d& a, const Vec3d& b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3d operator*(double s, const Vec3d& a) { return {s * a.x, s * a.y, s * a.z}; }
double dot(const Vec3d& a, const Vec3d& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
double lengthSquared(const Vec3d& a) { return dot(a, a); }
double length(const Vec3d& a) { return std::sqrt(lengthSquared(a)); }
Vec3d cross(const Vec3d& a, const Vec3d& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
          a.x * b.y - a.y * b.x};
}
Vec3d normalized(const Vec3d& a) {
  const double n = length(a);
  return n > kTiny ? (1.0 / n) * a : Vec3d{};
}

Vec3d positionOf(const Particles& p, std::size_t i) {
  return {p.px[i], p.py[i], p.pz[i]};
}
Vec3d velocityOf(const Particles& p, std::size_t i) {
  return {p.vx[i], p.vy[i], p.vz[i]};
}

double phaseMass(const std::vector<double>& mass, const Particles& p, std::size_t i) {
  const std::size_t phase = p.phase[i];
  return phase < mass.size() ? mass[phase] : 0.0;
}

double phaseViscosity(const std::vector<PhaseMaterial>& phases, const Particles& p,
                      std::size_t i) {
  const std::size_t phase = p.phase[i];
  return phase < phases.size() ? phases[phase].dynamicViscosity : 0.0;
}

void resizeScratch(std::size_t n, std::vector<float>& x, std::vector<float>& y,
                   std::vector<float>& z, std::vector<float>& vx, std::vector<float>& vy,
                   std::vector<float>& vz, std::vector<float>& fx, std::vector<float>& fy,
                   std::vector<float>& fz, std::vector<double>& error) {
  x.resize(n);
  y.resize(n);
  z.resize(n);
  vx.resize(n);
  vy.resize(n);
  vz.resize(n);
  fx.resize(n);
  fy.resize(n);
  fz.resize(n);
  error.resize(n);
}

void computeNumberDensity(Particles& p, const NeighbourGrid& grid, const VesselBoundary* boundary,
                          double support) {
  for (std::size_t i = 0; i < p.size(); ++i) {
    const Vec3d pi = positionOf(p, i);
    double delta = 0.0;
    grid.forEachCandidate(i, [&](std::size_t j) {
      const Vec3d r = pi - positionOf(p, j);
      const double radius = length(r);
      if (radius < support) delta += wendlandW(radius, support);
    });
    if (boundary != nullptr) {
      delta += boundary->boundaryDensity(boundary->query(pi.x, pi.y, pi.z).distance);
    }
    p.delta[i] = static_cast<float>(delta);
  }
}

void projectContact(double contactRadius, double wallFriction, const VesselBoundary& boundary,
                    Vec3d& position, Vec3d& velocity) {
  const SurfaceQuery wall = boundary.query(position.x, position.y, position.z);
  if (wall.distance <= -contactRadius) return;

  // VesselBoundary::query returns the glass-to-fluid normal. Moving by
  // (phi+r_c)n is therefore an inward projection to phi=-r_c.
  const Vec3d n{wall.nx, wall.ny, wall.nz};
  position = position + (wall.distance + contactRadius) * n;

  // Restitution zero removes only velocity into the glass. The remaining
  // tangential component is reduced by Coulomb-like wall loss; motion away
  // from the glass is retained so projection cannot make particles stick.
  const double vn = dot(velocity, n);
  const Vec3d tangent = velocity - vn * n;
  velocity = std::max(0.0, vn) * n + (1.0 - std::clamp(wallFriction, 0.0, 1.0)) * tangent;
}

std::vector<double> calibratePressureStiffness(const std::vector<double>& mass, double spacing,
                                               double support, double dt) {
  std::vector<double> result(mass.size(), 0.0);
  if (spacing <= 0.0 || support <= 0.0 || dt <= 0.0) return result;

  std::vector<Vec3d> lattice;
  const int extent = static_cast<int>(std::ceil(support / spacing));
  lattice.push_back({});
  for (int z = -extent; z <= extent; ++z) {
    for (int y = -extent; y <= extent; ++y) {
      for (int x = -extent; x <= extent; ++x) {
        if (x == 0 && y == 0 && z == 0) continue;
        const Vec3d q{spacing * x, spacing * y, spacing * z};
        if (length(q) < support) lattice.push_back(q);
      }
    }
  }
  const double delta0 = fluid::restNumberDensity(spacing, support);

  for (std::size_t phase = 0; phase < mass.size(); ++phase) {
    if (mass[phase] <= 0.0 || delta0 <= 0.0) continue;
    std::vector<Vec3d> force(lattice.size());
    // A one-pascal pressure correction is applied at the centre. The same
    // corrected density-contrast operator used by the solve moves the centre's
    // neighbours outward; the resulting centre-density deficit is the measured
    // response. Solenthaler and Pajarola 2009, eqs. 7-10 derive this scaling for
    // equal masses. Measuring the response with this phase's actual mass is the
    // variable-mass adaptation required by density-contrast SPH.
    for (std::size_t i = 0; i < lattice.size(); ++i) {
      Vec3d fi{};
      for (std::size_t j = 0; j < lattice.size(); ++j) {
        if (i == j) continue;
        const double pi = i == 0 ? 1.0 : 0.0;
        const double pj = j == 0 ? 1.0 : 0.0;
        if (pi + pj == 0.0) continue;
        const Vec3d rij = lattice[i] - lattice[j];
        const double r = length(rij);
        if (r <= kTiny || r >= support) continue;
        const Vec3d grad = (wendlandGradMagnitude(r, support) / r) * rij;
        fi = fi + (-(pi + pj) / (2.0 * delta0 * delta0)) * grad;
      }
      force[i] = fi;
    }

    double predictedDelta = 0.0;
    const Vec3d centre = lattice[0] + (dt * dt / mass[phase]) * force[0];
    for (std::size_t j = 0; j < lattice.size(); ++j) {
      const Vec3d predicted = lattice[j] + (dt * dt / mass[phase]) * force[j];
      predictedDelta += wendlandW(length(centre - predicted), support);
    }
    const double relativeDensityError = (predictedDelta - delta0) / delta0;
    result[phase] = 1.0 / std::max(1.0e-12, -relativeDensityError);
  }
  return result;
}

void computeColourGeometry(Particles& p, const NeighbourGrid& grid, double support,
                           std::vector<double>& curvature) {
  const std::size_t n = p.size();
  for (std::size_t i = 0; i < n; ++i) {
    const Vec3d pi = positionOf(p, i);
    double numerator = 0.0;
    double denominator = 0.0;
    grid.forEachCandidate(i, [&](std::size_t j) {
      const double r = length(pi - positionOf(p, j));
      if (r >= support || p.delta[j] <= 0.0f) return;
      const double volumeWeight = wendlandW(r, support) / static_cast<double>(p.delta[j]);
      denominator += volumeWeight;
      numerator += (p.phase[j] == 0 ? 0.0 : 1.0) * volumeWeight;
    });
    p.colour[i] = static_cast<float>(denominator > kTiny ? numerator / denominator
                                                         : (p.phase[i] == 0 ? 0.0 : 1.0));
  }

  // Solenthaler and Pajarola, SCA 2008, eqs. 22-24: use the normalised colour
  // field, its SPH gradient as the interface normal, then the divergence of the
  // unit normal as curvature. The 1/delta_j factor is particle volume in the
  // number-density formulation.
  for (std::size_t i = 0; i < n; ++i) {
    const Vec3d pi = positionOf(p, i);
    Vec3d normal{};
    grid.forEachCandidate(i, [&](std::size_t j) {
      if (i == j || p.delta[j] <= 0.0f) return;
      const Vec3d rij = pi - positionOf(p, j);
      const double r = length(rij);
      if (r <= kTiny || r >= support) return;
      const Vec3d grad = (wendlandGradMagnitude(r, support) / r) * rij;
      normal = normal + ((static_cast<double>(p.colour[j]) - p.colour[i]) /
                         static_cast<double>(p.delta[j])) *
                            grad;
    });
    p.nx[i] = static_cast<float>(normal.x);
    p.ny[i] = static_cast<float>(normal.y);
    p.nz[i] = static_cast<float>(normal.z);
  }

  curvature.assign(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    if (!(p.colour[i] > 0.05f && p.colour[i] < 0.95f)) continue;
    const Vec3d pi = positionOf(p, i);
    const Vec3d ni = normalized({p.nx[i], p.ny[i], p.nz[i]});
    double kappa = 0.0;
    grid.forEachCandidate(i, [&](std::size_t j) {
      if (i == j || p.delta[j] <= 0.0f) return;
      const Vec3d rij = pi - positionOf(p, j);
      const double r = length(rij);
      if (r <= kTiny || r >= support) return;
      const Vec3d grad = (wendlandGradMagnitude(r, support) / r) * rij;
      const Vec3d nj = normalized({p.nx[j], p.ny[j], p.nz[j]});
      kappa -= dot(nj - ni, grad) / static_cast<double>(p.delta[j]);
    });
    curvature[i] = kappa;
  }
}

void addSurfaceAcceleration(Particles& p, const NeighbourGrid& grid,
                            const std::vector<PhaseMaterial>& phases,
                            const std::vector<double>& mass, const InterfaceModel& interfaceModel,
                            double support, std::vector<double>& curvature) {
  if (!interfaceModel.calibrated || interfaceModel.cohesionGain <= 0.0 || phases.size() < 2) {
    std::fill(p.colour.begin(), p.colour.end(), 0.0f);
    std::fill(p.nx.begin(), p.nx.end(), 0.0f);
    std::fill(p.ny.begin(), p.ny.end(), 0.0f);
    std::fill(p.nz.begin(), p.nz.end(), 0.0f);
    return;
  }

  computeColourGeometry(p, grid, support, curvature);
  const std::size_t phaseCount = phases.size();
  for (std::size_t i = 0; i < p.size(); ++i) {
    const std::size_t phaseI = p.phase[i];
    if (phaseI >= phaseCount || p.delta[i] <= 0.0f) continue;
    const Vec3d pi = positionOf(p, i);
    const Vec3d ni = normalized({p.nx[i], p.ny[i], p.nz[i]});
    Vec3d acceleration{};
    grid.forEachCandidate(i, [&](std::size_t j) {
      if (i == j || p.phase[j] >= phaseCount || p.phase[j] == phaseI || p.delta[j] <= 0.0f)
        return;
      const std::size_t phaseJ = p.phase[j];
      const double sigma = interfaceModel.sigma[phaseI * phaseCount + phaseJ];
      // Akinci et al. 2013 explicitly note that their cohesion coefficient is
      // resolution dependent. It therefore enters only as sigma times the
      // Young-Laplace calibration gain, never as the measured N/m directly.
      const double coefficient = sigma * interfaceModel.cohesionGain;
      if (coefficient <= 0.0) return;
      const Vec3d rij = pi - positionOf(p, j);
      const double r = length(rij);
      if (r <= kTiny || r >= support) return;

      const double rhoi = mass[phaseI] * static_cast<double>(p.delta[i]);
      const double rhoj = mass[phaseJ] * static_cast<double>(p.delta[j]);
      const double restI = phases[phaseI].restDensity;
      const double restJ = phases[phaseJ].restDensity;
      const double densityCorrection =
          2.0 * std::sqrt(restI * restJ) / std::max(kTiny, rhoi + rhoj);
      const Vec3d nj = normalized({p.nx[j], p.ny[j], p.nz[j]});
      Vec3d pairShape = -cohesionC(r, support) * (1.0 / r) * rij;
      if ((p.colour[i] > 0.05f && p.colour[i] < 0.95f) ||
          (p.colour[j] > 0.05f && p.colour[j] < 0.95f)) {
        const Vec3d curvatureI = curvature[i] * ni;
        const Vec3d curvatureJ = curvature[j] * nj;
        const double volumeKernel =
            wendlandW(r, support) /
            std::sqrt(static_cast<double>(p.delta[i]) * static_cast<double>(p.delta[j]));
        pairShape = pairShape - volumeKernel * (curvatureI - curvatureJ);
      }

      // F_ij = coefficient*m_i*m_j*K_ij*pairShape. Every scalar is symmetric
      // and pairShape changes sign under i<->j, hence F_ji=-F_ij exactly at the
      // algebraic level. Gathering F_ij/m_i preserves momentum without a
      // schedule-dependent scatter.
      acceleration = acceleration +
                     (coefficient * mass[phaseJ] * densityCorrection) * pairShape;
    });
    p.ax[i] = static_cast<float>(static_cast<double>(p.ax[i]) + acceleration.x);
    p.ay[i] = static_cast<float>(static_cast<double>(p.ay[i]) + acceleration.y);
    p.az[i] = static_cast<float>(static_cast<double>(p.az[i]) + acceleration.z);
  }
}

void computePressureForce(const Particles& p, const NeighbourGrid& grid, double support,
                          const std::vector<float>& x, const std::vector<float>& y,
                          const std::vector<float>& z, double delta0,
                          const std::vector<double>& error, std::vector<float>& fx,
                          std::vector<float>& fy, std::vector<float>& fz) {
  std::fill(fx.begin(), fx.end(), 0.0f);
  std::fill(fy.begin(), fy.end(), 0.0f);
  std::fill(fz.begin(), fz.end(), 0.0f);
  for (std::size_t i = 0; i < p.size(); ++i) {
    const Vec3d pi{x[i], y[i], z[i]};
    const double deltai = std::max(kTiny, delta0 * (1.0 + error[i]));
    grid.forEachCandidate(i, [&](std::size_t j) {
      if (j <= i) return;
      const Vec3d rij = pi - Vec3d{x[j], y[j], z[j]};
      const double r = length(rij);
      if (r <= kTiny || r >= support) return;
      const double deltaj = std::max(kTiny, delta0 * (1.0 + error[j]));
      const Vec3d grad = (wendlandGradMagnitude(r, support) / r) * rij;
      const double scale = -(static_cast<double>(p.pressure[i]) + p.pressure[j]) /
                           (2.0 * deltai * deltaj);
      const float pairX = static_cast<float>(scale * grad.x);
      const float pairY = static_cast<float>(scale * grad.y);
      const float pairZ = static_cast<float>(scale * grad.z);
      // Evaluate each Solenthaler-Pajarola 2008 eq. 14 pair once and apply its
      // exact float negation. This prevents two gathers from rounding the same
      // delta product or kernel gradient differently.
      fx[i] += pairX;
      fy[i] += pairY;
      fz[i] += pairZ;
      fx[j] -= pairX;
      fy[j] -= pairY;
      fz[j] -= pairZ;
    });
  }
}

}  // namespace

void Solver::configure(const SolverConfig& config) {
  if (!(config.resolution.spacing > 0.0) || !(config.maxSubstepS > 0.0) ||
      !(config.cflNumber > 0.0) || !(config.accelerationSafety > 0.0) ||
      config.minPressureIterations < 1 ||
      config.maxPressureIterations < config.minPressureIterations ||
      !(config.densityTolerance > 0.0) || !(config.maxSpeed > 0.0) ||
      config.contactRadiusFactor < 0.0 || config.wallFriction < 0.0 ||
      config.wallFriction > 1.0 || config.xsphSmoothing < 0.0) {
    throw std::invalid_argument("invalid fluid solver configuration");
  }
  const bool resolutionChanged = config_.resolution.spacing != config.resolution.spacing;
  config_ = config;
  delta0_ = fluid::restNumberDensity(config_.resolution.spacing, config_.resolution.support());
  stiffness_ = calibratePressureStiffness(mass_, config_.resolution.spacing,
                                          config_.resolution.support(), config_.maxSubstepS);
  stats_ = {};
  if (resolutionChanged) {
    interface_.cohesionGain = 0.0;
    interface_.calibrated = false;
  }
  calibrationError_.clear();
}

void Solver::setPhases(const std::vector<PhaseMaterial>& phases,
                       const std::vector<double>& sigmaPairs) {
  if (phases.empty() || phases.size() > 255) {
    throw std::invalid_argument("fluid solver requires between one and 255 phases");
  }
  phases_ = phases;
  mass_.resize(phases_.size());
  for (std::size_t i = 0; i < phases_.size(); ++i) {
    if (!(phases_[i].restDensity > 0.0) || phases_[i].dynamicViscosity < 0.0) {
      throw std::invalid_argument("phase density must be positive and viscosity non-negative");
    }
    mass_[i] = phases_[i].restDensity * config_.resolution.particleVolume();
  }

  const std::size_t n = phases_.size();
  interface_.sigma.assign(n * n, 0.0);
  if (sigmaPairs.size() == n * n) {
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        const double a = sigmaPairs[i * n + j];
        const double b = sigmaPairs[j * n + i];
        interface_.sigma[i * n + j] = i == j ? 0.0 : std::max(0.0, 0.5 * (a + b));
      }
    }
  } else if (sigmaPairs.size() == n * (n - 1) / 2) {
    std::size_t pair = 0;
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = i + 1; j < n; ++j) {
        const double value = std::max(0.0, sigmaPairs[pair++]);
        interface_.sigma[i * n + j] = value;
        interface_.sigma[j * n + i] = value;
      }
    }
  } else if (!sigmaPairs.empty()) {
    throw std::invalid_argument("interfacial tension table has the wrong size");
  }

  stiffness_ = calibratePressureStiffness(mass_, config_.resolution.spacing,
                                          config_.resolution.support(), config_.maxSubstepS);
  stats_ = {};
  interface_.cohesionGain = 0.0;
  interface_.calibrated = false;
  calibrationError_.clear();
}

int Solver::advance(Particles& particles, const VesselBoundary& boundary,
                    const VesselMotion& motion, double timeS, double dt) {
  double stiffnessDt = stats_.substepS;
  stats_ = {};
  if (dt <= 0.0 || particles.empty()) return 0;
  if (phases_.empty()) throw std::logic_error("setPhases must precede advance");
  for (uint8_t phase : particles.phase) {
    if (phase >= phases_.size()) {
      throw std::invalid_argument("particle phase index is outside the material table");
    }
  }

  const double support = config_.resolution.support();
  const double contactRadius = config_.contactRadiusFactor * config_.resolution.spacing;
  const double displacementLimit = 0.25 * support;
  const double omega = 2.0 * kPi * std::abs(motion.shakeFrequencyHz);
  const double peakShakeAcceleration =
      motion.shaking && motion.shakeRemainingS > 0.0
          ? std::abs(motion.shakeAmplitudeM) * omega * omega
          : 0.0;
  resizeScratch(particles.size(), predictedX_, predictedY_, predictedZ_, predictedVX_,
                predictedVY_, predictedVZ_, forceX_, forceY_, forceZ_, densityError_);

  double pressureRelaxationLimit = 1.0;
  double elapsed = 0.0;
  double lastSubstep = 0.0;
  while (elapsed < dt) {
    const double remaining = dt - elapsed;
    if (!(remaining > std::numeric_limits<double>::epsilon() * std::max(1.0, dt))) break;

    double maxVelocity = 0.0;
    double maxAcceleration = 0.0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
      maxVelocity = std::max(maxVelocity, length(velocityOf(particles, i)));
      maxAcceleration =
          std::max(maxAcceleration, length({particles.ax[i], particles.ay[i], particles.az[i]}));
    }
    const FrameAcceleration limitingFrame = frameAcceleration(motion, timeS + elapsed);
    maxAcceleration =
        std::max(maxAcceleration,
                 length({limitingFrame.uniform[0], limitingFrame.uniform[1],
                         limitingFrame.uniform[2]}) +
                     peakShakeAcceleration);
    // Include A*(2*pi*f)^2 rather than sampling only the zero-acceleration
    // phase at the start of a long call. The positive root of
    // v*dt+a*dt^2=0.25H reserves displacement for velocity and acceleration
    // together; applying independent bounds let their sum cross 0.25 H.
    const double cflStep =
        config_.cflNumber * support / std::max(maxVelocity, 1.0e-9);
    const double accelerationStep =
        config_.accelerationSafety * std::sqrt(support / std::max(maxAcceleration, 1.0e-9));
    const double transportStep =
        1.9 * displacementLimit /
        (maxVelocity +
         std::sqrt(maxVelocity * maxVelocity + 4.0 * maxAcceleration * displacementLimit));
    double trialStep =
        std::min({remaining, config_.maxSubstepS, cflStep, accelerationStep, transportStep});

    bool accepted = false;
    double pressureRelaxation = pressureRelaxationLimit;
    int rejectionAttempts = 0;
    while (!accepted) {
      if (stiffnessDt != trialStep || stiffness_.size() != phases_.size()) {
        // PCISPH's pressure response scales with dt^2 (Solenthaler and Pajarola
        // 2009, eqs. 7-10), so every rejected or dynamically shortened trial
        // needs its own measured stiffness.
        stiffness_ =
            calibratePressureStiffness(mass_, config_.resolution.spacing, support, trialStep);
        stiffnessDt = trialStep;
      }

      grid_.build(particles, support);
      computeNumberDensity(particles, grid_, &boundary, support);
      const FrameAcceleration frame =
          frameAcceleration(motion, timeS + elapsed + 0.5 * trialStep);
      for (std::size_t i = 0; i < particles.size(); ++i) {
        Vec3d acceleration{frame.uniform[0], frame.uniform[1], frame.uniform[2]};
        if (frame.rotating) {
          const Vec3d angularVelocity{frame.omega[0], frame.omega[1], frame.omega[2]};
          const Vec3d angularAcceleration{frame.alpha[0], frame.alpha[1], frame.alpha[2]};
          const Vec3d q = positionOf(particles, i);
          const Vec3d u = velocityOf(particles, i);
          if (config_.enableCoriolis) {
            acceleration = acceleration - 2.0 * cross(angularVelocity, u);
          }
          acceleration = acceleration - cross(angularAcceleration, q) -
                         cross(angularVelocity, cross(angularVelocity, q));
        }
        particles.ax[i] = static_cast<float>(acceleration.x);
        particles.ay[i] = static_cast<float>(acceleration.y);
        particles.az[i] = static_cast<float>(acceleration.z);
      }

      // Monaghan, Rep. Prog. Phys. 68 (2005): the pair scalar below is
      // symmetric (harmonic mu, rho_i*rho_j, and r.grad W), while v_i-v_j
      // changes sign. Thus m_i*a_ij = -m_j*a_ji for unequal phase masses.
      const double regularizerSquared = std::pow(0.01 * support, 2.0);
      for (std::size_t i = 0; i < particles.size(); ++i) {
        const Vec3d pi = positionOf(particles, i);
        const Vec3d vi = velocityOf(particles, i);
        const double mi = phaseMass(mass_, particles, i);
        const double rhoi = mi * particles.delta[i];
        const double mui = phaseViscosity(phases_, particles, i);
        Vec3d viscous{};
        grid_.forEachCandidate(i, [&](std::size_t j) {
          if (i == j || particles.delta[j] <= 0.0f) return;
          const Vec3d rij = pi - positionOf(particles, j);
          const double r2 = lengthSquared(rij);
          const double r = std::sqrt(r2);
          if (r <= kTiny || r >= support) return;
          const double muj = phaseViscosity(phases_, particles, j);
          if (mui + muj <= kTiny) return;
          const double mj = phaseMass(mass_, particles, j);
          const double rhoj = mj * particles.delta[j];
          const double muij = 2.0 * mui * muj / (mui + muj);
          const Vec3d grad = (wendlandGradMagnitude(r, support) / r) * rij;
          const double pairForceScale =
              mi * mj * 2.0 * muij / std::max(kTiny, rhoi * rhoj) *
              dot(rij, grad) / (r2 + regularizerSquared);
          viscous = viscous + (pairForceScale / mi) * (vi - velocityOf(particles, j));
        });
        particles.ax[i] = static_cast<float>(static_cast<double>(particles.ax[i]) + viscous.x);
        particles.ay[i] = static_cast<float>(static_cast<double>(particles.ay[i]) + viscous.y);
        particles.az[i] = static_cast<float>(static_cast<double>(particles.az[i]) + viscous.z);
      }

      if (config_.enableSurfaceTension) {
        addSurfaceAcceleration(particles, grid_, phases_, mass_, interface_, support,
                               densityError_);
      } else {
        computeColourGeometry(particles, grid_, support, densityError_);
      }

      std::fill(particles.pressure.begin(), particles.pressure.end(), 0.0f);
      std::fill(forceX_.begin(), forceX_.end(), 0.0f);
      std::fill(forceY_.begin(), forceY_.end(), 0.0f);
      std::fill(forceZ_.begin(), forceZ_.end(), 0.0f);
      double finalError = std::numeric_limits<double>::infinity();
      int iterations = 0;
      for (; iterations < config_.maxPressureIterations; ++iterations) {
        for (std::size_t i = 0; i < particles.size(); ++i) {
          const double particleMass = phaseMass(mass_, particles, i);
          Vec3d predictedVelocity =
              velocityOf(particles, i) +
              trialStep *
                  Vec3d{static_cast<double>(particles.ax[i]) + forceX_[i] / particleMass,
                        static_cast<double>(particles.ay[i]) + forceY_[i] / particleMass,
                        static_cast<double>(particles.az[i]) + forceZ_[i] / particleMass};
          Vec3d predictedPosition = positionOf(particles, i) + trialStep * predictedVelocity;
          projectContact(contactRadius, config_.wallFriction, boundary, predictedPosition,
                         predictedVelocity);
          predictedX_[i] = static_cast<float>(predictedPosition.x);
          predictedY_[i] = static_cast<float>(predictedPosition.y);
          predictedZ_[i] = static_cast<float>(predictedPosition.z);
          predictedVX_[i] = static_cast<float>(predictedVelocity.x);
          predictedVY_[i] = static_cast<float>(predictedVelocity.y);
          predictedVZ_[i] = static_cast<float>(predictedVelocity.z);
        }

        finalError = 0.0;
        for (std::size_t i = 0; i < particles.size(); ++i) {
          const Vec3d pi{predictedX_[i], predictedY_[i], predictedZ_[i]};
          double delta = 0.0;
          grid_.forEachCandidate(i, [&](std::size_t j) {
            const double r =
                length(pi - Vec3d{predictedX_[j], predictedY_[j], predictedZ_[j]});
            if (r < support) delta += wendlandW(r, support);
          });
          delta += boundary.boundaryDensity(boundary.query(pi.x, pi.y, pi.z).distance);
          densityError_[i] = (delta - delta0_) / delta0_;
          finalError = std::max(finalError, std::abs(densityError_[i]));
        }

        for (std::size_t i = 0; i < particles.size(); ++i) {
          const std::size_t phase = particles.phase[i];
          const double pressure =
              static_cast<double>(particles.pressure[i]) +
              pressureRelaxation * stiffness_[phase] * densityError_[i];
          particles.pressure[i] = static_cast<float>(std::max(0.0, pressure));
        }
        computePressureForce(particles, grid_, support, predictedX_, predictedY_, predictedZ_,
                             delta0_, densityError_, forceX_, forceY_, forceZ_);
        if (iterations + 1 >= config_.minPressureIterations &&
            finalError <= config_.densityTolerance) {
          ++iterations;
          break;
        }
      }

      bool reject = false;
      for (std::size_t i = 0; i < particles.size(); ++i) {
        const double particleMass = phaseMass(mass_, particles, i);
        Vec3d velocity =
            velocityOf(particles, i) +
            trialStep *
                Vec3d{static_cast<double>(particles.ax[i]) + forceX_[i] / particleMass,
                      static_cast<double>(particles.ay[i]) + forceY_[i] / particleMass,
                      static_cast<double>(particles.az[i]) + forceZ_[i] / particleMass};
        Vec3d position = positionOf(particles, i) + trialStep * velocity;
        double displacement = length(position - positionOf(particles, i));
        const bool finite = std::isfinite(displacement) && std::isfinite(velocity.x) &&
                            std::isfinite(velocity.y) && std::isfinite(velocity.z);
        if ((!finite || displacement > displacementLimit) && rejectionAttempts < 12) {
          reject = true;
          break;
        }
        if (!finite) {
          velocity = {};
          position = positionOf(particles, i);
          ++stats_.clampedParticles;
        } else if (displacement > displacementLimit) {
          // Twelve rejected, dt-recalibrated trials are the last-resort point:
          // count the transport-speed clamp and advance exactly 0.25 H rather
          // than accepting an unbounded position or retrying forever.
          velocity = (displacementLimit / (trialStep * length(velocity))) * velocity;
          position = positionOf(particles, i) + trialStep * velocity;
          ++stats_.clampedParticles;
        }
        const double speed = length(velocity);
        if (speed > config_.maxSpeed) {
          velocity = (config_.maxSpeed / speed) * velocity;
          ++stats_.clampedParticles;
        }
        projectContact(contactRadius, config_.wallFriction, boundary, position, velocity);
        predictedX_[i] = static_cast<float>(position.x);
        predictedY_[i] = static_cast<float>(position.y);
        predictedZ_[i] = static_cast<float>(position.z);
        predictedVX_[i] = static_cast<float>(velocity.x);
        predictedVY_[i] = static_cast<float>(velocity.y);
        predictedVZ_[i] = static_cast<float>(velocity.z);
      }

      stats_.pressureIterations += iterations;
      if (reject) {
        ++stats_.rejectedSubsteps;
        ++rejectionAttempts;
        pressureRelaxation *= 0.5;
        pressureRelaxationLimit = pressureRelaxation;
        trialStep *= 0.5;
        continue;
      }

      particles.px = predictedX_;
      particles.py = predictedY_;
      particles.pz = predictedZ_;
      particles.vx = predictedVX_;
      particles.vy = predictedVY_;
      particles.vz = predictedVZ_;
      stats_.maxDensityError = finalError;
      accepted = true;
    }

    if (config_.xsphSmoothing > 0.0) {
      for (std::size_t i = 0; i < particles.size(); ++i) {
        const Vec3d pi = positionOf(particles, i);
        const Vec3d vi = velocityOf(particles, i);
        Vec3d correction{};
        grid_.forEachCandidate(i, [&](std::size_t j) {
          if (i == j || particles.delta[j] <= 0.0f) return;
          const double r = length(pi - positionOf(particles, j));
          if (r >= support) return;
          correction = correction +
                       (wendlandW(r, support) / static_cast<double>(particles.delta[j])) *
                           (velocityOf(particles, j) - vi);
        });
        const Vec3d smoothed = vi + config_.xsphSmoothing * correction;
        predictedVX_[i] = static_cast<float>(smoothed.x);
        predictedVY_[i] = static_cast<float>(smoothed.y);
        predictedVZ_[i] = static_cast<float>(smoothed.z);
      }
      particles.vx = predictedVX_;
      particles.vy = predictedVY_;
      particles.vz = predictedVZ_;
    }
    elapsed += trialStep;
    lastSubstep = trialStep;
    ++stats_.substeps;
  }

  grid_.build(particles, support);
  computeNumberDensity(particles, grid_, &boundary, support);
  computeColourGeometry(particles, grid_, support, densityError_);
  stats_.maxDensityError = 0.0;
  stats_.maxSpeed = 0.0;
  for (std::size_t i = 0; i < particles.size(); ++i) {
    stats_.maxDensityError =
        std::max(stats_.maxDensityError,
                 std::abs(static_cast<double>(particles.delta[i]) / delta0_ - 1.0));
    stats_.maxSpeed = std::max(stats_.maxSpeed, length(velocityOf(particles, i)));
  }
  stats_.substepS = lastSubstep;
  return stats_.substeps;
}

void Solver::calibrateInterface(const VesselBoundary& boundary) {
  (void)boundary;
  if (interface_.calibrated) return;
  interface_.cohesionGain = 0.0;
  interface_.calibrated = false;
  calibrationError_.clear();

  const bool hasSurfaceTension =
      std::any_of(interface_.sigma.begin(), interface_.sigma.end(),
                  [](double sigma) { return sigma > 0.0; });
  if (phases_.size() < 2 || !hasSurfaceTension) {
    // A one-phase or sigma=0 model has no interface to calibrate. Treat this as
    // a successful no-op so callers do not repeatedly perform meaningless work.
    interface_.calibrated = true;
    return;
  }

  try {
    const double dx = config_.resolution.spacing;
    const double support = config_.resolution.support();
    double radius = 3.0 * support;
    // Pressure samples must be farther than 2 H from the diffuse interface.
    // Grow the sphere until that resolved core contains a useful 3-D stencil.
    for (;;) {
      std::size_t coreCount = 0;
      const int radiusExtent = static_cast<int>(std::ceil(radius / dx));
      for (int z = -radiusExtent; z <= radiusExtent; ++z) {
        for (int y = -radiusExtent; y <= radiusExtent; ++y) {
          for (int x = -radiusExtent; x <= radiusExtent; ++x) {
            if (length({x * dx, y * dx, z * dx}) < radius - 2.0 * support) {
              ++coreCount;
            }
          }
        }
      }
      if (coreCount >= 27) break;
      radius += 0.5 * support;
    }
    const double outerRadius = radius + 3.5 * support;
    Particles droplet;
    const int extent = static_cast<int>(std::ceil(outerRadius / dx));
    for (int z = -extent; z <= extent; ++z) {
      for (int y = -extent; y <= extent; ++y) {
        for (int x = -extent; x <= extent; ++x) {
          const Vec3d q{x * dx, y * dx, z * dx};
          if (length(q) <= outerRadius) {
            droplet.add(static_cast<float>(q.x), static_cast<float>(q.y),
                        static_cast<float>(q.z),
                        static_cast<uint8_t>(length(q) <= radius ? 0 : 1));
          }
        }
      }
    }

    NeighbourGrid calibrationGrid;
    std::vector<double> curvature;
    std::vector<float> px(droplet.size()), py(droplet.size()), pz(droplet.size());
    std::vector<float> pvx(droplet.size()), pvy(droplet.size()), pvz(droplet.size());
    std::vector<float> fx(droplet.size()), fy(droplet.size()), fz(droplet.size());
    std::vector<double> error(droplet.size());
    InterfaceModel trialModel = interface_;
    trialModel.sigma.assign(phases_.size() * phases_.size(), 0.0);
    trialModel.sigma[1] = 1.0;
    trialModel.sigma[phases_.size()] = 1.0;
    // Akinci et al. 2013 make the cohesion coefficient resolution-dependent.
    // A small trial avoids distorting the calibration sphere; linear response
    // then maps it to the Young-Laplace coefficient.
    constexpr double kTrialCoefficient = 0.01;
    trialModel.cohesionGain = kTrialCoefficient;
    trialModel.calibrated = true;

    const double calibrationDt = std::min(config_.maxSubstepS, 0.05 * std::sqrt(support));
    const std::vector<double> calibrationStiffness =
        calibratePressureStiffness(mass_, dx, support, calibrationDt);
    constexpr int kRelaxationSteps = 96;
    const int calibrationIterations = std::min(config_.maxPressureIterations, 6);
    for (int step = 0; step < kRelaxationSteps; ++step) {
      calibrationGrid.build(droplet, support);
      computeNumberDensity(droplet, calibrationGrid, nullptr, support);
      std::fill(droplet.ax.begin(), droplet.ax.end(), 0.0f);
      std::fill(droplet.ay.begin(), droplet.ay.end(), 0.0f);
      std::fill(droplet.az.begin(), droplet.az.end(), 0.0f);
      addSurfaceAcceleration(droplet, calibrationGrid, phases_, mass_, trialModel, support,
                             curvature);
      std::fill(droplet.pressure.begin(), droplet.pressure.end(), 0.0f);
      std::fill(fx.begin(), fx.end(), 0.0f);
      std::fill(fy.begin(), fy.end(), 0.0f);
      std::fill(fz.begin(), fz.end(), 0.0f);

      for (int iteration = 0; iteration < calibrationIterations; ++iteration) {
        for (std::size_t i = 0; i < droplet.size(); ++i) {
          const double mi = phaseMass(mass_, droplet, i);
          const Vec3d v =
              velocityOf(droplet, i) +
              calibrationDt *
                  Vec3d{droplet.ax[i] + fx[i] / mi, droplet.ay[i] + fy[i] / mi,
                        droplet.az[i] + fz[i] / mi};
          const Vec3d q = positionOf(droplet, i) + calibrationDt * v;
          px[i] = static_cast<float>(q.x);
          py[i] = static_cast<float>(q.y);
          pz[i] = static_cast<float>(q.z);
          pvx[i] = static_cast<float>(v.x);
          pvy[i] = static_cast<float>(v.y);
          pvz[i] = static_cast<float>(v.z);
        }
        for (std::size_t i = 0; i < droplet.size(); ++i) {
          double delta = 0.0;
          const Vec3d qi{px[i], py[i], pz[i]};
          calibrationGrid.forEachCandidate(i, [&](std::size_t j) {
            const double r = length(qi - Vec3d{px[j], py[j], pz[j]});
            if (r < support) delta += wendlandW(r, support);
          });
          error[i] = (delta - delta0_) / delta0_;
          const std::size_t phase = droplet.phase[i];
          droplet.pressure[i] = static_cast<float>(
              std::max(0.0, static_cast<double>(droplet.pressure[i]) +
                                calibrationStiffness[phase] * error[i]));
        }
        computePressureForce(droplet, calibrationGrid, support, px, py, pz, delta0_, error, fx,
                             fy, fz);
      }

      for (std::size_t i = 0; i < droplet.size(); ++i) {
        Vec3d velocity{pvx[i], pvy[i], pvz[i]};
        const double speed = length(velocity);
        if (speed > config_.maxSpeed) velocity = (config_.maxSpeed / speed) * velocity;
        const Vec3d position = positionOf(droplet, i) + calibrationDt * velocity;
        droplet.px[i] = static_cast<float>(position.x);
        droplet.py[i] = static_cast<float>(position.y);
        droplet.pz[i] = static_cast<float>(position.z);
        droplet.vx[i] = static_cast<float>(0.95 * velocity.x);
        droplet.vy[i] = static_cast<float>(0.95 * velocity.y);
        droplet.vz[i] = static_cast<float>(0.95 * velocity.z);
      }
    }

    double insidePressure = 0.0;
    double outsidePressure = 0.0;
    std::size_t insideCount = 0;
    std::size_t outsideCount = 0;
    for (std::size_t i = 0; i < droplet.size(); ++i) {
      const double r = length(positionOf(droplet, i));
      if (droplet.phase[i] == 0 && r < radius - 2.0 * support) {
        insidePressure += droplet.pressure[i];
        ++insideCount;
      } else if (droplet.phase[i] == 1 && r > radius + 2.0 * support &&
                 r < outerRadius - support) {
        outsidePressure += droplet.pressure[i];
        ++outsideCount;
      }
    }
    if (insideCount > 0) insidePressure /= static_cast<double>(insideCount);
    if (outsideCount > 0) outsidePressure /= static_cast<double>(outsideCount);
    const double pressureJump = insidePressure - outsidePressure;
    const double sigmaEffective = 0.5 * pressureJump * radius;
    if (insideCount == 0 || outsideCount == 0 || !(sigmaEffective > 1.0e-12) ||
        !std::isfinite(sigmaEffective)) {
      calibrationError_ = "Young-Laplace calibration produced no resolved pressure jump";
      return;
    }
    // For the trial coefficient c_t, sigma_eff = dp*R/2. Therefore the
    // coefficient reproducing dp=2*sigma/R is c = sigma*c_t/sigma_eff.
    interface_.cohesionGain = kTrialCoefficient / sigmaEffective;
    interface_.calibrated = true;
  } catch (const std::exception& error) {
    calibrationError_ = error.what();
    interface_.cohesionGain = 0.0;
    interface_.calibrated = false;
  } catch (...) {
    calibrationError_ = "unknown Young-Laplace calibration failure";
    interface_.cohesionGain = 0.0;
    interface_.calibrated = false;
  }
}

}  // namespace chemcad::fluid

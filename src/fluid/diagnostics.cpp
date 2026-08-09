#include "fluid/diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <utility>

#include "fluid/kernels.hpp"

namespace chemcad::fluid {
namespace {

constexpr double kDisplayTimeConstantS = 0.15;

bool samePhaseNeighbour(const Particles& particles, std::size_t i, std::size_t j,
                        double supportSquared) {
  if (particles.phase[i] != particles.phase[j]) return false;
  const double dx = static_cast<double>(particles.px[i]) - particles.px[j];
  const double dy = static_cast<double>(particles.py[i]) - particles.py[j];
  const double dz = static_cast<double>(particles.pz[i]) - particles.pz[j];
  return dx * dx + dy * dy + dz * dz <= supportSquared;
}

}  // namespace

Diagnostics DiagnosticsEngine::compute(const Particles& particles, const NeighbourGrid& grid,
                                       const VesselBoundary& boundary,
                                       const std::vector<PhaseMaterial>& phases,
                                       double spacing) {
  Diagnostics result;
  result.phases.resize(phases.size());

  const std::size_t particleCount = particles.size();
  if (spacing <= 0.0 || grid.size() != particleCount) return result;

  // Solver particles have m_i=rho0_i*dx^3, so V_i=m_i/rho0_i=dx^3 for
  // every phase. Use that physical particle volume consistently for both the
  // neighbour quadrature V_j W_ij and the outer area quadrature V_i|grad c_i|;
  // 1/delta_i is a density estimate, not the particle's represented volume.
  const double particleVolume = spacing * spacing * spacing;
  const double support = grid.support();
  const double supportSquared = support * support;

  component_.assign(particleCount, -1);
  componentVol_.clear();
  componentPhase_.clear();
  componentCentroidZ_.clear();
  componentVol_.reserve(particleCount);
  componentPhase_.reserve(particleCount);
  componentCentroidZ_.reserve(particleCount);

  // The outer scan fixes component numbering. The explicit stack avoids a
  // data-dependent call stack and every candidate list is supplied in ascending
  // particle index by NeighbourGrid, making the traversal reproducible.
  std::vector<std::size_t> stack;
  stack.reserve(particleCount);
  for (std::size_t seed = 0; seed < particleCount; ++seed) {
    if (component_[seed] >= 0) continue;

    const int componentIndex = static_cast<int>(componentVol_.size());
    const uint8_t phase = particles.phase[seed];
    componentVol_.push_back(0.0);
    componentPhase_.push_back(phase);
    componentCentroidZ_.push_back(0.0);

    component_[seed] = componentIndex;
    stack.push_back(seed);

    while (!stack.empty()) {
      const std::size_t i = stack.back();
      stack.pop_back();


      grid.forEachCandidate(i, [&](std::size_t j) {
        if (j >= particleCount || component_[j] >= 0) return;
        if (!samePhaseNeighbour(particles, i, j, supportSquared)) return;
        component_[j] = componentIndex;
        stack.push_back(j);
      });
    }
  }
  std::vector<std::size_t> componentCount(componentVol_.size(), 0);
  for (std::size_t i = 0; i < particleCount; ++i) {
    const std::size_t componentIndex = static_cast<std::size_t>(component_[i]);
    ++componentCount[componentIndex];
    componentCentroidZ_[componentIndex] += static_cast<double>(particles.pz[i]);
  }
  for (std::size_t componentIndex = 0; componentIndex < componentVol_.size();
       ++componentIndex) {
    componentVol_[componentIndex] =
        static_cast<double>(componentCount[componentIndex]) * particleVolume;
    if (componentCount[componentIndex] != 0) {
      componentCentroidZ_[componentIndex] /=
          static_cast<double>(componentCount[componentIndex]);
    }
  }


  std::vector<int> largestComponent(phases.size(), -1);
  for (std::size_t componentIndex = 0; componentIndex < componentVol_.size();
       ++componentIndex) {
    const std::size_t phase = componentPhase_[componentIndex];
    if (phase >= phases.size()) continue;
    result.phases[phase].totalMl += componentVol_[componentIndex] * 1.0e6;
    const int current = largestComponent[phase];
    if (current < 0 || componentVol_[componentIndex] > componentVol_[current]) {
      largestComponent[phase] = static_cast<int>(componentIndex);
    }
  }

  std::vector<int> bulkComponent(phases.size(), -1);
  for (std::size_t phase = 0; phase < phases.size(); ++phase) {
    const int candidate = largestComponent[phase];
    const double totalVolume = result.phases[phase].totalMl * 1.0e-6;
    if (candidate < 0 || totalVolume <= 0.0) continue;

    // Five percent is the neutral decision. Once a bulk exists it is retained
    // down to four percent; after it disappears a component must exceed six
    // percent to reappear. This Schmitt trigger is applied to the previous
    // classification, not to a filtered volume, so labels cannot chatter.
    double threshold = 0.05;
    if (havePrevious_ && phase < previous_.phases.size()) {
      threshold = previous_.phases[phase].bulkResolved ? 0.04 : 0.06;
    }
    if (componentVol_[candidate] >= threshold * totalVolume) {
      bulkComponent[phase] = candidate;
      result.phases[phase].bulkResolved = true;
      result.phases[phase].bulkMl = componentVol_[candidate] * 1.0e6;
    }
    result.phases[phase].dispersedMl =
        result.phases[phase].totalMl - result.phases[phase].bulkMl;
  }

  double dispersedVolume = 0.0;
  double totalVolume = 0.0;
  double diameterCubedSum = 0.0;
  double diameterSquaredSum = 0.0;
  for (const PhaseDiagnostics& phase : result.phases) {
    dispersedVolume += phase.dispersedMl * 1.0e-6;
    totalVolume += phase.totalMl * 1.0e-6;
  }
  for (std::size_t componentIndex = 0; componentIndex < componentVol_.size();
       ++componentIndex) {
    const std::size_t phase = componentPhase_[componentIndex];
    if (phase >= phases.size() || bulkComponent[phase] == static_cast<int>(componentIndex)) {
      continue;
    }
    const double diameter = std::cbrt(6.0 * componentVol_[componentIndex] / kPi);
    diameterCubedSum += diameter * diameter * diameter;
    diameterSquaredSum += diameter * diameter;
    ++result.dispersedComponents;
  }
  // Zero denotes the absence of dispersed components. It is not a numerical
  // proxy for an infinitely fine emulsion.
  result.sauterDiameterM = diameterSquaredSum > 0.0
                               ? diameterCubedSum / diameterSquaredSum
                               : 0.0;
  result.dispersedFraction = totalVolume > 0.0 ? dispersedVolume / totalVolume : 0.0;

  // For phases a and b, particles of all other phases are omitted from both
  // numerator and denominator. Each unordered pair is integrated once: we do
  // not sum |grad c_a|+|grad c_b|, which would count the same transition
  // twice. Pairwise indicators generalise integral |grad chi| dV without
  // assigning an arbitrary scalar ordering to three or more liquids.
  for (std::size_t phaseA = 0; phaseA < phases.size(); ++phaseA) {
    for (std::size_t phaseB = phaseA + 1; phaseB < phases.size(); ++phaseB) {
      for (std::size_t i = 0; i < particleCount; ++i) {
        const std::size_t centrePhase = particles.phase[i];
        if (centrePhase != phaseA && centrePhase != phaseB) continue;

        double denominator = 0.0;
        double numerator = 0.0;
        double gradDenominator[3] = {0.0, 0.0, 0.0};
        double gradNumerator[3] = {0.0, 0.0, 0.0};
        grid.forEachCandidate(i, [&](std::size_t j) {
          if (j >= particleCount) return;
          const std::size_t neighbourPhase = particles.phase[j];
          if (neighbourPhase != phaseA && neighbourPhase != phaseB) return;

          const double dx = static_cast<double>(particles.px[i]) - particles.px[j];
          const double dy = static_cast<double>(particles.py[i]) - particles.py[j];
          const double dz = static_cast<double>(particles.pz[i]) - particles.pz[j];
          const double distanceSquared = dx * dx + dy * dy + dz * dz;
          if (distanceSquared >= supportSquared) return;
          const double distance = std::sqrt(distanceSquared);
          const double weight = particleVolume * wendlandW(distance, support);
          const double indicator = neighbourPhase == phaseB ? 1.0 : 0.0;
          denominator += weight;
          numerator += indicator * weight;
          if (distance <= 0.0) return;

          const double radialDerivative =
              particleVolume * wendlandGradMagnitude(distance, support) / distance;
          const double gx = radialDerivative * dx;
          const double gy = radialDerivative * dy;
          const double gz = radialDerivative * dz;
          gradDenominator[0] += gx;
          gradDenominator[1] += gy;
          gradDenominator[2] += gz;
          gradNumerator[0] += indicator * gx;
          gradNumerator[1] += indicator * gy;
          gradNumerator[2] += indicator * gz;
        });

        if (denominator <= 0.0) continue;
        // Solenthaler and Pajarola, SCA 2008, eqs. 22-24: for c=N/D,
        // grad(c)=(D grad(N)-N grad(D))/D^2. Keeping both denominator
        // derivatives is essential near a free surface or a third phase.
        const double inverseDenominatorSquared = 1.0 / (denominator * denominator);
        const double gx = (gradNumerator[0] * denominator -
                           numerator * gradDenominator[0]) *
                          inverseDenominatorSquared;
        const double gy = (gradNumerator[1] * denominator -
                           numerator * gradDenominator[1]) *
                          inverseDenominatorSquared;
        const double gz = (gradNumerator[2] * denominator -
                           numerator * gradDenominator[2]) *
                          inverseDenominatorSquared;
        result.interfacialAreaM2 +=
            particleVolume * std::sqrt(gx * gx + gy * gy + gz * gz);
      }
    }
  }

  std::vector<std::size_t> densityOrder(phases.size());
  std::iota(densityOrder.begin(), densityOrder.end(), std::size_t{0});
  std::stable_sort(densityOrder.begin(), densityOrder.end(),
                   [&](std::size_t a, std::size_t b) {
                     if (phases[a].restDensity != phases[b].restDensity) {
                       return phases[a].restDensity > phases[b].restDensity;
                     }
                     return a < b;
                   });
  double cumulativeBulkVolume = 0.0;
  for (std::size_t phase : densityOrder) {
    cumulativeBulkVolume += result.phases[phase].bulkMl * 1.0e-6;
    result.phases[phase].layerTopM = boundary.heightForVolume(cumulativeBulkVolume);
  }
  result.freeSurfaceM = boundary.heightForVolume(cumulativeBulkVolume);

  for (std::size_t i = 0; i < particleCount; ++i) {
    const std::size_t phase = particles.phase[i];
    if (phase >= phases.size()) continue;
    const double vx = particles.vx[i];
    const double vy = particles.vy[i];
    const double vz = particles.vz[i];
    const double mass = phases[phase].restDensity * particleVolume;
    result.kineticEnergyJ += 0.5 * mass * (vx * vx + vy * vy + vz * vz);
  }
  result.valid = true;

  // compute owns the hysteretic classification even when a caller chooses not
  // to request display smoothing. The other previous values remain the last
  // filtered values and are updated only by smooth().
  if (!havePrevious_ || previous_.phases.size() != result.phases.size()) {
    previous_ = result;
    havePrevious_ = true;
  } else {
    for (std::size_t phase = 0; phase < result.phases.size(); ++phase) {
      previous_.phases[phase].bulkResolved = result.phases[phase].bulkResolved;
    }
  }
  return result;
}

Diagnostics DiagnosticsEngine::smooth(const Diagnostics& raw, double dt) {
  if (!havePrevious_ || previous_.phases.size() != raw.phases.size()) {
    previous_ = raw;
    havePrevious_ = true;
    return raw;
  }

  Diagnostics filtered = raw;
  const double alpha = dt > 0.0 ? -std::expm1(-dt / kDisplayTimeConstantS) : 0.0;
  for (std::size_t phase = 0; phase < raw.phases.size(); ++phase) {
    filtered.phases[phase].layerTopM =
        previous_.phases[phase].layerTopM +
        alpha * (raw.phases[phase].layerTopM - previous_.phases[phase].layerTopM);
  }
  filtered.freeSurfaceM = previous_.freeSurfaceM +
                          alpha * (raw.freeSurfaceM - previous_.freeSurfaceM);
  filtered.dispersedFraction =
      previous_.dispersedFraction +
      alpha * (raw.dispersedFraction - previous_.dispersedFraction);
  filtered.sauterDiameterM =
      previous_.sauterDiameterM +
      alpha * (raw.sauterDiameterM - previous_.sauterDiameterM);

  previous_ = filtered;
  return filtered;
}

void DiagnosticsEngine::reset() {
  previous_ = {};
  havePrevious_ = false;
  component_.clear();
  componentVol_.clear();
  componentPhase_.clear();
  componentCentroidZ_.clear();
}

}  // namespace chemcad::fluid

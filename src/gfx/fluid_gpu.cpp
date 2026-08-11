#include "gfx/fluid_gpu.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "fluid/kernels.hpp"
#include "gfx/compute.hpp"

namespace chemcad::gfx {
namespace {

constexpr std::uint32_t kLocalSize = 128;
constexpr std::size_t kProfileSamples = 1025;
constexpr std::size_t kBoundarySamples = 4097;

struct alignas(16) GpuParticle {
  float positionPhase[4];
  float velocityDelta[4];
  float accelerationPressure[4];
  float colourNormal[4];
};

struct alignas(16) GpuScratch {
  float predictedPositionError[4];
  float predictedVelocity[4];
  float pressureForce[4];
  float wall[4];
};

struct alignas(16) GpuPhase {
  float values[4];  // mass, rest density, dynamic viscosity, pressure stiffness
};

struct GpuProfilePoint {
  float radius;
  float z;
};

const char* kCommonShader = R"GLSL(#version 430 core
layout(local_size_x = 128) in;

const float PI = 3.14159265358979323846;
const float TINY = 1.0e-12;

struct Particle {
  vec4 positionPhase;
  vec4 velocityDelta;
  vec4 accelerationPressure;
  vec4 colourNormal;
};
struct Scratch {
  vec4 predictedPositionError;
  vec4 predictedVelocity;
  vec4 pressureForce;
  vec4 wall;
};

layout(std430, binding = 0) buffer ParticleBuffer { Particle particles[]; };
layout(std430, binding = 1) buffer PhaseBuffer { vec4 phaseValues[]; };
// Binding 2 is deliberately vacant: the interfacial-tension table left the GPU
// with the pairwise cohesion term. Renumbering the rest would only churn.
layout(std430, binding = 3) buffer ProfileBuffer { vec2 profile[]; };
layout(std430, binding = 4) buffer BoundaryBuffer { float boundaryTable[]; };
layout(std430, binding = 5) buffer CountBuffer { uint cellCounts[]; };
layout(std430, binding = 6) buffer OffsetBuffer { uint cellOffsets[]; };
layout(std430, binding = 7) buffer CursorBuffer { uint cellCursors[]; };
layout(std430, binding = 8) buffer SortedBuffer { uint sortedIndices[]; };
layout(std430, binding = 9) buffer ScratchBuffer { Scratch scratchData[]; };
layout(std430, binding = 10) buffer ReductionBuffer { uint reductions[]; };
layout(std430, binding = 11) buffer CurvatureBuffer { float curvature[]; };

uniform uint particleCount;
uniform uint phaseCount;
uniform uint cellCount;
uniform uint profileCount;
uniform uint boundaryCount;
uniform ivec3 gridDims;
uniform vec3 gridMin;
uniform float support;
uniform float delta0;
uniform float stepS;
uniform float contactRadius;
uniform float wallFriction;
uniform float maxSpeed;
uniform float displacementLimit;
uniform float pressureRelaxation;
uniform float surfaceTension;
uniform vec3 frameUniform;
uniform vec3 frameOmega;
uniform vec3 frameAlpha;
uniform int frameRotating;
uniform int enableCoriolis;
uniform int surfaceEnabled;
uniform int rejectionAttempt;
uniform int reductionMode;

uint phaseOf(uint i) {
  return uint(particles[i].positionPhase.w + 0.5);
}

ivec3 cellOf(vec3 p) {
  return clamp(ivec3(floor((p - gridMin) / support)), ivec3(0), gridDims - 1);
}

uint flatCell(ivec3 cell) {
  return uint(cell.x + gridDims.x * (cell.y + gridDims.y * cell.z));
}

float wendlandW(float radius) {
  if (radius >= support || support <= 0.0) return 0.0;
  float q = radius / support;
  float oneMinusQ = 1.0 - q;
  float q4 = oneMinusQ * oneMinusQ * oneMinusQ * oneMinusQ;
  return 21.0 / (2.0 * PI * support * support * support) *
         q4 * (1.0 + 4.0 * q);
}

float wendlandGradOverR(float radius) {
  if (radius >= support || radius <= 0.0 || support <= 0.0) return 0.0;
  float q = radius / support;
  float oneMinusQ = 1.0 - q;
  float gradient = -210.0 / (PI * pow(support, 4.0)) * q *
                   oneMinusQ * oneMinusQ * oneMinusQ;
  return gradient / radius;
}

// Interfacial tension is a Continuum Surface Force on the colour field, so no
// cohesion kernel appears here; see the CPU solver for the derivation.

struct BoundaryResult {
  float squaredDistance;
  float normalS;
  float normalZ;
};

void considerSegment(vec2 point, vec2 a, vec2 b, vec2 normal,
                     inout BoundaryResult result) {
  vec2 edge = b - a;
  float edgeSquared = dot(edge, edge);
  float parameter = edgeSquared > 0.0
      ? clamp(dot(point - a, edge) / edgeSquared, 0.0, 1.0)
      : 0.0;
  vec2 error = point - (a + parameter * edge);
  float distanceSquared = dot(error, error);
  if (distanceSquared < result.squaredDistance) {
    result.squaredDistance = distanceSquared;
    float normalLength = length(normal);
    vec2 unitNormal = normalLength > 0.0 ? normal / normalLength : vec2(0.0, 1.0);
    result.normalS = unitNormal.x;
    result.normalZ = unitNormal.y;
  }
}

float profileRadius(float z) {
  if (profileCount < 2u) return 0.0;
  float height = profile[profileCount - 1u].y;
  if (z <= 0.0 || height <= 0.0) return profile[0].x;
  if (z >= height) return profile[profileCount - 1u].x;
  float scaled = z * float(profileCount - 1u) / height;
  uint lower = min(uint(scaled), profileCount - 2u);
  return mix(profile[lower].x, profile[lower + 1u].x, scaled - float(lower));
}

vec4 queryBoundary(vec3 position) {
  float radial = length(position.xy);
  vec2 point = vec2(radial, position.z);
  BoundaryResult result = BoundaryResult(3.402823466e+38, 0.0, 1.0);
  considerSegment(point, vec2(0.0), vec2(profile[0].x, 0.0),
                  vec2(0.0, 1.0), result);
  for (uint segment = 0u; segment + 1u < profileCount; ++segment) {
    vec2 a = profile[segment];
    vec2 b = profile[segment + 1u];
    vec2 edge = b - a;
    considerSegment(point, a, b, vec2(-edge.y, edge.x), result);
  }
  float height = profile[profileCount - 1u].y;
  considerSegment(point, vec2(0.0, height),
                  vec2(profile[profileCount - 1u].x, height),
                  vec2(0.0, -1.0), result);
  bool inside = position.z >= 0.0 && position.z <= height &&
                radial <= profileRadius(position.z);
  float distance = sqrt(max(0.0, result.squaredDistance));
  vec2 radialDirection = radial >= 1.0e-12
      ? position.xy / radial : vec2(1.0, 0.0);
  return vec4(inside ? -distance : distance,
              result.normalS * radialDirection.x,
              result.normalS * radialDirection.y, result.normalZ);
}

float boundaryDensity(float distance) {
  if (boundaryCount < 2u || support <= 0.0 || distance <= -support) return 0.0;
  if (distance >= support) return boundaryTable[boundaryCount - 1u];
  float scaled = (distance + support) * float(boundaryCount - 1u) /
                 (2.0 * support);
  uint lower = min(uint(scaled), boundaryCount - 2u);
  return mix(boundaryTable[lower], boundaryTable[lower + 1u],
             scaled - float(lower));
}

void eachCellBounds(ivec3 cell, out uint begin, out uint end) {
  // Not `flat`: that is an interpolation qualifier, and NVIDIA's compiler
  // rejects it as a variable name even in a compute shader.
  uint index = flatCell(cell);
  begin = cellOffsets[index];
  end = begin + cellCounts[index];
}
)GLSL";

const char* kGridCountShader = R"GLSL(
void main() {
  uint i = gl_GlobalInvocationID.x;
  if (i >= particleCount) return;
  atomicAdd(cellCounts[flatCell(cellOf(particles[i].positionPhase.xyz))], 1u);
}
)GLSL";

const char* kGridPrefixShader = R"GLSL(
void main() {
  if (gl_GlobalInvocationID.x != 0u) return;
  uint offset = 0u;
  for (uint cell = 0u; cell < cellCount; ++cell) {
    cellOffsets[cell] = offset;
    cellCursors[cell] = offset;
    offset += cellCounts[cell];
  }
}
)GLSL";

const char* kGridScatterShader = R"GLSL(
void main() {
  uint i = gl_GlobalInvocationID.x;
  if (i >= particleCount) return;
  uint cell = flatCell(cellOf(particles[i].positionPhase.xyz));
  sortedIndices[atomicAdd(cellCursors[cell], 1u)] = i;
}
)GLSL";

const char* kDensityShader = R"GLSL(
void main() {
  uint i = gl_GlobalInvocationID.x;
  if (i >= particleCount) return;
  vec3 position = particles[i].positionPhase.xyz;
  ivec3 centre = cellOf(position);
  float density = 0.0;
  for (int z = -1; z <= 1; ++z) {
    for (int y = -1; y <= 1; ++y) {
      for (int x = -1; x <= 1; ++x) {
        ivec3 cell = centre + ivec3(x, y, z);
        if (any(lessThan(cell, ivec3(0))) || any(greaterThanEqual(cell, gridDims))) continue;
        uint begin, end;
        eachCellBounds(cell, begin, end);
        for (uint entry = begin; entry < end; ++entry) {
          uint j = sortedIndices[entry];
          float radius = length(position - particles[j].positionPhase.xyz);
          if (radius < support) density += wendlandW(radius);
        }
      }
    }
  }
  vec4 wall = queryBoundary(position);
  density += boundaryDensity(wall.x);
  particles[i].velocityDelta.w = density;
  scratchData[i].wall = wall;
}
)GLSL";

const char* kColourShader = R"GLSL(
void main() {
  uint i = gl_GlobalInvocationID.x;
  if (i >= particleCount) return;
  vec3 position = particles[i].positionPhase.xyz;
  ivec3 centre = cellOf(position);
  float numerator = 0.0;
  float denominator = 0.0;
  for (int z = -1; z <= 1; ++z) for (int y = -1; y <= 1; ++y)
    for (int x = -1; x <= 1; ++x) {
      ivec3 cell = centre + ivec3(x, y, z);
      if (any(lessThan(cell, ivec3(0))) || any(greaterThanEqual(cell, gridDims))) continue;
      uint begin, end; eachCellBounds(cell, begin, end);
      for (uint entry = begin; entry < end; ++entry) {
        uint j = sortedIndices[entry];
        float neighbourDelta = particles[j].velocityDelta.w;
        if (neighbourDelta <= 0.0) continue;
        float radius = length(position - particles[j].positionPhase.xyz);
        if (radius >= support) continue;
        float weight = wendlandW(radius) / neighbourDelta;
        denominator += weight;
        numerator += (phaseOf(j) == 0u ? 0.0 : 1.0) * weight;
      }
    }
  particles[i].colourNormal.x = denominator > TINY
      ? numerator / denominator : (phaseOf(i) == 0u ? 0.0 : 1.0);
}
)GLSL";

// Mirrors fluid::computeSurfaceGeometry exactly, including the Bonet-Lok
// kernel-gradient correction. Any divergence between the two would show up as
// a different interfacial tension on machines that take the GPU path.
const char* kNormalShader = R"GLSL(
void main() {
  uint i = gl_GlobalInvocationID.x;
  if (i >= particleCount) return;
  vec3 position = particles[i].positionPhase.xyz;
  ivec3 centre = cellOf(position);
  vec3 gradient = vec3(0.0);
  for (int z = -1; z <= 1; ++z) for (int y = -1; y <= 1; ++y)
    for (int x = -1; x <= 1; ++x) {
      ivec3 cell = centre + ivec3(x, y, z);
      if (any(lessThan(cell, ivec3(0))) || any(greaterThanEqual(cell, gridDims))) continue;
      uint begin, end; eachCellBounds(cell, begin, end);
      for (uint entry = begin; entry < end; ++entry) {
        uint j = sortedIndices[entry];
        if (i == j || particles[j].velocityDelta.w <= 0.0) continue;
        vec3 difference = position - particles[j].positionPhase.xyz;
        float radius = length(difference);
        if (radius <= 0.0 || radius >= support) continue;
        float scale = (particles[j].colourNormal.x - particles[i].colourNormal.x) /
                      particles[j].velocityDelta.w * wendlandGradOverR(radius);
        gradient += scale * difference;
      }
    }

  // Correct only where there is a gradient to correct: building L for the bulk
  // would be most of the cost of this pass for none of its value.
  if (length(gradient) >= MIN_COLOUR_GRADIENT / support) {
    mat3 correction = mat3(0.0);
    for (int z = -1; z <= 1; ++z) for (int y = -1; y <= 1; ++y)
      for (int x = -1; x <= 1; ++x) {
        ivec3 cell = centre + ivec3(x, y, z);
        if (any(lessThan(cell, ivec3(0))) || any(greaterThanEqual(cell, gridDims))) continue;
        uint begin, end; eachCellBounds(cell, begin, end);
        for (uint entry = begin; entry < end; ++entry) {
          uint j = sortedIndices[entry];
          if (i == j || particles[j].velocityDelta.w <= 0.0) continue;
          vec3 difference = position - particles[j].positionPhase.xyz;
          float radius = length(difference);
          if (radius <= 0.0 || radius >= support) continue;
          float weight = -wendlandGradOverR(radius) / particles[j].velocityDelta.w;
          correction += weight * outerProduct(difference, difference);
        }
      }
    if (determinant(correction) > MIN_CORRECTION_DETERMINANT) {
      gradient = inverse(correction) * gradient;
    }
  }
  // The UNNORMALISED gradient: the Continuum Surface Force needs its magnitude
  // as the surface delta, and the curvature pass normalises where it needs a
  // direction. Storing only the unit vector would throw that magnitude away.
  particles[i].colourNormal.yzw = gradient;
}
)GLSL";

const char* kCurvatureShader = R"GLSL(
void main() {
  uint i = gl_GlobalInvocationID.x;
  if (i >= particleCount) return;
  curvature[i] = 0.0;
  float floorValue = MIN_COLOUR_GRADIENT / support;
  vec3 gradientI = particles[i].colourNormal.yzw;
  float lengthI = length(gradientI);
  if (lengthI < floorValue) return;
  vec3 ni = gradientI / lengthI;
  vec3 position = particles[i].positionPhase.xyz;
  ivec3 centre = cellOf(position);
  mat3 jacobian = mat3(0.0);
  mat3 correction = mat3(0.0);
  for (int z = -1; z <= 1; ++z) for (int y = -1; y <= 1; ++y)
    for (int x = -1; x <= 1; ++x) {
      ivec3 cell = centre + ivec3(x, y, z);
      if (any(lessThan(cell, ivec3(0))) || any(greaterThanEqual(cell, gridDims))) continue;
      uint begin, end; eachCellBounds(cell, begin, end);
      for (uint entry = begin; entry < end; ++entry) {
        uint j = sortedIndices[entry];
        float neighbourDelta = particles[j].velocityDelta.w;
        if (i == j || neighbourDelta <= 0.0) continue;
        // Neighbours that actually have a normal. A bulk particle has none, and
        // folding its zero vector in would report every flat patch as curved.
        vec3 gradientJ = particles[j].colourNormal.yzw;
        float lengthJ = length(gradientJ);
        if (lengthJ < floorValue) continue;
        vec3 difference = position - particles[j].positionPhase.xyz;
        float radius = length(difference);
        if (radius <= 0.0 || radius >= support) continue;
        float weight = -wendlandGradOverR(radius) / neighbourDelta;
        correction += weight * outerProduct(difference, difference);
        // Rows index the normal component, columns the derivative direction,
        // and grad W has the opposite sign to the accumulation above.
        jacobian -= weight * outerProduct(gradientJ / lengthJ - ni, difference);
      }
    }
  if (!(determinant(correction) > MIN_CORRECTION_DETERMINANT)) return;
  mat3 corrected = inverse(correction) * jacobian;
  float value = -(corrected[0][0] + corrected[1][1] + corrected[2][2]);
  // Same bound as the CPU: a curvature radius below the smoothing length is
  // discretisation noise, and acting on it detonates the vessel.
  float limit = MAX_INTERFACE_CURVATURE / support;
  curvature[i] = clamp(value, -limit, limit);
}
)GLSL";

const char* kForceShader = R"GLSL(
void main() {
  uint i = gl_GlobalInvocationID.x;
  if (i >= particleCount) return;
  vec3 position = particles[i].positionPhase.xyz;
  vec3 velocity = particles[i].velocityDelta.xyz;
  uint phaseI = phaseOf(i);
  vec3 acceleration = frameUniform;
  if (frameRotating != 0) {
    if (enableCoriolis != 0) acceleration -= 2.0 * cross(frameOmega, velocity);
    acceleration -= cross(frameAlpha, position) + cross(frameOmega, cross(frameOmega, position));
  }

  float massI = phaseValues[phaseI].x;
  float rhoI = massI * particles[i].velocityDelta.w;
  float muI = phaseValues[phaseI].z;
  float regularizerSquared = pow(0.01 * support, 2.0);
  vec3 viscous = vec3(0.0);
  vec3 surface = vec3(0.0);
  ivec3 centre = cellOf(position);
  for (int z = -1; z <= 1; ++z) for (int y = -1; y <= 1; ++y)
    for (int x = -1; x <= 1; ++x) {
      ivec3 cell = centre + ivec3(x, y, z);
      if (any(lessThan(cell, ivec3(0))) || any(greaterThanEqual(cell, gridDims))) continue;
      uint begin, end; eachCellBounds(cell, begin, end);
      for (uint entry = begin; entry < end; ++entry) {
        uint j = sortedIndices[entry];
        if (i == j || particles[j].velocityDelta.w <= 0.0) continue;
        vec3 difference = position - particles[j].positionPhase.xyz;
        float radiusSquared = dot(difference, difference);
        if (radiusSquared <= TINY || radiusSquared >= support * support) continue;
        float radius = sqrt(radiusSquared);
        uint phaseJ = phaseOf(j);
        float massJ = phaseValues[phaseJ].x;
        float rhoJ = massJ * particles[j].velocityDelta.w;
        float muJ = phaseValues[phaseJ].z;
        if (muI + muJ > TINY) {
          float harmonicMu = 2.0 * muI * muJ / (muI + muJ);
          float pairForceScale = massI * massJ * 2.0 * harmonicMu /
              max(TINY, rhoI * rhoJ) *
              (wendlandGradOverR(radius) * radiusSquared) /
              (radiusSquared + regularizerSquared);
          viscous += (pairForceScale / massI) *
                     (velocity - particles[j].velocityDelta.xyz);
        }
      }
    }
  // Continuum Surface Force: sigma * kappa * grad c / rho, a purely local term
  // once the colour gradient and its divergence are resident.
  float colourI = particles[i].colourNormal.x;
  if (surfaceEnabled != 0 && surfaceTension > 0.0 && rhoI > TINY &&
      colourI > 0.05 && colourI < 0.95) {
    surface = (surfaceTension * curvature[i] / rhoI) * particles[i].colourNormal.yzw;
  }
  particles[i].accelerationPressure.xyz = acceleration + viscous + surface;
}
)GLSL";

const char* kPressureResetShader = R"GLSL(
void main() {
  uint i = gl_GlobalInvocationID.x;
  if (i >= particleCount) return;
  particles[i].accelerationPressure.w = 0.0;
  scratchData[i].pressureForce = vec4(0.0);
}
)GLSL";

const char* kPressurePredictShader = R"GLSL(
void main() {
  uint i = gl_GlobalInvocationID.x;
  if (i >= particleCount) return;
  float mass = phaseValues[phaseOf(i)].x;
  vec3 velocity = particles[i].velocityDelta.xyz + stepS *
      (particles[i].accelerationPressure.xyz + scratchData[i].pressureForce.xyz / mass);
  vec3 position = particles[i].positionPhase.xyz + stepS * velocity;
  vec3 normal = scratchData[i].wall.yzw;
  float distance = scratchData[i].wall.x -
      dot(normal, position - particles[i].positionPhase.xyz);
  if (distance > -contactRadius) {
    position += (distance + contactRadius) * normal;
    float normalVelocity = dot(velocity, normal);
    vec3 tangent = velocity - normalVelocity * normal;
    velocity = max(0.0, normalVelocity) * normal + (1.0 - wallFriction) * tangent;
  }
  scratchData[i].predictedPositionError.xyz = position;
  scratchData[i].predictedVelocity.xyz = velocity;
}
)GLSL";

const char* kPressureDensityShader = R"GLSL(
void main() {
  uint i = gl_GlobalInvocationID.x;
  if (i >= particleCount) return;
  vec3 predicted = scratchData[i].predictedPositionError.xyz;
  vec3 original = particles[i].positionPhase.xyz;
  ivec3 centre = cellOf(original);
  float density = 0.0;
  for (int z = -1; z <= 1; ++z) for (int y = -1; y <= 1; ++y)
    for (int x = -1; x <= 1; ++x) {
      ivec3 cell = centre + ivec3(x, y, z);
      if (any(lessThan(cell, ivec3(0))) || any(greaterThanEqual(cell, gridDims))) continue;
      uint begin, end; eachCellBounds(cell, begin, end);
      for (uint entry = begin; entry < end; ++entry) {
        uint j = sortedIndices[entry];
        // PCISPH freezes the neighbour set for the whole correction loop and
        // its stiffness is derived on that assumption (Solenthaler & Pajarola
        // 2009). The 27-cell block is a SUPERSET of the committed neighbours,
        // so without this test the device counts pairs the reference solver
        // never sees, over-reads its own density error and under-corrects: it
        // measured 1.5-7.8% worst-frame compression against the CPU's 1.0%.
        if (distance(original, particles[j].positionPhase.xyz) >= support) continue;
        float radius = length(predicted - scratchData[j].predictedPositionError.xyz);
        if (radius < support) density += wendlandW(radius);
      }
    }
  float wallDistance = scratchData[i].wall.x -
      dot(scratchData[i].wall.yzw, predicted - original);
  density += boundaryDensity(wallDistance);
  float error = (density - delta0) / delta0;
  scratchData[i].predictedPositionError.w = error;
  atomicMax(reductions[0], floatBitsToUint(max(0.0, error)));
}
)GLSL";

const char* kPressureUpdateShader = R"GLSL(
void main() {
  uint i = gl_GlobalInvocationID.x;
  if (i >= particleCount) return;
  uint phase = phaseOf(i);
  float pressure = particles[i].accelerationPressure.w +
      pressureRelaxation * phaseValues[phase].w *
      scratchData[i].predictedPositionError.w;
  particles[i].accelerationPressure.w = max(0.0, pressure);
}
)GLSL";

const char* kPressureForceShader = R"GLSL(
void main() {
  uint i = gl_GlobalInvocationID.x;
  if (i >= particleCount) return;
  vec3 position = scratchData[i].predictedPositionError.xyz;
  vec3 original = particles[i].positionPhase.xyz;
  ivec3 centre = cellOf(original);
  float deltaI = max(TINY, delta0 * (1.0 + scratchData[i].predictedPositionError.w));
  vec3 force = vec3(0.0);
  for (int z = -1; z <= 1; ++z) for (int y = -1; y <= 1; ++y)
    for (int x = -1; x <= 1; ++x) {
      ivec3 cell = centre + ivec3(x, y, z);
      if (any(lessThan(cell, ivec3(0))) || any(greaterThanEqual(cell, gridDims))) continue;
      uint begin, end; eachCellBounds(cell, begin, end);
      for (uint entry = begin; entry < end; ++entry) {
        uint j = sortedIndices[entry];
        if (i == j) continue;
        // The same frozen set the density prediction uses: the force and the
        // error it corrects have to be taken over the same pairs or the
        // correction does not converge to the error it measured.
        if (distance(original, particles[j].positionPhase.xyz) >= support) continue;
        vec3 difference = position - scratchData[j].predictedPositionError.xyz;
        float radiusSquared = dot(difference, difference);
        if (radiusSquared <= TINY || radiusSquared >= support * support) continue;
        float radius = sqrt(radiusSquared);
        float deltaJ = max(TINY, delta0 *
            (1.0 + scratchData[j].predictedPositionError.w));
        float scale = -(particles[i].accelerationPressure.w +
                        particles[j].accelerationPressure.w) /
                      (2.0 * deltaI * deltaJ) * wendlandGradOverR(radius);
        force += scale * difference;
      }
    }
  scratchData[i].pressureForce.xyz = force;
}
)GLSL";

const char* kIntegrateShader = R"GLSL(
void main() {
  uint i = gl_GlobalInvocationID.x;
  if (i >= particleCount) return;
  vec3 oldPosition = particles[i].positionPhase.xyz;
  float mass = phaseValues[phaseOf(i)].x;
  vec3 velocity = particles[i].velocityDelta.xyz + stepS *
      (particles[i].accelerationPressure.xyz + scratchData[i].pressureForce.xyz / mass);
  vec3 position = oldPosition + stepS * velocity;
  float displacement = length(position - oldPosition);
  bool finiteState = !any(isnan(vec4(position, displacement))) &&
                     !any(isinf(vec4(position, displacement))) &&
                     !any(isnan(velocity)) && !any(isinf(velocity));
  if ((!finiteState || displacement > displacementLimit) && rejectionAttempt < 12) {
    atomicOr(reductions[1], 1u);
    return;
  }
  uint clamped = 0u;
  if (!finiteState) {
    velocity = vec3(0.0);
    position = oldPosition;
    clamped = 1u;
  } else if (displacement > displacementLimit) {
    velocity *= displacementLimit / max(TINY, stepS * length(velocity));
    position = oldPosition + stepS * velocity;
    clamped = 1u;
  }
  float speed = length(velocity);
  if (speed > maxSpeed) {
    velocity *= maxSpeed / speed;
    clamped = 1u;
  }
  if (scratchData[i].wall.x + displacement > -contactRadius) {
    for (int projection = 0; projection < 4; ++projection) {
      vec4 wall = queryBoundary(position);
      if (wall.x <= -contactRadius) break;
      vec3 normal = wall.yzw;
      position += (wall.x + contactRadius) * normal;
      float normalVelocity = dot(velocity, normal);
      vec3 tangent = velocity - normalVelocity * normal;
      velocity = max(0.0, normalVelocity) * normal +
                 (1.0 - wallFriction) * tangent;
    }
  }
  scratchData[i].predictedPositionError.xyz = position;
  scratchData[i].predictedVelocity.xyz = velocity;
  if (clamped != 0u) atomicAdd(reductions[2], 1u);
  atomicMax(reductions[3], floatBitsToUint(length(velocity)));
}
)GLSL";

const char* kCommitShader = R"GLSL(
void main() {
  uint i = gl_GlobalInvocationID.x;
  if (i >= particleCount) return;
  particles[i].positionPhase.xyz = scratchData[i].predictedPositionError.xyz;
  particles[i].velocityDelta.xyz = scratchData[i].predictedVelocity.xyz;
  particles[i].velocityDelta.w = delta0 *
      (1.0 + scratchData[i].predictedPositionError.w);
}
)GLSL";

const char* kReductionShader = R"GLSL(
void main() {
  uint i = gl_GlobalInvocationID.x;
  if (i >= particleCount) return;
  float speed = length(particles[i].velocityDelta.xyz);
  atomicMax(reductions[0], floatBitsToUint(speed));
  if (reductionMode == 0) {
    atomicMax(reductions[1], floatBitsToUint(length(particles[i].accelerationPressure.xyz)));
  } else {
    float error = particles[i].velocityDelta.w / delta0 - 1.0;
    atomicMax(reductions[1], floatBitsToUint(max(0.0, error)));
    atomicMax(reductions[2], floatBitsToUint(max(0.0, -error)));
  }
}
)GLSL";

std::string shaderSource(const char* body) {
  std::string source(kCommonShader);
  // Compiled in from the single definition in fluid/solver.hpp rather than
  // written out here, so the GPU cannot decide an interface is somewhere the
  // CPU does not.
  source += "const float MIN_COLOUR_GRADIENT = " +
            std::to_string(fluid::kInterfaceGradientFloor) + ";\n";
  source += "const float MIN_CORRECTION_DETERMINANT = " +
            std::to_string(fluid::kInterfaceCorrectionDeterminant) + ";\n";
  source += "const float MAX_INTERFACE_CURVATURE = " +
            std::to_string(fluid::kMaxInterfaceCurvature) + ";\n";
  source += body;
  return source;
}

std::uint32_t groupsFor(std::size_t count) {
  return static_cast<std::uint32_t>((count + kLocalSize - 1) / kLocalSize);
}

float uintAsFloat(std::uint32_t value) {
  float result = 0.0f;
  static_assert(sizeof(result) == sizeof(value));
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

double resolutionSubstepCeiling(const fluid::SolverConfig& config) {
  return std::min(config.maxSubstepS,
                  config.cflNumber * config.resolution.support());
}

double quantizedSubstepLimit(double resolutionCeiling, double dynamicLimit) {
  constexpr std::array<double, 10> fractions{
      1.0, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0.0625, 0.03125};
  for (double fraction : fractions) {
    const double candidate = resolutionCeiling * fraction;
    if (candidate <= dynamicLimit * (1.0 + 1.0e-12)) return candidate;
  }
  return dynamicLimit;
}

std::vector<double> pressureStiffness(const std::vector<double>& mass,
                                      double spacing, double support, double dt) {
  struct Vec3 { double x, y, z; };
  const auto subtract = [](const Vec3& a, const Vec3& b) {
    return Vec3{a.x - b.x, a.y - b.y, a.z - b.z};
  };
  const auto length = [](const Vec3& a) {
    return std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
  };
  std::vector<double> result(mass.size(), 0.0);
  if (spacing <= 0.0 || support <= 0.0 || dt <= 0.0) return result;
  std::vector<Vec3> lattice{{0.0, 0.0, 0.0}};
  const int extent = static_cast<int>(std::ceil(support / spacing));
  for (int z = -extent; z <= extent; ++z) for (int y = -extent; y <= extent; ++y)
    for (int x = -extent; x <= extent; ++x) {
      if (x == 0 && y == 0 && z == 0) continue;
      Vec3 q{spacing * x, spacing * y, spacing * z};
      if (length(q) < support) lattice.push_back(q);
    }
  const double restDelta = fluid::restNumberDensity(spacing, support);
  for (std::size_t phase = 0; phase < mass.size(); ++phase) {
    std::vector<Vec3> force(lattice.size(), {0.0, 0.0, 0.0});
    for (std::size_t i = 0; i < lattice.size(); ++i) {
      for (std::size_t j = 0; j < lattice.size(); ++j) {
        if (i == j) continue;
        const double pi = i == 0 ? 1.0 : 0.0;
        const double pj = j == 0 ? 1.0 : 0.0;
        if (pi + pj == 0.0) continue;
        const Vec3 difference = subtract(lattice[i], lattice[j]);
        const double radius = length(difference);
        if (radius <= 1.0e-12 || radius >= support) continue;
        const double scale = -(pi + pj) / (2.0 * restDelta * restDelta) *
                             fluid::wendlandGradMagnitude(radius, support) / radius;
        force[i].x += scale * difference.x;
        force[i].y += scale * difference.y;
        force[i].z += scale * difference.z;
      }
    }
    const auto predicted = [&](std::size_t i) {
      const double scale = dt * dt / mass[phase];
      return Vec3{lattice[i].x + scale * force[i].x,
                  lattice[i].y + scale * force[i].y,
                  lattice[i].z + scale * force[i].z};
    };
    const Vec3 centre = predicted(0);
    double predictedDelta = 0.0;
    for (std::size_t j = 0; j < lattice.size(); ++j) {
      predictedDelta += fluid::wendlandW(length(subtract(centre, predicted(j))), support);
    }
    const double relativeError = (predictedDelta - restDelta) / restDelta;
    result[phase] = 1.0 / std::max(1.0e-12, -relativeError);
  }
  return result;
}

}  // namespace

struct FluidGpuSolver::Impl {
  fluid::SolverConfig config;
  std::vector<fluid::PhaseMaterial> phases;
  std::vector<double> mass;
  fluid::Solver setupSolver;
  fluid::InterfaceModel interfaceModel;
  fluid::Solver::Stats stats;
  ComputeCapabilities capabilities;
  std::string error;
  bool programsReady = false;
  bool stateResident = false;
  std::size_t particleCount = 0;
  std::array<int, 3> gridDims{1, 1, 1};
  std::array<float, 3> gridMin{};
  std::size_t cellCount = 1;
  double uploadedBoundaryHeight = -1.0;
  double uploadedBoundaryRadius = -1.0;
  bool configurationSupported = true;

  ComputeProgram gridCount;
  ComputeProgram gridPrefix;
  ComputeProgram gridScatter;
  ComputeProgram density;
  ComputeProgram colour;
  ComputeProgram normal;
  ComputeProgram curvatureProgram;
  ComputeProgram force;
  ComputeProgram pressureReset;
  ComputeProgram pressurePredict;
  ComputeProgram pressureDensity;
  ComputeProgram pressureUpdate;
  ComputeProgram pressureForce;
  ComputeProgram integrate;
  ComputeProgram commit;
  ComputeProgram reduction;

  StorageBuffer particles;
  StorageBuffer phaseBuffer;
  StorageBuffer profileBuffer;
  StorageBuffer boundaryBuffer;
  StorageBuffer counts;
  StorageBuffer offsets;
  StorageBuffer cursors;
  StorageBuffer sorted;
  StorageBuffer scratch;
  StorageBuffer reductions;
  StorageBuffer curvatureBuffer;

  std::vector<GpuPhase> gpuPhases;
  std::vector<GpuProfilePoint> profile;
  std::vector<float> boundaryDensity;
  std::vector<GpuParticle> transfer;
  std::vector<std::uint32_t> ids;
  std::vector<double> activeStiffness;
  double stiffnessStep = -1.0;

  std::array<ComputeProgram*, 16> programs() {
    return {&gridCount, &gridPrefix, &gridScatter, &density, &colour, &normal,
            &curvatureProgram, &force, &pressureReset, &pressurePredict,
            &pressureDensity, &pressureUpdate, &pressureForce, &integrate,
            &commit, &reduction};
  }

  std::array<const char*, 16> shaderBodies() const {
    return {kGridCountShader, kGridPrefixShader, kGridScatterShader, kDensityShader,
            kColourShader, kNormalShader, kCurvatureShader, kForceShader,
            kPressureResetShader, kPressurePredictShader, kPressureDensityShader,
            kPressureUpdateShader, kPressureForceShader, kIntegrateShader,
            kCommitShader, kReductionShader};
  }

  bool initialise() {
    if (const char* requested = std::getenv("CHEMCAD_FLUID_GPU");
        requested != nullptr &&
        (std::string(requested) == "0" || std::string(requested) == "false" ||
         std::string(requested) == "off" || std::string(requested) == "cpu")) {
      error = "GPU fluid path disabled by CHEMCAD_FLUID_GPU";
      return false;
    }
    if (!configurationSupported) return false;
    if (programsReady) return true;
    capabilities = probeCompute();
    if (!capabilities.available) {
      error = capabilities.reason;
      return false;
    }
    auto shader = shaderBodies();
    auto target = programs();
    for (std::size_t i = 0; i < target.size(); ++i) {
      std::string log;
      const std::string source = shaderSource(shader[i]);
      if (!target[i]->compile(source, &log)) {
        error = "fluid compute shader " + std::to_string(i) + " failed: " + log;
        for (ComputeProgram* program : target) program->reset();
        return false;
      }
    }
    programsReady = true;
    error.clear();
    return true;
  }

  void bindBuffers() const {
    particles.bind(0);
    phaseBuffer.bind(1);
    profileBuffer.bind(3);
    boundaryBuffer.bind(4);
    counts.bind(5);
    offsets.bind(6);
    cursors.bind(7);
    sorted.bind(8);
    scratch.bind(9);
    reductions.bind(10);
    curvatureBuffer.bind(11);
  }

  void setCommon(ComputeProgram& program) const {
    program.use();
    program.setUInt("particleCount", static_cast<std::uint32_t>(particleCount));
    program.setUInt("phaseCount", static_cast<std::uint32_t>(phases.size()));
    program.setUInt("cellCount", static_cast<std::uint32_t>(cellCount));
    program.setUInt("profileCount", static_cast<std::uint32_t>(kProfileSamples));
    program.setUInt("boundaryCount", static_cast<std::uint32_t>(kBoundarySamples));
    program.setIVec3("gridDims", gridDims[0], gridDims[1], gridDims[2]);
    program.setVec3("gridMin", gridMin[0], gridMin[1], gridMin[2]);
    program.setFloat("support", static_cast<float>(config.resolution.support()));
    program.setFloat("delta0", static_cast<float>(
        fluid::restNumberDensity(config.resolution.spacing,
                                 config.resolution.support())));
  }

  void dispatchParticles(ComputeProgram& program) const {
    setCommon(program);
    program.dispatch(groupsFor(particleCount));
    ComputeProgram::memoryBarrier();
  }

  bool uploadBoundary(const fluid::VesselBoundary& boundary) {
    profile.resize(kProfileSamples);
    boundaryDensity.resize(kBoundarySamples);
    const double height = boundary.heightM();
    const double support = config.resolution.support();
    for (std::size_t i = 0; i < kProfileSamples; ++i) {
      const double z = height * static_cast<double>(i) /
                       static_cast<double>(kProfileSamples - 1);
      profile[i] = {static_cast<float>(boundary.radiusAt(z)), static_cast<float>(z)};
    }
    for (std::size_t i = 0; i < kBoundarySamples; ++i) {
      const double distance = -support + 2.0 * support * static_cast<double>(i) /
                                           static_cast<double>(kBoundarySamples - 1);
      boundaryDensity[i] = static_cast<float>(boundary.boundaryDensity(distance));
    }
    if (!profileBuffer.resize(profile.size() * sizeof(GpuProfilePoint)) ||
        !profileBuffer.upload(profile.data(), profile.size() * sizeof(GpuProfilePoint)) ||
        !boundaryBuffer.resize(boundaryDensity.size() * sizeof(float)) ||
        !boundaryBuffer.upload(boundaryDensity.data(), boundaryDensity.size() * sizeof(float))) {
      error = "failed to upload vessel SDF tables";
      stateResident = false;
      return false;
    }
    uploadedBoundaryHeight = boundary.heightM();
    uploadedBoundaryRadius = boundary.maxRadiusM();
    return true;
  }

  bool uploadPhaseData() {
    if (phases.empty()) return false;
    gpuPhases.resize(phases.size());
    for (std::size_t i = 0; i < phases.size(); ++i) {
      gpuPhases[i].values[0] = static_cast<float>(mass[i]);
      gpuPhases[i].values[1] = static_cast<float>(phases[i].restDensity);
      gpuPhases[i].values[2] = static_cast<float>(phases[i].dynamicViscosity);
      gpuPhases[i].values[3] = i < activeStiffness.size()
                                  ? static_cast<float>(activeStiffness[i]) : 0.0f;
    }
    return phaseBuffer.resize(gpuPhases.size() * sizeof(GpuPhase)) &&
           phaseBuffer.upload(gpuPhases.data(), gpuPhases.size() * sizeof(GpuPhase));
  }

  bool ensureStiffness(double step) {
    const double scale = std::max(std::abs(stiffnessStep), std::abs(step));
    if (stiffnessStep > 0.0 && std::abs(stiffnessStep - step) <= 1.0e-3 * scale) return true;
    activeStiffness = pressureStiffness(mass, config.resolution.spacing,
                                       config.resolution.support(), step);
    stiffnessStep = step;
    ++stats.pressureStiffnessCalibrations;
    return uploadPhaseData();
  }

  void buildGrid() {
    counts.clearUInt();
    ComputeProgram::memoryBarrier();
    dispatchParticles(gridCount);
    setCommon(gridPrefix);
    gridPrefix.dispatch(1);
    ComputeProgram::memoryBarrier();
    dispatchParticles(gridScatter);
  }

  std::array<std::uint32_t, 4> readReductions() const {
    std::array<std::uint32_t, 4> values{};
    reductions.download(values.data(), sizeof(values));
    return values;
  }

  void clearReductions() {
    reductions.clearUInt();
    ComputeProgram::memoryBarrier();
  }
};

FluidGpuSolver::FluidGpuSolver() : impl_(std::make_unique<Impl>()) {}
FluidGpuSolver::~FluidGpuSolver() = default;
FluidGpuSolver::FluidGpuSolver(FluidGpuSolver&&) noexcept = default;
FluidGpuSolver& FluidGpuSolver::operator=(FluidGpuSolver&&) noexcept = default;

bool FluidGpuSolver::initialise() { return impl_->initialise(); }

bool FluidGpuSolver::available() const {
  return impl_->configurationSupported && impl_->programsReady &&
         impl_->capabilities.available;
}

const std::string& FluidGpuSolver::unavailableReason() const { return impl_->error; }

void FluidGpuSolver::configure(const fluid::SolverConfig& config) {
  const bool resolutionChanged =
      impl_->config.resolution.spacing != config.resolution.spacing;
  impl_->config = config;
  impl_->setupSolver.configure(config);
  impl_->interfaceModel = impl_->setupSolver.interfaceModel();
  impl_->mass.resize(impl_->phases.size());
  for (std::size_t i = 0; i < impl_->phases.size(); ++i) {
    impl_->mass[i] = impl_->phases[i].restDensity *
                     config.resolution.particleVolume();
  }
  impl_->stateResident = false;
  impl_->stiffnessStep = -1.0;
  impl_->configurationSupported = config.xsphSmoothing <= 0.0;
  if (!impl_->configurationSupported) {
    impl_->error = "GPU fluid path does not support optional XSPH display smoothing";
  }
}

const fluid::SolverConfig& FluidGpuSolver::config() const { return impl_->config; }

void FluidGpuSolver::setPhases(
    const std::vector<fluid::PhaseMaterial>& phases,
    const std::vector<double>& sigmaPairs) {
  impl_->setupSolver.setPhases(phases, sigmaPairs);
  impl_->phases = phases;
  impl_->mass.resize(phases.size());
  for (std::size_t i = 0; i < phases.size(); ++i) {
    impl_->mass[i] = phases[i].restDensity * impl_->config.resolution.particleVolume();
  }
  // setupSolver owns the validation and symmetrisation of the tension table;
  // reading its model back keeps one definition of what sigma means.
  impl_->interfaceModel = impl_->setupSolver.interfaceModel();
  impl_->activeStiffness.assign(phases.size(), 0.0);
  impl_->stiffnessStep = -1.0;
  impl_->stateResident = false;
}

bool FluidGpuSolver::upload(const fluid::Particles& source,
                            const fluid::VesselBoundary& boundary) {
  if (!impl_->initialise() || source.empty() || impl_->phases.empty()) return false;
  for (std::uint8_t phase : source.phase) {
    if (phase >= impl_->phases.size()) {
      impl_->error = "particle phase index is outside the GPU material table";
      return false;
    }
  }
  impl_->particleCount = source.size();
  const double support = impl_->config.resolution.support();
  const double margin = support;
  impl_->gridMin = {static_cast<float>(-boundary.maxRadiusM() - margin),
                    static_cast<float>(-boundary.maxRadiusM() - margin),
                    static_cast<float>(-margin)};
  impl_->gridDims = {
      std::max(1, static_cast<int>(std::ceil((2.0 * boundary.maxRadiusM() + 2.0 * margin) / support))),
      std::max(1, static_cast<int>(std::ceil((2.0 * boundary.maxRadiusM() + 2.0 * margin) / support))),
      std::max(1, static_cast<int>(std::ceil((boundary.heightM() + 2.0 * margin) / support)))};
  impl_->cellCount = static_cast<std::size_t>(impl_->gridDims[0]) * impl_->gridDims[1] *
                     impl_->gridDims[2];

  impl_->transfer.resize(source.size());
  impl_->ids = source.id;
  for (std::size_t i = 0; i < source.size(); ++i) {
    GpuParticle& p = impl_->transfer[i];
    p.positionPhase[0] = source.px[i]; p.positionPhase[1] = source.py[i];
    p.positionPhase[2] = source.pz[i]; p.positionPhase[3] = static_cast<float>(source.phase[i]);
    p.velocityDelta[0] = source.vx[i]; p.velocityDelta[1] = source.vy[i];
    p.velocityDelta[2] = source.vz[i]; p.velocityDelta[3] = source.delta[i];
    p.accelerationPressure[0] = source.ax[i]; p.accelerationPressure[1] = source.ay[i];
    p.accelerationPressure[2] = source.az[i]; p.accelerationPressure[3] = source.pressure[i];
    p.colourNormal[0] = source.colour[i]; p.colourNormal[1] = source.nx[i];
    p.colourNormal[2] = source.ny[i]; p.colourNormal[3] = source.nz[i];
  }

  const bool buffersReady =
      impl_->particles.resize(impl_->transfer.size() * sizeof(GpuParticle)) &&
      impl_->particles.upload(impl_->transfer.data(), impl_->transfer.size() * sizeof(GpuParticle)) &&
      impl_->counts.resize(impl_->cellCount * sizeof(std::uint32_t)) &&
      impl_->offsets.resize(impl_->cellCount * sizeof(std::uint32_t)) &&
      impl_->cursors.resize(impl_->cellCount * sizeof(std::uint32_t)) &&
      impl_->sorted.resize(source.size() * sizeof(std::uint32_t)) &&
      impl_->scratch.resize(source.size() * sizeof(GpuScratch)) &&
      impl_->reductions.resize(4 * sizeof(std::uint32_t)) &&
      impl_->curvatureBuffer.resize(source.size() * sizeof(float)) &&
      impl_->uploadBoundary(boundary) && impl_->uploadPhaseData();
  if (!buffersReady) {
    impl_->error = "failed to allocate GPU fluid buffers";
    impl_->stateResident = false;
    return false;
  }
  impl_->bindBuffers();
  impl_->stateResident = true;
  impl_->error.clear();
  return true;
}

bool FluidGpuSolver::resident() const { return impl_->stateResident; }

int FluidGpuSolver::advance(const fluid::VesselBoundary& boundary,
                            const fluid::VesselMotion& motion,
                            double timeS, double dt) {
  impl_->stats = {};
  if (!available() || !impl_->stateResident || dt <= 0.0) return 0;
  if (boundary.heightM() != impl_->uploadedBoundaryHeight ||
      boundary.maxRadiusM() != impl_->uploadedBoundaryRadius) {
    if (!impl_->uploadBoundary(boundary)) return 0;
    impl_->bindBuffers();
  }

  const auto started = std::chrono::steady_clock::now();
  const double support = impl_->config.resolution.support();
  const double contactRadius = impl_->config.contactRadiusFactor * impl_->config.resolution.spacing;
  const double displacementLimit = 0.25 * support;
  const double shakeOmega = 2.0 * fluid::kPi * std::abs(motion.shakeFrequencyHz);
  const double peakShakeAcceleration = motion.shaking && motion.shakeRemainingS > 0.0
      ? std::abs(motion.shakeAmplitudeM) * shakeOmega * shakeOmega : 0.0;
  const bool legacySmallCharge = impl_->particleCount < 256;
  double relaxationLimit = legacySmallCharge ? 1.0 : 0.5;
  double elapsed = 0.0;
  double lastSubstep = 0.0;

  while (elapsed < dt) {
    const double remaining = dt - elapsed;
    if (!(remaining > std::numeric_limits<double>::epsilon() * std::max(1.0, dt))) break;
    impl_->clearReductions();
    impl_->setCommon(impl_->reduction);
    impl_->reduction.setInt("reductionMode", 0);
    impl_->reduction.dispatch(groupsFor(impl_->particleCount));
    ComputeProgram::memoryBarrier();
    const auto limiting = impl_->readReductions();
    const double maxVelocity = uintAsFloat(limiting[0]);
    double maxAcceleration = uintAsFloat(limiting[1]);
    const fluid::FrameAcceleration limitingFrame = fluid::frameAcceleration(motion, timeS + elapsed);
    const double uniformMagnitude = std::sqrt(
        limitingFrame.uniform[0] * limitingFrame.uniform[0] +
        limitingFrame.uniform[1] * limitingFrame.uniform[1] +
        limitingFrame.uniform[2] * limitingFrame.uniform[2]);
    maxAcceleration = std::max(maxAcceleration, uniformMagnitude + peakShakeAcceleration);
    const double cflStep = impl_->config.cflNumber * support / std::max(maxVelocity, 1.0e-9);
    const double accelerationStep = impl_->config.accelerationSafety *
        std::sqrt(support / std::max(maxAcceleration, 1.0e-9));
    const double transportStep = 1.9 * displacementLimit /
        (maxVelocity + std::sqrt(maxVelocity * maxVelocity +
                                 4.0 * maxAcceleration * displacementLimit));
    const double ceiling = resolutionSubstepCeiling(impl_->config);
    const double limit = quantizedSubstepLimit(
        ceiling, std::min({ceiling, cflStep, accelerationStep, transportStep}));
    const int planned = std::max(1, static_cast<int>(std::ceil(remaining / limit - 1.0e-12)));
    double trialStep = remaining / static_cast<double>(planned);
    double pressureRelaxation = legacySmallCharge ? relaxationLimit : 0.6;
    int rejectionAttempts = 0;
    bool accepted = false;

    while (!accepted) {
      if (!impl_->ensureStiffness(trialStep)) {
        impl_->stateResident = false;
        impl_->error = "failed to upload pressure stiffness";
        return impl_->stats.substeps;
      }
      impl_->bindBuffers();
      impl_->buildGrid();
      impl_->dispatchParticles(impl_->density);
      impl_->dispatchParticles(impl_->colour);
      const double surfaceTension = impl_->interfaceModel.interfacialTension();
      const bool surface = impl_->config.enableSurfaceTension && surfaceTension > 0.0 &&
                           impl_->phases.size() > 1;
      if (surface) {
        impl_->dispatchParticles(impl_->normal);
        impl_->dispatchParticles(impl_->curvatureProgram);
      }

      const fluid::FrameAcceleration frame =
          fluid::frameAcceleration(motion, timeS + elapsed + 0.5 * trialStep);
      impl_->setCommon(impl_->force);
      impl_->force.setVec3("frameUniform", static_cast<float>(frame.uniform[0]),
                           static_cast<float>(frame.uniform[1]), static_cast<float>(frame.uniform[2]));
      impl_->force.setVec3("frameOmega", static_cast<float>(frame.omega[0]),
                           static_cast<float>(frame.omega[1]), static_cast<float>(frame.omega[2]));
      impl_->force.setVec3("frameAlpha", static_cast<float>(frame.alpha[0]),
                           static_cast<float>(frame.alpha[1]), static_cast<float>(frame.alpha[2]));
      impl_->force.setInt("frameRotating", frame.rotating ? 1 : 0);
      impl_->force.setInt("enableCoriolis", impl_->config.enableCoriolis ? 1 : 0);
      impl_->force.setInt("surfaceEnabled", surface ? 1 : 0);
      impl_->force.setFloat("surfaceTension", static_cast<float>(surfaceTension));
      impl_->force.dispatch(groupsFor(impl_->particleCount));
      ComputeProgram::memoryBarrier();
      impl_->dispatchParticles(impl_->pressureReset);

      double finalCompression = std::numeric_limits<double>::infinity();
      double bestCompression = std::numeric_limits<double>::infinity();
      int withoutImprovement = 0;
      bool pressureStalled = false;
      int iterations = 0;
      for (; iterations < impl_->config.maxPressureIterations; ++iterations) {
        impl_->setCommon(impl_->pressurePredict);
        impl_->pressurePredict.setFloat("stepS", static_cast<float>(trialStep));
        impl_->pressurePredict.setFloat("contactRadius", static_cast<float>(contactRadius));
        impl_->pressurePredict.setFloat("wallFriction", static_cast<float>(impl_->config.wallFriction));
        impl_->pressurePredict.dispatch(groupsFor(impl_->particleCount));
        ComputeProgram::memoryBarrier();

        impl_->clearReductions();
        impl_->dispatchParticles(impl_->pressureDensity);
        finalCompression = uintAsFloat(impl_->readReductions()[0]);
        const double improvement = std::isfinite(bestCompression) ? 2.0e-3 * bestCompression : 0.0;
        if (finalCompression + improvement < bestCompression) {
          bestCompression = finalCompression;
          withoutImprovement = 0;
        } else {
          ++withoutImprovement;
        }

        impl_->setCommon(impl_->pressureUpdate);
        impl_->pressureUpdate.setFloat("pressureRelaxation", static_cast<float>(pressureRelaxation));
        impl_->pressureUpdate.dispatch(groupsFor(impl_->particleCount));
        ComputeProgram::memoryBarrier();
        impl_->dispatchParticles(impl_->pressureForce);
        if (iterations + 1 >= impl_->config.minPressureIterations &&
            finalCompression <= impl_->config.densityTolerance) {
          ++iterations;
          break;
        }
        if (iterations + 1 >= impl_->config.minPressureIterations && withoutImprovement >= 3) {
          pressureStalled = true;
          ++iterations;
          break;
        }
      }
      impl_->stats.pressureIterations += iterations;
      if (pressureStalled) ++impl_->stats.stalledPressureSubsteps;

      impl_->clearReductions();
      impl_->setCommon(impl_->integrate);
      impl_->integrate.setFloat("stepS", static_cast<float>(trialStep));
      impl_->integrate.setFloat("contactRadius", static_cast<float>(contactRadius));
      impl_->integrate.setFloat("wallFriction", static_cast<float>(impl_->config.wallFriction));
      impl_->integrate.setFloat("maxSpeed", static_cast<float>(impl_->config.maxSpeed));
      impl_->integrate.setFloat("displacementLimit", static_cast<float>(displacementLimit));
      impl_->integrate.setInt("rejectionAttempt", rejectionAttempts);
      impl_->integrate.dispatch(groupsFor(impl_->particleCount));
      ComputeProgram::memoryBarrier();
      const auto integration = impl_->readReductions();
      if (integration[1] != 0u) {
        ++impl_->stats.rejectedSubsteps;
        ++rejectionAttempts;
        pressureRelaxation *= 0.5;
        if (legacySmallCharge) relaxationLimit = pressureRelaxation;
        trialStep *= 0.5;
        continue;
      }
      impl_->stats.clampedParticles += static_cast<int>(integration[2]);
      impl_->stats.maxSpeed = std::max(impl_->stats.maxSpeed,
                                       static_cast<double>(uintAsFloat(integration[3])));
      impl_->dispatchParticles(impl_->commit);
      impl_->stats.maxDensityCompression = finalCompression;
      impl_->stats.maxDensityError = finalCompression;
      accepted = true;
    }
    elapsed += trialStep;
    lastSubstep = trialStep;
    ++impl_->stats.substeps;
  }

  impl_->dispatchParticles(impl_->colour);
  impl_->clearReductions();
  impl_->setCommon(impl_->reduction);
  impl_->reduction.setInt("reductionMode", 1);
  impl_->reduction.dispatch(groupsFor(impl_->particleCount));
  ComputeProgram::memoryBarrier();
  const auto finalValues = impl_->readReductions();
  impl_->stats.maxSpeed = uintAsFloat(finalValues[0]);
  impl_->stats.maxDensityCompression = uintAsFloat(finalValues[1]);
  impl_->stats.maxDensityDeficit = uintAsFloat(finalValues[2]);
  impl_->stats.maxDensityError = impl_->stats.maxDensityCompression;
  impl_->stats.substepS = lastSubstep;
  impl_->stats.workerCount = 1;
  impl_->stats.pressureStiffnessSubstepS = impl_->stiffnessStep;
  const double substeps = std::max(1, impl_->stats.substeps);
  impl_->stats.millisecondsPerSubstep =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - started).count() /
      substeps;
  return impl_->stats.substeps;
}

bool FluidGpuSolver::download(fluid::Particles& destination) const {
  if (!impl_->stateResident || impl_->particleCount == 0) return false;
  impl_->transfer.resize(impl_->particleCount);
  if (!impl_->particles.download(impl_->transfer.data(),
                                 impl_->transfer.size() * sizeof(GpuParticle))) return false;
  destination.resize(impl_->particleCount);
  for (std::size_t i = 0; i < impl_->particleCount; ++i) {
    const GpuParticle& p = impl_->transfer[i];
    destination.px[i] = p.positionPhase[0]; destination.py[i] = p.positionPhase[1];
    destination.pz[i] = p.positionPhase[2]; destination.phase[i] = static_cast<std::uint8_t>(p.positionPhase[3]);
    destination.vx[i] = p.velocityDelta[0]; destination.vy[i] = p.velocityDelta[1];
    destination.vz[i] = p.velocityDelta[2]; destination.delta[i] = p.velocityDelta[3];
    destination.ax[i] = p.accelerationPressure[0];
    destination.ay[i] = p.accelerationPressure[1];
    destination.az[i] = p.accelerationPressure[2];
    destination.pressure[i] = p.accelerationPressure[3];
    destination.colour[i] = p.colourNormal[0];
    destination.nx[i] = p.colourNormal[1];
    destination.ny[i] = p.colourNormal[2];
    destination.nz[i] = p.colourNormal[3];
    destination.id[i] =
        i < impl_->ids.size() ? impl_->ids[i] : static_cast<std::uint32_t>(i);
  }
  return true;
}

const fluid::Solver::Stats& FluidGpuSolver::stats() const { return impl_->stats; }

}  // namespace chemcad::gfx

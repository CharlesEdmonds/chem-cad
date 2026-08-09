#include "gfx/fluid_renderer.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <exception>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

#include "sol/funnel.hpp"

namespace chemcad::gfx {
namespace {
// Five-pass screen-space fluid pipeline:
//   1. sphere impostors write nearest eye depth, phase fields, and speed;
//   2. adaptive mean-curvature flow smooths that depth at a radius derived from
//      the projected particle size, without sampling the background;
//   3. additive sphere chords accumulate an independent optical path per phase;
//   4. inverse-projection derivatives reconstruct normals, then ordered
//      Beer-Lambert layers, nearest-phase lighting, and the meniscus are composited;
//   5. the analytic vessel shell is blended rear face first, then front face.
// The passes render off-screen and GlState restores every shared-context state
// touched here before the ImGui OpenGL backend draws.
// GL_RGBA32F is core in OpenGL 3.0; keep the token local because gl_api.hpp
// exposes only the subset shared by the rest of the application.
constexpr GLenum kRgba32f = 0x8814;

struct InstanceData {
  float x;
  float y;
  float z;
  float phase;
  float indicator;
  float speed;
};

struct GlassVertex {
  float x;
  float y;
  float z;
  float nx;
  float ny;
  float nz;
};

struct GlState {
  std::array<GLint, 4> viewport{};
  GLint drawFramebuffer = 0;
  GLint readFramebuffer = 0;
  GLint program = 0;
  GLint vao = 0;
  GLint arrayBuffer = 0;
  GLint renderbuffer = 0;
  GLint activeTexture = 0;
  std::array<GLint, 4> texture2d{};
  GLint blendSrcRgb = 0;
  GLint blendDstRgb = 0;
  GLint blendSrcAlpha = 0;
  GLint blendDstAlpha = 0;
  GLint blendEquationRgb = 0;
  GLint blendEquationAlpha = 0;
  GLint depthFunction = 0;
  GLint cullMode = 0;
  GLint frontFace = 0;
  GLboolean depthWrite = GL_TRUE;
  std::array<GLboolean, 4> colourWrite{GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE};
  std::array<GLfloat, 4> clearColour{};
  bool blend = false;
  bool depth = false;
  bool cull = false;
  bool scissor = false;
  bool stencil = false;
  bool framebufferSrgb = false;

  GlState() {
    glGetIntegerv(GL_VIEWPORT, viewport.data());
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &drawFramebuffer);
    glGetIntegerv(GL_READ_FRAMEBUFFER_BINDING, &readFramebuffer);
    glGetIntegerv(GL_CURRENT_PROGRAM, &program);
    glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &vao);
    glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &arrayBuffer);
    glGetIntegerv(GL_RENDERBUFFER_BINDING, &renderbuffer);
    glGetIntegerv(GL_ACTIVE_TEXTURE, &activeTexture);
    for (int unit = 0; unit < 4; ++unit) {
      glActiveTexture(GL_TEXTURE0 + static_cast<GLenum>(unit));
      glGetIntegerv(GL_TEXTURE_BINDING_2D, &texture2d[unit]);
    }
    glActiveTexture(static_cast<GLenum>(activeTexture));
    glGetIntegerv(GL_BLEND_SRC_RGB, &blendSrcRgb);
    glGetIntegerv(GL_BLEND_DST_RGB, &blendDstRgb);
    glGetIntegerv(GL_BLEND_SRC_ALPHA, &blendSrcAlpha);
    glGetIntegerv(GL_BLEND_DST_ALPHA, &blendDstAlpha);
    glGetIntegerv(GL_BLEND_EQUATION_RGB, &blendEquationRgb);
    glGetIntegerv(GL_BLEND_EQUATION_ALPHA, &blendEquationAlpha);
    glGetIntegerv(GL_DEPTH_FUNC, &depthFunction);
    glGetIntegerv(GL_CULL_FACE_MODE, &cullMode);
    glGetIntegerv(GL_FRONT_FACE, &frontFace);
    glGetBooleanv(GL_DEPTH_WRITEMASK, &depthWrite);
    glGetFloatv(GL_COLOR_CLEAR_VALUE, clearColour.data());
    glGetBooleanv(GL_COLOR_WRITEMASK, colourWrite.data());
    blend = glIsEnabled(GL_BLEND) == GL_TRUE;
    depth = glIsEnabled(GL_DEPTH_TEST) == GL_TRUE;
    cull = glIsEnabled(GL_CULL_FACE) == GL_TRUE;
    scissor = glIsEnabled(GL_SCISSOR_TEST) == GL_TRUE;
    stencil = glIsEnabled(GL_STENCIL_TEST) == GL_TRUE;
    framebufferSrgb = glIsEnabled(GL_FRAMEBUFFER_SRGB) == GL_TRUE;
  }

  ~GlState() {
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, static_cast<GLuint>(drawFramebuffer));
    glBindFramebuffer(GL_READ_FRAMEBUFFER, static_cast<GLuint>(readFramebuffer));
    glUseProgram(static_cast<GLuint>(program));
    glBindVertexArray(static_cast<GLuint>(vao));
    glBindBuffer(GL_ARRAY_BUFFER, static_cast<GLuint>(arrayBuffer));
    glBindRenderbuffer(GL_RENDERBUFFER, static_cast<GLuint>(renderbuffer));
    for (int unit = 0; unit < 4; ++unit) {
      glActiveTexture(GL_TEXTURE0 + static_cast<GLenum>(unit));
      glBindTexture(GL_TEXTURE_2D, static_cast<GLuint>(texture2d[unit]));
    }
    glActiveTexture(static_cast<GLenum>(activeTexture));
    glClearColor(clearColour[0], clearColour[1], clearColour[2], clearColour[3]);
    glBlendFuncSeparate(static_cast<GLenum>(blendSrcRgb), static_cast<GLenum>(blendDstRgb),
                        static_cast<GLenum>(blendSrcAlpha),
                        static_cast<GLenum>(blendDstAlpha));
    glBlendEquationSeparate(static_cast<GLenum>(blendEquationRgb),
                            static_cast<GLenum>(blendEquationAlpha));
    glDepthFunc(static_cast<GLenum>(depthFunction));
    glDepthMask(depthWrite);
    glCullFace(static_cast<GLenum>(cullMode));
    glFrontFace(static_cast<GLenum>(frontFace));
    glColorMask(colourWrite[0], colourWrite[1], colourWrite[2], colourWrite[3]);
    blend ? glEnable(GL_BLEND) : glDisable(GL_BLEND);
    depth ? glEnable(GL_DEPTH_TEST) : glDisable(GL_DEPTH_TEST);
    cull ? glEnable(GL_CULL_FACE) : glDisable(GL_CULL_FACE);
    scissor ? glEnable(GL_SCISSOR_TEST) : glDisable(GL_SCISSOR_TEST);
    stencil ? glEnable(GL_STENCIL_TEST) : glDisable(GL_STENCIL_TEST);
    framebufferSrgb ? glEnable(GL_FRAMEBUFFER_SRGB) : glDisable(GL_FRAMEBUFFER_SRGB);
    glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
  }
};

const char* kParticleVertex = R"GLSL(#version 330 core
layout(location = 0) in vec2 aCorner;
layout(location = 1) in vec3 aCentre;
layout(location = 2) in float aPhase;
layout(location = 3) in float aIndicator;
layout(location = 4) in float aSpeed;

uniform mat4 uModel;
uniform mat4 uView;
uniform mat4 uProjection;
uniform float uRadius;

out vec2 vCorner;
out vec3 vCentreEye;
out vec3 vCentreObject;
flat out float vPhase;
flat out float vIndicator;
flat out float vSpeed;

void main() {
  vec3 centreEye = (uView * uModel * vec4(aCentre, 1.0)).xyz;
  vec3 cornerEye = centreEye + vec3(aCorner * uRadius, 0.0);
  gl_Position = uProjection * vec4(cornerEye, 1.0);
  vCorner = aCorner;
  vCentreEye = centreEye;
  vCentreObject = aCentre;
  vPhase = aPhase;
  vIndicator = aIndicator;
  vSpeed = aSpeed;
}
)GLSL";

const char* kDepthFragment = R"GLSL(#version 330 core
in vec2 vCorner;
in vec3 vCentreEye;
in vec3 vCentreObject;
flat in float vPhase;
flat in float vIndicator;
flat in float vSpeed;

uniform mat4 uProjection;
uniform float uRadius;
// Vessel interior, sampled from the same analytic half-width profile that
// builds the glass shell, in vessel coordinates.
uniform mat3 uEyeToObject;
uniform float uVesselHeight;
uniform float uVesselMaxRadius;
uniform float uVesselProfile[72];

layout(location = 0) out vec4 outSurface;

void main() {
  float radius2 = dot(vCorner, vCorner);
  if (radius2 > 1.0) discard;

  // The quad is only a carrier. Reconstructing the positive-z hemisphere in
  // eye space makes each instance a true sphere impostor, including correct
  // per-fragment hardware depth for overlap between particles.
  vec3 sphereNormal = vec3(vCorner, sqrt(max(0.0, 1.0 - radius2)));
  vec3 surfaceEye = vCentreEye + uRadius * sphereNormal;

  // Rendering the surface at a full dx makes the union of spheres solid, but
  // it also lets a wall-adjacent particle bulge through the glass, and in the
  // stem -- narrower than dx -- it would swell into a string of beads. The
  // solver already confines the CENTRES, so clipping each fragment against the
  // vessel profile is what confines the SURFACE. Doing it here rather than in
  // the thickness pass is sufficient: shading reads thickness only where this
  // attachment carries a depth.
  vec3 surfaceObject = vCentreObject + uEyeToObject * (uRadius * sphereNormal);
  float level = clamp(surfaceObject.z / max(uVesselHeight, 1.0e-6), 0.0, 1.0) *
                float(uVesselProfile.length() - 1);
  int lower = int(floor(level));
  int upper = min(lower + 1, uVesselProfile.length() - 1);
  float wallRadius = uVesselMaxRadius *
                     mix(uVesselProfile[lower], uVesselProfile[upper], level - float(lower));
  if (dot(surfaceObject.xy, surfaceObject.xy) > wallRadius * wallRadius) discard;

  vec4 clip = uProjection * vec4(surfaceEye, 1.0);
  gl_FragDepth = 0.5 * (clip.z / clip.w) + 0.5;
  // Depth, nearest phase, solver-smoothed phase indicator, and surface speed
  // share one depth-tested attachment so shading cannot sample another particle.
  outSurface = vec4(-surfaceEye.z, step(0.5, vPhase),
                    clamp(vIndicator, 0.0, 1.0), max(vSpeed, 0.0));
}
)GLSL";

const char* kThicknessFragment = R"GLSL(#version 330 core
in vec2 vCorner;
flat in float vPhase;
// Metres. Calibrated on the CPU so one particle's footprint integrates to
// exactly the volume of liquid it represents; see the thickness amplitude in
// drawFrame.
uniform float uAmplitude;
layout(location = 0) out vec2 outThickness;

void main() {
  float radius2 = dot(vCorner, vCorner);
  if (radius2 > 1.0) discard;

  // Not the geometric sphere chord 2R sqrt(1-r^2). That profile has an
  // infinite slope at the rim, so every particle stamps a visible disc edge
  // and the summed field is mottled at particle scale even when the surface
  // above it is smooth. (1-r^2)^(3/2) lands on zero with zero slope, and
  // because the footprint spans a full dx the discs overlap enough for the sum
  // to read as one body of liquid. Additive blending keeps the two phases'
  // optical paths independent.
  float falloff = max(0.0, 1.0 - radius2);
  float weight = falloff * sqrt(falloff);
  float phaseB = step(0.5, vPhase);
  outThickness = uAmplitude * weight * vec2(1.0 - phaseB, phaseB);
}
)GLSL";

const char* kThicknessBlurFragment = R"GLSL(#version 330 core
in vec2 vUv;
uniform sampler2D uThickness;
uniform sampler2D uSurface;
uniform vec2 uViewport;
uniform vec2 uProjectionScale;
uniform float uParticleRadius;
uniform vec2 uDirection;
layout(location = 0) out vec2 outThickness;

void main() {
  // Optical path is an integral, so unlike depth it may be filtered directly:
  // a Gaussian removes the particle-scale granularity of the summed footprints
  // while preserving the total. Without this the surface above can be perfectly
  // smooth and the body still looks scaly, because Beer-Lambert turns a 17%
  // ripple in path length into a visible mottle.
  const float weights[7] =
      float[7](0.199681, 0.176221, 0.121107, 0.064836, 0.027019, 0.008766, 0.002216);

  ivec2 pixel = ivec2(gl_FragCoord.xy);
  vec2 centre = texelFetch(uThickness, pixel, 0).rg;
  float depth = texelFetch(uSurface, pixel, 0).x;
  if (depth <= 0.0) {
    outThickness = centre;
    return;
  }

  // Match the kernel to the projected particle radius so the filter follows
  // zoom instead of being a fixed pixel count. sigma is two taps.
  vec2 radiusPixels2 =
      0.5 * uParticleRadius * abs(uProjectionScale) * uViewport / depth;
  float projectedRadiusPixels = sqrt(max(1.0, radiusPixels2.x * radiusPixels2.y));
  vec2 stepUv = uDirection * max(1.0, 0.5 * projectedRadiusPixels) / uViewport;

  vec2 total = centre * weights[0];
  for (int i = 1; i < 7; ++i) {
    total += texture(uThickness, vUv + float(i) * stepUv).rg * weights[i];
    total += texture(uThickness, vUv - float(i) * stepUv).rg * weights[i];
  }
  outThickness = total;
}
)GLSL";

const char* kFullscreenVertex = R"GLSL(#version 330 core
out vec2 vUv;
void main() {
  // One oversized triangle covers the target without a diagonal interpolation seam.
  vec2 p = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
  vUv = p;
  gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
}
)GLSL";

const char* kBlurFragment = R"GLSL(#version 330 core
in vec2 vUv;
uniform sampler2D uSurface;
uniform mat4 uInverseProjection;
uniform vec2 uViewport;
uniform vec2 uProjectionScale;
uniform float uParticleRadius;
layout(location = 0) out vec4 outSurface;

vec3 eyePosition(vec2 uv, float linearDepth) {
  vec2 ndc = uv * 2.0 - 1.0;
  vec4 farPoint = uInverseProjection * vec4(ndc, 1.0, 1.0);
  vec3 ray = farPoint.xyz / farPoint.w;
  return ray * (linearDepth / max(1.0e-6, -ray.z));
}

vec4 guardedSample(ivec2 pixel, ivec2 pixelOffset, vec4 centre,
                   float maximumJump) {
  ivec2 limit = ivec2(uViewport) - ivec2(1);
  vec4 sampleSurface =
      texelFetch(uSurface, clamp(pixel + pixelOffset, ivec2(0), limit), 0);
  // A zero depth is the background, not a very distant piece of liquid. A
  // Neumann boundary (return the centre value) keeps both the silhouette and
  // surface attributes from being pulled toward that sentinel. texelFetch also
  // prevents bilinear interpolation from manufacturing a false foreground
  // depth between the liquid and its zero-depth background.
  if (sampleSurface.x <= 0.0 || abs(sampleSurface.x - centre.x) > maximumJump)
    return centre;
  return sampleSurface;
}

void main() {
  ivec2 pixel = ivec2(gl_FragCoord.xy);
  vec4 centre = texelFetch(uSurface, pixel, 0);
  if (centre.x <= 0.0) {
    outSurface = vec4(0.0);
    return;
  }

  // Perspective projects R metres to R*f/(2z) pixels. Sampling at 1.25
  // projected radii makes the curvature stencil cross particle centres at the
  // coarse dx=2R spacing instead of degenerating to a fixed few-pixel blur.
  vec2 radiusPixels2 =
      0.5 * uParticleRadius * abs(uProjectionScale) * uViewport / centre.x;
  float projectedRadiusPixels =
      sqrt(max(1.0, radiusPixels2.x * radiusPixels2.y));
  int stencilRadiusPixels =
      int(round(clamp(1.25 * projectedRadiusPixels, 1.0, 24.0)));
  vec2 stepUv = vec2(stencilRadiusPixels) / uViewport;

  // Reject a separate depth layer but retain the full slope of one particle.
  // This range guard only defines the curvature-flow domain; the smoothing
  // itself below is the mean-curvature PDE rather than a bilateral average.
  float maximumJump = max(2.5 * uParticleRadius, 0.012 * centre.x);
  vec4 left =
      guardedSample(pixel, ivec2(-stencilRadiusPixels, 0), centre, maximumJump);
  vec4 right =
      guardedSample(pixel, ivec2(stencilRadiusPixels, 0), centre, maximumJump);
  vec4 down =
      guardedSample(pixel, ivec2(0, -stencilRadiusPixels), centre, maximumJump);
  vec4 up =
      guardedSample(pixel, ivec2(0, stencilRadiusPixels), centre, maximumJump);
  vec4 leftDown = guardedSample(
      pixel, ivec2(-stencilRadiusPixels, -stencilRadiusPixels), centre, maximumJump);
  vec4 rightDown = guardedSample(
      pixel, ivec2(stencilRadiusPixels, -stencilRadiusPixels), centre, maximumJump);
  vec4 leftUp = guardedSample(
      pixel, ivec2(-stencilRadiusPixels, stencilRadiusPixels), centre, maximumJump);
  vec4 rightUp = guardedSample(
      pixel, ivec2(stencilRadiusPixels, stencilRadiusPixels), centre, maximumJump);

  // Express the depth graph in local eye-space metres before evaluating
  // z_t = H sqrt(1+|grad z|^2). This is the explicit mean-curvature-flow form
  // from van der Laan, Green & Sainz rather than an image-space Gaussian.
  vec3 flatLeft = eyePosition(vUv - vec2(stepUv.x, 0.0), centre.x);
  vec3 flatRight = eyePosition(vUv + vec2(stepUv.x, 0.0), centre.x);
  vec3 flatDown = eyePosition(vUv - vec2(0.0, stepUv.y), centre.x);
  vec3 flatUp = eyePosition(vUv + vec2(0.0, stepUv.y), centre.x);
  float hx = max(1.0e-6, 0.5 * length(flatRight - flatLeft));
  float hy = max(1.0e-6, 0.5 * length(flatUp - flatDown));

  float zx = (right.x - left.x) / (2.0 * hx);
  float zy = (up.x - down.x) / (2.0 * hy);
  float zxx = (right.x - 2.0 * centre.x + left.x) / (hx * hx);
  float zyy = (up.x - 2.0 * centre.x + down.x) / (hy * hy);
  float zxy =
      (rightUp.x - rightDown.x - leftUp.x + leftDown.x) / (4.0 * hx * hy);
  float gradient2 = zx * zx + zy * zy;
  float flow =
      ((1.0 + zy * zy) * zxx - 2.0 * zx * zy * zxy +
       (1.0 + zx * zx) * zyy) /
      (2.0 * (1.0 + gradient2));

  // 0.16*h^2 is below the 2-D explicit diffusion stability limit. The public
  // default requests four smoothing iterations; two stable PDE steps per
  // request give eight steps. Their sqrt(N)*1.25R influence spans neighbouring
  // centres at dx=2R while retaining an interactive, bounded pass count.
  float flowTime = 0.16 * min(hx * hx, hy * hy);
  float smoothedDepth = max(1.0e-5, centre.x + flowTime * flow);

  // Transport surface attributes with the same accepted stencil so the
  // meniscus and velocity cue cannot drift across a rejected depth layer.
  vec3 attributeMean =
      (4.0 * centre.yzw + left.yzw + right.yzw + down.yzw + up.yzw) / 8.0;
  vec3 smoothedAttributes = mix(centre.yzw, attributeMean, 0.35);
  outSurface = vec4(smoothedDepth, clamp(smoothedAttributes.xy, 0.0, 1.0),
                    max(smoothedAttributes.z, 0.0));
}
)GLSL";

const char* kShadeFragment = R"GLSL(#version 330 core
in vec2 vUv;
uniform sampler2D uSurface;
uniform sampler2D uThickness;
uniform mat4 uInverseProjection;
uniform vec2 uTexel;
uniform vec2 uViewport;
uniform vec2 uProjectionScale;
uniform float uParticleRadius;
uniform vec4 uPhaseA;
uniform vec4 uPhaseB;
uniform vec3 uAbsorptionA;
uniform vec3 uAbsorptionB;
uniform vec3 uBackground;
uniform float uAbsorptionScale;
uniform float uExposure;
uniform float uShowInterface;
layout(location = 0) out vec4 outColour;

vec3 eyePosition(vec2 uv, float linearDepth) {
  vec2 ndc = uv * 2.0 - 1.0;
  vec4 farPoint = uInverseProjection * vec4(ndc, 1.0, 1.0);
  vec3 ray = farPoint.xyz / farPoint.w;
  return ray * (linearDepth / max(1.0e-6, -ray.z));
}

float projectedParticleRadius(float linearDepth) {
  vec2 radiusPixels2 =
      0.5 * uParticleRadius * abs(uProjectionScale) * uViewport / linearDepth;
  return sqrt(max(1.0, radiusPixels2.x * radiusPixels2.y));
}

float meanTransmission(vec3 transmission) {
  return dot(transmission, vec3(1.0 / 3.0));
}

vec3 compositeLayer(vec3 behind, vec3 phaseTint, vec3 transmission,
                    float diffuse) {
  // This is the closed-form Beer-Lambert solution for a homogeneous layer
  // with phase-coloured in-scattered radiance. Applying it once per phase,
  // far layer before near layer, preserves both optical order and identity.
  float coverage = 1.0 - meanTransmission(transmission);
  float sourceStrength = coverage * (0.34 + 0.46 * diffuse);
  return behind * transmission + phaseTint * sourceStrength;
}

void main() {
  ivec2 pixel = ivec2(gl_FragCoord.xy);
  ivec2 limit = ivec2(uViewport) - ivec2(1);
  vec4 surface = texelFetch(uSurface, pixel, 0);
  float depth = surface.x;
  if (depth <= 0.0) {
    outColour = vec4(uBackground, 1.0);
    return;
  }

  vec2 offsetX = vec2(uTexel.x, 0.0);
  vec2 offsetY = vec2(0.0, uTexel.y);
  vec4 surfaceL =
      texelFetch(uSurface, clamp(pixel + ivec2(-1, 0), ivec2(0), limit), 0);
  vec4 surfaceR =
      texelFetch(uSurface, clamp(pixel + ivec2(1, 0), ivec2(0), limit), 0);
  vec4 surfaceD =
      texelFetch(uSurface, clamp(pixel + ivec2(0, -1), ivec2(0), limit), 0);
  vec4 surfaceU =
      texelFetch(uSurface, clamp(pixel + ivec2(0, 1), ivec2(0), limit), 0);
  bool haveL = surfaceL.x > 0.0;
  bool haveR = surfaceR.x > 0.0;
  bool haveD = surfaceD.x > 0.0;
  bool haveU = surfaceU.x > 0.0;

  vec3 positionEye = eyePosition(vUv, depth);
  vec3 tangentX;
  if (haveL && haveR) {
    tangentX = eyePosition(vUv + offsetX, surfaceR.x) -
               eyePosition(vUv - offsetX, surfaceL.x);
  } else if (haveR) {
    tangentX = eyePosition(vUv + offsetX, surfaceR.x) - positionEye;
  } else if (haveL) {
    tangentX = positionEye - eyePosition(vUv - offsetX, surfaceL.x);
  } else {
    tangentX = eyePosition(vUv + offsetX, depth) -
               eyePosition(vUv - offsetX, depth);
  }
  vec3 tangentY;
  if (haveD && haveU) {
    tangentY = eyePosition(vUv + offsetY, surfaceU.x) -
               eyePosition(vUv - offsetY, surfaceD.x);
  } else if (haveU) {
    tangentY = eyePosition(vUv + offsetY, surfaceU.x) - positionEye;
  } else if (haveD) {
    tangentY = positionEye - eyePosition(vUv - offsetY, surfaceD.x);
  } else {
    tangentY = eyePosition(vUv + offsetY, depth) -
               eyePosition(vUv - offsetY, depth);
  }
  // Each derivative uses exact inverse-projected eye positions. At a silhouette
  // it switches to a foreground one-sided derivative (or a same-depth tangent),
  // so a zero-depth background texel can never enter the reconstructed normal.
  vec3 normalEye = normalize(cross(tangentX, tangentY));
  if (normalEye.z < 0.0) normalEye = -normalEye;

  vec3 viewDirection = normalize(-positionEye);
  vec3 lightDirection = normalize(vec3(-0.35, 0.55, 0.76));
  vec3 halfDirection = normalize(lightDirection + viewDirection);
  float diffuse = max(dot(normalEye, lightDirection), 0.0);
  float specular = pow(max(dot(normalEye, halfDirection), 0.0), 88.0);
  float cosTheta = clamp(dot(normalEye, viewDirection), 0.0, 1.0);
  // Restrict the free-surface highlight to grazing normals; a camera-facing
  // normal has a zero rim term instead of receiving a broad Fresnel addition.
  float rim = pow(1.0 - cosTheta, 4.5);

  float surfacePhase = step(0.5, clamp(surface.y, 0.0, 1.0));
  vec3 surfaceTint = mix(uPhaseA.rgb, uPhaseB.rgb, surfacePhase);
  vec2 thicknessM = max(texture(uThickness, vUv).rg, vec2(0.0));
  float strength = max(uAbsorptionScale, 0.0);

  // The RG thickness attachment retains one optical path per phase. With the
  // CPU's 30 mm reference, a 30 mm path has transmission equal to the configured
  // material colour by construction.
  vec3 transmissionA = exp(-strength * uAbsorptionA * thicknessM.x);
  vec3 transmissionB = exp(-strength * uAbsorptionB * thicknessM.y);
  vec3 transmission = transmissionA * transmissionB;
  float presenceA = 1.0 - exp(-strength * thicknessM.x / 0.012);
  float presenceB = 1.0 - exp(-strength * thicknessM.y / 0.012);

  vec3 refractedBackground = uBackground * (0.92 + 0.18 * normalEye.x);
  vec3 body;
  if (surfacePhase < 0.5) {
    vec3 behindA =
        compositeLayer(refractedBackground, uPhaseB.rgb, transmissionB, diffuse);
    body = compositeLayer(behindA, uPhaseA.rgb, transmissionA, diffuse);
  } else {
    vec3 behindB =
        compositeLayer(refractedBackground, uPhaseA.rgb, transmissionA, diffuse);
    body = compositeLayer(behindB, uPhaseB.rgb, transmissionB, diffuse);
  }

  float surfacePresence = mix(presenceA, presenceB, surfacePhase);
  // Diffuse, specular, rim, and motion cues all use the depth-tested nearest
  // phase. The phase-coloured motion contribution is capped at 0.035.
  vec3 specularTint = mix(surfaceTint, vec3(0.94), 0.28);
  body += specularTint * surfacePresence * (0.08 * specular + 0.15 * rim);
  float velocityCue = 1.0 - exp(-max(surface.w, 0.0) / 0.25);
  body += surfaceTint * surfacePresence * (0.035 * velocityCue);
  body = max(body, vec3(0.012) + surfaceTint * (0.025 * surfacePresence));

  float indicator = clamp(surface.z, 0.0, 1.0);
  float indicatorL = haveL ? surfaceL.z : indicator;
  float indicatorR = haveR ? surfaceR.z : indicator;
  float indicatorD = haveD ? surfaceD.z : indicator;
  float indicatorU = haveU ? surfaceU.z : indicator;
  vec2 indicatorGradient =
      0.5 * vec2(indicatorR - indicatorL, indicatorU - indicatorD);
  float gradientStrength = length(indicatorGradient);
  float interfaceDistancePixels =
      abs(indicator - 0.5) / max(gradientStrength, 1.0e-4);
  // The meniscus is one pixel at minimum and otherwise 14% of the projected
  // particle radius. With no fixed upper clamp, its screen width follows zoom.
  float interfaceWidthPixels =
      max(1.0, 0.14 * projectedParticleRadius(depth));
  float interfaceBand =
      1.0 - smoothstep(0.35 * interfaceWidthPixels, interfaceWidthPixels,
                       interfaceDistancePixels);
  float interfaceCore =
      1.0 - smoothstep(0.0, 0.32 * interfaceWidthPixels,
                       interfaceDistancePixels);
  float showMeniscus = clamp(uShowInterface, 0.0, 1.0);
  body = mix(body, vec3(0.020, 0.026, 0.036),
             showMeniscus * 0.58 * interfaceBand);
  body += showMeniscus * interfaceCore * vec3(0.055, 0.065, 0.085);

  // Exposure is a true linear multiplier around 1.0. The fixed shoulder only
  // protects highlights, then explicit sRGB encoding prevents a linear-RGBA8
  // target from making mid-tones appear much darker than intended.
  vec3 exposed = max(uExposure, 0.0) * max(body, vec3(0.0));
  body = pow(exposed / (vec3(1.0) + exposed), vec3(1.0 / 2.2));

  float absorptionCoverage = 1.0 - meanTransmission(transmission);
  float maximumAlpha = max(uPhaseA.a, uPhaseB.a);
  float alpha =
      clamp(max(maximumAlpha * (0.24 + 0.76 * absorptionCoverage),
                0.10 + 0.20 * rim),
            0.0, 0.985);
  // Resolve the liquid once against the stage background and keep the texture
  // opaque. ImGui would otherwise apply alpha a second time.
  outColour = vec4(mix(uBackground, body, alpha), 1.0);
}
)GLSL";

const char* kGlassVertex = R"GLSL(#version 330 core
layout(location = 0) in vec3 aPosition;
layout(location = 1) in vec3 aNormal;
uniform mat4 uModel;
uniform mat4 uView;
uniform mat4 uProjection;
out vec3 vNormalEye;
out vec3 vPositionEye;
void main() {
  vec4 world = uModel * vec4(aPosition, 1.0);
  vec4 eye = uView * world;
  vPositionEye = eye.xyz;
  vNormalEye = mat3(uView) * mat3(uModel) * aNormal;
  gl_Position = uProjection * eye;
}
)GLSL";

const char* kGlassFragment = R"GLSL(#version 330 core
in vec3 vNormalEye;
in vec3 vPositionEye;
uniform vec3 uTint;
uniform float uBackFace;
layout(location = 0) out vec4 outColour;
void main() {
  vec3 n = normalize(vNormalEye);
  vec3 v = normalize(-vPositionEye);
  if (dot(n, v) < 0.0) n = -n;
  // Front and back surfaces are both visible in real thin glass. Their Schlick
  // rims add, with the rear wall deliberately dimmer to retain depth ordering.
  float facing = clamp(dot(n, v), 0.0, 1.0);
  float fresnel = 0.04 + 0.96 * pow(1.0 - facing, 5.0);
  vec3 lightDirection = normalize(vec3(-0.45, 0.62, 0.70));
  vec3 halfDirection = normalize(lightDirection + v);
  // A narrow Blinn-Phong streak supplies the perceptual cue that the otherwise
  // transparent revolved shell is glass rather than an outline.
  float streak = pow(max(dot(n, halfDirection), 0.0), 72.0);
  float side = mix(1.0, 0.52, uBackFace);
  float alpha = side * (0.040 + 0.24 * fresnel + 0.10 * streak);
  vec3 colour = mix(uTint, vec3(1.0), 0.55 * streak);
  outColour = vec4(colour, alpha);
}
)GLSL";

GLuint compileShader(GLenum type, const char* source, const char* label, std::string& error) {
  GLuint shader = glCreateShader(type);
  if (shader == 0) {
    error = std::string("OpenGL could not create the ") + label + " shader";
    return 0;
  }
  glShaderSource(shader, 1, &source, nullptr);
  glCompileShader(shader);
  GLint status = 0;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
  if (status == static_cast<GLint>(GL_TRUE)) return shader;

  GLint length = 0;
  glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);
  std::string log(static_cast<std::size_t>(std::max(length, 1)), '\0');
  GLsizei written = 0;
  glGetShaderInfoLog(shader, static_cast<GLsizei>(log.size()), &written, log.data());
  log.resize(static_cast<std::size_t>(std::max(written, 0)));
  error = std::string(label) + " shader compilation failed: " + log;
  glDeleteShader(shader);
  return 0;
}

GLuint linkProgram(const char* vertex, const char* fragment, const char* label,
                   std::string& error) {
  GLuint vertexShader = compileShader(GL_VERTEX_SHADER, vertex, label, error);
  if (vertexShader == 0) return 0;
  GLuint fragmentShader = compileShader(GL_FRAGMENT_SHADER, fragment, label, error);
  if (fragmentShader == 0) {
    glDeleteShader(vertexShader);
    return 0;
  }

  GLuint program = glCreateProgram();
  if (program == 0) {
    error = std::string("OpenGL could not create the ") + label + " program";
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    return 0;
  }
  glAttachShader(program, vertexShader);
  glAttachShader(program, fragmentShader);
  glLinkProgram(program);
  glDeleteShader(vertexShader);
  glDeleteShader(fragmentShader);

  GLint status = 0;
  glGetProgramiv(program, GL_LINK_STATUS, &status);
  if (status == static_cast<GLint>(GL_TRUE)) return program;

  GLint length = 0;
  glGetProgramiv(program, GL_INFO_LOG_LENGTH, &length);
  std::string log(static_cast<std::size_t>(std::max(length, 1)), '\0');
  GLsizei written = 0;
  glGetProgramInfoLog(program, static_cast<GLsizei>(log.size()), &written, log.data());
  log.resize(static_cast<std::size_t>(std::max(written, 0)));
  error = std::string(label) + " program link failed: " + log;
  glDeleteProgram(program);
  return 0;
}

Mat4 poseMatrix(const fluid::Pose& pose) {
  double w = pose.orientation[0];
  double x = pose.orientation[1];
  double y = pose.orientation[2];
  double z = pose.orientation[3];
  const double magnitude = std::sqrt(w * w + x * x + y * y + z * z);
  if (magnitude > 1.0e-12) {
    w /= magnitude;
    x /= magnitude;
    y /= magnitude;
    z /= magnitude;
  } else {
    w = 1.0;
    x = y = z = 0.0;
  }

  const float r00 = static_cast<float>(1.0 - 2.0 * (y * y + z * z));
  const float r01 = static_cast<float>(2.0 * (x * y - z * w));
  const float r02 = static_cast<float>(2.0 * (x * z + y * w));
  const float r10 = static_cast<float>(2.0 * (x * y + z * w));
  const float r11 = static_cast<float>(1.0 - 2.0 * (x * x + z * z));
  const float r12 = static_cast<float>(2.0 * (y * z - x * w));
  const float r20 = static_cast<float>(2.0 * (x * z - y * w));
  const float r21 = static_cast<float>(2.0 * (y * z + x * w));
  const float r22 = static_cast<float>(1.0 - 2.0 * (x * x + y * y));
  return {r00, r10, r20, 0.0f, r01, r11, r21, 0.0f,
          r02, r12, r22, 0.0f, static_cast<float>(pose.position[0]),
          static_cast<float>(pose.position[1]), static_cast<float>(pose.position[2]), 1.0f};
}

Mat4 inverse(const Mat4& matrix) {
  double augmented[4][8]{};
  for (int row = 0; row < 4; ++row) {
    for (int column = 0; column < 4; ++column)
      augmented[row][column] = matrix[column * 4 + row];
    augmented[row][row + 4] = 1.0;
  }

  for (int pivotColumn = 0; pivotColumn < 4; ++pivotColumn) {
    int pivotRow = pivotColumn;
    for (int row = pivotColumn + 1; row < 4; ++row) {
      if (std::abs(augmented[row][pivotColumn]) >
          std::abs(augmented[pivotRow][pivotColumn]))
        pivotRow = row;
    }
    if (std::abs(augmented[pivotRow][pivotColumn]) < 1.0e-14)
      return {1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
              0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f};
    if (pivotRow != pivotColumn)
      for (int column = 0; column < 8; ++column)
        std::swap(augmented[pivotRow][column], augmented[pivotColumn][column]);

    const double divisor = augmented[pivotColumn][pivotColumn];
    for (double& value : augmented[pivotColumn]) value /= divisor;
    for (int row = 0; row < 4; ++row) {
      if (row == pivotColumn) continue;
      const double factor = augmented[row][pivotColumn];
      for (int column = 0; column < 8; ++column)
        augmented[row][column] -= factor * augmented[pivotColumn][column];
    }
  }

  Mat4 result{};
  for (int row = 0; row < 4; ++row)
    for (int column = 0; column < 4; ++column)
      result[column * 4 + row] = static_cast<float>(augmented[row][column + 4]);
  return result;
}

std::array<float, 4> phaseColour(const fluid::Snapshot& snapshot, std::size_t index,
                                 const std::array<float, 4>& fallback) {
  if (index >= snapshot.phases.size()) return fallback;
  const auto& colour = snapshot.phases[index].colour;
  return {colour[0], colour[1], colour[2], colour[3]};
}

std::array<float, 3> absorptionFromColour(const std::array<float, 4>& colour) {
  // In the absence of measured extinction spectra, interpret the configured
  // display RGB as transmittance through a 30 mm path, matching the nominal
  // optical scale of one liquid layer rather than a metre-scale reference.
  constexpr float referencePathM = 0.030f;
  return {-std::log(std::clamp(colour[0], 0.03f, 0.98f)) / referencePathM,
          -std::log(std::clamp(colour[1], 0.03f, 0.98f)) / referencePathM,
          -std::log(std::clamp(colour[2], 0.03f, 0.98f)) / referencePathM};
}

std::string glErrorMessage(GLenum code) {
  std::ostringstream stream;
  stream << "OpenGL error 0x" << std::hex << code;
  return stream.str();
}

}  // namespace

struct FluidRenderer::Impl {
  bool valid = false;
  std::string failure;
  int targetWidth = 0;
  int targetHeight = 0;
  double frameMilliseconds = 0.0;

  GLuint depthProgram = 0;
  GLuint blurProgram = 0;
  GLuint thicknessProgram = 0;
  GLuint thicknessBlurProgram = 0;
  GLuint shadeProgram = 0;
  GLuint glassProgram = 0;

  GLuint particleVao = 0;
  GLuint quadVbo = 0;
  GLuint instanceVbo = 0;
  GLuint fullscreenVao = 0;
  GLuint glassVao = 0;
  GLuint glassVbo = 0;
  GLuint glassEbo = 0;
  GLsizei glassIndexCount = 0;

  GLuint depthFbo = 0;
  GLuint blurFbo[2]{};
  GLuint thicknessFbo = 0;
  GLuint thicknessBlurFbo = 0;
  GLuint compositeFbo = 0;
  GLuint rawSurfaceTexture = 0;
  GLuint blurTexture[2]{};
  GLuint thicknessTexture = 0;
  GLuint thicknessBlurTexture = 0;
  GLuint compositeTexture = 0;
  GLuint presentedTexture = 0;
  GLuint depthStencil = 0;

  int glassVessel = -1;
  double glassHeight = -1.0;
  double glassRadius = -1.0;
  std::vector<InstanceData> instances;

  struct ParticleUniforms {
    GLint model = -1;
    GLint view = -1;
    GLint projection = -1;
    GLint radius = -1;
  } depthUniforms, thicknessUniforms;
  GLint thicknessAmplitudeUniform = -1;
  // Vessel interior, uploaded to the depth program so the rendered surface can
  // never escape the glass. Only the depth pass clips; see kDepthFragment.
  static constexpr int kVesselLevels = 72;
  std::array<float, kVesselLevels> vesselProfile{};
  struct VesselClipUniforms {
    GLint eyeToObject = -1;
    GLint vesselHeight = -1;
    GLint vesselMaxRadius = -1;
    GLint vesselProfile = -1;
  } depthClipUniforms;
  struct BlurUniforms {
    GLint surface = -1;
    GLint inverseProjection = -1;
    GLint viewport = -1;
    GLint projectionScale = -1;
    GLint particleRadius = -1;
  } blurUniforms;
  struct ThicknessBlurUniforms {
    GLint thickness = -1;
    GLint surface = -1;
    GLint viewport = -1;
    GLint projectionScale = -1;
    GLint particleRadius = -1;
    GLint direction = -1;
  } thicknessBlurUniforms;
  struct ShadeUniforms {
    GLint surface = -1;
    GLint thickness = -1;
    GLint inverseProjection = -1;
    GLint texel = -1;
    GLint viewport = -1;
    GLint projectionScale = -1;
    GLint particleRadius = -1;
    GLint phaseA = -1;
    GLint phaseB = -1;
    GLint absorptionA = -1;
    GLint absorptionB = -1;
    GLint background = -1;
    GLint absorptionScale = -1;
    GLint exposure = -1;
    GLint showInterface = -1;
  } shadeUniforms;
  struct GlassUniforms {
    GLint model = -1;
    GLint view = -1;
    GLint projection = -1;
    GLint tint = -1;
    GLint backFace = -1;
  } glassUniforms;

  void setFailure(std::string message) {
    failure = std::move(message);
    valid = false;
  }

  void release() {
    if (!glLoaded()) {
      depthProgram = blurProgram = thicknessProgram = shadeProgram = glassProgram = 0;
      particleVao = fullscreenVao = glassVao = 0;
      quadVbo = instanceVbo = glassVbo = glassEbo = 0;
      depthFbo = thicknessFbo = thicknessBlurFbo = compositeFbo = 0;
      blurFbo[0] = blurFbo[1] = 0;
      rawSurfaceTexture = thicknessTexture = thicknessBlurTexture = compositeTexture = 0;
      blurTexture[0] = blurTexture[1] = 0;
      depthStencil = 0;
      valid = false;
      targetWidth = targetHeight = 0;
      presentedTexture = 0;
      return;
    }

    const GLuint programs[] = {depthProgram,        blurProgram,  thicknessProgram,
                               thicknessBlurProgram, shadeProgram, glassProgram};
    for (GLuint program : programs)
      if (program != 0) glDeleteProgram(program);
    const GLuint buffers[] = {quadVbo, instanceVbo, glassVbo, glassEbo};
    glDeleteBuffers(4, buffers);
    const GLuint vaos[] = {particleVao, fullscreenVao, glassVao};
    glDeleteVertexArrays(3, vaos);
    const GLuint textures[] = {rawSurfaceTexture, blurTexture[0],       blurTexture[1],
                               thicknessTexture,  thicknessBlurTexture, compositeTexture};
    glDeleteTextures(6, textures);
    const GLuint framebuffers[] = {depthFbo,         blurFbo[0],       blurFbo[1],
                                   thicknessFbo,     thicknessBlurFbo, compositeFbo};
    glDeleteFramebuffers(6, framebuffers);
    if (depthStencil != 0) glDeleteRenderbuffers(1, &depthStencil);

    depthProgram = blurProgram = thicknessProgram = thicknessBlurProgram = shadeProgram =
        glassProgram = 0;
    particleVao = fullscreenVao = glassVao = 0;
    quadVbo = instanceVbo = glassVbo = glassEbo = 0;
    depthFbo = thicknessFbo = thicknessBlurFbo = compositeFbo = 0;
    blurFbo[0] = blurFbo[1] = 0;
    rawSurfaceTexture = thicknessTexture = thicknessBlurTexture = compositeTexture = 0;
    blurTexture[0] = blurTexture[1] = 0;
    depthStencil = 0;
    glassIndexCount = 0;
    targetWidth = targetHeight = 0;
    presentedTexture = 0;
    glassVessel = -1;
    valid = false;
  }

  void queryUniforms() {
    depthUniforms = {glGetUniformLocation(depthProgram, "uModel"),
                     glGetUniformLocation(depthProgram, "uView"),
                     glGetUniformLocation(depthProgram, "uProjection"),
                     glGetUniformLocation(depthProgram, "uRadius")};
    depthClipUniforms = {
        glGetUniformLocation(depthProgram, "uEyeToObject"),
        glGetUniformLocation(depthProgram, "uVesselHeight"),
        glGetUniformLocation(depthProgram, "uVesselMaxRadius"),
        glGetUniformLocation(depthProgram, "uVesselProfile")};
    thicknessUniforms = {glGetUniformLocation(thicknessProgram, "uModel"),
                         glGetUniformLocation(thicknessProgram, "uView"),
                         glGetUniformLocation(thicknessProgram, "uProjection"),
                         glGetUniformLocation(thicknessProgram, "uRadius")};
    thicknessAmplitudeUniform = glGetUniformLocation(thicknessProgram, "uAmplitude");
    thicknessBlurUniforms = {
        glGetUniformLocation(thicknessBlurProgram, "uThickness"),
        glGetUniformLocation(thicknessBlurProgram, "uSurface"),
        glGetUniformLocation(thicknessBlurProgram, "uViewport"),
        glGetUniformLocation(thicknessBlurProgram, "uProjectionScale"),
        glGetUniformLocation(thicknessBlurProgram, "uParticleRadius"),
        glGetUniformLocation(thicknessBlurProgram, "uDirection")};
    blurUniforms = {glGetUniformLocation(blurProgram, "uSurface"),
                    glGetUniformLocation(blurProgram, "uInverseProjection"),
                    glGetUniformLocation(blurProgram, "uViewport"),
                    glGetUniformLocation(blurProgram, "uProjectionScale"),
                    glGetUniformLocation(blurProgram, "uParticleRadius")};
    shadeUniforms = {
        glGetUniformLocation(shadeProgram, "uSurface"),
        glGetUniformLocation(shadeProgram, "uThickness"),
        glGetUniformLocation(shadeProgram, "uInverseProjection"),
        glGetUniformLocation(shadeProgram, "uTexel"),
        glGetUniformLocation(shadeProgram, "uViewport"),
        glGetUniformLocation(shadeProgram, "uProjectionScale"),
        glGetUniformLocation(shadeProgram, "uParticleRadius"),
        glGetUniformLocation(shadeProgram, "uPhaseA"),
        glGetUniformLocation(shadeProgram, "uPhaseB"),
        glGetUniformLocation(shadeProgram, "uAbsorptionA"),
        glGetUniformLocation(shadeProgram, "uAbsorptionB"),
        glGetUniformLocation(shadeProgram, "uBackground"),
        glGetUniformLocation(shadeProgram, "uAbsorptionScale"),
        glGetUniformLocation(shadeProgram, "uExposure"),
        glGetUniformLocation(shadeProgram, "uShowInterface")};
    glassUniforms = {glGetUniformLocation(glassProgram, "uModel"),
                     glGetUniformLocation(glassProgram, "uView"),
                     glGetUniformLocation(glassProgram, "uProjection"),
                     glGetUniformLocation(glassProgram, "uTint"),
                     glGetUniformLocation(glassProgram, "uBackFace")};
  }

  bool createPrograms() {
    depthProgram = linkProgram(kParticleVertex, kDepthFragment, "fluid depth", failure);
    if (depthProgram == 0) return false;
    blurProgram = linkProgram(kFullscreenVertex, kBlurFragment, "fluid smoothing", failure);
    if (blurProgram == 0) return false;
    thicknessProgram =
        linkProgram(kParticleVertex, kThicknessFragment, "fluid thickness", failure);
    if (thicknessProgram == 0) return false;
    thicknessBlurProgram = linkProgram(kFullscreenVertex, kThicknessBlurFragment,
                                       "fluid thickness smoothing", failure);
    if (thicknessBlurProgram == 0) return false;
    shadeProgram = linkProgram(kFullscreenVertex, kShadeFragment, "fluid shading", failure);
    if (shadeProgram == 0) return false;
    glassProgram = linkProgram(kGlassVertex, kGlassFragment, "vessel glass", failure);
    if (glassProgram == 0) return false;
    queryUniforms();
    return true;
  }

  void createGeometry() {
    static constexpr float quad[] = {-1.0f, -1.0f, 1.0f, -1.0f, 1.0f, 1.0f,
                                     -1.0f, -1.0f, 1.0f, 1.0f, -1.0f, 1.0f};
    glGenVertexArrays(1, &particleVao);
    glGenBuffers(1, &quadVbo);
    glGenBuffers(1, &instanceVbo);
    glBindVertexArray(particleVao);
    glBindBuffer(GL_ARRAY_BUFFER, quadVbo);
    glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(sizeof(quad)), quad, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * static_cast<GLsizei>(sizeof(float)),
                          nullptr);
    glBindBuffer(GL_ARRAY_BUFFER, instanceVbo);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, static_cast<GLsizei>(sizeof(InstanceData)),
                          reinterpret_cast<const void*>(offsetof(InstanceData, x)));
    glVertexAttribDivisor(1, 1);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, static_cast<GLsizei>(sizeof(InstanceData)),
                          reinterpret_cast<const void*>(offsetof(InstanceData, phase)));
    glVertexAttribDivisor(2, 1);
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, static_cast<GLsizei>(sizeof(InstanceData)),
                          reinterpret_cast<const void*>(offsetof(InstanceData, indicator)));
    glVertexAttribDivisor(3, 1);
    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, static_cast<GLsizei>(sizeof(InstanceData)),
                          reinterpret_cast<const void*>(offsetof(InstanceData, speed)));
    glVertexAttribDivisor(4, 1);

    glGenVertexArrays(1, &fullscreenVao);
    glGenVertexArrays(1, &glassVao);
    glGenBuffers(1, &glassVbo);
    glGenBuffers(1, &glassEbo);
  }

  void createTargets() {
    glGenFramebuffers(1, &depthFbo);
    glGenFramebuffers(2, blurFbo);
    glGenFramebuffers(1, &thicknessFbo);
    glGenFramebuffers(1, &thicknessBlurFbo);
    glGenFramebuffers(1, &compositeFbo);
    glGenTextures(1, &rawSurfaceTexture);
    glGenTextures(2, blurTexture);
    glGenTextures(1, &thicknessTexture);
    glGenTextures(1, &thicknessBlurTexture);
    glGenTextures(1, &compositeTexture);
    glGenRenderbuffers(1, &depthStencil);
  }

  void allocateTexture(GLuint texture, GLint internalFormat, GLenum format, GLenum type,
                       GLint filtering, int width, int height) {
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, filtering);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, filtering);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, width, height, 0, format, type, nullptr);
  }

  bool checkFramebuffer(const char* label) {
    const GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    if (status == GL_FRAMEBUFFER_COMPLETE) return true;
    std::ostringstream stream;
    stream << label << " framebuffer is incomplete (0x" << std::hex << status << ')';
    setFailure(stream.str());
    return false;
  }

  bool resizeTargets(int width, int height) {
    if (width == targetWidth && height == targetHeight) return true;

    for (int i = 0; i < 32 && glGetError() != GL_NO_ERROR; ++i) {}
    allocateTexture(rawSurfaceTexture, static_cast<GLint>(kRgba32f), GL_RGBA, GL_FLOAT,
                    GL_NEAREST, width, height);
    // Depth/phase fields use exact texel fetches: interpolation across the
    // zero-depth silhouette would create geometry that no particle contributed.
    allocateTexture(blurTexture[0], static_cast<GLint>(kRgba32f), GL_RGBA, GL_FLOAT,
                    GL_NEAREST, width, height);
    allocateTexture(blurTexture[1], static_cast<GLint>(kRgba32f), GL_RGBA, GL_FLOAT,
                    GL_NEAREST, width, height);
    allocateTexture(thicknessTexture, static_cast<GLint>(GL_RG16F), GL_RG, GL_FLOAT, GL_LINEAR,
                    width, height);
    allocateTexture(thicknessBlurTexture, static_cast<GLint>(GL_RG16F), GL_RG, GL_FLOAT,
                    GL_LINEAR, width, height);
    allocateTexture(compositeTexture, static_cast<GLint>(GL_RGBA8), GL_RGBA, GL_UNSIGNED_BYTE,
                    GL_LINEAR, width, height);

    glBindRenderbuffer(GL_RENDERBUFFER, depthStencil);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, width, height);

    const GLenum colourAttachment = GL_COLOR_ATTACHMENT0;
    glBindFramebuffer(GL_FRAMEBUFFER, depthFbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                           rawSurfaceTexture, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER,
                              depthStencil);
    glDrawBuffers(1, &colourAttachment);
    if (!checkFramebuffer("depth")) return false;

    for (int i = 0; i < 2; ++i) {
      glBindFramebuffer(GL_FRAMEBUFFER, blurFbo[i]);
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                             blurTexture[i], 0);
      glDrawBuffers(1, &colourAttachment);
      if (!checkFramebuffer("smoothing")) return false;
    }

    glBindFramebuffer(GL_FRAMEBUFFER, thicknessFbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, thicknessTexture,
                           0);
    glDrawBuffers(1, &colourAttachment);
    if (!checkFramebuffer("thickness")) return false;

    glBindFramebuffer(GL_FRAMEBUFFER, thicknessBlurFbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                           thicknessBlurTexture, 0);
    glDrawBuffers(1, &colourAttachment);
    if (!checkFramebuffer("thickness smoothing")) return false;

    glBindFramebuffer(GL_FRAMEBUFFER, compositeFbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, compositeTexture,
                           0);
    glDrawBuffers(1, &colourAttachment);
    if (!checkFramebuffer("composite")) return false;

    const GLenum errorCode = glGetError();
    if (errorCode != GL_NO_ERROR) {
      setFailure("Render target allocation failed: " + glErrorMessage(errorCode));
      return false;
    }
    presentedTexture = 0;
    targetWidth = width;
    targetHeight = height;
    return true;
  }

  void rebuildGlass(const fluid::Snapshot& snapshot) {
    const int vessel = static_cast<int>(snapshot.vessel);
    if (vessel == glassVessel && snapshot.vesselHeightM == glassHeight &&
        snapshot.maxRadiusM == glassRadius)
      return;

    constexpr int levelCount = kVesselLevels;
    constexpr int sliceCount = 96;
    constexpr double pi = 3.14159265358979323846;
    std::array<double, levelCount> radii{};
    for (int level = 0; level < levelCount; ++level) {
      const double fraction = static_cast<double>(level) / static_cast<double>(levelCount - 1);
      const double width = sol::vesselWidthAt(snapshot.vessel, fraction);
      // The depth pass clips against the same profile the glass is revolved
      // from, kept normalised so it survives a vessel-height change untouched.
      vesselProfile[level] = static_cast<float>(width);
      radii[level] = snapshot.maxRadiusM * width;
    }

    std::vector<GlassVertex> vertices;
    vertices.reserve(levelCount * sliceCount);
    for (int level = 0; level < levelCount; ++level) {
      const double fraction = static_cast<double>(level) / static_cast<double>(levelCount - 1);
      const double z = snapshot.vesselHeightM * fraction;
      const int lower = std::max(0, level - 1);
      const int upper = std::min(levelCount - 1, level + 1);
      const double dz = snapshot.vesselHeightM * static_cast<double>(upper - lower) /
                        static_cast<double>(levelCount - 1);
      const double slope = dz > 0.0 ? (radii[upper] - radii[lower]) / dz : 0.0;
      const double normalScale = 1.0 / std::sqrt(1.0 + slope * slope);
      for (int slice = 0; slice < sliceCount; ++slice) {
        const double angle = 2.0 * pi * static_cast<double>(slice) /
                             static_cast<double>(sliceCount);
        const double cosine = std::cos(angle);
        const double sine = std::sin(angle);
        vertices.push_back({static_cast<float>(radii[level] * cosine),
                            static_cast<float>(radii[level] * sine), static_cast<float>(z),
                            static_cast<float>(cosine * normalScale),
                            static_cast<float>(sine * normalScale),
                            static_cast<float>(-slope * normalScale)});
      }
    }

    std::vector<GLuint> indices;
    indices.reserve((levelCount - 1) * sliceCount * 6);
    for (int level = 0; level + 1 < levelCount; ++level) {
      for (int slice = 0; slice < sliceCount; ++slice) {
        const int next = (slice + 1) % sliceCount;
        const GLuint a = static_cast<GLuint>(level * sliceCount + slice);
        const GLuint b = static_cast<GLuint>(level * sliceCount + next);
        const GLuint c = static_cast<GLuint>((level + 1) * sliceCount + slice);
        const GLuint d = static_cast<GLuint>((level + 1) * sliceCount + next);
        indices.insert(indices.end(), {a, b, c, b, d, c});
      }
    }

    glBindVertexArray(glassVao);
    glBindBuffer(GL_ARRAY_BUFFER, glassVbo);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(vertices.size() * sizeof(GlassVertex)), vertices.data(),
                 GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, glassEbo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(indices.size() * sizeof(GLuint)), indices.data(),
                 GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, static_cast<GLsizei>(sizeof(GlassVertex)),
                          reinterpret_cast<const void*>(offsetof(GlassVertex, x)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, static_cast<GLsizei>(sizeof(GlassVertex)),
                          reinterpret_cast<const void*>(offsetof(GlassVertex, nx)));
    glassIndexCount = static_cast<GLsizei>(indices.size());
    glassVessel = vessel;
    glassHeight = snapshot.vesselHeightM;
    glassRadius = snapshot.maxRadiusM;
  }

  void uploadInstances(const fluid::Snapshot& snapshot) {
    const std::size_t count = std::min(
        {snapshot.px.size(), snapshot.py.size(), snapshot.pz.size(), snapshot.phase.size()});
    instances.resize(count);
    const bool haveIndicator = snapshot.colour.size() >= count;
    const bool haveSpeed = snapshot.speed.size() >= count;
    for (std::size_t i = 0; i < count; ++i) {
      const float phase = snapshot.phase[i] == 0 ? 0.0f : 1.0f;
      const float speed =
          haveSpeed ? std::max(0.0f, static_cast<float>(snapshot.speed[i])) : 0.0f;
      instances[i] = {snapshot.px[i], snapshot.py[i], snapshot.pz[i], phase,
                      haveIndicator ? snapshot.colour[i] : phase, speed};
    }
    glBindBuffer(GL_ARRAY_BUFFER, instanceVbo);
    glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(instances.size() * sizeof(InstanceData)),
                 instances.empty() ? nullptr : instances.data(), GL_STREAM_DRAW);
  }

  void setParticleUniforms(const ParticleUniforms& uniforms, const Mat4& model,
                           const Mat4& view, const Mat4& projection, float radius) {
    glUniformMatrix4fv(uniforms.model, 1, GL_FALSE, model.data());
    glUniformMatrix4fv(uniforms.view, 1, GL_FALSE, view.data());
    glUniformMatrix4fv(uniforms.projection, 1, GL_FALSE, projection.data());
    glUniform1f(uniforms.radius, radius);
  }

  GLuint drawFrame(const fluid::Snapshot& snapshot, const Camera3D& camera, int width, int height,
                   const FluidRenderSettings& settings) {
    const float aspect = static_cast<float>(width) / static_cast<float>(height);
    const Mat4 model = poseMatrix(snapshot.pose);
    const Mat4 view = camera.view();
    const Mat4 projection = camera.projection(aspect);
    const Mat4 inverseProjection = inverse(projection);
    const float radius = std::max(1.0e-5f, static_cast<float>(snapshot.particleRadiusM));
    // The thickness pass measures optical path, so its footprint is calibrated
    // by volume rather than by silhouette. Each particle stands for one lattice
    // cell, dx^3 of liquid, and its (1-r^2)^(3/2) footprint of radius R
    // integrates to A * pi * R^2 / 2.5. Setting that equal to dx^3 gives the
    // amplitude below, so the summed field is a true path length in metres no
    // matter what radius the surface is drawn at. Getting this wrong is
    // visible: the old dx/2 sphere chord under-reported the path by 1/0.524
    // and made 100 mL of dichloromethane look like haze.
    const double spacing = std::max(1.0e-6, snapshot.particleSpacingM);
    constexpr double pi = 3.14159265358979323846;
    const float thicknessAmplitude = static_cast<float>(
        2.5 * spacing * spacing * spacing / (pi * static_cast<double>(radius) * radius));
    const GLsizei particleCount = static_cast<GLsizei>(instances.size());
    const GLenum colourAttachment = GL_COLOR_ATTACHMENT0;
    const GLfloat zero[4] = {0.0f, 0.0f, 0.0f, 0.0f};

    glViewport(0, 0, width, height);

    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glDisable(GL_SCISSOR_TEST);
    glDisable(GL_STENCIL_TEST);
    glDisable(GL_FRAMEBUFFER_SRGB);
    // Pass 1: sphere-impostor surface depth and phase indicator.
    glBindFramebuffer(GL_FRAMEBUFFER, depthFbo);
    glDrawBuffers(1, &colourAttachment);
    glClearBufferfv(GL_COLOR, 0, zero);
    const GLfloat farDepth = 1.0f;
    glClearBufferfv(GL_DEPTH, 0, &farDepth);
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glDepthMask(GL_TRUE);
    glDisable(GL_CULL_FACE);
    glUseProgram(depthProgram);
    setParticleUniforms(depthUniforms, model, view, projection, radius);
    // The clip needs the eye-space impostor offset back in vessel coordinates.
    // model and view are both rigid (a normalised quaternion and a look-at), so
    // the inverse of the 3x3 of view*model is simply its transpose. Mat4 is
    // column-major, so element (row, col) lives at [col * 4 + row] and the
    // transpose is written out directly.
    float eyeToObject[9];
    for (int row = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col) {
        float product = 0.0f;
        for (int k = 0; k < 3; ++k) product += view[k * 4 + row] * model[col * 4 + k];
        eyeToObject[row * 3 + col] = product;  // transposed on write
      }
    }
    glUniformMatrix3fv(depthClipUniforms.eyeToObject, 1, GL_FALSE, eyeToObject);
    glUniform1f(depthClipUniforms.vesselHeight,
                static_cast<float>(snapshot.vesselHeightM));
    glUniform1f(depthClipUniforms.vesselMaxRadius,
                static_cast<float>(snapshot.maxRadiusM));
    glUniform1fv(depthClipUniforms.vesselProfile, kVesselLevels, vesselProfile.data());
    glBindVertexArray(particleVao);
    glDrawArraysInstanced(GL_TRIANGLES, 0, 6, particleCount);

    // Pass 2: adaptive mean-curvature flow, ping-ponged once per stable PDE step.
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    glDisable(GL_BLEND);
    glUseProgram(blurProgram);
    glUniform1i(blurUniforms.surface, 0);
    glUniformMatrix4fv(blurUniforms.inverseProjection, 1, GL_FALSE,
                       inverseProjection.data());
    glUniform2f(blurUniforms.viewport, static_cast<float>(width),
                static_cast<float>(height));
    glUniform2f(blurUniforms.projectionScale, projection[0], projection[5]);
    glUniform1f(blurUniforms.particleRadius, radius);
    glBindVertexArray(fullscreenVao);
    GLuint smoothedSurface = rawSurfaceTexture;
    const int requestedSmoothingIterations =
        settings.showParticles
            ? 0
            : std::clamp(
                  static_cast<int>(std::lround(settings.smoothingIterations)), 0, 12);
    const int curvatureSteps = 2 * requestedSmoothingIterations;
    for (int step = 0; step < curvatureSteps; ++step) {
      const int destination = step & 1;
      glBindFramebuffer(GL_FRAMEBUFFER, blurFbo[destination]);
      glDrawBuffers(1, &colourAttachment);
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, smoothedSurface);
      glDrawArrays(GL_TRIANGLES, 0, 3);
      smoothedSurface = blurTexture[destination];
    }

    // Pass 3: per-phase Beer-Lambert optical path, accumulated additively.
    glBindFramebuffer(GL_FRAMEBUFFER, thicknessFbo);
    glDrawBuffers(1, &colourAttachment);
    glClearBufferfv(GL_COLOR, 0, zero);
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    glEnable(GL_BLEND);
    glBlendEquationSeparate(GL_FUNC_ADD, GL_FUNC_ADD);
    glBlendFuncSeparate(GL_ONE, GL_ONE, GL_ONE, GL_ONE);
    glUseProgram(thicknessProgram);
    // The footprint spans the same radius as the surface so the discs overlap;
    // the amplitude, not the radius, carries the volume calibration.
    setParticleUniforms(thicknessUniforms, model, view, projection, radius);
    glUniform1f(thicknessAmplitudeUniform, thicknessAmplitude);
    glBindVertexArray(particleVao);
    glDrawArraysInstanced(GL_TRIANGLES, 0, 6, particleCount);

    // Pass 3b: separable Gaussian over the optical path. The surface can be
    // perfectly smooth and the body still look scaly, because the summed
    // footprints carry a particle-scale ripple that Beer-Lambert makes visible.
    glDisable(GL_BLEND);
    glUseProgram(thicknessBlurProgram);
    glUniform1i(thicknessBlurUniforms.thickness, 0);
    glUniform1i(thicknessBlurUniforms.surface, 1);
    glUniform2f(thicknessBlurUniforms.viewport, static_cast<float>(width),
                static_cast<float>(height));
    glUniform2f(thicknessBlurUniforms.projectionScale, projection[0], projection[5]);
    glUniform1f(thicknessBlurUniforms.particleRadius, radius);
    glActiveTexture(GL_TEXTURE0 + 1);
    glBindTexture(GL_TEXTURE_2D, smoothedSurface);
    glBindVertexArray(fullscreenVao);
    const GLuint thicknessPingPong[2] = {thicknessBlurFbo, thicknessFbo};
    const GLuint thicknessSource[2] = {thicknessTexture, thicknessBlurTexture};
    GLuint filteredThickness = thicknessTexture;
    for (int axis = 0; axis < 2; ++axis) {
      glBindFramebuffer(GL_FRAMEBUFFER, thicknessPingPong[axis]);
      glDrawBuffers(1, &colourAttachment);
      glUniform2f(thicknessBlurUniforms.direction, axis == 0 ? 1.0f : 0.0f,
                  axis == 0 ? 0.0f : 1.0f);
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, thicknessSource[axis]);
      glDrawArrays(GL_TRIANGLES, 0, 3);
      filteredThickness = axis == 0 ? thicknessBlurTexture : thicknessTexture;
    }

    // Pass 4: reconstruct the smooth surface, shade, and write straight RGBA.
    glBindFramebuffer(GL_FRAMEBUFFER, compositeFbo);
    glDrawBuffers(1, &colourAttachment);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glDisable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    glUseProgram(shadeProgram);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, smoothedSurface);
    glActiveTexture(GL_TEXTURE0 + 1);
    glBindTexture(GL_TEXTURE_2D, filteredThickness);
    glUniform1i(shadeUniforms.surface, 0);
    glUniform1i(shadeUniforms.thickness, 1);
    glUniformMatrix4fv(shadeUniforms.inverseProjection, 1, GL_FALSE, inverseProjection.data());
    glUniform2f(shadeUniforms.texel, 1.0f / static_cast<float>(width),
                1.0f / static_cast<float>(height));
    glUniform2f(shadeUniforms.viewport, static_cast<float>(width),
                static_cast<float>(height));
    glUniform2f(shadeUniforms.projectionScale, projection[0], projection[5]);
    glUniform1f(shadeUniforms.particleRadius, radius);

    const auto phaseA = phaseColour(snapshot, 0, {0.25f, 0.55f, 0.90f, 0.82f});
    const auto phaseB = phaseColour(snapshot, 1, {0.92f, 0.58f, 0.18f, 0.82f});
    const auto absorptionA = absorptionFromColour(phaseA);
    const auto absorptionB = absorptionFromColour(phaseB);
    glUniform4f(shadeUniforms.phaseA, phaseA[0], phaseA[1], phaseA[2], phaseA[3]);
    glUniform4f(shadeUniforms.phaseB, phaseB[0], phaseB[1], phaseB[2], phaseB[3]);
    glUniform3f(shadeUniforms.absorptionA, absorptionA[0], absorptionA[1], absorptionA[2]);
    glUniform3f(shadeUniforms.absorptionB, absorptionB[0], absorptionB[1], absorptionB[2]);
    glUniform3f(shadeUniforms.background, settings.backgroundTint[0],
                settings.backgroundTint[1], settings.backgroundTint[2]);
    glUniform1f(shadeUniforms.absorptionScale, std::max(0.0f, settings.absorptionScale));
    glUniform1f(shadeUniforms.exposure, std::max(0.0f, settings.exposure));
    glUniform1f(shadeUniforms.showInterface, settings.showInterface ? 1.0f : 0.0f);
    glBindVertexArray(fullscreenVao);
    glDrawArrays(GL_TRIANGLES, 0, 3);

    // Pass 5: thin glass shell, rear wall first and front wall second.
    if (settings.showGlass && glassIndexCount > 0) {
      glEnable(GL_BLEND);
      glBlendEquationSeparate(GL_FUNC_ADD, GL_FUNC_ADD);
      // The composite target stores straight colour because ImGui supplies the
      // final source-alpha blend; use the matching over operator for glass too.
      glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE,
                          GL_ONE_MINUS_SRC_ALPHA);
      glDisable(GL_DEPTH_TEST);
      glDepthMask(GL_FALSE);
      glEnable(GL_CULL_FACE);
      glFrontFace(GL_CCW);
      glUseProgram(glassProgram);
      glUniformMatrix4fv(glassUniforms.model, 1, GL_FALSE, model.data());
      glUniformMatrix4fv(glassUniforms.view, 1, GL_FALSE, view.data());
      glUniformMatrix4fv(glassUniforms.projection, 1, GL_FALSE, projection.data());
      glUniform3f(glassUniforms.tint, 0.72f, 0.88f, 1.0f);
      glBindVertexArray(glassVao);
      glCullFace(GL_FRONT);
      glUniform1f(glassUniforms.backFace, 1.0f);
      glDrawElements(GL_TRIANGLES, glassIndexCount, GL_UNSIGNED_INT, nullptr);
      glCullFace(GL_BACK);
      glUniform1f(glassUniforms.backFace, 0.0f);
      glDrawElements(GL_TRIANGLES, glassIndexCount, GL_UNSIGNED_INT, nullptr);
    }

    return compositeTexture;
  }
};

FluidRenderer::FluidRenderer() : impl_(new Impl) {}

FluidRenderer::~FluidRenderer() {
  shutdown();
  delete impl_;
}

bool FluidRenderer::initialise() {
  if (!glLoaded()) {
    impl_->setFailure("OpenGL 3.3 entry points have not been loaded");
    return false;
  }

  GlState savedState;
  impl_->release();
  impl_->failure.clear();
  if (!impl_->createPrograms()) {
    impl_->release();
    return false;
  }
  impl_->createGeometry();
  impl_->createTargets();
  const GLenum errorCode = glGetError();
  if (errorCode != GL_NO_ERROR) {
    impl_->setFailure("OpenGL renderer initialisation failed: " + glErrorMessage(errorCode));
    impl_->release();
    return false;
  }
  impl_->valid = true;
  return true;
}

bool FluidRenderer::ready() const { return impl_->valid; }

const std::string& FluidRenderer::error() const { return impl_->failure; }

void FluidRenderer::shutdown() { impl_->release(); }

std::uint32_t FluidRenderer::render(const fluid::Snapshot& snapshot, const Camera3D& camera,
                                    int width, int height,
                                    const FluidRenderSettings& settings) {
  if (!impl_->valid || width <= 0 || height <= 0) return 0;
  const auto started = std::chrono::steady_clock::now();
  const auto finishFailure = [&]() {
    impl_->frameMilliseconds =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - started)
            .count();
    return std::uint32_t{0};
  };

  try {
    for (int i = 0; i < 32 && glGetError() != GL_NO_ERROR; ++i) {}
    GlState savedState;
    if (!impl_->resizeTargets(width, height)) return finishFailure();

    impl_->uploadInstances(snapshot);
    impl_->rebuildGlass(snapshot);
    const GLuint texture = impl_->drawFrame(snapshot, camera, width, height, settings);
    const GLenum errorCode = glGetError();
    if (errorCode != GL_NO_ERROR) {
      impl_->setFailure("Fluid rendering failed: " + glErrorMessage(errorCode));
      return finishFailure();
    }

    impl_->frameMilliseconds =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - started)
            .count();
    impl_->presentedTexture = texture;
    return texture;
  } catch (const std::exception& exception) {
    // Snapshot upload and glass tessellation allocate; converting any failure to
    // renderer state keeps exceptions from escaping the application's UI seam.
    impl_->setFailure(std::string("Fluid rendering failed: ") + exception.what());
    return finishFailure();
  } catch (...) {
    impl_->setFailure("Fluid rendering failed with an unknown error");
    return finishFailure();
  }
}

std::uint32_t FluidRenderer::colourTexture() const {
  return impl_->valid ? impl_->presentedTexture : 0;
}

int FluidRenderer::width() const { return impl_->targetWidth; }

int FluidRenderer::height() const { return impl_->targetHeight; }

double FluidRenderer::lastFrameMs() const { return impl_->frameMilliseconds; }

}  // namespace chemcad::gfx

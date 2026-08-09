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
// GL_RGBA32F is core in OpenGL 3.0; keep the token local because gl_api.hpp
// exposes only the subset shared by the rest of the application.
constexpr GLenum kRgba32f = 0x8814;

struct InstanceData {
  float x;
  float y;
  float z;
  float phase;
  float indicator;
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

uniform mat4 uModel;
uniform mat4 uView;
uniform mat4 uProjection;
uniform float uRadius;

out vec2 vCorner;
out vec3 vCentreEye;
flat out float vPhase;
flat out float vIndicator;

void main() {
  vec3 centreEye = (uView * uModel * vec4(aCentre, 1.0)).xyz;
  vec3 cornerEye = centreEye + vec3(aCorner * uRadius, 0.0);
  gl_Position = uProjection * vec4(cornerEye, 1.0);
  vCorner = aCorner;
  vCentreEye = centreEye;
  vPhase = aPhase;
  vIndicator = aIndicator;
}
)GLSL";

const char* kDepthFragment = R"GLSL(#version 330 core
in vec2 vCorner;
in vec3 vCentreEye;
flat in float vPhase;
flat in float vIndicator;

uniform mat4 uProjection;
uniform float uRadius;

layout(location = 0) out vec3 outSurface;

void main() {
  float radius2 = dot(vCorner, vCorner);
  if (radius2 > 1.0) discard;

  // The quad is only a carrier. Reconstructing the positive-z hemisphere in
  // eye space makes each instance a true sphere impostor, including correct
  // per-fragment hardware depth for overlap between particles.
  vec3 sphereNormal = vec3(vCorner, sqrt(max(0.0, 1.0 - radius2)));
  vec3 surfaceEye = vCentreEye + uRadius * sphereNormal;
  vec4 clip = uProjection * vec4(surfaceEye, 1.0);
  gl_FragDepth = 0.5 * (clip.z / clip.w) + 0.5;
  // Depth, discrete particle phase, and the solver-smoothed colour indicator
  // must undergo identical filtering or the A/B boundary drifts off the surface.
  outSurface = vec3(-surfaceEye.z, step(0.5, vPhase),
                    clamp(vIndicator, 0.0, 1.0));
}
)GLSL";

const char* kThicknessFragment = R"GLSL(#version 330 core
in vec2 vCorner;
flat in float vPhase;
uniform float uRadius;
layout(location = 0) out vec2 outThickness;

void main() {
  float radius2 = dot(vCorner, vCorner);
  if (radius2 > 1.0) discard;

  // A ray through a sphere at normalised radius r travels 2R sqrt(1-r^2).
  // Additive blending sums that optical path independently for the two phases.
  float chord = 2.0 * uRadius * sqrt(max(0.0, 1.0 - radius2));
  float phaseB = step(0.5, vPhase);
  outThickness = chord * vec2(1.0 - phaseB, phaseB);
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
uniform vec2 uDirection;
uniform vec2 uTexel;
layout(location = 0) out vec3 outSurface;

void main() {
  vec3 centre = texture(uSurface, vUv).rgb;
  if (centre.x <= 0.0) {
    outSurface = vec3(0.0);
    return;
  }

  // Fixed support and range rejection are the separable bilateral form of the
  // screen-space smoothing in van der Laan, Green & Sainz, I3D 2009. Applying
  // the same weights to depth and colour keeps the interface on that surface.
  const float spatial[5] = float[](1.0, 0.8825, 0.6065, 0.3247, 0.1353);
  float threshold = max(0.0025, 0.015 * centre.x);
  vec3 sum = centre * spatial[0];
  float weight = spatial[0];
  for (int i = 1; i <= 4; ++i) {
    for (int signIndex = 0; signIndex < 2; ++signIndex) {
      float directionSign = signIndex == 0 ? -1.0 : 1.0;
      vec3 sampleSurface =
          texture(uSurface, vUv + directionSign * float(i) * uDirection * uTexel).rgb;
      float difference = abs(sampleSurface.x - centre.x);
      if (sampleSurface.x > 0.0 && difference < threshold) {
        float rangeWeight =
            exp(-difference * difference / max(1.0e-7, 0.18 * threshold * threshold));
        float w = spatial[i] * rangeWeight;
        sum += sampleSurface * w;
        weight += w;
      }
    }
  }
  outSurface = sum / weight;
}
)GLSL";

const char* kShadeFragment = R"GLSL(#version 330 core
in vec2 vUv;
uniform sampler2D uSurface;
uniform sampler2D uThickness;
uniform mat4 uInverseProjection;
uniform vec2 uTexel;
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

void main() {
  vec3 surface = texture(uSurface, vUv).rgb;
  float depth = surface.x;
  if (depth <= 0.0) {
    outColour = vec4(uBackground, 1.0);
    return;
  }

  vec3 surfaceL = texture(uSurface, vUv - vec2(uTexel.x, 0.0)).rgb;
  vec3 surfaceR = texture(uSurface, vUv + vec2(uTexel.x, 0.0)).rgb;
  vec3 surfaceD = texture(uSurface, vUv - vec2(0.0, uTexel.y)).rgb;
  vec3 surfaceU = texture(uSurface, vUv + vec2(0.0, uTexel.y)).rgb;
  surfaceL = surfaceL.x > 0.0 ? surfaceL : surface;
  surfaceR = surfaceR.x > 0.0 ? surfaceR : surface;
  surfaceD = surfaceD.x > 0.0 ? surfaceD : surface;
  surfaceU = surfaceU.x > 0.0 ? surfaceU : surface;

  vec3 positionEye = eyePosition(vUv, depth);
  vec3 tangentX = eyePosition(vUv + vec2(uTexel.x, 0.0), surfaceR.x) -
                  eyePosition(vUv - vec2(uTexel.x, 0.0), surfaceL.x);
  vec3 tangentY = eyePosition(vUv + vec2(0.0, uTexel.y), surfaceU.x) -
                  eyePosition(vUv - vec2(0.0, uTexel.y), surfaceD.x);
  vec3 normalEye = normalize(cross(tangentX, tangentY));
  // Screen-space winding can flip after projection; face the reconstructed
  // normal toward the eye so diffuse and Fresnel terms cannot black out a surface.
  if (normalEye.z < 0.0) normalEye = -normalEye;

  vec3 viewDirection = normalize(-positionEye);
  vec3 lightDirection = normalize(vec3(-0.35, 0.55, 0.76));
  vec3 halfDirection = normalize(lightDirection + viewDirection);
  float diffuse = max(dot(normalEye, lightDirection), 0.0);
  float specular = pow(max(dot(normalEye, halfDirection), 0.0), 72.0);

  // Schlick's approximation preserves the grazing-angle rise of a dielectric
  // interface without pretending that the display colours are optical IOR data.
  float cosTheta = clamp(dot(normalEye, viewDirection), 0.0, 1.0);
  float fresnel = 0.020 + 0.980 * pow(1.0 - cosTheta, 5.0);
  float rim = pow(1.0 - cosTheta, 2.4);

  // Base hue follows the actual surface particle phase; the continuous solver
  // indicator is reserved for locating the interface rather than muddying phases.
  float phase = step(0.5, clamp(surface.y, 0.0, 1.0));
  vec4 liquid = mix(uPhaseA, uPhaseB, phase);
  vec2 thicknessM = max(texture(uThickness, vUv).rg, vec2(0.0));
  // Beer-Lambert T = exp(-mu L): thickness is already in metres and mu is 1/m.
  // Keeping that unit contract explicit avoids the 1000x millimetre black bias.
  vec3 transmission =
      exp(-uAbsorptionScale *
          (uAbsorptionA * thicknessM.x + uAbsorptionB * thicknessM.y));

  vec3 refractedBackground = uBackground * (0.92 + 0.18 * normalEye.x);
  float absorptionCoverage = 1.0 - dot(transmission, vec3(1.0 / 3.0));
  // Transmission remains spectrally per-phase, while a small tinted scattering
  // term preserves configured liquid hue against the dark laboratory backdrop.
  // Its coverage dependence keeps a thin film pale and a 100 mL layer saturated.
  vec3 body = refractedBackground * transmission +
              liquid.rgb * (0.08 + 0.72 * absorptionCoverage);
  body *= 0.34 + 0.66 * diffuse;
  // A diffuse/ambient floor preserves the phase hue even where the key light
  // misses; a readable laboratory liquid must not collapse to silhouette black.
  body = max(body, vec3(0.035) + 0.12 * liquid.rgb);
  body += vec3(1.0) * (0.72 * specular + 0.24 * fresnel);
  body += liquid.rgb * (0.18 * rim);

  float indicator = clamp(surface.z, 0.0, 1.0);
  vec2 indicatorGradient =
      0.5 * vec2(surfaceR.z - surfaceL.z, surfaceU.z - surfaceD.z);
  float gradientStrength = length(indicatorGradient);
  // Dividing distance from 0.5 by the local indicator slope makes a roughly
  // one-pixel band exactly where the smoothed A/B field crosses its interface.
  float interfaceDistance = abs(indicator - 0.5) / max(gradientStrength, 1.0e-4);
  float interfaceBand = (1.0 - smoothstep(0.25, 1.25, interfaceDistance)) *
                        smoothstep(0.004, 0.06, gradientStrength);
  body += uShowInterface * interfaceBand * vec3(0.18, 0.22, 0.28);

  // Exposure is a true linear multiplier around 1.0. The fixed shoulder only
  // protects highlights, then explicit sRGB encoding prevents a linear-RGBA8
  // target from making mid-tones appear much darker than intended.
  vec3 exposed = max(uExposure, 0.0) * max(body, vec3(0.0));
  body = pow(exposed / (vec3(1.0) + exposed), vec3(1.0 / 2.2));

  float alpha = clamp(max(liquid.a * (0.18 + 0.82 * absorptionCoverage),
                          0.08 + 0.30 * fresnel + 0.14 * rim),
                      0.0, 0.985);
  // Resolve the liquid once against the stage background and keep the texture
  // opaque. ImGui would otherwise apply alpha a second time, recreating the
  // near-black bias and making transparent glass rims disappear.
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
  // display RGB as transmittance through a 25 mm reference cuvette. This keeps
  // Beer-Lambert dimensionally correct (mu in 1/m) without inventing chemistry.
  constexpr float referencePathM = 0.025f;
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
  GLuint compositeFbo = 0;
  GLuint rawSurfaceTexture = 0;
  GLuint blurTexture[2]{};
  GLuint thicknessTexture = 0;
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
  struct BlurUniforms {
    GLint surface = -1;
    GLint direction = -1;
    GLint texel = -1;
  } blurUniforms;
  struct ShadeUniforms {
    GLint surface = -1;
    GLint thickness = -1;
    GLint inverseProjection = -1;
    GLint texel = -1;
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
      depthFbo = thicknessFbo = compositeFbo = 0;
      blurFbo[0] = blurFbo[1] = 0;
      rawSurfaceTexture = thicknessTexture = compositeTexture = 0;
      blurTexture[0] = blurTexture[1] = 0;
      depthStencil = 0;
      valid = false;
      targetWidth = targetHeight = 0;
      presentedTexture = 0;
      return;
    }

    const GLuint programs[] = {depthProgram, blurProgram, thicknessProgram, shadeProgram,
                               glassProgram};
    for (GLuint program : programs)
      if (program != 0) glDeleteProgram(program);
    const GLuint buffers[] = {quadVbo, instanceVbo, glassVbo, glassEbo};
    glDeleteBuffers(4, buffers);
    const GLuint vaos[] = {particleVao, fullscreenVao, glassVao};
    glDeleteVertexArrays(3, vaos);
    const GLuint textures[] = {rawSurfaceTexture, blurTexture[0], blurTexture[1],
                               thicknessTexture, compositeTexture};
    glDeleteTextures(5, textures);
    const GLuint framebuffers[] = {depthFbo, blurFbo[0], blurFbo[1], thicknessFbo,
                                   compositeFbo};
    glDeleteFramebuffers(5, framebuffers);
    if (depthStencil != 0) glDeleteRenderbuffers(1, &depthStencil);

    depthProgram = blurProgram = thicknessProgram = shadeProgram = glassProgram = 0;
    particleVao = fullscreenVao = glassVao = 0;
    quadVbo = instanceVbo = glassVbo = glassEbo = 0;
    depthFbo = thicknessFbo = compositeFbo = 0;
    blurFbo[0] = blurFbo[1] = 0;
    rawSurfaceTexture = thicknessTexture = compositeTexture = 0;
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
    thicknessUniforms = {glGetUniformLocation(thicknessProgram, "uModel"),
                         glGetUniformLocation(thicknessProgram, "uView"),
                         glGetUniformLocation(thicknessProgram, "uProjection"),
                         glGetUniformLocation(thicknessProgram, "uRadius")};
    blurUniforms = {glGetUniformLocation(blurProgram, "uSurface"),
                    glGetUniformLocation(blurProgram, "uDirection"),
                    glGetUniformLocation(blurProgram, "uTexel")};
    shadeUniforms = {
        glGetUniformLocation(shadeProgram, "uSurface"),
        glGetUniformLocation(shadeProgram, "uThickness"),
        glGetUniformLocation(shadeProgram, "uInverseProjection"),
        glGetUniformLocation(shadeProgram, "uTexel"),
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

    glGenVertexArrays(1, &fullscreenVao);
    glGenVertexArrays(1, &glassVao);
    glGenBuffers(1, &glassVbo);
    glGenBuffers(1, &glassEbo);
  }

  void createTargets() {
    glGenFramebuffers(1, &depthFbo);
    glGenFramebuffers(2, blurFbo);
    glGenFramebuffers(1, &thicknessFbo);
    glGenFramebuffers(1, &compositeFbo);
    glGenTextures(1, &rawSurfaceTexture);
    glGenTextures(2, blurTexture);
    glGenTextures(1, &thicknessTexture);
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
    allocateTexture(blurTexture[0], static_cast<GLint>(kRgba32f), GL_RGBA, GL_FLOAT, GL_LINEAR,
                    width, height);
    allocateTexture(blurTexture[1], static_cast<GLint>(kRgba32f), GL_RGBA, GL_FLOAT, GL_LINEAR,
                    width, height);
    allocateTexture(thicknessTexture, static_cast<GLint>(GL_RG16F), GL_RG, GL_FLOAT, GL_LINEAR,
                    width, height);
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

    constexpr int levelCount = 72;
    constexpr int sliceCount = 96;
    constexpr double pi = 3.14159265358979323846;
    std::array<double, levelCount> radii{};
    for (int level = 0; level < levelCount; ++level) {
      const double fraction = static_cast<double>(level) / static_cast<double>(levelCount - 1);
      radii[level] = snapshot.maxRadiusM * sol::vesselWidthAt(snapshot.vessel, fraction);
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
    for (std::size_t i = 0; i < count; ++i) {
      const float phase = snapshot.phase[i] == 0 ? 0.0f : 1.0f;
      instances[i] = {snapshot.px[i], snapshot.py[i], snapshot.pz[i], phase,
                      haveIndicator ? snapshot.colour[i] : phase};
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
    glBindVertexArray(particleVao);
    glDrawArraysInstanced(GL_TRIANGLES, 0, 6, particleCount);

    // Pass 2: separable bilateral smoothing, ping-ponged for each iteration.
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    glDisable(GL_BLEND);
    glUseProgram(blurProgram);
    glUniform1i(blurUniforms.surface, 0);
    glUniform2f(blurUniforms.texel, 1.0f / static_cast<float>(width),
                1.0f / static_cast<float>(height));
    glBindVertexArray(fullscreenVao);
    GLuint smoothedSurface = rawSurfaceTexture;
    const int smoothingPasses = settings.showParticles
                                    ? 0
                                    : std::clamp(static_cast<int>(std::lround(
                                                     settings.smoothingIterations)),
                                                 0, 12);
    for (int iteration = 0; iteration < smoothingPasses; ++iteration) {
      glBindFramebuffer(GL_FRAMEBUFFER, blurFbo[0]);
      glDrawBuffers(1, &colourAttachment);
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, smoothedSurface);
      glUniform2f(blurUniforms.direction, 1.0f, 0.0f);
      glDrawArrays(GL_TRIANGLES, 0, 3);

      glBindFramebuffer(GL_FRAMEBUFFER, blurFbo[1]);
      glBindTexture(GL_TEXTURE_2D, blurTexture[0]);
      glUniform2f(blurUniforms.direction, 0.0f, 1.0f);
      glDrawArrays(GL_TRIANGLES, 0, 3);
      smoothedSurface = blurTexture[1];
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
    setParticleUniforms(thicknessUniforms, model, view, projection, radius);
    glBindVertexArray(particleVao);
    glDrawArraysInstanced(GL_TRIANGLES, 0, 6, particleCount);

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
    glBindTexture(GL_TEXTURE_2D, thicknessTexture);
    glUniform1i(shadeUniforms.surface, 0);
    glUniform1i(shadeUniforms.thickness, 1);
    glUniformMatrix4fv(shadeUniforms.inverseProjection, 1, GL_FALSE, inverseProjection.data());
    glUniform2f(shadeUniforms.texel, 1.0f / static_cast<float>(width),
                1.0f / static_cast<float>(height));

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

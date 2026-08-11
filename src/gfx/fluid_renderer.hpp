#pragma once
// Real-time 3D renderer for the multiphase fluid, drawn with OpenGL 3.3 core
// into an off-screen target that the UI composites with ImGui::Image.
//
// Why off-screen: the panel is inside a docked ImGui window, so the fluid has
// to obey that window's rectangle, its clipping and its z-order. Rendering into
// a framebuffer object and handing ImGui a texture achieves that exactly, keeps
// every GL call on the main thread inside the app's render seam, and leaves the
// headless panel tests (which have no GL context at all) able to run the same
// panel code with the renderer absent.
//
// Technique: screen-space fluid rendering (van der Laan, Green & Sainz 2009,
// "Screen space fluid rendering with curvature flow"), which is the standard
// way to get a liquid surface out of particles without meshing them:
//
//   1. sphere-impostor depth pass -- each particle is a camera-facing quad,
//      its fragments discarded outside the disc and their depth pushed to the
//      sphere's surface, so the depth buffer holds the particle surface;
//   2. smoothing pass -- separable bilateral / curvature-flow blur of that
//      depth, which turns a bumpy particle field into a fluid surface without
//      bleeding across silhouettes;
//   3. thickness pass -- additive, no depth test, giving optical path length
//      per phase for Beer-Lambert absorption (this is what makes a 100 mL
//      layer of dichloromethane read as a dense liquid rather than a cloud);
//   4. shading pass -- normals reconstructed from smoothed depth, Blinn-Phong
//      specular, Fresnel reflectance, per-phase absorption from the thickness,
//      and refraction of the background;
//   5. glass pass -- the vessel as a revolved shell of the analytic profile,
//      drawn with front/back faces and a Fresnel rim so the fluid is seen
//      through glassware rather than floating in space.
//
// The renderer owns every GL object it creates and releases them in shutdown();
// it never touches ImGui state, and it restores the GL state it changes so the
// ImGui backend's own draw call is unaffected.

#include <cstdint>
#include <string>

#include "fluid/simulation.hpp"
#include "gfx/camera3d.hpp"
#include "gfx/gl_api.hpp"

namespace chemcad::gfx {

struct FluidRenderSettings {
  bool showGlass = true;
  bool showParticles = false;   // debug: draw raw impostors instead of the surface
  bool showInterface = true;    // tint the A/B interface where the colour field turns
  float exposure = 1.0f;
  float smoothingIterations = 4.0f;
  float absorptionScale = 1.0f;
  float backgroundTint[3] = {0.05f, 0.06f, 0.08f};
  // Opacity of the stage backdrop, i.e. of every texel the apparatus does not
  // cover. 1 is the docked stage, which wants a bench behind the glass. 0 makes
  // the backdrop vanish so only the vessel and its contents are composited,
  // which is what lets the funnel be carried across the application without
  // the stage rectangle blacking out everything it passes over. The liquid
  // still refracts `backgroundTint` either way -- a transparent backdrop is not
  // a transparent liquid.
  float backgroundAlpha = 1.0f;
};

class FluidRenderer {
 public:
  FluidRenderer();
  ~FluidRenderer();
  FluidRenderer(const FluidRenderer&) = delete;
  FluidRenderer& operator=(const FluidRenderer&) = delete;

  // Compiles shaders and allocates targets. Requires a current GL context and
  // a successful gfx::loadGl. Returns false and sets `error()` on failure; the
  // caller must then fall back to the 2D schematic view.
  bool initialise();
  bool ready() const;
  const std::string& error() const;
  void shutdown();

  // Renders one frame of `snapshot` seen through `camera` into an internal
  // target of `width` x `height` pixels, reallocating it when the size
  // changes. Returns the colour texture, or 0 when not ready.
  //
  // `pose` is where to DRAW the vessel, and is supplied by the caller rather
  // than read from the snapshot: the snapshot's pose is only as fresh as the
  // last completed physics step, whereas the panel knows the hand-driven part
  // of the vessel's motion every frame. Both the glass and the particles are
  // placed by this one matrix, so they always move as one rigid body.
  std::uint32_t render(const fluid::Snapshot&, const fluid::Pose& pose, const Camera3D&,
                       int width, int height, const FluidRenderSettings&);

  // Texture of the most recent successful render, 0 when there is none.
  std::uint32_t colourTexture() const;
  int width() const;
  int height() const;

  // Milliseconds spent in the last render, for the diagnostics line.
  double lastFrameMs() const;

 private:
  struct Impl;
  Impl* impl_;
};

}  // namespace chemcad::gfx

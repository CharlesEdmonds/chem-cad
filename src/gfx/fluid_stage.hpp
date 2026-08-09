#pragma once
// The hand-off between the fluid panel (which runs during ImGui submission and
// must not touch GL) and the app's render seam (which owns the GL context).
//
// Sequence per frame:
//   1. the panel sizes its viewport, updates the camera from input, publishes a
//      request here, and draws ImGui::Image with the texture from the PREVIOUS
//      frame (a one-frame latency, invisible at 60 Hz and the standard way to
//      composite an FBO into an ImGui window);
//   2. after ImGui::Render and before the ImGui backend's draw call, the app
//      renders the request into the FBO and stores the resulting texture here.
//
// The struct deliberately holds no GL types, so it can live in AppState and be
// compiled into the headless panel tests. When `texture` stays 0 -- headless, or
// a renderer that failed to initialise -- the panel falls back to the 2D
// schematic view.

#include <cstdint>
#include <memory>
#include <string>

#include "fluid/simulation.hpp"
#include "gfx/camera3d.hpp"
#include "gfx/fluid_renderer.hpp"

namespace chemcad::gfx {

struct FluidStage {
  // ---- written by the panel ------------------------------------------
  bool requested = false;                        // panel wants a frame this tick
  int width = 0;                                 // requested pixel size
  int height = 0;
  Camera3D camera;
  FluidRenderSettings settings;
  std::shared_ptr<const fluid::Snapshot> snapshot;  // state to draw
  // Where to draw the vessel this frame. The snapshot's own pose is only as
  // fresh as the last completed physics step; the panel knows the hand-driven
  // part of the motion every frame and writes the combined result here, so a
  // dragged vessel tracks the pointer instead of the solver's publish rate.
  fluid::Pose pose;

  // ---- written by the app render seam --------------------------------
  std::uint32_t texture = 0;   // colour texture of the last completed render
  int textureWidth = 0;
  int textureHeight = 0;
  bool available = false;      // a renderer exists and initialised
  std::string status;          // driver / timing / failure reason for the UI
  double lastFrameMs = 0.0;
};

}  // namespace chemcad::gfx

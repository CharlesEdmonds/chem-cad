#pragma once
// Adapts the OpenGL compute fluid backend to fluid::Accelerator, so
// fluid::Simulation can hand its integration step to the GPU without knowing
// that a GPU, or OpenGL, exists.
//
// The context callbacks are supplied by whoever owns the window system. gfx
// deliberately does not link GLFW: every GL entry point here comes through
// gfx::loadGl, and the one thing this module cannot do for itself is create a
// second context and make it current on the physics worker thread.

#include <cstddef>
#include <functional>
#include <memory>

#include "fluid/accelerator.hpp"

namespace chemcad::gfx {

// `makeCurrent` runs on the physics worker thread and must bind a context that
// is NOT the renderer's -- a context belongs to one thread at a time.
// `clearCurrent` runs on the same thread when the worker stops.
//
// `minimumParticles` is the charge below which the CPU solver is simply faster;
// see docs/gpu-acceleration.md for the measurement behind the default.
std::shared_ptr<fluid::Accelerator> makeFluidAccelerator(
    std::function<bool()> makeCurrent, std::function<void()> clearCurrent,
    std::size_t minimumParticles = 600);

}  // namespace chemcad::gfx

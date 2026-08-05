#pragma once
// Framebuffer capture for "Export PNG". Must be called with a current GL
// context, after rendering and before the buffer swap.

#include <string>

namespace chemcad::app {

// Grabs the rectangle (x, y, w, h) in *window* pixel coordinates with the
// origin at the top-left, flips it, and writes a PNG.
// Returns false and sets `error` on failure.
bool capturePng(int x, int y, int w, int h, int windowHeight,
                const std::string& path, std::string& error);

}  // namespace chemcad::app

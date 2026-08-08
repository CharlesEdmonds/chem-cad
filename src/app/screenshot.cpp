#include "app/screenshot.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

// <GL/gl.h> on Windows is a bare extension of <windows.h> and does not compile
// without the WINGDIAPI/APIENTRY macros it defines.
#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif
#include <GL/gl.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#include "third_party/stb_image_write.h"
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace chemcad::app {

bool capturePng(int x, int y, int w, int h, int windowHeight,
                const std::string& path, std::string& error) {
  error.clear();
  if (w <= 0 || h <= 0) {
    error = "nothing to capture";
    return false;
  }

  constexpr std::size_t kChannels = 4;
  const std::size_t rowBytes = static_cast<std::size_t>(w) * kChannels;
  std::vector<unsigned char> framebuffer(rowBytes * static_cast<std::size_t>(h));

  GLint previousPackAlignment = 0;
  glGetIntegerv(GL_PACK_ALIGNMENT, &previousPackAlignment);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(x, windowHeight - y - h, w, h, GL_RGBA, GL_UNSIGNED_BYTE,
               framebuffer.data());
  glPixelStorei(GL_PACK_ALIGNMENT, previousPackAlignment);

  std::vector<unsigned char> pngPixels(framebuffer.size());
  for (int row = 0; row < h; ++row) {
    const std::size_t sourceOffset = static_cast<std::size_t>(h - row - 1) * rowBytes;
    const std::size_t destinationOffset = static_cast<std::size_t>(row) * rowBytes;
    std::copy_n(framebuffer.data() + sourceOffset, rowBytes,
                pngPixels.data() + destinationOffset);
  }

  if (stbi_write_png(path.c_str(), w, h, static_cast<int>(kChannels), pngPixels.data(),
                     static_cast<int>(rowBytes)) == 0) {
    error = "cannot write PNG " + path;
    return false;
  }
  return true;
}

}  // namespace chemcad::app

#include "core/paths.hpp"

#include <cstdlib>

#ifndef CHEMCAD_DATA_DIR
#define CHEMCAD_DATA_DIR "data"
#endif

namespace chemcad::core {
namespace fs = std::filesystem;

fs::path dataDir() {
  if (const char* env = std::getenv("CHEMCAD_DATA_DIR"); env && *env) return fs::path(env);
  return fs::path(CHEMCAD_DATA_DIR);
}

fs::path cacheDir() {
  fs::path base;
  if (const char* xdg = std::getenv("XDG_CACHE_HOME"); xdg && *xdg) {
    base = fs::path(xdg);
  } else if (const char* home = std::getenv("HOME"); home && *home) {
    base = fs::path(home) / ".cache";
  } else {
    base = fs::temp_directory_path();
  }
  const fs::path dir = base / "chemcad";
  std::error_code ec;
  fs::create_directories(dir, ec);
  return dir;
}

}  // namespace chemcad::core

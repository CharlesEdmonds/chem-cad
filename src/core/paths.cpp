#include "core/paths.hpp"

#include <cstdlib>
#include <system_error>

#ifndef CHEMCAD_DATA_DIR
#define CHEMCAD_DATA_DIR "data"
#endif
#ifndef CHEMCAD_ASSETS_DIR
#define CHEMCAD_ASSETS_DIR "assets"
#endif

namespace chemcad::core {
namespace fs = std::filesystem;

fs::path dataDir() {
  if (const char* env = std::getenv("CHEMCAD_DATA_DIR"); env && *env) return fs::path(env);
  return fs::path(CHEMCAD_DATA_DIR);
}

fs::path assetsDir() {
  if (const char* env = std::getenv("CHEMCAD_ASSETS_DIR"); env && *env) return fs::path(env);
  return fs::path(CHEMCAD_ASSETS_DIR);
}

fs::path cacheDir() {
  fs::path base;
#ifdef _WIN32
  // %LOCALAPPDATA% is the Windows equivalent of $XDG_CACHE_HOME; fall back to
  // the profile root the same way the XDG spec falls back to $HOME.
  if (const char* local = std::getenv("LOCALAPPDATA"); local && *local) {
    base = fs::path(local);
  } else if (const char* profile = std::getenv("USERPROFILE"); profile && *profile) {
    base = fs::path(profile) / "AppData" / "Local";
  }
#else
  if (const char* xdg = std::getenv("XDG_CACHE_HOME"); xdg && *xdg) {
    base = fs::path(xdg);
  } else if (const char* home = std::getenv("HOME"); home && *home) {
    base = fs::path(home) / ".cache";
  }
#endif
  if (base.empty()) {
    std::error_code ec;
    base = fs::temp_directory_path(ec);
  }
  const fs::path dir = base / "chemcad";
  std::error_code ec;
  fs::create_directories(dir, ec);
  return dir;
}

}  // namespace chemcad::core

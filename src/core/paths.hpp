#pragma once
// Filesystem locations. Kept in core so every layer (naming, rxn, ui) can find
// data and cache directories without depending on the GUI/app layer.

#include <filesystem>

namespace chemcad::core {

// Bundled read-only data (data/ptable.json, data/reactions/*.json).
// Compile-time source path, overridable at runtime with $CHEMCAD_DATA_DIR.
std::filesystem::path dataDir();

// Bundled read-only assets (assets/fonts/*.ttf).
// Compile-time source path, overridable at runtime with $CHEMCAD_ASSETS_DIR.
std::filesystem::path assetsDir();

// Writable per-user cache, created on demand: %LOCALAPPDATA%\chemcad on
// Windows, $XDG_CACHE_HOME/chemcad (or ~/.cache/chemcad) elsewhere.
std::filesystem::path cacheDir();

}  // namespace chemcad::core

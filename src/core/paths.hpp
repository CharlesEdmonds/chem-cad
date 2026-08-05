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

// Writable per-user cache (~/.cache/chemcad), created on demand.
std::filesystem::path cacheDir();

}  // namespace chemcad::core

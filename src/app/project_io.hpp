#pragma once
// .chemcad project files.
// Schema (version 1):
//   { "version": 1,
//     "molecules": [ { "atoms": [ {"z","x","y","q","iso","h"} ],
//                      "bonds": [ {"a","b","o","s"} ] } ],
//     "planner": { "starts": ["SMILES", ...], "target": "SMILES" } }
// Bond "a"/"b" are POSITIONAL indices into that molecule's atoms array, not
// runtime AtomIds, so files stay stable across sessions.

#include <string>
#include <vector>

#include "core/model.hpp"

namespace chemcad::app {

struct Project {
  core::Document doc;
  std::vector<std::string> plannerStarts;
  std::string plannerTarget;
};

// Both throw std::runtime_error with a human-readable message on failure.
void saveProject(const Project&, const std::string& path);
Project loadProject(const std::string& path);

// Serialization without touching the filesystem (used by the roundtrip test).
std::string serializeProject(const Project&);
Project deserializeProject(const std::string& json);

}  // namespace chemcad::app

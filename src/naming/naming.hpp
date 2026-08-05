#pragma once
// Structure <-> name services.
//   name  -> structure : OPSIN (local jar, offline), PubChem as fallback
//   structure -> name  : PubChem PUG REST (no good offline IUPAC namer exists)
// Both directions are cached on disk at ~/.cache/chemcad/namecache.json.
//
// Every entry point BLOCKS on subprocess/network work. The UI must call these
// through chemcad::app::TaskRunner, never on the frame thread.

#include <string>

namespace chemcad::naming {

struct Result {
  bool ok = false;
  std::string value;
  std::string error;  // human-readable: "offline", "not found", ...
};

Result nameToSmiles(const std::string& name);
Result smilesToName(const std::string& smiles);

// Test/diagnostic hooks.
void setCachePath(const std::string& path);  // must be called before first use
std::string cachePath();
bool opsinAvailable();  // jar present AND java on PATH

}  // namespace chemcad::naming

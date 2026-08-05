// STUB -- implemented by the naming workstream.
#include "naming/naming.hpp"

namespace chemcad::naming {

Result nameToSmiles(const std::string&) { return {false, "", "unimplemented"}; }
Result smilesToName(const std::string&) { return {false, "", "unimplemented"}; }
void setCachePath(const std::string&) {}
std::string cachePath() { return {}; }
bool opsinAvailable() { return false; }

}  // namespace chemcad::naming

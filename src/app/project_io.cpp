// STUB -- implemented by the panels/export workstream.
#include "app/project_io.hpp"

#include <stdexcept>

namespace chemcad::app {

void saveProject(const Project&, const std::string&) {
  throw std::runtime_error("unimplemented");
}
Project loadProject(const std::string&) { throw std::runtime_error("unimplemented"); }
std::string serializeProject(const Project&) { throw std::runtime_error("unimplemented"); }
Project deserializeProject(const std::string&) { throw std::runtime_error("unimplemented"); }

}  // namespace chemcad::app

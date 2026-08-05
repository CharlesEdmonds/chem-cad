// STUB -- implemented by the reaction-engine workstream.
#include "rxn/llm.hpp"

namespace chemcad::rxn::llm {

bool configured() { return false; }
std::string baseUrl() { return "https://api.openai.com"; }
std::string model() { return "gpt-4o-mini"; }
std::vector<Route> proposeRoutes(const std::vector<std::string>&, const std::string&, int) {
  return {};
}
void enrich(std::vector<Route>&, const std::string&) {}

}  // namespace chemcad::rxn::llm

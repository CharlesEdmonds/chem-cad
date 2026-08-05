#pragma once
// OpenAI-compatible chat client used as the reaction planner's fallback lane.
//   CHEMCAD_LLM_API_KEY   required; absent => LLM lane disabled entirely
//   CHEMCAD_LLM_BASE_URL  default https://api.openai.com
//   CHEMCAD_LLM_MODEL     default gpt-4o-mini

#include <string>
#include <vector>

#include "rxn/engine.hpp"

namespace chemcad::rxn::llm {

bool configured();
std::string baseUrl();
std::string model();

// Ask the model for whole routes. Returns {} on any failure (network, bad key,
// unparseable JSON) -- the planner must degrade to KB-only, never throw.
std::vector<Route> proposeRoutes(const std::vector<std::string>& startSmiles,
                                 const std::string& targetSmiles, int maxRoutes);

// Enrich KB routes in place with practical notes. No-op on failure.
void enrich(std::vector<Route>& routes, const std::string& targetSmiles);

}  // namespace chemcad::rxn::llm

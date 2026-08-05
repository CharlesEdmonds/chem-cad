#pragma once
// Reaction planning: given starting materials and a target, propose routes with
// reagents, conditions and predicted side products.
//
// Hybrid strategy:
//   1. curated reaction-SMARTS knowledge base (data/reactions/*.json) - offline,
//      deterministic, explainable
//   2. an OpenAI-compatible LLM, used only to fill gaps the KB cannot reach and
//      to enrich KB routes with practical notes
//
// suggestRoutes() BLOCKS (RDKit search + optional network). Call it through
// chemcad::app::TaskRunner.

#include <string>
#include <vector>

namespace chemcad::rxn {

struct Step {
  std::string reactionName;                    // "Fischer esterification"
  std::vector<std::string> reactantSmiles;     // inputs consumed by this step
  std::vector<std::string> reagents;           // "H2SO4 (cat.)"
  std::string conditions;                      // "reflux"
  std::string productSmiles;                   // main product
  std::vector<std::string> sideProductSmiles;  // predicted co-products
  std::string notes;
  enum class Source { KB, LLM } source = Source::KB;
};

struct Route {
  std::vector<Step> steps;
  double score = 0;  // higher is better
  bool usesLlm() const {
    for (const Step& s : steps)
      if (s.source == Step::Source::LLM) return true;
    return false;
  }
};

struct Request {
  std::vector<std::string> startSmiles;
  std::string targetSmiles;
  int maxDepth = 3;
  int maxRoutes = 5;
  bool allowLlm = true;
};

std::vector<Route> suggestRoutes(const Request&);

// True when CHEMCAD_LLM_API_KEY is set in the environment.
bool llmAvailable();

// Number of reaction templates loaded from data/reactions/*.json.
size_t knowledgeBaseSize();

}  // namespace chemcad::rxn

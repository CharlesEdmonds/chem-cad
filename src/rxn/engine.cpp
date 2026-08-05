#include "rxn/engine.hpp"

#include <algorithm>
#include <cstddef>
#include <deque>
#include <iterator>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include "chem/bridge.hpp"
#include "rxn/kb.hpp"
#include "rxn/llm.hpp"

namespace chemcad::rxn {
namespace {

constexpr size_t kMaxApplications = 500;
constexpr double kKbScoreBase = 1.0e12;
constexpr double kDepthScoreUnit = 1.0e6;

struct SearchState {
  std::set<std::string> available;
  std::vector<Step> steps;
  std::vector<std::string> routeSignature;
  int prioritySum = 0;
};

std::string setKey(const std::set<std::string>& values) {
  std::string key;
  for (const std::string& value : values) {
    key += std::to_string(value.size());
    key.push_back(':');
    key += value;
    key.push_back(';');
  }
  return key;
}

std::string routeKey(const std::vector<std::string>& signature) {
  std::string key;
  for (const std::string& value : signature) {
    key += std::to_string(value.size());
    key.push_back(':');
    key += value;
    key.push_back(';');
  }
  return key;
}

std::shared_ptr<RDKit::ChemicalReaction> compiledReaction(const ReactionTemplate& reaction) {
  struct Cache {
    std::mutex mutex;
    std::unordered_map<std::string, std::shared_ptr<RDKit::ChemicalReaction>> reactions;
  };
  static Cache cache;

  std::lock_guard lock(cache.mutex);
  if (auto found = cache.reactions.find(reaction.id); found != cache.reactions.end()) {
    return found->second;
  }

  std::shared_ptr<RDKit::ChemicalReaction> compiled;
  try {
    compiled.reset(RDKit::RxnSmartsToChemicalReaction(reaction.smarts));
    if (compiled) compiled->initReactantMatchers(true);
  } catch (...) {
    compiled.reset();
  }
  cache.reactions.emplace(reaction.id, compiled);
  return compiled;
}

std::optional<std::string> sanitizedSmiles(const RDKit::ROMol& molecule) {
  try {
    RDKit::RWMol sanitized(molecule);
    RDKit::MolOps::sanitizeMol(sanitized);
    return RDKit::MolToSmiles(sanitized, true, false, -1, true);
  } catch (...) {
    return std::nullopt;
  }
}

void appendUnique(std::vector<std::string>& values, std::unordered_set<std::string>& seen,
                  std::string value) {
  if (seen.insert(value).second) values.push_back(std::move(value));
}

double kbScore(size_t stepCount, int prioritySum) {
  // The application cap also bounds route length, so this preserves the requested
  // lexicographic order: KB lane, then fewer steps, then summed priority.
  return kKbScoreBase + (static_cast<double>(kMaxApplications + 1) - stepCount) *
                            kDepthScoreUnit +
         prioritySum;
}

void assignLlmScores(std::vector<Route>& routes) {
  for (Route& route : routes) {
    route.score = (static_cast<double>(kMaxApplications + 1) - route.steps.size()) *
                  kDepthScoreUnit;
  }
}

void sortAndTruncate(std::vector<Route>& routes, int maxRoutes) {
  std::stable_sort(routes.begin(), routes.end(), [](const Route& lhs, const Route& rhs) {
    return lhs.score > rhs.score;
  });
  if (maxRoutes >= 0 && routes.size() > static_cast<size_t>(maxRoutes)) {
    routes.resize(static_cast<size_t>(maxRoutes));
  }
}

}  // namespace

size_t knowledgeBaseSize() { return knowledgeBase().size(); }

bool llmAvailable() { return llm::configured(); }

std::vector<Route> suggestRoutes(const Request& request) {
  std::vector<Route> routes;
  if (request.startSmiles.empty() || request.targetSmiles.empty() || request.maxRoutes <= 0) {
    return routes;
  }

  std::vector<std::string> starts;
  std::string target;
  try {
    starts.reserve(request.startSmiles.size());
    for (const std::string& smiles : request.startSmiles) {
      starts.push_back(chem::canonicalize(smiles));
    }
    target = chem::canonicalize(request.targetSmiles);
  } catch (...) {
    return {};
  }

  try {
    SearchState initial;
    initial.available.insert(starts.begin(), starts.end());
    std::deque<SearchState> frontier;
    frontier.push_back(std::move(initial));

    std::unordered_set<std::string> expandedStates;
    std::unordered_set<std::string> completedRoutes;
    std::unordered_map<std::string, RDKit::ROMOL_SPTR> molecules;
    size_t applications = 0;
    bool stopSearch = request.maxDepth <= 0;

    const auto moleculeFor = [&molecules](const std::string& smiles) -> RDKit::ROMOL_SPTR {
      if (auto found = molecules.find(smiles); found != molecules.end()) return found->second;
      RDKit::ROMOL_SPTR molecule(RDKit::SmilesToMol(smiles));
      molecules.emplace(smiles, molecule);
      return molecule;
    };

    const std::vector<ReactionTemplate>& templates = knowledgeBase();
    while (!frontier.empty() && !stopSearch) {
      SearchState state = std::move(frontier.front());
      frontier.pop_front();
      if (state.steps.size() >= static_cast<size_t>(request.maxDepth)) continue;
      if (!expandedStates.insert(setKey(state.available)).second) continue;

      const std::vector<std::string> available(state.available.begin(), state.available.end());
      for (const ReactionTemplate& reactionTemplate : templates) {
        const std::shared_ptr<RDKit::ChemicalReaction> reaction =
            compiledReaction(reactionTemplate);
        if (!reaction || reaction->getNumReactantTemplates() !=
                             static_cast<unsigned int>(reactionTemplate.arity)) {
          continue;
        }

        const size_t combinationCount = reactionTemplate.arity == 1
                                            ? available.size()
                                            : available.size() * available.size();
        for (size_t combination = 0; combination < combinationCount; ++combination) {
          if (applications >= kMaxApplications) {
            stopSearch = true;
            break;
          }
          ++applications;

          std::vector<std::string> consumed;
          if (reactionTemplate.arity == 1) {
            consumed.push_back(available[combination]);
          } else {
            consumed.push_back(available[combination / available.size()]);
            consumed.push_back(available[combination % available.size()]);
          }

          RDKit::MOL_SPTR_VECT reactants;
          bool validReactants = true;
          for (const std::string& smiles : consumed) {
            RDKit::ROMOL_SPTR molecule = moleculeFor(smiles);
            if (!molecule) {
              validReactants = false;
              break;
            }
            reactants.push_back(std::move(molecule));
          }
          if (!validReactants) continue;

          std::vector<RDKit::MOL_SPTR_VECT> productSets;
          try {
            productSets = reaction->runReactants(reactants);
          } catch (...) {
            continue;
          }

          for (const RDKit::MOL_SPTR_VECT& productSet : productSets) {
            if (productSet.empty() || !productSet.front()) continue;
            const std::optional<std::string> mainProduct =
                sanitizedSmiles(*productSet.front());
            if (!mainProduct) continue;

            Step step;
            step.reactionName = reactionTemplate.name;
            step.reactantSmiles = consumed;
            step.reagents = reactionTemplate.reagents;
            step.conditions = reactionTemplate.conditions;
            step.productSmiles = *mainProduct;
            step.notes = reactionTemplate.notes;
            step.source = Step::Source::KB;

            std::unordered_set<std::string> sideSeen;
            for (size_t index = 1; index < productSet.size(); ++index) {
              if (!productSet[index]) continue;
              if (auto side = sanitizedSmiles(*productSet[index])) {
                appendUnique(step.sideProductSmiles, sideSeen, std::move(*side));
              }
            }
            for (const std::string& byproduct : reactionTemplate.byproducts) {
              try {
                appendUnique(step.sideProductSmiles, sideSeen,
                             chem::canonicalize(byproduct));
              } catch (...) {
                // A bad optional byproduct annotation must not invalidate a usable transformation.
              }
            }

            SearchState next = state;
            next.steps.push_back(std::move(step));
            next.routeSignature.push_back(reactionTemplate.id + "\n" + *mainProduct);
            next.prioritySum += reactionTemplate.priority;

            if (*mainProduct == target) {
              const std::string signature = routeKey(next.routeSignature);
              if (completedRoutes.insert(signature).second) {
                Route route;
                route.steps = next.steps;
                route.score = kbScore(route.steps.size(), next.prioritySum);
                routes.push_back(std::move(route));
                if (routes.size() >= static_cast<size_t>(request.maxRoutes)) {
                  stopSearch = true;
                  break;
                }
              }
            }

            if (next.steps.size() < static_cast<size_t>(request.maxDepth) &&
                next.available.insert(*mainProduct).second) {
              frontier.push_back(std::move(next));
            }
          }
          if (stopSearch) break;
        }
        if (stopSearch) break;
      }
    }
  } catch (...) {
    // A malformed rule or unexpected RDKit failure degrades to the routes already found.
  }

  if (!routes.empty() && llm::configured()) {
    try {
      llm::enrich(routes, target);
    } catch (...) {
    }
  }
  if (request.allowLlm && llm::configured() &&
      routes.size() < static_cast<size_t>(request.maxRoutes)) {
    try {
      auto proposed = llm::proposeRoutes(starts, target, request.maxRoutes - routes.size());
      assignLlmScores(proposed);
      routes.insert(routes.end(), std::make_move_iterator(proposed.begin()),
                    std::make_move_iterator(proposed.end()));
    } catch (...) {
    }
  }

  try {
    sortAndTruncate(routes, request.maxRoutes);
  } catch (...) {
  }
  return routes;
}

}  // namespace chemcad::rxn

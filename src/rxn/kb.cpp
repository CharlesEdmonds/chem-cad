// STUB -- implemented by the reaction-engine workstream.
#include "rxn/kb.hpp"

namespace chemcad::rxn {

const std::vector<ReactionTemplate>& knowledgeBase() {
  static const std::vector<ReactionTemplate> empty;
  return empty;
}

std::vector<ReactionTemplate> loadReactionFile(const std::string&) { return {}; }

}  // namespace chemcad::rxn

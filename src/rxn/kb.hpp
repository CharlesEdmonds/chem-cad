#pragma once
// Curated named-reaction knowledge base loaded from data/reactions/*.json.
// Schema (one JSON array of these objects per file):
//   { "id", "name", "smarts", "arity", "reagents":[], "conditions",
//     "byproducts":[], "priority", "notes", "tags":[] }

#include <string>
#include <vector>

namespace chemcad::rxn {

struct ReactionTemplate {
  std::string id;
  std::string name;
  std::string smarts;                    // reaction SMARTS, main product first
  int arity = 1;                         // 1 or 2 reactants
  std::vector<std::string> reagents;
  std::string conditions;
  std::vector<std::string> byproducts;   // SMILES not expressible in the mapping
  int priority = 5;                      // 1..10, higher preferred
  std::string notes;
  std::vector<std::string> tags;
  std::string sourceFile;                // for diagnostics
};

// Loads and caches every data/reactions/*.json. Throws std::runtime_error when
// the directory is missing or a file is malformed.
const std::vector<ReactionTemplate>& knowledgeBase();

// Parse a single file; used by the KB validation test.
std::vector<ReactionTemplate> loadReactionFile(const std::string& path);

}  // namespace chemcad::rxn

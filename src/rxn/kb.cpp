#include "rxn/kb.hpp"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "core/paths.hpp"

namespace chemcad::rxn {
namespace {

std::runtime_error entryError(const std::string& path, size_t index, const std::string& detail) {
  return std::runtime_error("Reaction file " + path + ", entry " + std::to_string(index) +
                            ": " + detail);
}

template <typename T>
T optionalValue(const nlohmann::json& entry, const char* key, T fallback,
                const std::string& path, size_t index) {
  if (!entry.contains(key)) return fallback;
  try {
    return entry.at(key).get<T>();
  } catch (const std::exception& e) {
    throw entryError(path, index, std::string("invalid '") + key + "': " + e.what());
  }
}

}  // namespace

std::vector<ReactionTemplate> loadReactionFile(const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("Could not open reaction file " + path);

  nlohmann::json root;
  try {
    input >> root;
  } catch (const std::exception& e) {
    throw std::runtime_error("Could not parse reaction file " + path + ": " + e.what());
  }
  if (!root.is_array()) {
    throw std::runtime_error("Reaction file " + path + " must contain a JSON array");
  }

  std::vector<ReactionTemplate> templates;
  templates.reserve(root.size());
  for (size_t index = 0; index < root.size(); ++index) {
    const nlohmann::json& entry = root[index];
    if (!entry.is_object()) throw entryError(path, index, "expected an object");

    for (const char* key : {"id", "name", "smarts", "arity"}) {
      if (!entry.contains(key)) {
        throw entryError(path, index, std::string("missing required key '") + key + "'");
      }
    }

    ReactionTemplate reaction;
    try {
      reaction.id = entry.at("id").get<std::string>();
      reaction.name = entry.at("name").get<std::string>();
      reaction.smarts = entry.at("smarts").get<std::string>();
      reaction.arity = entry.at("arity").get<int>();
    } catch (const std::exception& e) {
      throw entryError(path, index, std::string("invalid required field: ") + e.what());
    }
    if (reaction.arity != 1 && reaction.arity != 2) {
      throw entryError(path, index, "'arity' must be 1 or 2");
    }

    reaction.reagents = optionalValue<std::vector<std::string>>(
        entry, "reagents", {}, path, index);
    reaction.conditions = optionalValue<std::string>(entry, "conditions", {}, path, index);
    reaction.byproducts = optionalValue<std::vector<std::string>>(
        entry, "byproducts", {}, path, index);
    reaction.priority = optionalValue<int>(entry, "priority", reaction.priority, path, index);
    reaction.notes = optionalValue<std::string>(entry, "notes", {}, path, index);
    reaction.tags = optionalValue<std::vector<std::string>>(entry, "tags", {}, path, index);
    reaction.sourceFile = path;
    templates.push_back(std::move(reaction));
  }
  return templates;
}

const std::vector<ReactionTemplate>& knowledgeBase() {
  static std::once_flag once;
  static std::vector<ReactionTemplate> templates;
  std::call_once(once, [] {
    const char* overrideDir = std::getenv("CHEMCAD_REACTIONS_DIR");
    const std::filesystem::path directory = overrideDir && *overrideDir
                                                ? std::filesystem::path(overrideDir)
                                                : core::dataDir() / "reactions";
    if (!std::filesystem::is_directory(directory)) {
      throw std::runtime_error("Reaction directory does not exist: " + directory.string());
    }

    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
      if (entry.is_regular_file() && entry.path().extension() == ".json") {
        files.push_back(entry.path());
      }
    }
    std::sort(files.begin(), files.end(), [](const auto& lhs, const auto& rhs) {
      return lhs.filename().string() < rhs.filename().string();
    });

    for (const auto& file : files) {
      auto loaded = loadReactionFile(file.string());
      templates.insert(templates.end(), std::make_move_iterator(loaded.begin()),
                       std::make_move_iterator(loaded.end()));
    }
  });
  return templates;
}

}  // namespace chemcad::rxn

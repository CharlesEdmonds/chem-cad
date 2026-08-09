#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <GraphMol/Atom.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <memory>
#include <regex>
#include <set>
#include <string>
#include <vector>

#include "core/paths.hpp"
#include "rxn/kb.hpp"

namespace {

using Json = nlohmann::json;

Json readJson(const std::filesystem::path& path) {
  std::ifstream input(path);
  REQUIRE_MESSAGE(input.good(), "Could not open " << path.string());
  Json value;
  REQUIRE_NOTHROW(input >> value);
  return value;
}

int reactantComponentCount(const std::string& smarts) {
  const std::size_t arrow = smarts.find(">>");
  if (arrow == std::string::npos || arrow == 0) return 0;

  int components = 1;
  int bracketDepth = 0;
  int parenDepth = 0;
  for (std::size_t i = 0; i < arrow; ++i) {
    switch (smarts[i]) {
      case '[': ++bracketDepth; break;
      case ']': --bracketDepth; break;
      case '(':
        if (bracketDepth == 0) ++parenDepth;
        break;
      case ')':
        if (bracketDepth == 0) --parenDepth;
        break;
      case '.':
        if (bracketDepth == 0 && parenDepth == 0) ++components;
        break;
      default: break;
    }
    if (bracketDepth < 0 || parenDepth < 0) return 0;
  }
  return bracketDepth == 0 && parenDepth == 0 ? components : 0;
}

void checkStringArray(const Json& entry, const char* key) {
  INFO("key: " << key);
  REQUIRE(entry.contains(key));
  REQUIRE(entry.at(key).is_array());
  for (const auto& item : entry.at(key)) CHECK(item.is_string());
}

}  // namespace

TEST_CASE("reaction knowledge base schema and SMARTS are valid") {
  const auto reactionDir = chemcad::core::dataDir() / "reactions";
  REQUIRE(std::filesystem::is_directory(reactionDir));

  std::vector<std::filesystem::path> files;
  for (const auto& item : std::filesystem::directory_iterator(reactionDir)) {
    if (item.is_regular_file() && item.path().extension() == ".json") files.push_back(item.path());
  }
  std::sort(files.begin(), files.end());
  REQUIRE_FALSE(files.empty());

  const std::array<const char*, 13> keys = {
      "id",         "name",       "smarts",     "arity",    "reagents",
      "conditions", "substrate",  "outcome",    "source",   "byproducts",
      "priority",   "notes",      "tags"};
  const std::regex idPattern("^[a-z][a-z0-9]*(?:_[a-z0-9]+)*$");
  const std::regex substratePattern("^[a-z][a-z0-9 -]*s$");
  std::set<std::string> ids;
  std::size_t templateCount = 0;

  for (const auto& file : files) {
    INFO("file: " << file.string());
    const Json root = readJson(file);
    REQUIRE(root.is_array());
    CHECK(root.size() >= 12);

    const auto loaded = chemcad::rxn::loadReactionFile(file.string());
    REQUIRE(loaded.size() == root.size());
    REQUIRE_FALSE(loaded.empty());
    CHECK(chemcad::rxn::reactionType(loaded.front()) == file.stem().string());

    for (std::size_t entryIndex = 0; entryIndex < root.size(); ++entryIndex) {
      const Json& entry = root.at(entryIndex);
      const chemcad::rxn::ReactionTemplate& loadedReaction = loaded.at(entryIndex);
      ++templateCount;
      REQUIRE(entry.is_object());
      const std::string id = entry.contains("id") && entry.at("id").is_string()
                                 ? entry.at("id").get<std::string>()
                                 : "<missing-or-invalid-id>";
      INFO("reaction id: " << id);
      INFO("file: " << file.filename().string());

      CHECK(entry.size() == keys.size());
      for (const char* key : keys) CHECK(entry.contains(key));

      REQUIRE(entry.contains("id"));
      REQUIRE(entry.at("id").is_string());
      CHECK(std::regex_match(id, idPattern));
      CHECK(ids.insert(id).second);

      for (const char* key : {"name", "smarts", "conditions", "substrate", "outcome",
                              "source", "notes"}) {
        INFO("key: " << key);
        REQUIRE(entry.contains(key));
        REQUIRE(entry.at(key).is_string());
      }
      CHECK_FALSE(id.empty());
      CHECK_FALSE(entry.at("name").get<std::string>().empty());
      CHECK_FALSE(entry.at("smarts").get<std::string>().empty());
      CHECK_FALSE(entry.at("conditions").get<std::string>().empty());
      const std::string substrate = entry.at("substrate").get<std::string>();
      CHECK_FALSE(substrate.empty());
      CHECK(std::regex_match(substrate, substratePattern));
      CHECK_FALSE(entry.at("source").get<std::string>().empty());
      CHECK(loadedReaction.id == id);
      CHECK(loadedReaction.substrate == substrate);
      CHECK(loadedReaction.outcome == entry.at("outcome").get<std::string>());
      CHECK(loadedReaction.source == entry.at("source").get<std::string>());
      checkStringArray(entry, "reagents");
      checkStringArray(entry, "byproducts");
      checkStringArray(entry, "tags");

      REQUIRE(entry.contains("arity"));
      REQUIRE(entry.at("arity").is_number_integer());
      const int arity = entry.at("arity").get<int>();
      CHECK((arity == 1 || arity == 2));

      REQUIRE(entry.contains("priority"));
      REQUIRE(entry.at("priority").is_number_integer());
      const int priority = entry.at("priority").get<int>();
      CHECK(priority >= 1);
      CHECK(priority <= 10);

      const std::string smarts = entry.at("smarts").get<std::string>();
      CHECK(smarts.find(">>") != std::string::npos);
      CHECK(reactantComponentCount(smarts) == arity);

      std::unique_ptr<RDKit::ChemicalReaction> reaction;
      try {
        reaction.reset(RDKit::RxnSmartsToChemicalReaction(smarts));
      } catch (const std::exception& error) {
        FAIL_CHECK("SMARTS parse threw: " << error.what());
      }
      REQUIRE(reaction != nullptr);
      CHECK(reaction->getNumReactantTemplates() >= 1);
      CHECK(reaction->getNumProductTemplates() >= 1);
      CHECK(reaction->getNumReactantTemplates() == static_cast<unsigned int>(arity));
      CHECK_NOTHROW(reaction->initReactantMatchers());

      for (const auto& byproductValue : entry.at("byproducts")) {
        const std::string smiles = byproductValue.get<std::string>();
        INFO("byproduct: " << smiles);
        std::unique_ptr<RDKit::ROMol> molecule;
        CHECK_NOTHROW(molecule.reset(RDKit::SmilesToMol(smiles)));
        REQUIRE(molecule != nullptr);
        int formalCharge = 0;
        for (const RDKit::Atom* atom : molecule->atoms()) formalCharge += atom->getFormalCharge();
        CHECK(formalCharge == 0);
      }
    }
  }

  CHECK(templateCount >= 132);
}

TEST_CASE("periodic table contains all 118 elements") {
  const auto path = chemcad::core::dataDir() / "ptable.json";
  const Json root = readJson(path);
  REQUIRE(root.is_array());
  REQUIRE(root.size() == 118);

  const std::set<std::string> categories = {
      "alkali-metal",   "alkaline-earth", "transition-metal", "post-transition",
      "metalloid",      "nonmetal",       "halogen",          "noble-gas",
      "lanthanide",     "actinide"};
  const std::array<const char*, 7> keys = {"z", "symbol", "name", "mass",
                                           "group", "period", "category"};
  std::set<std::string> symbols;

  for (std::size_t i = 0; i < root.size(); ++i) {
    const Json& element = root.at(i);
    INFO("periodic table entry: " << (i + 1));
    REQUIRE(element.is_object());
    CHECK(element.size() == keys.size());
    for (const char* key : keys) CHECK(element.contains(key));

    REQUIRE(element.at("z").is_number_integer());
    CHECK(element.at("z").get<int>() == static_cast<int>(i + 1));
    REQUIRE(element.at("symbol").is_string());
    CHECK(symbols.insert(element.at("symbol").get<std::string>()).second);
    REQUIRE(element.at("name").is_string());
    CHECK_FALSE(element.at("name").get<std::string>().empty());
    REQUIRE(element.at("mass").is_number());
    CHECK(element.at("mass").get<double>() > 0.0);
    REQUIRE(element.at("period").is_number_integer());
    CHECK(element.at("period").get<int>() >= 1);
    CHECK(element.at("period").get<int>() <= 7);
    REQUIRE(element.at("group").is_number_integer());
    const int z = element.at("z").get<int>();
    const int group = element.at("group").get<int>();
    if ((z >= 57 && z <= 71) || (z >= 89 && z <= 103)) {
      CHECK(group == 0);
    } else {
      CHECK(group >= 1);
      CHECK(group <= 18);
    }
    REQUIRE(element.at("category").is_string());
    CHECK(categories.contains(element.at("category").get<std::string>()));
  }
}

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "chem/bridge.hpp"
#include "rxn/engine.hpp"

namespace {

// setenv/unsetenv are POSIX; MSVC only has _putenv_s, where an empty value
// removes the variable outright.
void setEnvironment(const char* name, const char* value) {
#ifdef _WIN32
  _putenv_s(name, value ? value : "");
#else
  if (value) {
    setenv(name, value, 1);
  } else {
    unsetenv(name);
  }
#endif
}

struct TestEnvironment {
  TestEnvironment() {
    const std::filesystem::path fixtures =
        std::filesystem::path(__FILE__).parent_path() / "fixtures" / "reactions";
    setEnvironment("CHEMCAD_REACTIONS_DIR", fixtures.string().c_str());
    setEnvironment("CHEMCAD_LLM_API_KEY", nullptr);
  }
};

TestEnvironment testEnvironment;

bool containsSmiles(const std::vector<std::string>& values, const std::string& expected) {
  return std::find(values.begin(), values.end(), expected) != values.end();
}

}  // namespace

TEST_CASE("Fischer esterification is found in one step with water") {
  chemcad::rxn::Request request;
  request.startSmiles = {"CC(=O)O", "CCO"};
  request.targetSmiles = "CCOC(C)=O";
  request.maxDepth = 3;
  request.maxRoutes = 5;
  request.allowLlm = false;

  const auto routes = chemcad::rxn::suggestRoutes(request);
  std::cout << "one-step routes: " << routes.size() << '\n';
  REQUIRE_FALSE(routes.empty());
  REQUIRE(routes.front().steps.size() == 1);

  const auto& step = routes.front().steps.front();
  CHECK(chemcad::chem::canonicalize(step.productSmiles) ==
        chemcad::chem::canonicalize(request.targetSmiles));
  CHECK(containsSmiles(step.sideProductSmiles, chemcad::chem::canonicalize("O")));
}

TEST_CASE("SN2 followed by oxidation reaches a carboxylic acid in two steps") {
  chemcad::rxn::Request request;
  request.startSmiles = {"CCCCBr"};
  request.targetSmiles = "CCCC(=O)O";
  request.maxDepth = 2;
  request.maxRoutes = 5;
  request.allowLlm = false;

  const auto routes = chemcad::rxn::suggestRoutes(request);
  std::cout << "two-step routes: " << routes.size() << '\n';
  const auto route = std::find_if(routes.begin(), routes.end(), [](const auto& candidate) {
    return candidate.steps.size() == 2;
  });
  REQUIRE(route != routes.end());

  const std::string alcohol = chemcad::chem::canonicalize("CCCCO");
  CHECK(chemcad::chem::canonicalize(route->steps[0].productSmiles) == alcohol);
  CHECK(containsSmiles(route->steps[1].reactantSmiles, alcohol));
}

TEST_CASE("maximum depth prevents a required second step") {
  chemcad::rxn::Request request;
  request.startSmiles = {"CCCCBr"};
  request.targetSmiles = "CCCC(=O)O";
  request.maxDepth = 1;
  request.maxRoutes = 5;
  request.allowLlm = false;

  const auto routes = chemcad::rxn::suggestRoutes(request);
  std::cout << "depth-one routes: " << routes.size() << '\n';
  CHECK(routes.empty());
}

TEST_CASE("an unreachable target returns no routes without throwing") {
  chemcad::rxn::Request request;
  request.startSmiles = {"C"};
  request.targetSmiles = "c1ccccc1";
  request.maxDepth = 3;
  request.maxRoutes = 5;
  request.allowLlm = false;

  std::vector<chemcad::rxn::Route> routes;
  CHECK_NOTHROW(routes = chemcad::rxn::suggestRoutes(request));
  std::cout << "unreachable routes: " << routes.size() << '\n';
  CHECK(routes.empty());
}

TEST_CASE("invalid requests return no routes without throwing") {
  chemcad::rxn::Request request;
  std::vector<chemcad::rxn::Route> routes;

  request.targetSmiles = "CCO";
  CHECK_NOTHROW(routes = chemcad::rxn::suggestRoutes(request));
  CHECK(routes.empty());

  request.startSmiles = {"CCO"};
  request.targetSmiles.clear();
  CHECK_NOTHROW(routes = chemcad::rxn::suggestRoutes(request));
  CHECK(routes.empty());

  request.startSmiles = {"this-is-not-smiles"};
  request.targetSmiles = "CCO";
  CHECK_NOTHROW(routes = chemcad::rxn::suggestRoutes(request));
  CHECK(routes.empty());
}

TEST_CASE("LLM availability follows the API key environment variable") {
  setEnvironment("CHEMCAD_LLM_API_KEY", nullptr);
  CHECK_FALSE(chemcad::rxn::llmAvailable());
}

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "naming/naming.hpp"

#include <nlohmann/json.hpp>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <random>

namespace {

namespace fs = std::filesystem;
using json = nlohmann::json;

fs::path uniqueCachePath(const char* label) {
  const auto tick = std::chrono::steady_clock::now().time_since_epoch().count();
  static std::random_device entropy;
  return fs::temp_directory_path() /
         ("chemcad-namecache-" + std::string(label) + "-" + std::to_string(entropy()) + "-" +
          std::to_string(tick) + ".json");
}

struct RemoveOnExit {
  fs::path path;
  ~RemoveOnExit() {
    std::error_code ec;
    fs::remove(path, ec);
  }
};

void writeJson(const fs::path& path, const json& value) {
  std::ofstream output(path, std::ios::trunc);
  REQUIRE(output.good());
  output << value.dump(2) << '\n';
  REQUIRE(output.good());
}

}  // namespace

TEST_CASE("naming cache hits are hermetic and name keys ignore case") {
  const fs::path path = uniqueCachePath("hits");
  RemoveOnExit cleanup{path};
  std::error_code ec;
  fs::remove(path, ec);

  const std::string syntheticSmiles = "synthetic-smiles://cache-only";
  writeJson(path, {{"n2s", {{"testonlyunobtainium", syntheticSmiles}}},
                   {"s2n", {{syntheticSmiles, "cache-only systematic name"}}}});
  chemcad::naming::setCachePath(path.string());

  const auto lower = chemcad::naming::nameToSmiles("testonlyunobtainium");
  CHECK(lower.ok);
  CHECK(lower.value == syntheticSmiles);
  CHECK(lower.error.empty());

  const auto mixedCase = chemcad::naming::nameToSmiles("  TestOnlyUnobtainium  ");
  CHECK(mixedCase.ok);
  CHECK(mixedCase.value == syntheticSmiles);

  const auto reverse = chemcad::naming::smilesToName(syntheticSmiles);
  CHECK(reverse.ok);
  CHECK(reverse.value == "cache-only systematic name");
  CHECK(reverse.error.empty());
}

TEST_CASE("empty input and corrupt cache are harmless") {
  const fs::path path = uniqueCachePath("corrupt");
  RemoveOnExit cleanup{path};
  std::error_code ec;
  fs::remove(path, ec);
  {
    std::ofstream output(path, std::ios::trunc);
    REQUIRE(output.good());
    output << "{ this is not json";
  }
  chemcad::naming::setCachePath(path.string());
  CHECK(chemcad::naming::cachePath() == path.string());

  chemcad::naming::Result emptyName;
  CHECK_NOTHROW(emptyName = chemcad::naming::nameToSmiles(" \t\n "));
  CHECK_FALSE(emptyName.ok);
  CHECK(emptyName.value.empty());
  CHECK_FALSE(emptyName.error.empty());

  chemcad::naming::Result emptySmiles;
  CHECK_NOTHROW(emptySmiles = chemcad::naming::smilesToName(" \t\n "));
  CHECK_FALSE(emptySmiles.ok);
  CHECK(emptySmiles.value.empty());
}

TEST_CASE("optional live OPSIN subprocess") {
  const char* enabled = std::getenv("CHEMCAD_TEST_OPSIN");
  if (!enabled || std::string(enabled) != "1") return;

  REQUIRE(chemcad::naming::opsinAvailable());
  const fs::path path = uniqueCachePath("opsin");
  RemoveOnExit cleanup{path};
  std::error_code ec;
  fs::remove(path, ec);
  chemcad::naming::setCachePath(path.string());

  const auto result = chemcad::naming::nameToSmiles("ethanol");
  REQUIRE(result.ok);
  std::cout << "OPSIN subprocess ethanol: " << result.value << '\n';
}

TEST_CASE("optional live PubChem lookup") {
  const char* enabled = std::getenv("CHEMCAD_TEST_NETWORK");
  if (!enabled || std::string(enabled) != "1") return;

  const fs::path path = uniqueCachePath("network");
  RemoveOnExit cleanup{path};
  std::error_code ec;
  fs::remove(path, ec);
  chemcad::naming::setCachePath(path.string());

  const auto result = chemcad::naming::smilesToName("CCO");
  CHECK(result.ok);
  CHECK_FALSE(result.value.empty());
}

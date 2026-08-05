#include "rxn/llm.hpp"

#include <algorithm>
#include <cstdlib>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include <curl/curl.h>
#include <nlohmann/json.hpp>

#include "chem/bridge.hpp"

namespace chemcad::rxn::llm {
namespace {

size_t receiveBody(char* data, size_t size, size_t count, void* destination) {
  const size_t bytes = size * count;
  static_cast<std::string*>(destination)->append(data, bytes);
  return bytes;
}

std::string endpointUrl() {
  std::string url = baseUrl();
  while (!url.empty() && url.back() == '/') url.pop_back();
  return url + "/v1/chat/completions";
}

std::optional<std::string> postJson(const nlohmann::json& request) {
  try {
    static std::once_flag curlOnce;
    static CURLcode curlInit = CURLE_FAILED_INIT;
    std::call_once(curlOnce, [] { curlInit = curl_global_init(CURL_GLOBAL_DEFAULT); });
    const char* apiKey = std::getenv("CHEMCAD_LLM_API_KEY");
    if (curlInit != CURLE_OK || !apiKey || !*apiKey) return std::nullopt;

    CURL* curl = curl_easy_init();
    if (!curl) return std::nullopt;

    std::string response;
    const std::string body = request.dump();
    const std::string authorization = std::string("Authorization: Bearer ") + apiKey;
    curl_slist* headers = nullptr;
    headers = curl_slist_append(headers, "Content-Type: application/json");
    headers = curl_slist_append(headers, authorization.c_str());

    const std::string url = endpointUrl();
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
    curl_easy_setopt(curl, CURLOPT_POST, 1L);
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, body.c_str());
    curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, static_cast<long>(body.size()));
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 30L);
    curl_easy_setopt(curl, CURLOPT_CONNECTTIMEOUT, 10L);
    curl_easy_setopt(curl, CURLOPT_NOSIGNAL, 1L);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, receiveBody);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);

    const CURLcode result = curl_easy_perform(curl);
    long status = 0;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &status);
    curl_slist_free_all(headers);
    curl_easy_cleanup(curl);

    if (result != CURLE_OK || status < 200 || status >= 300) return std::nullopt;
    return response;
  } catch (...) {
    return std::nullopt;
  }
}

std::string trim(std::string value) {
  const auto whitespace = [](unsigned char c) {
    return c == ' ' || c == '\t' || c == '\r' || c == '\n';
  };
  auto first = std::find_if_not(value.begin(), value.end(), whitespace);
  auto last = std::find_if_not(value.rbegin(), value.rend(), whitespace).base();
  if (first >= last) return {};
  return std::string(first, last);
}

std::string stripJsonFence(std::string content) {
  content = trim(std::move(content));
  if (!content.starts_with("```")) return content;
  const size_t lineEnd = content.find('\n');
  if (lineEnd == std::string::npos) return content;
  const size_t fenceEnd = content.rfind("```");
  if (fenceEnd <= lineEnd) return content;
  return trim(content.substr(lineEnd + 1, fenceEnd - lineEnd - 1));
}

std::optional<nlohmann::json> chatCompletion(const std::string& systemPrompt,
                                             const std::string& userPrompt) {
  try {
    nlohmann::json request = {
        {"model", model()},
        {"temperature", 0.2},
        {"response_format", {{"type", "json_object"}}},
        {"messages",
         {{{"role", "system"}, {"content", systemPrompt}},
          {{"role", "user"}, {"content", userPrompt}}}}};
    const auto responseBody = postJson(request);
    if (!responseBody) return std::nullopt;

    const nlohmann::json response = nlohmann::json::parse(*responseBody);
    const std::string content =
        response.at("choices").at(0).at("message").at("content").get<std::string>();
    return nlohmann::json::parse(stripJsonFence(content));
  } catch (...) {
    return std::nullopt;
  }
}

std::optional<std::vector<std::string>> canonicalInputs(
    const std::vector<std::string>& smiles) {
  try {
    std::vector<std::string> canonical;
    canonical.reserve(smiles.size());
    for (const std::string& value : smiles) canonical.push_back(chem::canonicalize(value));
    return canonical;
  } catch (...) {
    return std::nullopt;
  }
}

}  // namespace

bool configured() {
  const char* key = std::getenv("CHEMCAD_LLM_API_KEY");
  return key && *key;
}

std::string baseUrl() {
  const char* value = std::getenv("CHEMCAD_LLM_BASE_URL");
  return value && *value ? value : "https://api.openai.com";
}

std::string model() {
  const char* value = std::getenv("CHEMCAD_LLM_MODEL");
  return value && *value ? value : "gpt-4o-mini";
}

std::vector<Route> proposeRoutes(const std::vector<std::string>& startSmiles,
                                 const std::string& targetSmiles, int maxRoutes) {
  try {
    if (!configured() || startSmiles.empty() || targetSmiles.empty() || maxRoutes <= 0) return {};
    const auto starts = canonicalInputs(startSmiles);
    if (!starts) return {};
    std::string target;
    try {
      target = chem::canonicalize(targetSmiles);
    } catch (...) {
      return {};
    }

    const std::string systemPrompt =
        "You are an expert synthetic organic chemist. Answer ONLY with valid JSON and no prose.";
    nlohmann::json promptInputs = *starts;
    const std::string userPrompt =
        "Starting material SMILES: " + promptInputs.dump() + "\nTarget SMILES: " + target +
        "\nPropose at most " + std::to_string(maxRoutes) +
        " practical routes. Return exactly this shape: "
        R"({"routes":[{"steps":[{"reaction_name":"","reagents":[],"conditions":"","product_smiles":"","side_product_smiles":[],"notes":""}]}]})";

    const auto answer = chatCompletion(systemPrompt, userPrompt);
    if (!answer || !answer->contains("routes") || !answer->at("routes").is_array()) return {};

    std::vector<Route> routes;
    for (const nlohmann::json& routeJson : answer->at("routes")) {
      if (routes.size() >= static_cast<size_t>(maxRoutes)) break;
      if (!routeJson.is_object() || !routeJson.contains("steps") ||
          !routeJson.at("steps").is_array()) {
        return {};
      }

      Route route;
      std::vector<std::string> currentInputs = *starts;
      for (const nlohmann::json& stepJson : routeJson.at("steps")) {
        if (!stepJson.is_object()) return {};
        for (const char* key : {"reaction_name", "reagents", "conditions", "product_smiles",
                                "side_product_smiles", "notes"}) {
          if (!stepJson.contains(key)) return {};
        }

        Step step;
        try {
          step.reactionName = stepJson.at("reaction_name").get<std::string>();
          step.reagents = stepJson.at("reagents").get<std::vector<std::string>>();
          step.conditions = stepJson.at("conditions").get<std::string>();
          step.notes = stepJson.at("notes").get<std::string>();
          step.reactantSmiles = currentInputs;
          step.productSmiles =
              chem::canonicalize(stepJson.at("product_smiles").get<std::string>());
          step.source = Step::Source::LLM;
        } catch (const chem::ChemError&) {
          continue;
        } catch (...) {
          return {};
        }

        if (!stepJson.at("side_product_smiles").is_array()) return {};
        std::unordered_set<std::string> sideSeen;
        for (const nlohmann::json& sideJson : stepJson.at("side_product_smiles")) {
          if (!sideJson.is_string()) return {};
          try {
            std::string side = chem::canonicalize(sideJson.get<std::string>());
            if (sideSeen.insert(side).second) step.sideProductSmiles.push_back(std::move(side));
          } catch (...) {
          }
        }

        currentInputs = {step.productSmiles};
        route.steps.push_back(std::move(step));
      }
      if (!route.steps.empty()) routes.push_back(std::move(route));
    }
    return routes;
  } catch (...) {
    return {};
  }
}

void enrich(std::vector<Route>& routes, const std::string& targetSmiles) {
  try {
    if (!configured() || routes.empty()) return;
    const std::string target = chem::canonicalize(targetSmiles);

    nlohmann::json summary = nlohmann::json::array();
    for (const Route& route : routes) {
      nlohmann::json steps = nlohmann::json::array();
      for (const Step& step : route.steps) {
        steps.push_back({{"reaction_name", step.reactionName},
                         {"reactants", step.reactantSmiles},
                         {"reagents", step.reagents},
                         {"conditions", step.conditions},
                         {"product", step.productSmiles}});
      }
      summary.push_back(std::move(steps));
    }

    const std::string systemPrompt =
        "You are an expert synthetic organic chemist. Answer ONLY with valid JSON and no prose.";
    const std::string userPrompt =
        "Target SMILES: " + target + "\nThese are curated knowledge-base routes: " +
        summary.dump() +
        "\nGive one short practical note per step, preserving route and step order. Return exactly "
        R"({"notes":[["...","..."]]})";
    const auto answer = chatCompletion(systemPrompt, userPrompt);
    if (!answer || !answer->contains("notes") || !answer->at("notes").is_array()) return;

    const auto& routeNotes = answer->at("notes");
    for (size_t routeIndex = 0; routeIndex < routes.size() && routeIndex < routeNotes.size();
         ++routeIndex) {
      if (!routeNotes[routeIndex].is_array()) continue;
      for (size_t stepIndex = 0;
           stepIndex < routes[routeIndex].steps.size() &&
           stepIndex < routeNotes[routeIndex].size();
           ++stepIndex) {
        if (!routeNotes[routeIndex][stepIndex].is_string()) continue;
        const std::string note = routeNotes[routeIndex][stepIndex].get<std::string>();
        std::string& existing = routes[routeIndex].steps[stepIndex].notes;
        if (!note.empty() && (existing.empty() || existing.size() < note.size())) existing = note;
      }
    }
  } catch (...) {
  }
}

}  // namespace chemcad::rxn::llm

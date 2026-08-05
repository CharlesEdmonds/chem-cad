#include "naming/naming.hpp"

#include "core/paths.hpp"

#include <curl/curl.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cerrno>
#include <csignal>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <mutex>
#include <optional>
#include <poll.h>
#include <shared_mutex>
#include <string_view>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
#include <unordered_map>

namespace chemcad::naming {
namespace {

namespace fs = std::filesystem;
using json = nlohmann::json;

struct Cache {
  fs::path path;
  bool pathOverridden = false;
  bool loaded = false;
  std::unordered_map<std::string, std::string> n2s;
  std::unordered_map<std::string, std::string> s2n;
};

Cache gCache;
std::mutex gCacheMutex;
std::shared_mutex gPathGuard;
std::once_flag gCurlInitFlag;

std::string trim(std::string_view value) {
  auto first = std::find_if_not(value.begin(), value.end(), [](unsigned char c) {
    return std::isspace(c) != 0;
  });
  auto last = std::find_if_not(value.rbegin(), value.rend(), [](unsigned char c) {
                return std::isspace(c) != 0;
              }).base();
  if (first >= last) return {};
  return std::string(first, last);
}

std::string nameKey(std::string_view name) {
  std::string key = trim(name);
  std::transform(key.begin(), key.end(), key.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return key;
}

fs::path effectiveCachePathLocked() {
  if (gCache.pathOverridden) return gCache.path;
  return core::cacheDir() / "namecache.json";
}

void loadStringMap(const json& root, const char* key,
                   std::unordered_map<std::string, std::string>& output) {
  const auto section = root.find(key);
  if (section == root.end() || !section->is_object()) return;
  for (auto it = section->begin(); it != section->end(); ++it) {
    if (it.value().is_string()) output.emplace(it.key(), it.value().get<std::string>());
  }
}

void ensureCacheLoadedLocked() {
  if (gCache.loaded) return;
  gCache.loaded = true;
  gCache.n2s.clear();
  gCache.s2n.clear();

  try {
    std::ifstream input(effectiveCachePathLocked());
    if (!input) return;
    const json root = json::parse(input);
    loadStringMap(root, "n2s", gCache.n2s);
    loadStringMap(root, "s2n", gCache.s2n);
  } catch (...) {
    gCache.n2s.clear();
    gCache.s2n.clear();
  }
}

void saveCacheLocked() noexcept {
  try {
    const fs::path path = effectiveCachePathLocked();
    if (!path.parent_path().empty()) {
      std::error_code ec;
      fs::create_directories(path.parent_path(), ec);
    }
    json root = {{"n2s", gCache.n2s}, {"s2n", gCache.s2n}};
    std::ofstream output(path, std::ios::trunc);
    if (output) output << root.dump(2) << '\n';
  } catch (...) {
    // A read-only cache directory must not make naming itself fail.
  }
}

std::optional<std::string> cachedValue(bool namesToSmiles, const std::string& key) {
  std::lock_guard lock(gCacheMutex);
  ensureCacheLoadedLocked();
  const auto& entries = namesToSmiles ? gCache.n2s : gCache.s2n;
  const auto found = entries.find(key);
  if (found == entries.end()) return std::nullopt;
  return found->second;
}

void cacheValue(bool namesToSmiles, const std::string& key, const std::string& value) {
  std::lock_guard lock(gCacheMutex);
  ensureCacheLoadedLocked();
  auto& entries = namesToSmiles ? gCache.n2s : gCache.s2n;
  entries[key] = value;
  saveCacheLocked();
}

#ifdef CHEMCAD_OPSIN_JAR
fs::path opsinJarPath() { return fs::path(CHEMCAD_OPSIN_JAR); }
#endif

std::optional<fs::path> javaPath() {
  const char* pathValue = std::getenv("PATH");
  if (!pathValue) return std::nullopt;

  std::string_view path(pathValue);
  std::size_t start = 0;
  while (start <= path.size()) {
    const std::size_t end = path.find(':', start);
    const std::string_view part = path.substr(start, end - start);
    const fs::path candidate = (part.empty() ? fs::path(".") : fs::path(part)) / "java";
    std::error_code ec;
    if (fs::is_regular_file(candidate, ec) && ::access(candidate.c_str(), X_OK) == 0) {
      return candidate;
    }
    if (end == std::string_view::npos) break;
    start = end + 1;
  }
  return std::nullopt;
}

struct OpsinResult {
  bool attempted = false;
  std::optional<std::string> smiles;
};

bool sendAllNoSignal(int fd, std::string_view data) {
  while (!data.empty()) {
    const ssize_t written = ::send(fd, data.data(), data.size(), MSG_NOSIGNAL);
    if (written > 0) {
      data.remove_prefix(static_cast<std::size_t>(written));
      continue;
    }
    if (written < 0 && errno == EINTR) continue;
    return false;
  }
  return true;
}

OpsinResult runOpsin(const std::string& name) {
#ifndef CHEMCAD_OPSIN_JAR
  (void)name;
  return {};
#else
  const fs::path jar = opsinJarPath();
  const auto java = javaPath();
  std::error_code ec;
  if (!java || !fs::is_regular_file(jar, ec)) return {};

  int inputSockets[2] = {-1, -1};
  int outputPipe[2] = {-1, -1};
  if (::socketpair(AF_UNIX, SOCK_STREAM, 0, inputSockets) != 0) return {};
  if (::pipe(outputPipe) != 0) {
    ::close(inputSockets[0]);
    ::close(inputSockets[1]);
    return {};
  }

  const std::string javaString = java->string();
  const std::string jarString = jar.string();
  const pid_t child = ::fork();
  if (child == 0) {
    ::dup2(inputSockets[1], STDIN_FILENO);
    ::dup2(outputPipe[1], STDOUT_FILENO);
    ::close(inputSockets[0]);
    ::close(inputSockets[1]);
    ::close(outputPipe[0]);
    ::close(outputPipe[1]);
    ::execl(javaString.c_str(), javaString.c_str(), "-jar", jarString.c_str(), "-o", "smi",
            static_cast<char*>(nullptr));
    ::_exit(127);
  }

  ::close(inputSockets[1]);
  ::close(outputPipe[1]);
  if (child < 0) {
    ::close(inputSockets[0]);
    ::close(outputPipe[0]);
    return {};
  }

  OpsinResult result;
  result.attempted = true;
  const std::string input = name + '\n';
  sendAllNoSignal(inputSockets[0], input);
  ::close(inputSockets[0]);

  std::string output;
  bool pipeOpen = true;
  bool childDone = false;
  bool timedOut = false;
  int status = 0;
  const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(10);

  while (pipeOpen || !childDone) {
    if (!childDone) {
      pid_t waited;
      do {
        waited = ::waitpid(child, &status, WNOHANG);
      } while (waited < 0 && errno == EINTR);
      childDone = waited == child || (waited < 0 && errno == ECHILD);
    }

    if (pipeOpen) {
      pollfd descriptor{outputPipe[0], POLLIN | POLLHUP, 0};
      const auto remaining = std::chrono::duration_cast<std::chrono::milliseconds>(
          deadline - std::chrono::steady_clock::now());
      const int waitMs = static_cast<int>(std::clamp<long long>(remaining.count(), 0, 100));
      int ready;
      do {
        ready = ::poll(&descriptor, 1, waitMs);
      } while (ready < 0 && errno == EINTR);

      if (ready > 0 && (descriptor.revents & (POLLIN | POLLHUP | POLLERR))) {
        char buffer[4096];
        const ssize_t count = ::read(outputPipe[0], buffer, sizeof(buffer));
        if (count > 0) {
          output.append(buffer, static_cast<std::size_t>(count));
        } else if (count == 0 || (count < 0 && errno != EINTR)) {
          ::close(outputPipe[0]);
          pipeOpen = false;
        }
      }
    } else if (!childDone) {
      ::poll(nullptr, 0, 10);
    }

    if (std::chrono::steady_clock::now() >= deadline && (!childDone || pipeOpen)) {
      timedOut = true;
      break;
    }
  }

  if (pipeOpen) ::close(outputPipe[0]);
  if (!childDone) {
    if (timedOut) ::kill(child, SIGKILL);
    pid_t waited;
    do {
      waited = ::waitpid(child, &status, 0);
    } while (waited < 0 && errno == EINTR);
  }

  if (timedOut || !WIFEXITED(status) || WEXITSTATUS(status) != 0) return result;

  std::size_t lineStart = 0;
  while (lineStart <= output.size()) {
    const std::size_t lineEnd = output.find('\n', lineStart);
    const std::string line = trim(std::string_view(output).substr(lineStart, lineEnd - lineStart));
    if (!line.empty()) {
      result.smiles = line;
      break;
    }
    if (lineEnd == std::string::npos) break;
    lineStart = lineEnd + 1;
  }
  return result;
#endif
}

struct CurlDeleter {
  void operator()(CURL* curl) const {
    if (curl) curl_easy_cleanup(curl);
  }
};
using CurlHandle = std::unique_ptr<CURL, CurlDeleter>;

void initializeCurl() { std::call_once(gCurlInitFlag, [] { curl_global_init(CURL_GLOBAL_DEFAULT); }); }

size_t appendResponse(char* data, size_t size, size_t count, void* destination) noexcept {
  const size_t bytes = size * count;
  try {
    static_cast<std::string*>(destination)->append(data, bytes);
    return bytes;
  } catch (...) {
    return 0;
  }
}

std::optional<std::string> urlEncode(std::string_view value) {
  initializeCurl();
  CurlHandle curl(curl_easy_init());
  if (!curl) return std::nullopt;
  char* encoded = curl_easy_escape(curl.get(), value.data(), static_cast<int>(value.size()));
  if (!encoded) return std::nullopt;
  std::string result(encoded);
  curl_free(encoded);
  return result;
}

struct HttpResponse {
  CURLcode code = CURLE_FAILED_INIT;
  long status = 0;
  std::string body;
};

HttpResponse request(const std::string& url, const std::optional<std::string>& postBody = std::nullopt) {
  initializeCurl();
  HttpResponse response;
  for (int attempt = 0; attempt < 2; ++attempt) {
    CurlHandle curl(curl_easy_init());
    if (!curl) {
      response.code = CURLE_FAILED_INIT;
      continue;
    }

    response.body.clear();
    response.status = 0;
    curl_easy_setopt(curl.get(), CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl.get(), CURLOPT_TIMEOUT, 10L);
    curl_easy_setopt(curl.get(), CURLOPT_CONNECTTIMEOUT, 5L);
    curl_easy_setopt(curl.get(), CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl.get(), CURLOPT_USERAGENT, "chemcad/0.1");
    curl_easy_setopt(curl.get(), CURLOPT_NOSIGNAL, 1L);
    curl_easy_setopt(curl.get(), CURLOPT_WRITEFUNCTION, appendResponse);
    curl_easy_setopt(curl.get(), CURLOPT_WRITEDATA, &response.body);
    if (postBody) {
      curl_easy_setopt(curl.get(), CURLOPT_POST, 1L);
      curl_easy_setopt(curl.get(), CURLOPT_POSTFIELDS, postBody->c_str());
      curl_easy_setopt(curl.get(), CURLOPT_POSTFIELDSIZE,
                       static_cast<long>(postBody->size()));
    }

    response.code = curl_easy_perform(curl.get());
    if (response.code == CURLE_OK) {
      curl_easy_getinfo(curl.get(), CURLINFO_RESPONSE_CODE, &response.status);
      return response;
    }
  }
  return response;
}

std::optional<std::string> propertyValue(const std::string& body, const char* property) {
  try {
    const json root = json::parse(body);
    const auto table = root.find("PropertyTable");
    if (table == root.end() || !table->is_object()) return std::nullopt;
    const auto properties = table->find("Properties");
    if (properties == table->end() || !properties->is_array() || properties->empty()) {
      return std::nullopt;
    }
    const auto value = (*properties)[0].find(property);
    if (value == (*properties)[0].end() || !value->is_string()) return std::nullopt;
    const std::string result = trim(value->get<std::string>());
    if (result.empty()) return std::nullopt;
    return result;
  } catch (...) {
    return std::nullopt;
  }
}

enum class LookupStatus { Success, NotFound, Offline };
struct LookupResult {
  LookupStatus status = LookupStatus::NotFound;
  std::string value;
};

LookupResult pubChemNameToSmiles(const std::string& name) {
  const auto encoded = urlEncode(name);
  if (!encoded) return {LookupStatus::Offline, {}};
  const std::string prefix =
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + *encoded + "/property/";

  HttpResponse response = request(prefix + "CanonicalSMILES/JSON");
  if (response.code != CURLE_OK) return {LookupStatus::Offline, {}};
  if (response.status == 404) return {LookupStatus::NotFound, {}};
  if (response.status != 400) {
    if (auto value = propertyValue(response.body, "CanonicalSMILES")) {
      return {LookupStatus::Success, std::move(*value)};
    }
  }

  response = request(prefix + "SMILES/JSON");
  if (response.code != CURLE_OK) return {LookupStatus::Offline, {}};
  if (response.status == 404) return {LookupStatus::NotFound, {}};
  if (auto value = propertyValue(response.body, "SMILES")) {
    return {LookupStatus::Success, std::move(*value)};
  }
  return {LookupStatus::NotFound, {}};
}

LookupResult pubChemSmilesToName(const std::string& smiles) {
  const auto encoded = urlEncode(smiles);
  if (!encoded) return {LookupStatus::Offline, {}};
  const std::string url =
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/property/IUPACName/JSON";
  const HttpResponse response = request(url, "smiles=" + *encoded);
  if (response.code != CURLE_OK) return {LookupStatus::Offline, {}};
  if (response.status == 404) return {LookupStatus::NotFound, {}};
  if (auto value = propertyValue(response.body, "IUPACName")) {
    return {LookupStatus::Success, std::move(*value)};
  }
  return {LookupStatus::NotFound, {}};
}

}  // namespace

Result nameToSmiles(const std::string& name) {
  std::shared_lock pathLock(gPathGuard);
  const std::string key = nameKey(name);
  if (key.empty()) {
    std::lock_guard cacheLock(gCacheMutex);
    ensureCacheLoadedLocked();
    return {false, "", "enter a name"};
  }
  if (const auto cached = cachedValue(true, key)) return {true, *cached, ""};

  const OpsinResult opsin = runOpsin(trim(name));
  if (opsin.smiles) {
    cacheValue(true, key, *opsin.smiles);
    return {true, *opsin.smiles, ""};
  }

  const LookupResult pubChem = pubChemNameToSmiles(trim(name));
  if (pubChem.status == LookupStatus::Success) {
    cacheValue(true, key, pubChem.value);
    return {true, pubChem.value, ""};
  }
  if (pubChem.status == LookupStatus::Offline) return {false, "", "offline"};
  return {false, "", opsin.attempted ? "could not parse that name" : "not found"};
}

Result smilesToName(const std::string& smiles) {
  std::shared_lock pathLock(gPathGuard);
  const std::string key = trim(smiles);
  if (key.empty()) {
    std::lock_guard cacheLock(gCacheMutex);
    ensureCacheLoadedLocked();
    return {false, "", ""};
  }
  if (const auto cached = cachedValue(false, key)) return {true, *cached, ""};

  const LookupResult pubChem = pubChemSmilesToName(key);
  if (pubChem.status == LookupStatus::Success) {
    cacheValue(false, key, pubChem.value);
    return {true, pubChem.value, ""};
  }
  if (pubChem.status == LookupStatus::Offline) return {false, "", "offline"};
  return {false, "", "not found"};
}

void setCachePath(const std::string& path) {
  std::unique_lock pathLock(gPathGuard);
  std::lock_guard cacheLock(gCacheMutex);
  gCache.path = fs::path(path);
  gCache.pathOverridden = true;
  gCache.loaded = false;
  gCache.n2s.clear();
  gCache.s2n.clear();
}

std::string cachePath() {
  std::shared_lock pathLock(gPathGuard);
  std::lock_guard cacheLock(gCacheMutex);
  return effectiveCachePathLocked().string();
}

bool opsinAvailable() {
#ifdef CHEMCAD_OPSIN_JAR
  std::error_code ec;
  return fs::is_regular_file(opsinJarPath(), ec) && javaPath().has_value();
#else
  return false;
#endif
}

}  // namespace chemcad::naming

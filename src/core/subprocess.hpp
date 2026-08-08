#pragma once
// Child-process execution with captured stdout. Kept in core so the naming
// layer (and anything else that shells out) does not have to carry a POSIX or
// Win32 code path of its own.

#include <chrono>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace chemcad::core {

struct ProcessResult {
  bool launched = false;  // false when the child could not be spawned at all
  bool timedOut = false;  // the child outlived the deadline and was killed
  int exitCode = -1;      // only meaningful when launched && !timedOut
  std::string output;     // everything the child wrote to stdout
};

// Resolves `name` against $PATH. On Windows every $PATHEXT suffix is tried, so
// "java" finds java.exe. Returns nullopt when nothing executable matches.
std::optional<std::filesystem::path> findExecutable(std::string_view name);

// Runs `executable` with `arguments`, feeds `input` to its stdin, closes stdin
// and captures stdout until the child exits or `timeout` elapses (in which case
// the child is killed and timedOut is set). Never throws; stderr is inherited.
ProcessResult runCaptured(const std::filesystem::path& executable,
                          const std::vector<std::string>& arguments,
                          std::string_view input,
                          std::chrono::milliseconds timeout);

}  // namespace chemcad::core

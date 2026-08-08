#include "core/subprocess.hpp"

#include <algorithm>
#include <cstdlib>
#include <system_error>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

#include <future>
#include <type_traits>
#else
#include <cerrno>
#include <csignal>
#include <poll.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace chemcad::core {
namespace fs = std::filesystem;
namespace {

#ifdef _WIN32
constexpr char kPathSeparator = ';';
#else
constexpr char kPathSeparator = ':';
#endif

std::vector<std::string> splitList(std::string_view value, char separator) {
  std::vector<std::string> parts;
  std::size_t start = 0;
  while (start <= value.size()) {
    const std::size_t end = value.find(separator, start);
    std::string_view part = value.substr(start, end - start);
    if (part.size() >= 2 && part.front() == '"' && part.back() == '"') {
      part = part.substr(1, part.size() - 2);
    }
    if (!part.empty()) parts.emplace_back(part);
    if (end == std::string_view::npos) break;
    start = end + 1;
  }
  return parts;
}

std::vector<fs::path> searchPath() {
  std::vector<fs::path> directories;
  const char* value = std::getenv("PATH");
  if (!value) return directories;
  for (std::string& part : splitList(value, kPathSeparator)) {
    directories.emplace_back(std::move(part));
  }
  return directories;
}

#ifdef _WIN32

std::wstring widen(std::string_view value) {
  if (value.empty()) return {};
  const int length = ::MultiByteToWideChar(CP_UTF8, 0, value.data(),
                                           static_cast<int>(value.size()), nullptr, 0);
  if (length <= 0) return {};
  std::wstring result(static_cast<std::size_t>(length), L'\0');
  ::MultiByteToWideChar(CP_UTF8, 0, value.data(), static_cast<int>(value.size()), result.data(),
                        length);
  return result;
}

// CommandLineToArgvW quoting: a run of backslashes is only special when it
// precedes a quote, and then it has to be doubled.
std::wstring quoteArgument(const std::wstring& argument) {
  if (!argument.empty() && argument.find_first_of(L" \t\n\v\"") == std::wstring::npos) {
    return argument;
  }
  std::wstring quoted(1, L'"');
  for (auto it = argument.begin();; ++it) {
    std::size_t backslashes = 0;
    while (it != argument.end() && *it == L'\\') {
      ++it;
      ++backslashes;
    }
    if (it == argument.end()) {
      quoted.append(backslashes * 2, L'\\');
      break;
    }
    quoted.append(*it == L'"' ? backslashes * 2 + 1 : backslashes, L'\\');
    quoted.push_back(*it);
  }
  quoted.push_back(L'"');
  return quoted;
}

class Handle {
 public:
  Handle() = default;
  explicit Handle(HANDLE handle) : handle_(handle) {}
  Handle(const Handle&) = delete;
  Handle& operator=(const Handle&) = delete;
  ~Handle() { reset(); }

  void reset(HANDLE handle = nullptr) {
    if (handle_ && handle_ != INVALID_HANDLE_VALUE) ::CloseHandle(handle_);
    handle_ = handle;
  }
  [[nodiscard]] HANDLE get() const { return handle_; }
  [[nodiscard]] bool valid() const { return handle_ && handle_ != INVALID_HANDLE_VALUE; }
  HANDLE* address() { return &handle_; }

 private:
  HANDLE handle_ = nullptr;
};

#else

bool writeAll(int descriptor, std::string_view data) {
  while (!data.empty()) {
    const ssize_t written = ::send(descriptor, data.data(), data.size(), MSG_NOSIGNAL);
    if (written > 0) {
      data.remove_prefix(static_cast<std::size_t>(written));
      continue;
    }
    if (written < 0 && errno == EINTR) continue;
    return false;
  }
  return true;
}

#endif

}  // namespace

std::optional<fs::path> findExecutable(std::string_view name) {
  if (name.empty()) return std::nullopt;
  std::error_code ec;

#ifdef _WIN32
  std::vector<std::string> suffixes;
  if (const char* pathExt = std::getenv("PATHEXT"); pathExt && *pathExt) {
    suffixes = splitList(pathExt, ';');
  }
  if (suffixes.empty()) suffixes = {".COM", ".EXE", ".BAT", ".CMD"};
  const bool alreadySuffixed = fs::path(name).has_extension();

  for (const fs::path& directory : searchPath()) {
    if (alreadySuffixed && fs::is_regular_file(directory / name, ec)) return directory / name;
    for (const std::string& suffix : suffixes) {
      const fs::path candidate = directory / (std::string(name) + suffix);
      if (fs::is_regular_file(candidate, ec)) return candidate;
    }
  }
#else
  for (const fs::path& directory : searchPath()) {
    const fs::path candidate = directory / name;
    if (fs::is_regular_file(candidate, ec) && ::access(candidate.c_str(), X_OK) == 0) {
      return candidate;
    }
  }
#endif
  return std::nullopt;
}

#ifdef _WIN32

ProcessResult runCaptured(const fs::path& executable, const std::vector<std::string>& arguments,
                          std::string_view input, std::chrono::milliseconds timeout) {
  SECURITY_ATTRIBUTES inheritable{sizeof(SECURITY_ATTRIBUTES), nullptr, TRUE};
  Handle childInput, parentInput, parentOutput, childOutput;
  if (!::CreatePipe(childInput.address(), parentInput.address(), &inheritable, 0)) return {};
  if (!::CreatePipe(parentOutput.address(), childOutput.address(), &inheritable, 0)) return {};
  // Only the child ends may cross the process boundary.
  ::SetHandleInformation(parentInput.get(), HANDLE_FLAG_INHERIT, 0);
  ::SetHandleInformation(parentOutput.get(), HANDLE_FLAG_INHERIT, 0);

  // A GUI process has no stderr, and STARTF_USESTDHANDLES with a null handle
  // hands the child a broken descriptor. NUL always works and never fills up.
  Handle nullDevice(::CreateFileW(L"NUL", GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE,
                                  &inheritable, OPEN_EXISTING, 0, nullptr));

  std::wstring commandLine = quoteArgument(executable.wstring());
  for (const std::string& argument : arguments) {
    commandLine.push_back(L' ');
    commandLine += quoteArgument(widen(argument));
  }
  std::vector<wchar_t> mutableCommandLine(commandLine.begin(), commandLine.end());
  mutableCommandLine.push_back(L'\0');

  STARTUPINFOW startup{};
  startup.cb = sizeof(startup);
  startup.dwFlags = STARTF_USESTDHANDLES;
  startup.hStdInput = childInput.get();
  startup.hStdOutput = childOutput.get();
  startup.hStdError = nullDevice.get();

  PROCESS_INFORMATION created{};
  const std::wstring executablePath = executable.wstring();
  if (!::CreateProcessW(executablePath.c_str(), mutableCommandLine.data(), nullptr, nullptr, TRUE,
                        CREATE_NO_WINDOW, nullptr, nullptr, &startup, &created)) {
    return {};
  }
  Handle process(created.hProcess);
  Handle thread(created.hThread);
  childInput.reset();
  childOutput.reset();
  nullDevice.reset();

  // Callers feed a single short line, so a blocking write cannot deadlock
  // against a child that has not started reading yet.
  while (!input.empty()) {
    DWORD written = 0;
    const DWORD chunk = static_cast<DWORD>(std::min<std::size_t>(input.size(), 1u << 16));
    if (!::WriteFile(parentInput.get(), input.data(), chunk, &written, nullptr) || written == 0) {
      break;
    }
    input.remove_prefix(written);
  }
  parentInput.reset();

  const HANDLE readEnd = parentOutput.get();
  std::future<std::string> pending = std::async(std::launch::async, [readEnd] {
    std::string collected;
    char buffer[4096];
    DWORD count = 0;
    while (::ReadFile(readEnd, buffer, sizeof(buffer), &count, nullptr) && count > 0) {
      collected.append(buffer, count);
    }
    return collected;
  });

  ProcessResult result;
  result.launched = true;
  if (pending.wait_for(timeout) == std::future_status::timeout) {
    // Killing the child drops its pipe end, which unblocks the reader.
    ::TerminateProcess(process.get(), 1);
    result.timedOut = true;
  }
  result.output = pending.get();

  ::WaitForSingleObject(process.get(), INFINITE);
  DWORD exitCode = 0;
  if (::GetExitCodeProcess(process.get(), &exitCode)) result.exitCode = static_cast<int>(exitCode);
  return result;
}

#else

ProcessResult runCaptured(const fs::path& executable, const std::vector<std::string>& arguments,
                          std::string_view input, std::chrono::milliseconds timeout) {
  int inputSockets[2] = {-1, -1};
  int outputPipe[2] = {-1, -1};
  if (::socketpair(AF_UNIX, SOCK_STREAM, 0, inputSockets) != 0) return {};
  if (::pipe(outputPipe) != 0) {
    ::close(inputSockets[0]);
    ::close(inputSockets[1]);
    return {};
  }

  // Everything the child needs must exist before the fork: only async-signal
  // safe calls are legal between fork() and exec() in a threaded process.
  const std::string executableString = executable.string();
  std::vector<char*> argv;
  argv.reserve(arguments.size() + 2);
  argv.push_back(const_cast<char*>(executableString.c_str()));
  for (const std::string& argument : arguments) argv.push_back(const_cast<char*>(argument.c_str()));
  argv.push_back(nullptr);

  const pid_t child = ::fork();
  if (child == 0) {
    ::dup2(inputSockets[1], STDIN_FILENO);
    ::dup2(outputPipe[1], STDOUT_FILENO);
    ::close(inputSockets[0]);
    ::close(inputSockets[1]);
    ::close(outputPipe[0]);
    ::close(outputPipe[1]);
    ::execv(executableString.c_str(), argv.data());
    ::_exit(127);
  }

  ::close(inputSockets[1]);
  ::close(outputPipe[1]);
  if (child < 0) {
    ::close(inputSockets[0]);
    ::close(outputPipe[0]);
    return {};
  }

  ProcessResult result;
  result.launched = true;
  writeAll(inputSockets[0], input);
  ::close(inputSockets[0]);

  bool pipeOpen = true;
  bool childDone = false;
  int status = 0;
  const auto deadline = std::chrono::steady_clock::now() + timeout;

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
          result.output.append(buffer, static_cast<std::size_t>(count));
        } else if (count == 0 || (count < 0 && errno != EINTR)) {
          ::close(outputPipe[0]);
          pipeOpen = false;
        }
      }
    } else if (!childDone) {
      ::poll(nullptr, 0, 10);
    }

    if (std::chrono::steady_clock::now() >= deadline && (!childDone || pipeOpen)) {
      result.timedOut = true;
      break;
    }
  }

  if (pipeOpen) ::close(outputPipe[0]);
  if (!childDone) {
    if (result.timedOut) ::kill(child, SIGKILL);
    pid_t waited;
    do {
      waited = ::waitpid(child, &status, 0);
    } while (waited < 0 && errno == EINTR);
  }
  if (!result.timedOut && WIFEXITED(status)) result.exitCode = WEXITSTATUS(status);
  return result;
}

#endif

}  // namespace chemcad::core

#pragma once
// Tiny worker pool that keeps blocking work (PubChem lookups, OPSIN, route
// search, LLM calls) off the frame thread.
//
// Usage from the UI:
//   state.tasks.run<naming::Result>(
//       [smiles] { return naming::smilesToName(smiles); },
//       [&state](naming::Result r) { state.props.name = r.value; });
//
// The completion runs on the UI thread inside pump(), so it may touch AppState
// freely. R must be default-constructible: if the worker throws, the completion
// receives a default-constructed R (task bodies are expected to return their own
// error representation rather than throw).

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <functional>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

namespace chemcad::app {

class TaskRunner {
 public:
  explicit TaskRunner(unsigned workers = 3);
  ~TaskRunner();
  TaskRunner(const TaskRunner&) = delete;
  TaskRunner& operator=(const TaskRunner&) = delete;

  template <class R>
  void run(std::function<R()> fn, std::function<void(R)> done) {
    inflight_.fetch_add(1, std::memory_order_relaxed);
    submit([this, fn = std::move(fn), done = std::move(done)]() mutable {
      R result{};
      try {
        result = fn();
      } catch (...) {
        result = R{};
      }
      {
        std::lock_guard<std::mutex> lk(completionMu_);
        completions_.push_back(
            [done = std::move(done), result = std::move(result)]() mutable {
              done(std::move(result));
            });
      }
      inflight_.fetch_sub(1, std::memory_order_release);
    });
  }

  // Call once per frame on the UI thread.
  void pump();

  bool busy() const { return inflight_.load(std::memory_order_acquire) > 0; }
  size_t inflight() const {
    return static_cast<size_t>(inflight_.load(std::memory_order_acquire));
  }

 private:
  void submit(std::function<void()> job);
  void workerLoop();

  std::vector<std::thread> workers_;
  std::vector<std::function<void()>> queue_;
  std::mutex queueMu_;
  std::condition_variable queueCv_;
  bool stopping_ = false;

  std::vector<std::function<void()>> completions_;
  std::mutex completionMu_;
  std::atomic<int> inflight_{0};
};

}  // namespace chemcad::app

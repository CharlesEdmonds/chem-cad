#include "app/task_runner.hpp"

namespace chemcad::app {

TaskRunner::TaskRunner(unsigned workers) {
  if (workers == 0) workers = 1;
  workers_.reserve(workers);
  for (unsigned i = 0; i < workers; ++i)
    workers_.emplace_back([this] { workerLoop(); });
}

TaskRunner::~TaskRunner() {
  {
    std::lock_guard<std::mutex> lk(queueMu_);
    stopping_ = true;
  }
  queueCv_.notify_all();
  for (std::thread& t : workers_)
    if (t.joinable()) t.join();
}

void TaskRunner::submit(std::function<void()> job) {
  {
    std::lock_guard<std::mutex> lk(queueMu_);
    queue_.push_back(std::move(job));
  }
  queueCv_.notify_one();
}

void TaskRunner::workerLoop() {
  for (;;) {
    std::function<void()> job;
    {
      std::unique_lock<std::mutex> lk(queueMu_);
      queueCv_.wait(lk, [this] { return stopping_ || !queue_.empty(); });
      if (stopping_ && queue_.empty()) return;
      job = std::move(queue_.front());
      queue_.erase(queue_.begin());
    }
    job();
  }
}

void TaskRunner::pump() {
  std::vector<std::function<void()>> ready;
  {
    std::lock_guard<std::mutex> lk(completionMu_);
    ready.swap(completions_);
  }
  for (std::function<void()>& c : ready) c();
}

}  // namespace chemcad::app

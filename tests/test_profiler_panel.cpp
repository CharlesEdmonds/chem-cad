// The profiler, exercised rather than trusted.
//
// It shipped wired to a View menu item and F2 and had never been opened, which
// is the same standing the GPU backend had when its shader turned out not to
// compile. This drives the core through nested zones on two threads and then
// draws the real panel headlessly, in both its empty and its populated state,
// across the display shapes the layout contract already covers.
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "imgui.h"

#include "core/profiler.hpp"
#include "ui/app_state.hpp"
#include "ui/ui.hpp"

using namespace chemcad;

namespace {

struct HeadlessImGui {
  HeadlessImGui(float width, float height, float fontScale) {
    IMGUI_CHECKVERSION();
    ctx = ImGui::CreateContext();
    ImGui::SetCurrentContext(ctx);
    ImGuiIO& io = ImGui::GetIO();
    io.DisplaySize = ImVec2(width, height);
    io.DeltaTime = 1.0f / 60.0f;
    io.IniFilename = nullptr;
    io.LogFilename = nullptr;
    io.FontGlobalScale = fontScale;
    unsigned char* pixels = nullptr;
    int w = 0;
    int h = 0;
    io.Fonts->GetTexDataAsRGBA32(&pixels, &w, &h);
    io.Fonts->SetTexID(static_cast<ImTextureID>(1));
  }
  ~HeadlessImGui() { ImGui::DestroyContext(ctx); }
  HeadlessImGui(const HeadlessImGui&) = delete;
  HeadlessImGui& operator=(const HeadlessImGui&) = delete;
  ImGuiContext* ctx = nullptr;
};

// Busy-waits rather than sleeping: the profiler measures CPU time inside a
// zone, and a sleeping thread would report a zone that took no work at all.
void burn(double milliseconds) {
  const auto deadline =
      std::chrono::steady_clock::now() + std::chrono::duration<double, std::milli>(milliseconds);
  volatile double sink = 0.0;
  while (std::chrono::steady_clock::now() < deadline) sink += 1.0;
  (void)sink;
}

const core::ProfileZone* findZone(const core::ProfileFrame& frame, const char* name) {
  for (const core::ProfileZone& zone : frame.zones) {
    if (zone.name == name) return &zone;
  }
  return nullptr;
}

float drawOnce(ui::AppState& state) {
  float overflow = 0.0f;
  // Two frames, matching the layout suite: docking and animated values settle
  // on the second, and a panel that only overflows once settled must not pass.
  for (int frame = 0; frame < 2; ++frame) {
    ImGui::NewFrame();
    const ImGuiIO& io = ImGui::GetIO();
    ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f));
    ImGui::SetNextWindowSize(io.DisplaySize);
    if (ImGui::Begin("Performance", nullptr,
                     ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                         ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings)) {
      ui::drawProfiler(state);
      overflow = ImGui::GetScrollMaxY();
    }
    ImGui::End();
    ImGui::Render();
  }
  return overflow;
}

// Records `count` frames of a small, deliberately nested workload.
void recordFrames(int count) {
  core::Profiler& performance = core::profiler();
  performance.setEnabled(true);
  for (int frame = 0; frame < count; ++frame) {
    performance.beginFrame();
    {
      CHEMCAD_PROFILE_ZONE("outer");
      burn(0.4);
      {
        CHEMCAD_PROFILE_ZONE("inner");
        burn(0.6);
      }
    }
    performance.counter("particles", 512.0 + frame);
    performance.endFrame();
  }
}

}  // namespace

TEST_CASE("the profiler separates inclusive from exclusive time and keeps its ring") {
  core::Profiler& performance = core::profiler();
  performance.reset();
  recordFrames(4);

  std::vector<core::ProfileFrame> frames = performance.frames();
  REQUIRE(frames.size() == 4);
  const core::ProfileFrame& last = frames.back();
  const core::ProfileZone* outer = findZone(last, "outer");
  const core::ProfileZone* inner = findZone(last, "inner");
  REQUIRE(outer != nullptr);
  REQUIRE(inner != nullptr);

  CHECK(outer->calls == 1);
  CHECK(inner->calls == 1);
  // The distinction the panel's table exists to show: a parent's inclusive time
  // contains its children, its exclusive time does not.
  CHECK(inner->inclusiveMs > 0.0);
  CHECK(outer->inclusiveMs >= inner->inclusiveMs);
  CHECK(outer->exclusiveMs < outer->inclusiveMs);
  CHECK(outer->exclusiveMs + inner->inclusiveMs ==
        doctest::Approx(outer->inclusiveMs).epsilon(0.02));
  CHECK(inner->exclusiveMs == doctest::Approx(inner->inclusiveMs).epsilon(0.02));
  CHECK(last.frameMs >= outer->inclusiveMs);
  REQUIRE(last.counters.size() == 1);
  CHECK(last.counters[0].name == "particles");
  CHECK(last.counters[0].value == doctest::Approx(515.0));

  // Sequence numbers are a lifetime frame count, not a ring index: reset()
  // drops the retained samples but keeps counting, so the frame-time plot's x
  // axis still says which frame of the session it is showing.
  CHECK(frames.front().sequence + 3 == frames.back().sequence);
  const std::uint64_t beforeReset = frames.back().sequence;
  performance.reset();
  CHECK(performance.frames().empty());

  recordFrames(static_cast<int>(core::Profiler::frameCapacity) + 12);
  frames = performance.frames();
  CHECK(frames.size() == core::Profiler::frameCapacity);
  CHECK(frames.front().sequence > beforeReset);
  // Contiguous after the wrap, which is the property the plot depends on.
  CHECK(frames.back().sequence - frames.front().sequence ==
        core::Profiler::frameCapacity - 1);
  for (std::size_t i = 1; i < frames.size(); ++i) {
    REQUIRE(frames[i].sequence == frames[i - 1].sequence + 1);
  }

  // Disabled means disabled: no zone recorded and no frame retained.
  performance.reset();
  performance.setEnabled(false);
  performance.beginFrame();
  {
    CHEMCAD_PROFILE_ZONE("ignored");
    burn(0.2);
  }
  performance.counter("ignored", 1.0);
  performance.endFrame();
  CHECK(performance.frames().empty());
  performance.setEnabled(true);
}

TEST_CASE("the profiler CSV round-trips every zone of every frame") {
  core::Profiler& performance = core::profiler();
  performance.reset();
  recordFrames(3);

  const std::string csv = performance.exportCsv();
  const std::filesystem::path destination =
      std::filesystem::temp_directory_path() / "chemcad-profile-test.csv";
  {
    std::ofstream output(destination, std::ios::binary | std::ios::trunc);
    REQUIRE(output.good());
    output << csv;
    REQUIRE(output.good());
  }

  std::ifstream input(destination, std::ios::binary);
  REQUIRE(input.good());
  std::string header;
  REQUIRE(std::getline(input, header));
  CHECK(header == "frame,frame_ms,cpu_ms,zone,parent,calls,inclusive_ms,exclusive_ms,share_percent");

  std::size_t rows = 0;
  std::size_t innerRows = 0;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    ++rows;
    // Nine fields, and the parent column is what makes the nesting recoverable
    // from the file rather than only on screen.
    CHECK(std::count(line.begin(), line.end(), ',') == 8);
    if (line.find(",inner,outer,") != std::string::npos) ++innerRows;
  }
  input.close();
  std::filesystem::remove(destination);

  CHECK(rows == 6);  // three frames, two zones each
  CHECK(innerRows == 3);
}

TEST_CASE("the profiler panel draws and fits its page, empty and populated") {
  const struct {
    const char* name;
    float width;
    float height;
  } surfaces[] = {
      {"1366x768", 1366.0f, 768.0f},
      {"1920x1080", 1920.0f, 1080.0f},
      {"2560x1080", 2560.0f, 1080.0f},
      {"1280x1024", 1280.0f, 1024.0f},
  };
  const float scales[] = {1.0f, 1.25f, 1.75f};

  for (const auto& surface : surfaces) {
    for (const float scale : scales) {
      // Empty first: the panel has to survive being opened before anything has
      // been recorded, which is exactly how a user meets it.
      core::profiler().reset();
      {
        HeadlessImGui gui(surface.width, surface.height, scale);
        ui::AppState state;
        const float overflow = drawOnce(state);
        INFO("empty profiler at " << surface.name << " scale " << scale << " overflowed by "
                                  << overflow << " px");
        CHECK(overflow == 0.0f);
      }

      recordFrames(24);
      {
        HeadlessImGui gui(surface.width, surface.height, scale);
        ui::AppState state;
        const float overflow = drawOnce(state);
        INFO("populated profiler at " << surface.name << " scale " << scale << " overflowed by "
                                      << overflow << " px");
        CHECK(overflow == 0.0f);
      }
    }
  }
  core::profiler().reset();
}

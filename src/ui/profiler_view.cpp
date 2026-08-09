#include "ui/ui.hpp"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "core/paths.hpp"
#include "core/profiler.hpp"
#include "ui/charts.hpp"
#include "ui/icons.hpp"
#include "ui/layout.hpp"
#include "ui/theme.hpp"
#include "ui/widgets.hpp"

namespace chemcad::ui {
namespace {

struct CounterHistory {
  explicit CounterHistory(std::string counterName)
      : name(std::move(counterName)), trace(core::Profiler::frameCapacity) {}

  std::string name;
  charts::Trace trace;
  double latest = 0.0;
};

std::vector<CounterHistory> buildCounterHistories(
    const std::vector<core::ProfileFrame>& frames) {
  std::vector<CounterHistory> histories;
  for (const core::ProfileFrame& frame : frames) {
    for (const core::ProfileCounter& counter : frame.counters) {
      const auto found = std::find_if(histories.begin(), histories.end(),
                                      [&](const CounterHistory& history) {
                                        return history.name == counter.name;
                                      });
      if (found == histories.end()) histories.emplace_back(counter.name);
    }
  }

  for (const core::ProfileFrame& frame : frames) {
    for (CounterHistory& history : histories) {
      const auto found = std::find_if(frame.counters.begin(), frame.counters.end(),
                                      [&](const core::ProfileCounter& counter) {
                                        return counter.name == history.name;
                                      });
      if (found != frame.counters.end()) {
        history.trace.push(found->value);
        history.latest = found->value;
      }
    }
  }
  return histories;
}

void drawZoneRow(const core::ProfileFrame& frame,
                 const std::vector<std::vector<int>>& children, int zoneIndex, int depth,
                 std::uint64_t selectedZone) {
  const core::ProfileZone& zone = frame.zones[static_cast<std::size_t>(zoneIndex)];
  if (zone.id == selectedZone) {
    widgets::dataRow(style::col::Accent);
  } else {
    widgets::dataRow();
  }

  bool open = true;
  if (children[static_cast<std::size_t>(zoneIndex)].empty()) {
    std::string indented(static_cast<std::size_t>(depth) * 2, ' ');
    indented += zone.name;
    widgets::dataCell(indented.c_str());
  } else {
    ImGui::TableNextColumn();
    const float indent = style::metrics().gap * static_cast<float>(depth);
    if (indent > 0.0f) ImGui::Indent(indent);
    const std::string id = "##profiler_zone_" + std::to_string(zone.id);
    open = widgets::disclosure(id.c_str(), zone.name.c_str(), "", depth == 0,
                               icons::Icon::None, style::col::Data);
    if (indent > 0.0f) ImGui::Unindent(indent);
  }

  widgets::dataCellf("%u", zone.calls);
  widgets::dataCellf("%.3f", zone.inclusiveMs);
  widgets::dataCellf("%.3f", zone.exclusiveMs);
  const double share = frame.cpuMs > 0.0 ? zone.exclusiveMs * 100.0 / frame.cpuMs : 0.0;
  widgets::dataCellf("%.1f", share);

  if (!open) return;
  for (const int child : children[static_cast<std::size_t>(zoneIndex)]) {
    drawZoneRow(frame, children, child, depth + 1, selectedZone);
  }
}

void drawZoneTable(const core::ProfileFrame& frame, const layout::Frame& panel,
                   std::uint64_t selectedZone) {
  std::vector<std::vector<int>> children(frame.zones.size());
  std::vector<int> roots;
  for (std::size_t i = 0; i < frame.zones.size(); ++i) {
    const int parent = frame.zones[i].parent;
    if (parent >= 0 && parent < static_cast<int>(frame.zones.size()) &&
        parent != static_cast<int>(i)) {
      children[static_cast<std::size_t>(parent)].push_back(static_cast<int>(i));
    } else {
      roots.push_back(static_cast<int>(i));
    }
  }

  const auto byExclusive = [&](int left, int right) {
    return frame.zones[static_cast<std::size_t>(left)].exclusiveMs >
           frame.zones[static_cast<std::size_t>(right)].exclusiveMs;
  };
  std::sort(roots.begin(), roots.end(), byExclusive);
  for (std::vector<int>& siblings : children) std::sort(siblings.begin(), siblings.end(), byExclusive);

  const widgets::Column columns[] = {
      {"Zone", false, true, nullptr, 12.0f},
      {"Calls", true, false, nullptr, 4.0f},
      {"Inclusive", true, false, "ms", 6.0f},
      {"Exclusive", true, false, "ms", 6.0f},
      {"Share", true, false, "%", 5.0f},
  };
  const float visibleRows = static_cast<float>(std::min<std::size_t>(frame.zones.size() + 1, 14));
  if (widgets::beginDataTable("##profiler_zones", columns, 5,
                              ImVec2(ImGui::GetContentRegionAvail().x,
                                     panel.row * visibleRows))) {
    for (const int root : roots) drawZoneRow(frame, children, root, 0, selectedZone);
    widgets::endDataTable();
  }
}

}  // namespace

void drawProfiler(AppState& state) {
  (void)state;
  core::Profiler& performance = core::profiler();
  static std::string exportMessage;
  static bool exportFailed = false;
  static std::uint64_t selectedZone = 0;

  const layout::Frame panel = layout::measure();
  widgets::sectionHeader("PERFORMANCE PROFILER", style::col::Data);

  widgets::beginToolbar("##profiler_toolbar");
  const bool recording = performance.enabled();
  const icons::Icon transportIcon = recording ? icons::Icon::Pause : icons::Icon::Play;
  const char* transportLabel = recording ? "Pause" : "Record";
  if (widgets::actionButton("##profiler_record", transportIcon, transportLabel,
                            ImVec2(0.0f, panel.control), true,
                            recording ? "Pause performance capture"
                                      : "Record performance capture")) {
    performance.setEnabled(!recording);
  }
  ImGui::SameLine();
  widgets::toolbarSeparator();
  ImGui::SameLine();
  if (widgets::actionButton("##profiler_reset", icons::Icon::Rewind, "Reset",
                            ImVec2(0.0f, panel.control), false,
                            "Clear retained profiler frames")) {
    performance.reset();
    selectedZone = 0;
    exportMessage.clear();
  }
  ImGui::SameLine();
  if (widgets::actionButton("##profiler_export", icons::Icon::Save, "Export CSV",
                            ImVec2(0.0f, panel.control), false,
                            "Write the retained zone samples to the cache directory")) {
    const std::filesystem::path destination = core::cacheDir() / "chemcad-profile.csv";
    std::ofstream output(destination, std::ios::binary | std::ios::trunc);
    output << performance.exportCsv();
    exportFailed = !output.good();
    exportMessage = exportFailed ? "Could not write profiler CSV: " + destination.string()
                                 : "Profiler CSV: " + destination.string();
  }
  widgets::endToolbar();

  if (!exportMessage.empty()) {
    widgets::notice(exportFailed ? icons::Icon::Warning : icons::Icon::Save,
                    exportMessage.c_str(),
                    exportFailed ? style::col::Danger : style::col::Data);
  }

  const std::vector<core::ProfileFrame> frames = performance.frames();
  if (frames.empty()) {
    widgets::emptyState(icons::Icon::Gauge, "No performance samples",
                        "Start recording to capture frame time, scoped CPU zones and counters.");
    return;
  }

  std::vector<double> frameX(frames.size());
  std::vector<double> frameTimes(frames.size());
  for (std::size_t i = 0; i < frames.size(); ++i) {
    frameX[i] = static_cast<double>(i);
    frameTimes[i] = frames[i].frameMs;
  }
  const charts::Series frameSeries = {
      "Frame", frameTimes.data(), frameX.data(), static_cast<int>(frameTimes.size()),
      style::col::Data, true, false, false};
  charts::PlotStyle framePlot;
  framePlot.xLabel = "frame time (ms)";
  framePlot.yLabel = "retained frame";
  framePlot.legend = false;
  framePlot.grid = true;
  framePlot.hasCursor = true;
  framePlot.cursorX = 16.7;
  framePlot.xMin = 0.0;
  framePlot.xMax = std::max(16.7, *std::max_element(frameTimes.begin(), frameTimes.end())) * 1.1;

  if (widgets::beginCard("##profiler_timeline", ImVec2(panel.size.x, 0.0f),
                         style::col::BgSurface)) {
    widgets::cardHeader(icons::Icon::ChartLine, "Frame time",
                        "Frames right of the 16.7 ms cursor exceed budget", style::col::Data);
    charts::linePlot("##profiler_frame_plot", &frameSeries, 1,
                     ImVec2(ImGui::GetContentRegionAvail().x, panel.row * 9.0f), framePlot);
    widgets::endCard();
  }

  const core::ProfileFrame& latest = frames.back();
  std::vector<int> rankedIndices(latest.zones.size());
  for (std::size_t i = 0; i < latest.zones.size(); ++i) rankedIndices[i] = static_cast<int>(i);
  std::sort(rankedIndices.begin(), rankedIndices.end(), [&](int left, int right) {
    return latest.zones[static_cast<std::size_t>(left)].exclusiveMs >
           latest.zones[static_cast<std::size_t>(right)].exclusiveMs;
  });
  rankedIndices.erase(
      std::remove_if(rankedIndices.begin(), rankedIndices.end(), [&](int index) {
        return latest.zones[static_cast<std::size_t>(index)].calls == 0;
      }),
      rankedIndices.end());
  if (rankedIndices.size() > 8) rankedIndices.resize(8);

  std::vector<std::string> annotations(rankedIndices.size());
  std::vector<charts::BarRow> bars;
  bars.reserve(rankedIndices.size());
  for (std::size_t i = 0; i < rankedIndices.size(); ++i) {
    const core::ProfileZone& zone = latest.zones[static_cast<std::size_t>(rankedIndices[i])];
    char value[48];
    std::snprintf(value, sizeof(value), "%.3f ms", zone.exclusiveMs);
    annotations[i] = value;
    bars.push_back(charts::BarRow{zone.name.c_str(), zone.exclusiveMs,
                                  annotations[i].c_str(), style::col::Data,
                                  zone.id == selectedZone});
  }

  if (!bars.empty() &&
      widgets::beginCard("##profiler_hotspots", ImVec2(panel.size.x, 0.0f),
                         style::col::BgSurface)) {
    widgets::cardHeader(icons::Icon::ChartBars, "CPU hot spots",
                        "Latest frame, ranked by exclusive time", style::col::Data);
    const int selected = charts::rankedBars(
        "##profiler_ranked", bars.data(), static_cast<int>(bars.size()),
        ImVec2(ImGui::GetContentRegionAvail().x,
               panel.row * static_cast<float>(bars.size())));
    if (selected >= 0) {
      selectedZone = latest.zones[static_cast<std::size_t>(rankedIndices[static_cast<std::size_t>(selected)])].id;
    }
    widgets::endCard();
  }

  if (widgets::beginCard("##profiler_zone_tree", ImVec2(panel.size.x, 0.0f),
                         style::col::BgSurface)) {
    widgets::cardHeader(icons::Icon::Table, "CPU zone tree",
                        "Inclusive and direct-child-exclusive timings", style::col::Data);
    drawZoneTable(latest, panel, selectedZone);
    widgets::endCard();
  }

  std::vector<CounterHistory> histories = buildCounterHistories(frames);
  if (!histories.empty()) {
    widgets::sectionHeader("COUNTERS", style::col::Data);
    const layout::Frame counterFrame = layout::measure(ImGui::GetContentRegionAvail());
    const int columns = std::min<int>(layout::columnsThatFit(counterFrame, 15.0f),
                                      static_cast<int>(histories.size()));
    const float tileWidth = layout::columnWidth(counterFrame, columns);
    charts::SparklineStyle sparkStyle;
    sparkStyle.accent = style::col::Data;
    sparkStyle.fill = true;
    sparkStyle.showLatest = true;
    for (std::size_t i = 0; i < histories.size(); ++i) {
      if (i % static_cast<std::size_t>(columns) != 0) ImGui::SameLine();
      CounterHistory& history = histories[i];
      char value[64];
      std::snprintf(value, sizeof(value), "%.4g", history.latest);
      const std::string id = "##profiler_counter_" + history.name;
      charts::instrument(id.c_str(), history.name.c_str(), value, nullptr, history.trace,
                         ImVec2(tileWidth, panel.row * 5.0f), sparkStyle);
    }
  }
}

}  // namespace chemcad::ui

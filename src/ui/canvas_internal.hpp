#pragma once

#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

#include "imgui.h"

#include "ui/app_state.hpp"

namespace chemcad::ui::canvas {

struct CanvasRect {
  core::Vec2 origin;
  core::Vec2 size;
  ImVec2 min;
  ImVec2 max;
};

enum class Gesture {
  None,
  Bond,
  Chain,
  Marquee,
  Move,
};

struct MoveOrigin {
  AtomRef ref;
  core::Vec2 pos;
};

struct RingGeometry {
  std::vector<core::Vec2> positions;
  std::vector<std::pair<int, int>> edges;
  std::vector<core::BondOrder> orders;
  std::vector<AtomRef> reuse;
};

struct Runtime {
  Gesture gesture = Gesture::None;
  core::Vec2 downWorld;
  core::Vec2 currentWorld;
  ImVec2 downScreen{};
  ImVec2 currentScreen{};
  AtomRef downAtom;
  BondRef downBond;
  bool dragged = false;
  bool shiftAtStart = false;
  std::vector<MoveOrigin> moveOrigins;

  uint64_t hydrogenRevision = ~0ull;
  std::unordered_map<uint64_t, int> hydrogenCounts;
};

void hitTest(AppState& st, const CanvasRect& rect, bool canvasHovered);
void handleInput(AppState& st, Runtime& rt, const CanvasRect& rect, bool canvasHovered,
                 bool canvasActive);
void render(AppState& st, Runtime& rt, const CanvasRect& rect);

RingGeometry makeRingGeometry(const AppState& st, core::Vec2 cursor);
std::vector<core::Vec2> makeChainPreview(const AppState& st, const Runtime& rt);

}  // namespace chemcad::ui::canvas

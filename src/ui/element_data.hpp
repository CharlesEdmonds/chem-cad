#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace chemcad::ui {

struct ElementData {
  uint8_t z = 0;
  std::string symbol;
  std::string name;
  double mass = 0;
  int group = 0;
  int period = 0;
  std::string category;
};

// The returned table is loaded once. A built-in subset is returned when the
// external table is unavailable so element picking never becomes unusable.
const std::vector<ElementData>& elementTable();
const std::string& elementTableLoadError();
const ElementData* findElement(uint8_t z);

}  // namespace chemcad::ui

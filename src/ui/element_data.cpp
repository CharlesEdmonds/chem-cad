#include "ui/element_data.hpp"

#include <algorithm>
#include <array>
#include <fstream>
#include <stdexcept>
#include <unordered_set>

#include <nlohmann/json.hpp>

#include "core/paths.hpp"

namespace chemcad::ui {
namespace {

struct LoadedElements {
  std::vector<ElementData> elements;
  std::string error;
};

std::vector<ElementData> fallbackElements() {
  return {
      {1, "H", "Hydrogen", 1.008, 1, 1, "nonmetal"},
      {5, "B", "Boron", 10.81, 13, 2, "metalloid"},
      {6, "C", "Carbon", 12.011, 14, 2, "nonmetal"},
      {7, "N", "Nitrogen", 14.007, 15, 2, "nonmetal"},
      {8, "O", "Oxygen", 15.999, 16, 2, "nonmetal"},
      {9, "F", "Fluorine", 18.998, 17, 2, "halogen"},
      {10, "Ne", "Neon", 20.180, 18, 2, "noble-gas"},
      {11, "Na", "Sodium", 22.990, 1, 3, "alkali-metal"},
      {12, "Mg", "Magnesium", 24.305, 2, 3, "alkaline-earth"},
      {13, "Al", "Aluminium", 26.982, 13, 3, "post-transition"},
      {14, "Si", "Silicon", 28.085, 14, 3, "metalloid"},
      {15, "P", "Phosphorus", 30.974, 15, 3, "nonmetal"},
      {16, "S", "Sulfur", 32.06, 16, 3, "nonmetal"},
      {17, "Cl", "Chlorine", 35.45, 17, 3, "halogen"},
      {19, "K", "Potassium", 39.098, 1, 4, "alkali-metal"},
      {20, "Ca", "Calcium", 40.078, 2, 4, "alkaline-earth"},
      {26, "Fe", "Iron", 55.845, 8, 4, "transition-metal"},
      {29, "Cu", "Copper", 63.546, 11, 4, "transition-metal"},
      {30, "Zn", "Zinc", 65.38, 12, 4, "transition-metal"},
      {35, "Br", "Bromine", 79.904, 17, 4, "halogen"},
      {53, "I", "Iodine", 126.904, 17, 5, "halogen"},
  };
}

bool validCategory(const std::string& category) {
  static constexpr std::array<const char*, 10> categories{
      "alkali-metal", "alkaline-earth", "transition-metal", "post-transition", "metalloid",
      "nonmetal", "halogen", "noble-gas", "lanthanide", "actinide"};
  return std::find(categories.begin(), categories.end(), category) != categories.end();
}

LoadedElements loadElements() {
  LoadedElements loaded;
  const auto path = core::dataDir() / "ptable.json";
  try {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("could not open " + path.string());

    const nlohmann::json root = nlohmann::json::parse(input);
    if (!root.is_array()) throw std::runtime_error("periodic table root is not an array");

    std::unordered_set<int> seen;
    loaded.elements.reserve(root.size());
    for (const auto& item : root) {
      ElementData e;
      const int z = item.at("z").get<int>();
      e.symbol = item.at("symbol").get<std::string>();
      e.name = item.at("name").get<std::string>();
      e.mass = item.at("mass").get<double>();
      e.group = item.at("group").get<int>();
      e.period = item.at("period").get<int>();
      e.category = item.at("category").get<std::string>();
      const bool detached = (z >= 57 && z <= 71) || (z >= 89 && z <= 103);
      if (z < 1 || z > 118 || e.symbol.empty() || e.name.empty() || e.mass <= 0 ||
          e.group < 0 || e.group > 18 || e.period < 1 || e.period > 7 ||
          (detached ? e.group != 0 : e.group == 0) || !validCategory(e.category) ||
          !seen.insert(z).second) {
        throw std::runtime_error("periodic table contains an invalid element entry");
      }
      e.z = static_cast<uint8_t>(z);
      loaded.elements.push_back(std::move(e));
    }
    if (loaded.elements.size() != 118) {
      throw std::runtime_error("periodic table does not contain all 118 elements");
    }
  } catch (const std::exception& e) {
    loaded.elements = fallbackElements();
    loaded.error = "Periodic table data unavailable; using common elements (";
    loaded.error += e.what();
    loaded.error += ")";
  }
  return loaded;
}

const LoadedElements& loadedElements() {
  static const LoadedElements loaded = loadElements();
  return loaded;
}

}  // namespace

const std::vector<ElementData>& elementTable() { return loadedElements().elements; }

const std::string& elementTableLoadError() { return loadedElements().error; }

const ElementData* findElement(uint8_t z) {
  for (const auto& e : elementTable()) {
    if (e.z == z) return &e;
  }
  return nullptr;
}

}  // namespace chemcad::ui

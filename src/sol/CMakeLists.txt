# Solubility modelling and the liquid-liquid separation simulation.
# Sits above chem (it needs descriptors for the solute) and below ui.
add_library(chemcad_sol STATIC
  solvent.cpp
  solubility.cpp
  anchors.cpp
  kirkwood_buff.cpp
  funnel.cpp)
target_include_directories(chemcad_sol PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(chemcad_sol
  PUBLIC chemcad_core
  PRIVATE chemcad_chem nlohmann_json::nlohmann_json)

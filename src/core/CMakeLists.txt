# Pure data model + geometry + paths. STL only: no RDKit, no ImGui, no network.
add_library(chemcad_core STATIC
  model.cpp
  sprout.cpp
  paths.cpp)
target_include_directories(chemcad_core PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_compile_definitions(chemcad_core PRIVATE
  CHEMCAD_DATA_DIR="${CHEMCAD_DATA_DIR}")

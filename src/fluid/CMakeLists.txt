# 3D multiphase fluid solver (PCISPH + density-contrast SPH).
# Pure computation: no GL, no ImGui, no RDKit, so it links into the test
# binaries and can run on a worker thread. It depends on sol only for the
# analytic vessel profile the boundary SDF revolves.
add_library(chemcad_fluid STATIC
  kernels.cpp
  grid.cpp
  vessel_sdf.cpp
  frame.cpp
  solver.cpp
  diagnostics.cpp
  simulation.cpp)
target_include_directories(chemcad_fluid PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(chemcad_fluid
  PUBLIC chemcad_core chemcad_sol Threads::Threads)

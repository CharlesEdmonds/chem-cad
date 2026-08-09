# OpenGL renderer plus an optional GL 4.3 compute fluid backend. Every GL call
# that is not ImGui's own lives here, behind gfx::loadGl, so the rest of the
# code and headless tests never see a GL symbol.
add_library(chemcad_gfx STATIC
  gl_api.cpp
  compute.cpp
  camera3d.cpp
  fluid_renderer.cpp
  fluid_gpu.cpp)
target_include_directories(chemcad_gfx PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(chemcad_gfx
  PUBLIC chemcad_core chemcad_fluid
  PRIVATE OpenGL::GL)

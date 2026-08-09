# OpenGL 3.3 renderer for the fluid stage. Every GL call in the project that
# is not ImGui's own lives here, behind gfx::loadGl, so the rest of the code
# (and the headless tests) never sees a GL symbol.
add_library(chemcad_gfx STATIC
  gl_api.cpp
  camera3d.cpp
  fluid_renderer.cpp)
target_include_directories(chemcad_gfx PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(chemcad_gfx
  PUBLIC chemcad_core chemcad_fluid
  PRIVATE OpenGL::GL)

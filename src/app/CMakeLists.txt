# Application services: async plumbing, project files, framebuffer capture.
# main.cpp is NOT here -- it belongs to the executable target.
file(GLOB CHEMCAD_APP_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
list(REMOVE_ITEM CHEMCAD_APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")

add_library(chemcad_app STATIC ${CHEMCAD_APP_SOURCES})
target_include_directories(chemcad_app PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(chemcad_app
  PUBLIC chemcad_core Threads::Threads
  PRIVATE chemcad_chem nlohmann_json::nlohmann_json OpenGL::GL)

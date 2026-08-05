# Globbed so parallel UI work can add translation units without contending on
# this file. CONFIGURE_DEPENDS re-globs when the directory changes.
file(GLOB CHEMCAD_UI_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_library(chemcad_ui STATIC ${CHEMCAD_UI_SOURCES})
target_include_directories(chemcad_ui PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(chemcad_ui
  PUBLIC chemcad_core chemcad_app imgui_lib
  PRIVATE chemcad_chem chemcad_naming chemcad_rxn nlohmann_json::nlohmann_json)

# One doctest binary per module keeps parallel slices from contending here.
# Globbed: dropping a tests/test_<name>.cpp file registers a new test target.
file(GLOB CHEMCAD_TEST_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/test_*.cpp")

foreach(test_src IN LISTS CHEMCAD_TEST_SOURCES)
  get_filename_component(test_name "${test_src}" NAME_WE)
  # ui/camera.cpp is pure geometry with no ImGui dependency, so tests compile it
  # directly rather than linking all of chemcad_ui (which would make every suite
  # unbuildable whenever a panel is mid-edit).
  add_executable(${test_name} "${test_src}" "${CMAKE_SOURCE_DIR}/src/ui/camera.cpp")
  # chemcad_chem keeps RDKit PRIVATE on purpose (only bridge.cpp may see it), so
  # suites that validate reaction SMARTS directly need their own RDKit link.
  target_link_libraries(${test_name} PRIVATE
    doctest::doctest
    chemcad_core
    chemcad_chem
    chemcad_naming
    chemcad_rxn
    chemcad_sol
    chemcad_fluid
    chemcad_gfx
    chemcad_app
    nlohmann_json::nlohmann_json
    ${CHEMCAD_RDKIT_LIBS})
  # The headless panel suites drive the real panels through a null ImGui
  # backend, so they need the UI layer; other suites stay lightweight. Any
  # test_*_panel.cpp is a headless panel suite by convention.
  if(test_name STREQUAL "test_ui_interaction" OR test_name MATCHES "_panel$")
    target_link_libraries(${test_name} PRIVATE chemcad_ui imgui_lib)
  endif()
  # The GPU fluid suite creates its own hidden OpenGL 4.3 context so it can
  # exercise the compute backend on real hardware; everything else stays
  # headless and window-system free.
  if(test_name STREQUAL "test_fluid_gpu")
    target_link_libraries(${test_name} PRIVATE glfw)
  endif()
  target_include_directories(${test_name} PRIVATE ${CMAKE_SOURCE_DIR}/src)
  chemcad_copy_runtime_dlls(${test_name})
  add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()

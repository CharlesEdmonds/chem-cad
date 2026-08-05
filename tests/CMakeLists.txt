# One doctest binary per module keeps parallel slices from contending here.
# Globbed: dropping a tests/test_<name>.cpp file registers a new test target.
file(GLOB CHEMCAD_TEST_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/test_*.cpp")

foreach(test_src IN LISTS CHEMCAD_TEST_SOURCES)
  get_filename_component(test_name "${test_src}" NAME_WE)
  add_executable(${test_name} "${test_src}")
  target_link_libraries(${test_name} PRIVATE
    doctest::doctest
    chemcad_core
    chemcad_chem
    chemcad_naming
    chemcad_rxn
    chemcad_app
    nlohmann_json::nlohmann_json)
  target_include_directories(${test_name} PRIVATE ${CMAKE_SOURCE_DIR}/src)
  add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()

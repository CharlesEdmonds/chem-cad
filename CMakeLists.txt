cmake_minimum_required(VERSION 3.24)
project(chemcad VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

if(MSVC)
  # /utf-8: the sources carry UTF-8 punctuation in UI strings.
  # /bigobj: the RDKit-facing translation units blow past the 64k section limit.
  add_compile_options(/utf-8 /bigobj /EHsc)
  # WIN32 (which MSVC does not predefine) is half of the condition that turns
  # RDKit's export macros into __declspec(dllimport); chemcad_rdkit supplies the
  # other half. Without both, MSVC resolves functions through the import library
  # but leaves exported *data* symbols such as rdErrorLog unresolved.
  # BOOST_ALL_NO_LIB turns off Boost's #pragma auto-linking, which guesses
  # library names that no packaged Boost actually uses.
  add_compile_definitions(NOMINMAX WIN32_LEAN_AND_MEAN _USE_MATH_DEFINES
                          _CRT_SECURE_NO_WARNINGS WIN32 BOOST_ALL_NO_LIB)
endif()

option(CHEMCAD_BUILD_TESTS "Build the doctest suites" ON)

include(cmake/Dependencies.cmake)

# Data files (reaction KB, periodic table) are loaded at runtime.
# Dev builds read them straight from the source tree; the env var
# CHEMCAD_DATA_DIR overrides at runtime for installed/relocated copies.
set(CHEMCAD_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")
set(CHEMCAD_ASSETS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/assets")

add_subdirectory(src/core)
add_subdirectory(src/chem)
add_subdirectory(src/naming)
add_subdirectory(src/rxn)
add_subdirectory(src/sol)
add_subdirectory(src/app)
add_subdirectory(src/ui)

add_executable(chemcad src/app/main.cpp)
target_link_libraries(chemcad PRIVATE chemcad_ui chemcad_app)
chemcad_copy_runtime_dlls(chemcad)

if(CHEMCAD_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

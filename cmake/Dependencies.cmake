# Third-party dependencies.
#   - Fetched at configure time: imgui, glfw, nlohmann_json, doctest
#   - Provided by scripts/setup_deps.sh: RDKit, OPSIN jar
#   - System packages: libcurl, OpenGL, Threads

include(FetchContent)
set(FETCHCONTENT_QUIET OFF)

find_package(OpenGL REQUIRED)
find_package(CURL REQUIRED)
find_package(Threads REQUIRED)

# ------------------------------------------------------------------ RDKit
if(WIN32)
  set(CHEMCAD_DEPS_PREFIX "$ENV{LOCALAPPDATA}/chemcad-deps" CACHE PATH
      "Prefix where scripts/setup_deps.ps1 installed RDKit and OPSIN")
else()
  set(CHEMCAD_DEPS_PREFIX "$ENV{HOME}/.local/share/chemcad-deps" CACHE PATH
      "Prefix where scripts/setup_deps.sh installed RDKit and OPSIN")
endif()

# Every packaged RDKit (Debian's librdkit-dev, conda-forge's librdkit-dev) names
# Cairo::Cairo in MolDraw2D's link interface without exporting a find_package
# for it, so rdkit-targets.cmake hard-errors unless the target already exists.
if(NOT TARGET Cairo::Cairo)
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    pkg_check_modules(CHEMCAD_CAIRO QUIET IMPORTED_TARGET cairo)
  endif()
  if(CHEMCAD_CAIRO_FOUND)
    add_library(Cairo::Cairo ALIAS PkgConfig::CHEMCAD_CAIRO)
  else()
    find_path(CHEMCAD_CAIRO_INCLUDE_DIR cairo.h PATH_SUFFIXES cairo)
    find_library(CHEMCAD_CAIRO_LIBRARY NAMES cairo)
    if(CHEMCAD_CAIRO_INCLUDE_DIR AND CHEMCAD_CAIRO_LIBRARY)
      add_library(Cairo::Cairo UNKNOWN IMPORTED)
      set_target_properties(Cairo::Cairo PROPERTIES
        IMPORTED_LOCATION "${CHEMCAD_CAIRO_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${CHEMCAD_CAIRO_INCLUDE_DIR}")
    endif()
  endif()
endif()

find_package(RDKit CONFIG QUIET
  HINTS "${CHEMCAD_DEPS_PREFIX}/rdkit" "${CHEMCAD_DEPS_PREFIX}"
  PATH_SUFFIXES lib64/cmake/rdkit lib/cmake/rdkit share/RDKit/cmake
                Library/lib/cmake/rdkit)
if(NOT RDKit_FOUND)
  if(WIN32)
    message(FATAL_ERROR
      "RDKit not found. Run scripts/setup_deps.ps1 first, then configure with "
      "-DCMAKE_PREFIX_PATH=<prefix>/rdkit/Library "
      "(or pass -DCHEMCAD_DEPS_PREFIX=<prefix>).")
  else()
    message(FATAL_ERROR
      "RDKit not found. Run ./scripts/setup_deps.sh first "
      "(or pass -DCHEMCAD_DEPS_PREFIX=/path/to/prefix).")
  endif()
endif()
message(STATUS "Found RDKit: ${RDKit_DIR}")

# RDKit's shared libraries pull in siblings (libexpatpp, libcoordgen, ...) via
# plain DT_NEEDED entries, so the loader needs the RDKit lib dir on the RPATH.
# Windows has no RPATH; chemcad_copy_runtime_dlls() below stages the DLLs next
# to each executable instead.
get_filename_component(CHEMCAD_RDKIT_LIBDIR "${RDKit_DIR}/../.." ABSOLUTE)
if(UNIX)
  list(APPEND CMAKE_BUILD_RPATH "${CHEMCAD_RDKIT_LIBDIR}")
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
endif()

# Packaged RDKit keeps its import libraries in lib/ and its DLLs in the sibling
# bin/, and pulls in Boost, cairo, freetype and friends through plain PE
# imports that no CMake target models.
get_filename_component(CHEMCAD_RUNTIME_BINDIR "${CHEMCAD_RDKIT_LIBDIR}/../bin" ABSOLUTE)

# Stages the transitive DLL closure of an executable into its output directory.
# A no-op on platforms with a real runtime loader search path.
function(chemcad_copy_runtime_dlls target)
  if(NOT WIN32)
    return()
  endif()
  add_custom_command(TARGET ${target} POST_BUILD
    COMMAND ${CMAKE_COMMAND}
            -D "EXE=$<TARGET_FILE:${target}>"
            -D "DEST=$<TARGET_FILE_DIR:${target}>"
            -D "SEARCH=${CHEMCAD_RUNTIME_BINDIR}"
            -P "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/StageRuntimeDlls.cmake"
    VERBATIM)
endfunction()

# RDKit's exported targets are namespaced RDKit::<Component>. They are wrapped
# in one INTERFACE target so the consumer-side compile definitions travel with
# the link line.
add_library(chemcad_rdkit INTERFACE)
target_link_libraries(chemcad_rdkit INTERFACE
  RDKit::GraphMol
  RDKit::SmilesParse
  RDKit::FileParsers
  RDKit::Depictor
  RDKit::ChemReactions
  RDKit::Descriptors
  RDKit::MolDraw2D
  RDKit::SubstructMatch
  RDKit::DistGeomHelpers
  RDKit::ForceFieldHelpers
  RDKit::ForceField
  RDKit::RDGeneral
  RDKit::DataStructs)
if(WIN32)
  # RDKit ships DLLs here. Without RDKIT_DYN_LINK its headers declare exported
  # data (rdErrorLog, RDDepict::preferCoordGen) with no __declspec(dllimport),
  # and MSVC cannot auto-import data symbols the way it can functions.
  target_compile_definitions(chemcad_rdkit INTERFACE RDKIT_DYN_LINK)
endif()
set(CHEMCAD_RDKIT_LIBS chemcad_rdkit)

# ------------------------------------------------------------------ OPSIN
set(CHEMCAD_OPSIN_JAR "${CHEMCAD_DEPS_PREFIX}/share/opsin/opsin.jar"
    CACHE FILEPATH "Path to the OPSIN cli jar (name -> structure)")
if(NOT EXISTS "${CHEMCAD_OPSIN_JAR}")
  message(WARNING
    "OPSIN jar not found at ${CHEMCAD_OPSIN_JAR}. "
    "name->structure will fall back to PubChem only. Run ./scripts/setup_deps.sh.")
endif()

# ------------------------------------------------------------------ GLFW
FetchContent_Declare(glfw
  GIT_REPOSITORY https://github.com/glfw/glfw.git
  GIT_TAG 3.4
  GIT_SHALLOW TRUE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL        OFF CACHE BOOL "" FORCE)
if(UNIX AND NOT APPLE)
  set(GLFW_BUILD_WAYLAND ON CACHE BOOL "" FORCE)
  set(GLFW_BUILD_X11     ON CACHE BOOL "" FORCE)
endif()
FetchContent_MakeAvailable(glfw)

# ------------------------------------------------------------------ JSON
FetchContent_Declare(nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG v3.11.3
  GIT_SHALLOW TRUE)
set(JSON_BuildTests OFF CACHE INTERNAL "")
FetchContent_MakeAvailable(nlohmann_json)

# ------------------------------------------------------------------ doctest
# doctest 2.4.11 still declares cmake_minimum_required(VERSION 3.0), which
# CMake 4 rejects outright. Scope the compatibility floor to this subproject.
FetchContent_Declare(doctest
  GIT_REPOSITORY https://github.com/doctest/doctest.git
  GIT_TAG v2.4.11
  GIT_SHALLOW TRUE)
set(DOCTEST_WITH_TESTS OFF CACHE INTERNAL "")
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
FetchContent_MakeAvailable(doctest)
unset(CMAKE_POLICY_VERSION_MINIMUM)

# ------------------------------------------------------------------ Dear ImGui
FetchContent_Declare(imgui
  GIT_REPOSITORY https://github.com/ocornut/imgui.git
  GIT_TAG v1.92.9b-docking
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(imgui)

add_library(imgui_lib STATIC
  ${imgui_SOURCE_DIR}/imgui.cpp
  ${imgui_SOURCE_DIR}/imgui_draw.cpp
  ${imgui_SOURCE_DIR}/imgui_tables.cpp
  ${imgui_SOURCE_DIR}/imgui_widgets.cpp
  ${imgui_SOURCE_DIR}/imgui_demo.cpp
  ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
  ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp)
target_include_directories(imgui_lib PUBLIC
  ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends)
target_link_libraries(imgui_lib PUBLIC glfw OpenGL::GL)
target_compile_definitions(imgui_lib PUBLIC IMGUI_DEFINE_MATH_OPERATORS)

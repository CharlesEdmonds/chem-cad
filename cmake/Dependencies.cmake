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
set(CHEMCAD_DEPS_PREFIX "$ENV{HOME}/.local/share/chemcad-deps" CACHE PATH
    "Prefix where scripts/setup_deps.sh installed RDKit and OPSIN")

# Debian/Ubuntu ship librdkit-dev, whose rdkit-targets.cmake names Cairo::Cairo
# in MolDraw2D's link interface without exporting a find_package for it. Define
# the target from pkg-config first so the config file resolves.
if(NOT TARGET Cairo::Cairo)
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    pkg_check_modules(CHEMCAD_CAIRO QUIET IMPORTED_TARGET cairo)
    if(CHEMCAD_CAIRO_FOUND)
      add_library(Cairo::Cairo ALIAS PkgConfig::CHEMCAD_CAIRO)
    endif()
  endif()
endif()

find_package(RDKit CONFIG QUIET
  HINTS "${CHEMCAD_DEPS_PREFIX}/rdkit"
  PATH_SUFFIXES lib64/cmake/rdkit lib/cmake/rdkit share/RDKit/cmake)
if(NOT RDKit_FOUND)
  message(FATAL_ERROR
    "RDKit not found. Run ./scripts/setup_deps.sh first "
    "(or pass -DCHEMCAD_DEPS_PREFIX=/path/to/prefix).")
endif()
message(STATUS "Found RDKit: ${RDKit_DIR}")

# RDKit's shared libraries pull in siblings (libexpatpp, libcoordgen, ...) via
# plain DT_NEEDED entries, so the loader needs the RDKit lib dir on the RPATH.
get_filename_component(CHEMCAD_RDKIT_LIBDIR "${RDKit_DIR}/../.." ABSOLUTE)
list(APPEND CMAKE_BUILD_RPATH "${CHEMCAD_RDKIT_LIBDIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

# RDKit's exported targets are namespaced RDKit::<Component>.
set(CHEMCAD_RDKIT_LIBS
  RDKit::GraphMol
  RDKit::SmilesParse
  RDKit::FileParsers
  RDKit::Depictor
  RDKit::ChemReactions
  RDKit::Descriptors
  RDKit::MolDraw2D
  RDKit::SubstructMatch
  RDKit::RDGeneral
  RDKit::DataStructs)

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
set(GLFW_BUILD_WAYLAND  ON  CACHE BOOL "" FORCE)
set(GLFW_BUILD_X11      ON  CACHE BOOL "" FORCE)
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

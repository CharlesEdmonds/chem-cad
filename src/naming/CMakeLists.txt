add_library(chemcad_naming STATIC
  naming.cpp)
target_include_directories(chemcad_naming PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(chemcad_naming
  PUBLIC chemcad_core
  PRIVATE CURL::libcurl nlohmann_json::nlohmann_json)
target_compile_definitions(chemcad_naming PRIVATE
  CHEMCAD_OPSIN_JAR="${CHEMCAD_OPSIN_JAR}")

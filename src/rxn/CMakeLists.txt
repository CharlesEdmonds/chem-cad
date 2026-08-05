add_library(chemcad_rxn STATIC
  engine.cpp
  kb.cpp
  llm.cpp)
target_include_directories(chemcad_rxn PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(chemcad_rxn
  PUBLIC chemcad_core
  PRIVATE chemcad_chem CURL::libcurl nlohmann_json::nlohmann_json
          ${CHEMCAD_RDKIT_LIBS})

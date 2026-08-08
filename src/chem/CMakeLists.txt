# RDKit boundary layer. The only target that sees RDKit headers.
add_library(chemcad_chem STATIC
  bridge.cpp
  embed3d.cpp)
target_include_directories(chemcad_chem PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(chemcad_chem
  PUBLIC chemcad_core
  PRIVATE ${CHEMCAD_RDKIT_LIBS})

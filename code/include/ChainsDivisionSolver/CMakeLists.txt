set(CDS_SOURCES ChainsDivisionSolver.hpp)

add_library(CDS OBJECT ${CDS_SOURCES})
set_target_properties(CDS PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(CDS PRIVATE DAWG blossom5)

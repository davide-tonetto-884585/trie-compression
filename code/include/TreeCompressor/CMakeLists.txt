set(TreeCompressorSources TreeCompressor.hpp)

add_library(TreeCompressor OBJECT ${TreeCompressorSources})
set_target_properties(TreeCompressor PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(TreeCompressor PUBLIC DAWG Triplet)

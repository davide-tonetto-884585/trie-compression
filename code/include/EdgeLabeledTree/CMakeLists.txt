set(EDGE_LABELED_TREE_SOURCES EdgeLabeledTree.hpp)

add_library(EdgeLabeledTree OBJECT ${EDGE_LABELED_TREE_SOURCES})
set_target_properties(EdgeLabeledTree PROPERTIES LINKER_LANGUAGE CXX)

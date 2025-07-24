set(LABELED_TREE_SOURCES LabeledTree.hpp)

add_library(LabeledTree OBJECT ${LABELED_TREE_SOURCES})
set_target_properties(LabeledTree PROPERTIES LINKER_LANGUAGE CXX)

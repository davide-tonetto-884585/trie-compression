add_subdirectory(Triplet)
target_include_directories(Triplet INTERFACE Triplet)

add_subdirectory(XBWT)
target_include_directories(XBWT PUBLIC XBWT)

add_subdirectory(LabeledTree)
target_include_directories(LabeledTree PUBLIC LabeledTree)

add_subdirectory(DAWG)
target_include_directories(DAWG PUBLIC DAWG)

add_subdirectory(ChainsDivisionSolver)
target_include_directories(CDS PUBLIC ChainsDivisionSolver)

add_subdirectory(TreeCompressor)
target_include_directories(TreeCompressor PUBLIC TreeCompressor)

add_subdirectory(EdgeLabeledTree)
target_include_directories(EdgeLabeledTree PUBLIC EdgeLabeledTree)
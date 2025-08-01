add_subdirectory(XBWT)
target_include_directories(XBWT PUBLIC XBWT)

add_subdirectory(LabeledTree)
target_include_directories(LabeledTree PUBLIC LabeledTree)

add_subdirectory(DAWG)
target_include_directories(DAWG PUBLIC DAWG)

add_subdirectory(ChainsDivisionSolver)
target_include_directories(CDS PUBLIC ChainsDivisionSolver)

set(XBWT_SOURCES XBWT.hpp)

add_library(XBWT OBJECT ${XBWT_SOURCES})
set_target_properties(XBWT PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(XBWT PUBLIC sdsl Triplet)

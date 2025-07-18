set(XBWT_SOURCES XBWT.hpp)

add_library(XBWT OBJECT ${XBWT_SOURCES})
target_link_libraries(XBWT PUBLIC sdsl)

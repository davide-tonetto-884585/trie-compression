set(TripletSources Triplet.hpp)

add_library(Triplet INTERFACE)
target_include_directories(Triplet INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
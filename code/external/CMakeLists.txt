include(FetchContent)

##########################################################################

## sdsl-lite
FetchContent_Declare(
  sdsl
  GIT_REPOSITORY https://github.com/simongog/sdsl-lite.git
  GIT_TAG v2.1.1
)

FetchContent_GetProperties(sdsl)
if(NOT sdsl_POPULATED)
  FetchContent_Populate(sdsl)
  
  # Apply the fix only when using Clang compiler
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    execute_process(
      COMMAND sed -i.bak "s/tree\\.m_select1/tree.m_bv_select1/g; s/tree\\.m_select0/tree.m_bv_select0/g" 
              "${sdsl_SOURCE_DIR}/include/sdsl/louds_tree.hpp"
      RESULT_VARIABLE sed_result
    )
    
    if(sed_result EQUAL 0)
      message(STATUS "Successfully applied SDSL louds_tree.hpp fix for Clang")
    else()
      message(WARNING "Failed to apply SDSL louds_tree.hpp fix for Clang")
    endif()
  else()
    message(STATUS "Skipping SDSL louds_tree.hpp fix (not using Clang compiler)")
  endif()

  set(GENERATE_DOC OFF CACHE BOOL "Do not generate doxygen for sdsl-lite")
  
  add_subdirectory(${sdsl_SOURCE_DIR} ${sdsl_BINARY_DIR} EXCLUDE_FROM_ALL)
  
  # Explicitly add divsufsort include directories to sdsl target
  target_include_directories(sdsl PUBLIC 
    "${sdsl_BINARY_DIR}/external/libdivsufsort/include"
    "${sdsl_BINARY_DIR}/include"
  )
endif()

##########################################################################



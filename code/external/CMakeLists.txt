include(FetchContent)

##########################################################################

## sdsl-lite
FetchContent_Declare(
  sdsl
  GIT_REPOSITORY https://github.com/simongog/sdsl-lite.git
  GIT_TAG v2.1.1
  GIT_SHALLOW TRUE
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
# Blossom V - Minimum Cost Perfect Matching

FetchContent_Declare(
  blossom5
  URL https://pub.ista.ac.at/~vnk/software/blossom5-v2.05.src.tar.gz
)

FetchContent_GetProperties(blossom5)
if(NOT blossom5_POPULATED)
  FetchContent_Populate(blossom5)
  
  # Fix pointer type issue for 64-bit systems
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)  # 64-bit system
    file(READ ${blossom5_SOURCE_DIR}/PMimplementation.h PMIMPL_CONTENT)
    string(REPLACE "typedef unsigned long POINTER_TYPE;" 
                   "typedef uintptr_t POINTER_TYPE;" 
                   PMIMPL_CONTENT "${PMIMPL_CONTENT}")
    string(REPLACE "#include <stdio.h>" 
                   "#include <stdio.h>\n#include <stdint.h>" 
                   PMIMPL_CONTENT "${PMIMPL_CONTENT}")
    file(WRITE ${blossom5_SOURCE_DIR}/PMimplementation.h "${PMIMPL_CONTENT}")
    message(STATUS "Applied 64-bit pointer fix to Blossom V")
  endif()
  
  # Create CMakeLists.txt for Blossom V
  file(WRITE ${blossom5_SOURCE_DIR}/CMakeLists.txt "
    cmake_minimum_required(VERSION 3.13)
    project(blossom5)
    
    set(BLOSSOM5_SOURCES
        PMduals.cpp
        PMexpand.cpp
        PMinit.cpp
        PMinterface.cpp
        PMmain.cpp
        PMrepair.cpp
        PMshrink.cpp
        MinCost/MinCost.cpp
    )
    
    add_library(blossom5 STATIC \${BLOSSOM5_SOURCES})
    target_include_directories(blossom5 PUBLIC .)
    set_target_properties(blossom5 PROPERTIES CXX_STANDARD 11)
    
    # Disable problematic warnings
    if(MSVC)
        target_compile_options(blossom5 PRIVATE /wd4267 /wd4244)
    else()
        target_compile_options(blossom5 PRIVATE -Wno-pointer-to-int-cast -Wno-int-to-pointer-cast)
    endif()
  ")
  
  add_subdirectory(${blossom5_SOURCE_DIR} ${blossom5_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

##########################################################################

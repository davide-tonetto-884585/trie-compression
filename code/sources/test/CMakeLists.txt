# Build the test executable
add_executable(dawg_test test.cpp)
target_link_libraries(dawg_test PUBLIC DAWG CDS XBWT)

# Add debug symbols for better Valgrind output
target_compile_options(dawg_test PRIVATE -g -O0)

# Basic test
add_test(NAME basic_test COMMAND test)

# Valgrind tests (only if Valgrind is found)
if(VALGRIND_FOUND)
    # Memory leak detection test
    add_test(
        NAME memcheck_test
        COMMAND ${VALGRIND_EXECUTABLE}
            --tool=memcheck
            --leak-check=full
            --show-leak-kinds=all
            --track-origins=yes
            --verbose
            --error-exitcode=1
            --suppressions=${CMAKE_CURRENT_SOURCE_DIR}/valgrind.supp
            $<TARGET_FILE:test>
    )
    
    # Cache profiling test
    add_test(
        NAME cachegrind_test
        COMMAND ${VALGRIND_EXECUTABLE}
            --tool=cachegrind
            --cache-sim=yes
            --branch-sim=yes
            --cachegrind-out-file=${CMAKE_BINARY_DIR}/cachegrind.out.%p
            $<TARGET_FILE:test>
    )
    
    # Heap profiling test
    add_test(
        NAME massif_test
        COMMAND ${VALGRIND_EXECUTABLE}
            --tool=massif
            --heap=yes
            --stacks=yes
            --massif-out-file=${CMAKE_BINARY_DIR}/massif.out.%p
            $<TARGET_FILE:test>
    )
    
    # Call graph profiling test
    add_test(
        NAME callgrind_test
        COMMAND ${VALGRIND_EXECUTABLE}
            --tool=callgrind
            --callgrind-out-file=${CMAKE_BINARY_DIR}/callgrind.out.%p
            --collect-jumps=yes
            --collect-systime=yes
            $<TARGET_FILE:test>
    )
    
    # Set test properties
    set_tests_properties(memcheck_test PROPERTIES
        TIMEOUT 300
        LABELS "memory;valgrind"
    )
    
    set_tests_properties(cachegrind_test PROPERTIES
        TIMEOUT 300
        LABELS "performance;cache;valgrind"
    )
    
    set_tests_properties(massif_test PROPERTIES
        TIMEOUT 300
        LABELS "memory;heap;valgrind"
    )
    
    set_tests_properties(callgrind_test PROPERTIES
        TIMEOUT 300
        LABELS "performance;profiling;valgrind"
    )
    
    # Custom targets for easier execution
    add_custom_target(memcheck
        COMMAND ${VALGRIND_EXECUTABLE}
            --tool=memcheck
            --leak-check=full
            --show-leak-kinds=all
            --track-origins=yes
            --verbose
            --suppressions=${CMAKE_CURRENT_SOURCE_DIR}/valgrind.supp
            $<TARGET_FILE:test>
        DEPENDS test
        COMMENT "Running memory check with Valgrind"
    )
    
    add_custom_target(cachegrind
        COMMAND ${VALGRIND_EXECUTABLE}
            --tool=cachegrind
            --cache-sim=yes
            --branch-sim=yes
            --cachegrind-out-file=${CMAKE_BINARY_DIR}/cachegrind.out
            $<TARGET_FILE:test>
        DEPENDS test
        COMMENT "Running cache profiling with Cachegrind"
    )
    
    add_custom_target(massif
        COMMAND ${VALGRIND_EXECUTABLE}
            --tool=massif
            --heap=yes
            --stacks=yes
            --massif-out-file=${CMAKE_BINARY_DIR}/massif.out
            $<TARGET_FILE:test>
        DEPENDS test
        COMMENT "Running heap profiling with Massif"
    )
    
    add_custom_target(callgrind
        COMMAND ${VALGRIND_EXECUTABLE}
            --tool=callgrind
            --callgrind-out-file=${CMAKE_BINARY_DIR}/callgrind.out
            --collect-jumps=yes
            --collect-systime=yes
            $<TARGET_FILE:test>
        DEPENDS test
        COMMENT "Running call graph profiling with Callgrind"
    )
    
else()
    message(STATUS "Valgrind not found - memory tests disabled")
endif()

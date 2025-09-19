# Build the test executable
add_executable(dawg_experiments time_experiments.cpp)
target_link_libraries(dawg_experiments PUBLIC DAWG CDS XBWT TreeCompressor)

# Add debug symbols for better Valgrind output
target_compile_options(dawg_experiments PRIVATE -g -O0)

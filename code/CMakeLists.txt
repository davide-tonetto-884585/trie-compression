cmake_minimum_required(VERSION 3.10)
project(DAWGMinimizer LANGUAGES CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Debug build type
set(CMAKE_BUILD_TYPE Debug)

# Debug flags for macOS
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g3 -Wall -Wextra -fsanitize=address -fno-omit-frame-pointer")

# Set optimization flags for Release mode
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Wall")

add_executable(dawg_minimizer main.cpp DAWG.hpp)

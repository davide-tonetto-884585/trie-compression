cmake_minimum_required(VERSION 3.10)
project(DAWGMinimizer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(dawg_minimizer main.cpp DAWG.hpp)

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Stampa un messaggio dopo la configurazione
message(STATUS "CMake configuration complete. Use 'make' or your chosen build system to build ${PROJECT_NAME}.")

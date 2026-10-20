# Specify the minimum version for CMake.
cmake_minimum_required(VERSION 3.8.0)
# Project's name.
project(jeanbaptiste.test VERSION 0.1.0)

enable_testing()

# Set the output folder where your program will be created.
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)
# Set the output folder for executable files.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# Use C++23.
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the main include directory.
include_directories(include)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

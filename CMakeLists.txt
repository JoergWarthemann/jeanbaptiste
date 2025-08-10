cmake_minimum_required(VERSION 3.15)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

# Only set the cxx_standard if it is not set by someone else
if (NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 23)
endif()

# strongly encouraged to enable this globally to avoid conflicts between
# -Wpedantic being enabled and -std=c++20 and -std=gnu++20 for example
# when compiling with PCH enabled
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Disable cmake searching the system for local Boost installation.  
set(BOOST_ROOT "" CACHE PATH "Boost library path")
set(Boost_NO_SYSTEM_PATHS on CACHE BOOL "Do not search system for Boost")

# Set the project name and language
project(
    jeanbaptiste
    VERSION 0.0.1
    DESCRIPTION "Implements FFT algorithms"
    LANGUAGES CXX C)

set(ROOT_SOURCE_DIR "${PROJECT_SOURCE_DIR}")

add_custom_target(build-jean-baptiste-tests)

enable_testing()

add_subdirectory(jean-baptiste-app)
add_subdirectory(jean-baptiste-lib)

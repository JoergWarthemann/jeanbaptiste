# Specify the minimum version for CMake.
cmake_minimum_required(VERSION 3.8.0)
# Project's name.
project(jeanbaptiste.test VERSION 0.1.0)

enable_testing()

# Set the output folder where your program will be created.
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)
# Set the output folder for executable files.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# Use C++17.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the main include directory.
include_directories(include)

# Use boost.
set(Boost_INCLUDE_DIRS /usr/local/include/)
set(Boost_LIBRARY_DIR /usr/local/lib/)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_DEBUG_LIBS OFF)
set(Boost_USE_RELEASE_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

#find_package(Boost 1.68.0 COMPONENTS)
find_package(Boost 1.68.0 REQUIRED COMPONENTS unit_test_framework filesystem)
if(Boost_FOUND)   
    # Set boost include directory.
    include_directories(${Boost_INCLUDE_DIRS})

    set(SOURCE test.cpp)
    add_executable(jeanbaptiste.test src/test.cpp)
    target_include_directories(jeanbaptiste.test PRIVATE ${Boost_INCLUDE_DIRS})
    target_link_libraries(jeanbaptiste.test
        ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY})
    add_test(NAME jbt COMMAND jeanbaptiste.test)
endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

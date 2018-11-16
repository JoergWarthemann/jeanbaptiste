# Specify the minimum version for CMake.
cmake_minimum_required(VERSION 3.0.0)

option(BUILD_JEANBAPTISTE "Enable building jeanbaptiste" ON)
option(BUILD_JEANBAPTISTE_TEST "Enable building jeanbaptiste.test" ON)

if(BUILD_JEANBAPTISTE)
    add_subdirectory(JeanBaptiste)
endif()
if(BUILD_JEANBAPTISTE_TEST)
    add_subdirectory(JeanBaptiste.Test)
endif()
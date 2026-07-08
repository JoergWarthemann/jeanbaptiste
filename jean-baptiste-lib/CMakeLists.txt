cmake_minimum_required(VERSION 3.15)

set(Boost_INCLUDE_DIR /usr/lib/boost_1_88_0/include)
set(Boost_LIBRARY_DIR /usr/lib/boost_1_88_0/lib)

set(LIBRARY_NAME jean-baptiste-lib)
set(JEAN_BAPTISTE_LIB_SOURCES
    include/Algorithm.hpp
    include/AlgorithmFactory.hpp
    include/IExecutableAlgorithm.hpp
    include/Options.hpp
    include/core/Radix2.hpp
    include/core/Radix4.hpp
    include/core/RadixSplit24.hpp
    include/core/SplitRadix.hpp
    include/normalization/DivisionByLengthNormalization.hpp
    include/normalization/NoNormalization.hpp
    include/normalization/SquareRootNormalization.hpp
    include/tools/Abs.hpp
    include/tools/BitReversalIndexSwapping.hpp
    include/tools/Ceil.hpp
    include/tools/Floor.hpp
    include/tools/HeronSquareRoot.hpp
    include/tools/RealComplexConversion.hpp
    include/tools/Round.hpp
    include/tools/SineCosine.hpp
    include/tools/SubTask.hpp
    include/windows/BartlettWindow.hpp
    include/windows/BlackmanHarrisWindow.hpp
    include/windows/BlackmanWindow.hpp
    include/windows/CosineWindow.hpp
    include/windows/ExecuteWindowOnComplexData.hpp
    include/windows/FlatTopWindow.hpp
    include/windows/HammingWindow.hpp
    include/windows/NoWindow.hpp
    include/windows/VonHannWindow.hpp
    include/windows/WelchWindow.hpp
    )

find_package(Boost COMPONENTS system REQUIRED)
if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost library was not found")
endif()

add_library(${LIBRARY_NAME} INTERFACE
    ${JEAN_BAPTISTE_LIB_SOURCES}
    )

target_include_directories(${LIBRARY_NAME}
    INTERFACE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(${LIBRARY_NAME}
    INTERFACE
        Boost::boost Boost::system)

###############################################

set(TEST_WORKING_DIR "${CMAKE_CURRENT_BINARY_DIR}/utest")

find_package(GTest REQUIRED)
find_library(GOOGLEMOCK_LIB libgmock.a)

add_executable(${LIBRARY_NAME}-googletest
    ${JEAN_BAPTISTE_LIB_SOURCES}
    testsrc/src/BitReversalTest.cpp
    testsrc/src/Radix2Test.cpp
    testsrc/src/Radix4Test.cpp
    testsrc/src/SinCosTest.cpp

    testsrc/include/AlgorithmFixture.hpp
    testsrc/include/AlgorithmResultAnalysis.hpp
    testsrc/include/StringConversion.hpp
    testsrc/include/TestCaseLoader.hpp
    testsrc/include/WindowingAnalysis.hpp
    testsrc/src/main.cpp
    testsrc/src/TestCaseLoader.cpp
    )

target_include_directories(${LIBRARY_NAME}-googletest
    PRIVATE
        testsrc/include
    )

target_link_libraries(${LIBRARY_NAME}-googletest
    PRIVATE
        ${LIBRARY_NAME}
        ${GOOGLEMOCK_LIB}
        ${GTEST_LIBRARIES}
    )

target_compile_definitions(${LIBRARY_NAME}-googletest
    PRIVATE
        TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/testsrc/test-data"
    )

add_custom_command(TARGET ${LIBRARY_NAME}-googletest POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
         "${CMAKE_CURRENT_SOURCE_DIR}/testsrc/test-data" "${TEST_WORKING_DIR}/test-data"
    )

add_dependencies(build-jean-baptiste-tests
    ${LIBRARY_NAME}-googletest
    )
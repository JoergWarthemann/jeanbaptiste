cmake_minimum_required(VERSION 3.22)
project(inspect VERSION 0.1.0 LANGUAGES C CXX)

set(COMPONENT_NAME inspect)

set(JUCE_SOURCE_DIR "" CACHE PATH "JUCE source directory")
if (NOT JUCE_SOURCE_DIR AND DEFINED ENV{JUCE_SOURCE_DIR})
  set(JUCE_SOURCE_DIR "$ENV{JUCE_SOURCE_DIR}" CACHE PATH "JUCE source directory" FORCE)
endif()
if (NOT JUCE_SOURCE_DIR AND DEFINED ENV{JUCE_PATH})
  set(JUCE_SOURCE_DIR "$ENV{JUCE_PATH}" CACHE PATH "JUCE source directory" FORCE)
endif()
if (NOT JUCE_SOURCE_DIR)
  message(FATAL_ERROR "JUCE_SOURCE_DIR is not set. Configure with --preset linux_x86_Debug_local, set the JUCE_PATH environment variable, or pass -DJUCE_SOURCE_DIR=/path/to/JUCE.")
endif()
if (NOT EXISTS "${JUCE_SOURCE_DIR}/CMakeLists.txt")
  message(FATAL_ERROR "JUCE_SOURCE_DIR does not point to a JUCE source tree: ${JUCE_SOURCE_DIR}")
endif()
add_subdirectory("${JUCE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/JUCE")

juce_add_gui_app(inspect
    # VERSION ...                       # Set this if the app version is different to the project version
    # ICON_BIG ...                      # ICON_* arguments specify a path to an image file to use as an icon
    # ICON_SMALL ...
    # DOCUMENT_EXTENSIONS ...           # Specify file extensions that should be associated with this app
    # COMPANY_NAME ...                  # Specify the name of the app's author
    PRODUCT_NAME ${COMPONENT_NAME}      # The name of the final executable, which can differ from the target name
    NEEDS_WEB_BROWSER FALSE             # Set to TRUE if you need to use webviews
    NEEDS_CURL FALSE)                   # Set to TRUE if you need to use curl

juce_generate_juce_header(${COMPONENT_NAME})

target_sources(${COMPONENT_NAME}
    PRIVATE
        include/MainComponent.hpp
        src/Main.cpp
        src/MainComponent.cpp)

target_include_directories(${COMPONENT_NAME}
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_compile_definitions(${COMPONENT_NAME}
    PRIVATE
        JUCE_WEB_BROWSER=0
        JUCE_USE_CURL=0
        JUCE_APPLICATION_NAME_STRING="$<TARGET_PROPERTY:inspect,JUCE_PRODUCT_NAME>"
        JUCE_APPLICATION_VERSION_STRING="$<TARGET_PROPERTY:inspect,JUCE_VERSION>")

target_link_libraries(${COMPONENT_NAME}
    PRIVATE
        jean-baptiste-lib
        juce::juce_audio_utils
        juce::juce_gui_extra
    PUBLIC
        juce::juce_recommended_config_flags
        juce::juce_recommended_lto_flags
        juce::juce_recommended_warning_flags)

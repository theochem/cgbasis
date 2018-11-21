# src-level CMakeLists

cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# Specify the project, its version and the used language
project(gbasis
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

# Set some compiler optimizations (https://stackoverflow.com/a/41361741/7930415)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Set CMake variables
include(${CMAKE_SOURCE_DIR}/cmake/SetCMakeEnvironment.cmake)

message(STATUS "Building ${LIBRARY_NAME} for ${CMAKE_BUILD_TYPE}")

# Find all required packages. Update CMake's MODULE_PATH to let CMake know that we have supplied our own FindXXX.cmake files in the cmake directory
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
include(${CMAKE_SOURCE_DIR}/cmake/FindPackages.cmake)

# Make a library target
add_library(${LIBRARY_NAME} ${LIBRARY_TYPE})

target_include_directories(${LIBRARY_NAME} PRIVATE ${Libint2_INCLUDE_DIRS})
target_link_libraries(${LIBRARY_NAME} PRIVATE ${Libint2_LIBRARIES})
include(CMakeListsVersion.txt.in)
set_target_properties(${LIBRARY_NAME} PROPERTIES SOVERSION ${GIT_TAG_SOVERSION})

# Add the sources for the library
include(src/CMakeLists.txt)
include(src/ones/CMakeLists.txt)
include(src/twos/CMakeLists.txt)
include(src/grids/CMakeLists.txt)
include(src/sparse/CMakeLists.txt)

# Install the headers
install(DIRECTORY "${CMAKE_SOURCE_DIR}/src/" # source directory
        DESTINATION "include/${LIBRARY_NAME}" # target directory
        FILES_MATCHING # install only matched files
        PATTERN "*.h" # select header files
)


# Make the install target for this library
include(${CMAKE_SOURCE_DIR}/cmake/InstallLibrary.cmake)

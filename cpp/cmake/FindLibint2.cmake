# FindLibint2.cmake. Find<package>.cmake-template from (https://cmake.org/Wiki/CMake:How_To_Find_Libraries#Writing_find_modules)
# Try to find Libint2

# The following installations of Libint2 are supported:
#   - /usr/local/libint
#   - any directory, but with an environment variable ${LIBINTROOT}, which contains the include, lib and share folders

# If found, this will define
#  Libint2_FOUND            Libint2 is available on the system
#  Libint2_INCLUDE_DIRS     the Libint2 include directories + dependency include directories
#  Libint2_LIBRARIES        the Libint2 library + dependency libraries


# We have to find the version! Luckily, Libint2 -when installed defaultly- provides a directory /usr/local/libint/x.y.z
# When the user has set ${LIBINTROOT} in the enviroment, this path can also be used
find_path(LIBINT_PREFIX libint2.h HINTS /usr/local/libint/*/ ENV{LIBINTROOT} ${LIBINT_ROOT} PATH_SUFFIXES libint2)

if("${LIBINT_PREFIX}" STREQUAL "LIBINT_PREFIX-NOTFOUND")
    message(FATAL_ERROR "Libint2 was not found in the default location /usr/local/libint/x.y.z or through the environment variables")
else()
    # Set FOUND
    set(Libint2_FOUND TRUE)

    # Set the INCLUDE_DIRS
    set(Libint2_INCLUDE_DIRS "${Libint2_INCLUDE_DIRS};${LIBINT_PREFIX}")

    # Set the LIBRARIES
    find_library(Libint2_LIBRARIES
            NAMES int2
            PATHS ${Libint2_LIBRARIES} ${LIBINT_PREFIX}/lib)

    message(STATUS "Libint2 was found at ${LIBINT_PREFIX}")
endif()

# In this CMake file, all CMake variables will be set


# Parse the project name into uppercase and lowercase
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPERCASE)  # ppercase is needed in version.hpp.in
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWERCASE)

# The name of the library should be equal to the project name
if(NOT LIBRARY_NAME)
    set(LIBRARY_NAME ${PROJECT_NAME})
endif()

# We want to make a shared library
set(LIBRARY_TYPE SHARED)
set(EXPORT_TYPE LIBRARY)

# Find the source folder
set(PROJECT_SOURCE_FOLDER ${CMAKE_SOURCE_DIR}/src)

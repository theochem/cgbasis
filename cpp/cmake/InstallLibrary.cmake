# In this CMake file, we will provide the necessary commands to make the install target for this library


# The target of this project is a library called "cpputil" (${PROJECT_NAME})
# To specify that this target should also be exported, we add the EXPORT option. This is used in conjuction with the
# install(EXPORT) command below
install(TARGETS ${LIBRARY_NAME}
        EXPORT ${LIBRARY_NAME} ${EXPORT_TYPE}
        DESTINATION lib) 

# Install the header files (not including version.hpp.in)
install(FILES ${PROJECT_INCLUDE_FILES} 
        DESTINATION include)


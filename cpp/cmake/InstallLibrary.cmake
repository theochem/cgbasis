# In this CMake file, we will provide the necessary commands to make the install target for this library


# The target of this project is a library called "gbasis" (${PROJECT_NAME})
install(TARGETS ${LIBRARY_NAME}
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include
        )

target_sources(${LIBRARY_NAME}
        PUBLIC
            "${CMAKE_CURRENT_LIST_DIR}/cholesky.cpp"
            "${CMAKE_CURRENT_LIST_DIR}/cholesky.h"
        )

target_include_directories(${LIBRARY_NAME} PUBLIC "${CMAKE_CURRENT_LIST_DIR}")

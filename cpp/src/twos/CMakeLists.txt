target_sources(${LIBRARY_NAME}
        PUBLIC
            "${CMAKE_CURRENT_LIST_DIR}/boys.cpp"
            "${CMAKE_CURRENT_LIST_DIR}/boys.h"
            "${CMAKE_CURRENT_LIST_DIR}/boys_inc.cpp"
            "${CMAKE_CURRENT_LIST_DIR}/gbw.cpp"
            "${CMAKE_CURRENT_LIST_DIR}/gbw.h"
            "${CMAKE_CURRENT_LIST_DIR}/ints4.cpp"
            "${CMAKE_CURRENT_LIST_DIR}/ints4.h"
            "${CMAKE_CURRENT_LIST_DIR}/iter_gb4.cpp"
            "${CMAKE_CURRENT_LIST_DIR}/iter_gb4.h"
        )

target_include_directories(${LIBRARY_NAME} PUBLIC "${CMAKE_CURRENT_LIST_DIR}")

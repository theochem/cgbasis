target_sources(${LIBRARY_NAME}
        PUBLIC
            "${CMAKE_CURRENT_LIST_DIR}/ints2.cpp"
            "${CMAKE_CURRENT_LIST_DIR}/iter_gb2.cpp"
            "${CMAKE_CURRENT_LIST_DIR}/iter_pow2.cpp"
            "${CMAKE_CURRENT_LIST_DIR}/nucpot.cpp"
            "${CMAKE_CURRENT_LIST_DIR}/ints2.h"
            "${CMAKE_CURRENT_LIST_DIR}/iter_gb2.h"
            "${CMAKE_CURRENT_LIST_DIR}/iter_pow2.h"
            "${CMAKE_CURRENT_LIST_DIR}/nucpot.h"
        )

target_include_directories(${LIBRARY_NAME} PUBLIC "${CMAKE_CURRENT_LIST_DIR}")

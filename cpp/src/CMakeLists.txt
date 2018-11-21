# src-level CMakeLists.txt

target_sources(${LIBRARY_NAME}
        PUBLIC
            "${CMAKE_CURRENT_LIST_DIR}/calc.cpp"
            "${CMAKE_CURRENT_LIST_DIR}/cartpure.cpp"
            "${CMAKE_CURRENT_LIST_DIR}/common.cpp"
            "${CMAKE_CURRENT_LIST_DIR}/gbasis.cpp"
            "${CMAKE_CURRENT_LIST_DIR}/calc.h"
            "${CMAKE_CURRENT_LIST_DIR}/cartpure.h"
            "${CMAKE_CURRENT_LIST_DIR}/common.h"
            "${CMAKE_CURRENT_LIST_DIR}/gbasis.h"
        )

target_include_directories(${LIBRARY_NAME} PUBLIC "${CMAKE_CURRENT_LIST_DIR}")

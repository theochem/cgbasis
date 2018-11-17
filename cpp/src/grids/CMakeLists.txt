target_sources(${LIBRARY_NAME}
        PUBLIC
            "${CMAKE_CURRENT_LIST_DIR}/fns.cpp"
            "${CMAKE_CURRENT_LIST_DIR}/iter_pow1.cpp"
            "${CMAKE_CURRENT_LIST_DIR}/nucpot.cpp"
	    "${CMAKE_CURRENT_LIST_DIR}/iter_gb1.cpp"
            "${CMAKE_CURRENT_LIST_DIR}/fns.h"
            "${CMAKE_CURRENT_LIST_DIR}/iter_pow1.h"
            "${CMAKE_CURRENT_LIST_DIR}/nucpot.h"
	    "${CMAKE_CURRENT_LIST_DIR}/iter_gb1.h"
        )

target_include_directories(${LIBRARY_NAME} PUBLIC "${CMAKE_CURRENT_LIST_DIR}")

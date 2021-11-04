set(SOURCE_HEADERS
        ${CMAKE_CURRENT_LIST_DIR}/include/HubbardKanamori.h
        ${CMAKE_CURRENT_LIST_DIR}/include/HubbardKanamoriAll.h
        ${CMAKE_CURRENT_LIST_DIR}/include/HubbardKanamoriMeasureFixSDSD.h
        ${CMAKE_CURRENT_LIST_DIR}/include/HubbardKanamoriMeasureObserveSDSD.h
        ${CMAKE_CURRENT_LIST_DIR}/include/HubbardKanamoriSDOperation.h
        )

add_subdirectory(include)
add_subdirectory(source)
add_subdirectory(test)
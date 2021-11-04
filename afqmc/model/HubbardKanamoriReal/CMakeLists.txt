set(SOURCE_HEADERS
        ${CMAKE_CURRENT_LIST_DIR}/include/HubbardKanamoriReal.h
        ${CMAKE_CURRENT_LIST_DIR}/include/HubbardKanamoriRealAll.h
        ${CMAKE_CURRENT_LIST_DIR}/include/HubbardKanamoriRealMeasureFixSDSD.h
        ${CMAKE_CURRENT_LIST_DIR}/include/HubbardKanamoriRealMeasureObserveSDSD.h
        ${CMAKE_CURRENT_LIST_DIR}/include/HubbardKanamoriRealSDOperation.h
        )

add_subdirectory(include)
add_subdirectory(source)
add_subdirectory(test)
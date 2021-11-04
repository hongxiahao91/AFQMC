set(SOURCE_HEADERS
        ${CMAKE_CURRENT_LIST_DIR}/include/KanamoriInteractReal.h
        ${CMAKE_CURRENT_LIST_DIR}/include/KanamoriInteractRealAux.h
        ${CMAKE_CURRENT_LIST_DIR}/include/KanamoriInteractRealForce.h
        ${CMAKE_CURRENT_LIST_DIR}/include/KanamoriInteractRealSample.h
        )

add_subdirectory(include)
add_subdirectory(source)
add_subdirectory(test)
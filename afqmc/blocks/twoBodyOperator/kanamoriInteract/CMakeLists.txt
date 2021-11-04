set(SOURCE_HEADERS
        ${CMAKE_CURRENT_LIST_DIR}/include/KanamoriInteract.h
        ${CMAKE_CURRENT_LIST_DIR}/include/KanamoriInteractAux.h
        ${CMAKE_CURRENT_LIST_DIR}/include/KanamoriInteractForce.h
        ${CMAKE_CURRENT_LIST_DIR}/include/KanamoriInteractSample.h
        )

add_subdirectory(include)
add_subdirectory(source)
add_subdirectory(test)
set(SOURCE_FILES
        KanamoriInteract.cpp
        KanamoriInteractAux.cpp
        KanamoriInteractForce.cpp)

add_library(afqmcKanamoriInteractHaoObject OBJECT ${SOURCE_HEADERS} ${SOURCE_FILES})
add_library(afqmcKanamoriInteractHao STATIC $<TARGET_OBJECTS:afqmcKanamoriInteractHaoObject>)
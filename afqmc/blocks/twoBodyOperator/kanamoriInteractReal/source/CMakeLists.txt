set(SOURCE_FILES
        KanamoriInteractReal.cpp
        KanamoriInteractRealAux.cpp
        KanamoriInteractRealForce.cpp)

add_library(afqmcKanamoriInteractRealHaoObject OBJECT ${SOURCE_HEADERS} ${SOURCE_FILES})
add_library(afqmcKanamoriInteractRealHao STATIC $<TARGET_OBJECTS:afqmcKanamoriInteractRealHaoObject>)
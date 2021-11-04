set(SOURCE_FILES
        HubbardKanamori.cpp
        HubbardKanamoriMeasureFixSDSD.cpp
        HubbardKanamoriMeasureObserveSDSD.cpp
        HubbardKanamoriSDOperation.cpp)

add_library(afqmcHubbardKanamoriHaoObject OBJECT ${SOURCE_HEADERS} ${SOURCE_FILES})
add_library(afqmcHubbardKanamoriHao STATIC $<TARGET_OBJECTS:afqmcHubbardKanamoriHaoObject>)
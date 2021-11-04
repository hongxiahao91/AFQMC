set(SOURCE_FILES
        HubbardKanamoriReal.cpp
        HubbardKanamoriRealMeasureFixSDSD.cpp
        HubbardKanamoriRealMeasureObserveSDSD.cpp
        HubbardKanamoriRealSDOperation.cpp)

add_library(afqmcHubbardKanamoriRealHaoObject OBJECT ${SOURCE_HEADERS} ${SOURCE_FILES})
add_library(afqmcHubbardKanamoriRealHao STATIC $<TARGET_OBJECTS:afqmcHubbardKanamoriRealHaoObject>)
set(SOURCE_FILES
        logHopSDOperation.cpp)

add_library(afqmcLogHopSDOperationHaoObject OBJECT ${SOURCE_HEADERS} ${SOURCE_FILES})
add_library(afqmcLogHopSDOperationHao STATIC $<TARGET_OBJECTS:afqmcLogHopSDOperationHaoObject>)

set(SOURCE_FILES
        logHop.cpp)

add_library(afqmcLogHopHaoObject OBJECT ${SOURCE_HEADERS} ${SOURCE_FILES})
add_library(afqmcLogHopHao STATIC $<TARGET_OBJECTS:afqmcLogHopHaoObject>)

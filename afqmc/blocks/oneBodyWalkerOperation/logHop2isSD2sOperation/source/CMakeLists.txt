set(SOURCE_FILES
        logHop2isSD2sOperation.cpp)

add_library(afqmcLogHop2isSD2sOperationHaoObject OBJECT ${SOURCE_HEADERS} ${SOURCE_FILES})
add_library(afqmcLogHop2isSD2sOperationHao STATIC $<TARGET_OBJECTS:afqmcLogHop2isSD2sOperationHaoObject>)

set(SOURCE_FILES
        hop2isSD2sOperation.cpp)

add_library(afqmcHop2isSD2sOperationHaoObject OBJECT ${SOURCE_HEADERS} ${SOURCE_FILES})
add_library(afqmcHop2isSD2sOperationHao STATIC $<TARGET_OBJECTS:afqmcHop2isSD2sOperationHaoObject>)

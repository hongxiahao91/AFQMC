set(SOURCE_FILES
        SD2sSD2sOperation.cpp)

add_library(afqmcSD2sSD2sOperationHaoObject OBJECT ${SOURCE_HEADERS} ${SOURCE_FILES})
add_library(afqmcSD2sSD2sOperationHao STATIC $<TARGET_OBJECTS:afqmcSD2sSD2sOperationHaoObject>)

set(SOURCE_FILES
        casBasis.cpp)

add_library(afqmcCasBasisHaoObject OBJECT ${SOURCE_HEADERS} ${SOURCE_FILES})
add_library(afqmcCasBasisHao STATIC $<TARGET_OBJECTS:afqmcCasBasisHaoObject>)

set(SOURCE_FILES
        MDCas2sSD2isOperation.cpp)

add_library(afqmcMDCas2sSD2isOperationHaoObject OBJECT ${SOURCE_HEADERS} ${SOURCE_FILES})
add_library(afqmcMDCas2sSD2isOperationHao STATIC $<TARGET_OBJECTS:afqmcMDCas2sSD2isOperationHaoObject>)

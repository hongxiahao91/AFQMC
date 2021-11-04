set(SOURCE_FILES
        MDCas2sSD2sOperation.cpp)

add_library(afqmcMDCas2sSD2sOperationHaoObject OBJECT ${SOURCE_HEADERS} ${SOURCE_FILES})
add_library(afqmcMDCas2sSD2sOperationHao STATIC $<TARGET_OBJECTS:afqmcMDCas2sSD2sOperationHaoObject>)

set(SOURCE_FILES
        MDCas2s.cpp)

add_library(afqmcMDCas2sHaoObject OBJECT ${SOURCE_HEADERS} ${SOURCE_FILES})
add_library(afqmcMDCas2sHao STATIC $<TARGET_OBJECTS:afqmcMDCas2sHaoObject>)

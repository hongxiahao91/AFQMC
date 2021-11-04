set(SOURCE_FILES
        realMaterialMolecule.cpp
        realMaterialMoleculeMeasureFixedMDCas2sSD2is.cpp
        realMaterialMoleculeMeasureFixedMDCas2sSD2s.cpp
        realMaterialMoleculeMeasureFixedMDCas2sSD2sGreen.cpp
        realMaterialMoleculeMeasureFixedSD2sSD2is.cpp
        realMaterialMoleculeMeasureFixedSD2sSD2s.cpp
        realMaterialMoleculeMeasureObserveSD2sSD2s.cpp
        realMaterialMoleculeSD2isOperation.cpp
        realMaterialMoleculeSD2sOperation.cpp
        realMaterialMoleculeMDCas2sOperation.cpp)

add_library(afqmcRealMaterialMoleculeHaoObject OBJECT ${SOURCE_HEADERS} ${SOURCE_FILES})
add_library(afqmcRealMaterialMoleculeHao STATIC $<TARGET_OBJECTS:afqmcRealMaterialMoleculeHaoObject>)

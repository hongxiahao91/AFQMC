set(SOURCE_TESTS
        MDCas2sSD2sOperationTest.cpp)

add_executable(afqmcMDCas2sSD2sOperationTest ${SOURCE_HEADERS} ${SOURCE_TESTS} )

target_link_libraries(afqmcMDCas2sSD2sOperationTest afqmcMDCas2sSD2sOperationHao afqmcMDCas2sHao afqmcCasBasisHao afqmcSD2sHao
        commonHao gtest_mainHao ${EXT_LIBRARIES})

add_my_test(afqmcMDCas2sSD2sOperationTest "$<TARGET_FILE:afqmcMDCas2sSD2sOperationTest>")
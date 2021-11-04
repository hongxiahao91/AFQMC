set(SOURCE_TESTS
        MDCas2sSD2isOperationTest.cpp)

add_executable(afqmcMDCas2sSD2isOperationTest ${SOURCE_HEADERS} ${SOURCE_TESTS} )

target_link_libraries(afqmcMDCas2sSD2isOperationTest afqmcMDCas2sSD2isOperationHao afqmcMDCas2sHao afqmcCasBasisHao afqmcSD2isHao
        commonHao gtest_mainHao ${EXT_LIBRARIES})

add_my_test(afqmcMDCas2sSD2isOperationTest "$<TARGET_FILE:afqmcMDCas2sSD2isOperationTest>")
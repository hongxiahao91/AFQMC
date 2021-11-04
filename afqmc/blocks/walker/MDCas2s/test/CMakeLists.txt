set(SOURCE_TESTS
        MDCas2sTest.cpp)

add_executable(afqmcMDCas2sTest ${SOURCE_HEADERS} ${SOURCE_TESTS} )

target_link_libraries(afqmcMDCas2sTest afqmcMDCas2sHao afqmcCasBasisHao commonHao gtest_mainHao ${EXT_LIBRARIES})

add_my_test(afqmcMDCas2sTest "$<TARGET_FILE:afqmcMDCas2sTest>")

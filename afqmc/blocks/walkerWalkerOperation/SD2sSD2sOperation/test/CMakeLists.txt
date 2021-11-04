set(SOURCE_TESTS
        SD2sSD2sOperationTest.cpp)

add_executable(afqmcSD2sSD2sOperationTest ${SOURCE_HEADERS} ${SOURCE_TESTS} )

target_link_libraries(afqmcSD2sSD2sOperationTest afqmcSD2sSD2sOperationHao afqmcSD2sHao
        commonHao gtest_mainHao ${EXT_LIBRARIES})

add_my_test(afqmcSD2sSD2sOperationTest "$<TARGET_FILE:afqmcSD2sSD2sOperationTest>")
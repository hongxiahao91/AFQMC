set(SOURCE_TESTS
        logHopSDOperationTest.cpp)

add_executable(afqmcLogHopSDOperationTest ${SOURCE_HEADERS} ${SOURCE_TESTS} )

target_link_libraries(afqmcLogHopSDOperationTest afqmcLogHopSDOperationHao
        afqmcLogHopHao afqmcSDHao commonHao gtest_mainHao ${EXT_LIBRARIES})

add_my_test(afqmcLogHopSDOperationTest "$<TARGET_FILE:afqmcLogHopSDOperationTest>")
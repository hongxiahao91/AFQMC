set(SOURCE_TESTS
        logHop2isSD2sOperationTest.cpp)

add_executable(afqmcLogHop2isSD2sOperationTest ${SOURCE_HEADERS} ${SOURCE_TESTS} )

target_link_libraries(afqmcLogHop2isSD2sOperationTest afqmcLogHop2isSD2sOperationHao
        afqmcLogHop2isHao afqmcSD2sHao commonHao gtest_mainHao ${EXT_LIBRARIES})

add_my_test(afqmcLogHop2isSD2sOperationTest "$<TARGET_FILE:afqmcLogHop2isSD2sOperationTest>")
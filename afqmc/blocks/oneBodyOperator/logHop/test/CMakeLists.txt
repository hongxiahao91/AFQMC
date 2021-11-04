set(SOURCE_TESTS
        logHopTest.cpp)

add_executable(afqmcLogHopTest ${SOURCE_HEADERS} ${SOURCE_TESTS})

target_link_libraries(afqmcLogHopTest afqmcLogHopHao commonHao gtest_mainHao ${EXT_LIBRARIES})

add_my_test(afqmcLogHopTest "$<TARGET_FILE:afqmcLogHopTest>")

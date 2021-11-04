set(SOURCE_TESTS
        HubbardKanamoriTest.cpp)

add_executable(afqmcHubbardKanamoriTest ${SOURCE_HEADERS} ${SOURCE_TESTS})

target_link_libraries(afqmcHubbardKanamoriTest afqmcHubbardKanamoriHao afqmcUtilitiesHao afqmcBlocksHao
        commonHao gtest_mainHao ${EXT_LIBRARIES})

add_my_test(afqmcHubbardKanamoriTest "$<TARGET_FILE:afqmcHubbardKanamoriTest>")
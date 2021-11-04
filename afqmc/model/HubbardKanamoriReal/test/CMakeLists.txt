set(SOURCE_TESTS
        HubbardKanamoriRealTest.cpp)

add_executable(afqmcHubbardKanamoriRealTest ${SOURCE_HEADERS} ${SOURCE_TESTS})

target_link_libraries(afqmcHubbardKanamoriRealTest afqmcHubbardKanamoriRealHao afqmcUtilitiesHao afqmcBlocksHao
        commonHao gtest_mainHao ${EXT_LIBRARIES})

add_my_test(afqmcHubbardKanamoriRealTest "$<TARGET_FILE:afqmcHubbardKanamoriRealTest>")
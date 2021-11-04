set(SOURCE_TESTS
        KanamoriInteractRealTest.cpp
        )

add_executable(afqmcKanamoriInteractRealTest ${SOURCE_HEADERS} ${SOURCE_TESTS} )

target_link_libraries(afqmcKanamoriInteractRealTest afqmcKanamoriInteractRealHao afqmcLogHopHao commonHao gtest_mainHao ${EXT_LIBRARIES})

add_my_test(afqmcKanamoriInteractRealTest "$<TARGET_FILE:afqmcKanamoriInteractRealTest>")
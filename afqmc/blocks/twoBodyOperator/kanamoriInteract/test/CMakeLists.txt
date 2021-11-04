set(SOURCE_TESTS
        KanamoriInteractTest.cpp
        )

add_executable(afqmcKanamoriInteractTest ${SOURCE_HEADERS} ${SOURCE_TESTS} )

target_link_libraries(afqmcKanamoriInteractTest afqmcKanamoriInteractHao afqmcLogHopHao commonHao gtest_mainHao ${EXT_LIBRARIES})

add_my_test(afqmcKanamoriInteractTest "$<TARGET_FILE:afqmcKanamoriInteractTest>")
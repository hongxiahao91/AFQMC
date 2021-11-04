set(SOURCE_TESTS
        hop2isSD2sOperationTest.cpp)

add_executable(afqmcHop2isSD2sOperationTest ${SOURCE_HEADERS} ${SOURCE_TESTS} )

target_link_libraries(afqmcHop2isSD2sOperationTest afqmcHop2isSD2sOperationHao
        afqmcHop2isHao afqmcSD2sHao commonHao gtest_mainHao ${EXT_LIBRARIES})

add_my_test(afqmcHop2isSD2sOperationTest "$<TARGET_FILE:afqmcHop2isSD2sOperationTest>")

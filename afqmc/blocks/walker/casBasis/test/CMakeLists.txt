set(SOURCE_TESTS
        casBasisTest.cpp)

add_executable(afqmcCasBasisTest ${SOURCE_HEADERS} ${SOURCE_TESTS} )

target_link_libraries(afqmcCasBasisTest afqmcCasBasisHao commonHao gtest_mainHao ${EXT_LIBRARIES})

add_my_test(afqmcCasBasisTest "$<TARGET_FILE:afqmcCasBasisTest>")

add_subdirectory(NiNj)
add_subdirectory(NiupNidn)
add_subdirectory(choleskyReal)
add_subdirectory(kanamoriInteract)
add_subdirectory(kanamoriInteractReal)

file(RELATIVE_PATH relativeDirectorPath ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_LIST_DIR})
install(FILES twoBodyOperator.h DESTINATION "include/${relativeDirectorPath}")
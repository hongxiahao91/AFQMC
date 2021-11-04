add_subdirectory(SDSDOperation)
add_subdirectory(SD2sSD2isOperation)
add_subdirectory(SD2sSD2sOperation)
add_subdirectory(MDCas2sSD2isOperation)
add_subdirectory(MDCas2sSD2sOperation)

file(RELATIVE_PATH relativeDirectorPath ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_LIST_DIR})
install(FILES walkerWalkerOperation.h DESTINATION "include/${relativeDirectorPath}")
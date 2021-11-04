add_subdirectory(SD)
add_subdirectory(SD2s)
add_subdirectory(SD2is)
add_subdirectory(casBasis)
add_subdirectory(MDCas2s)

file(RELATIVE_PATH relativeDirectorPath ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_LIST_DIR})
install(FILES walker.h DESTINATION "include/${relativeDirectorPath}")

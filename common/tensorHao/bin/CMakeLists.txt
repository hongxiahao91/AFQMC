#For bin:diag_cpu
add_executable(diag_cpu ${SOURCE_HEADERS} diag_cpu.cpp)
target_link_libraries(diag_cpu tensorHao mpiHao ${EXT_LIBRARIES})
set_target_properties(diag_cpu PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install(TARGETS diag_cpu DESTINATION bin)

#For bin:diag_magma
add_executable(diag_magma ${SOURCE_HEADERS} diag_magma.cpp)
target_link_libraries(diag_magma tensorHao mpiHao ${EXT_LIBRARIES})
set_target_properties(diag_magma PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install(TARGETS diag_magma DESTINATION bin)
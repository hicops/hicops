
add_executable(hicops ${_EXCLUDE}
    ${CMAKE_CURRENT_LIST_DIR}/hicops.cpp)

target_include_directories(hicops PRIVATE ${CMAKE_CURRENT_LIST_DIR}/../core/include)

target_link_libraries(hicops PRIVATE hicops-core)
target_link_libraries(hicops PRIVATE ${MPI_LIBRARIES})

if(USE_MPI AND MPI_FOUND)
    target_compile_options(hicops PUBLIC -DUSE_MPI)
endif()

set_target_properties(hicops PROPERTIES INSTALL_RPATH_USE_LINK_PATH ON)

# installation
install(TARGETS hicops DESTINATION ${CMAKE_INSTALL_LIBDIR})

file(GLOB c_headers                      ${CMAKE_CURRENT_LIST_DIR}/include/*.h)
file(GLOB cxx_headers                    ${CMAKE_CURRENT_LIST_DIR}/include/*.hpp)
file(GLOB c_sources                 ${CMAKE_CURRENT_LIST_DIR}/*.c)
file(GLOB cxx_sources               ${CMAKE_CURRENT_LIST_DIR}/*.cpp)

add_executable(counter ${cxx_sources} ${c_sources})

target_include_directories(counter PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include)
set_target_properties(counter PROPERTIES INSTALL_RPATH_USE_LINK_PATH ON)

# installation
install(TARGETS counter DESTINATION ${CMAKE_INSTALL_LIBDIR})
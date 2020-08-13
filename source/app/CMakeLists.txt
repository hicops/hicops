
project(hicops LANGUAGES C CXX)

# do we need timemory
if (USE_TIMEMORY)
    set(timemory_FIND_COMPONENTS_INTERFACE timemory-hicops)
    set(COMPONENTS headers compile-options analysis-tools cpu-roofline OPTIONAL_COMPONENTS cxx)

    find_package(timemory REQUIRED COMPONENTS ${COMPONENTS} OPTIONAL_COMPONENTS mpi)

    add_library(cpu-fp-single INTERFACE)
    add_library(cpu-fp-double INTERFACE)

    target_compile_definitions(cpu-fp-single INTERFACE ROOFLINE_FP_BYTES=4)
    target_compile_definitions(cpu-fp-double INTERFACE ROOFLINE_FP_BYTES=8)
endif()

add_executable(hicops ${_EXCLUDE}
    ${CMAKE_CURRENT_LIST_DIR}/hicops.cpp)

if(USE_MPI AND MPI_FOUND)
    target_compile_options(hicops PUBLIC -DUSE_MPI)
endif()

target_include_directories(hicops PRIVATE ${CMAKE_CURRENT_LIST_DIR}/../core/include)
target_link_libraries(hicops hicops-core ${MPI_LIBRARIES} timemory-hicops)

set_target_properties(hicops-core
    PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
        INSTALL_RPATH_USE_LINK_PATH ON
)

# installation
install(TARGETS hicops DESTINATION ${CMAKE_INSTALL_LIBDIR})
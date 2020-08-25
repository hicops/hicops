
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

include_directories(BEFORE ${CMAKE_CURRENT_LIST_DIR}/include)

set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME library)

#----------------------------------------------------------------------------------------#
# Locate sources and headers for this project
# - headers are included so they will show up in IDEs
#
file(GLOB c_headers                 ${CMAKE_CURRENT_LIST_DIR}/include/*.h)
file(GLOB cxx_headers               ${CMAKE_CURRENT_LIST_DIR}/include/*.hpp)
file(GLOB c_sources                 ${CMAKE_CURRENT_LIST_DIR}/*.c)
file(GLOB cxx_sources               ${CMAKE_CURRENT_LIST_DIR}/*.cpp)

# all the headers
set(hicops_glob_headers
    ${c_headers} ${cxx_headers})

if (USE_TIMEMORY)
    set(timemory_FIND_COMPONENTS_INTERFACE timemory-hicops-core)
    set(COMPONENTS headers compile-options analysis-tools papi
        #cpu-roofline
        OPTIONAL_COMPONENTS cxx)

#    add_library(cpu-fp-single INTERFACE)
#    add_library(cpu-fp-double INTERFACE)

#    target_compile_definitions(cpu-fp-single INTERFACE ROOFLINE_FP_BYTES=4)
#    target_compile_definitions(cpu-fp-double INTERFACE ROOFLINE_FP_BYTES=8)

    find_package(timemory REQUIRED COMPONENTS ${COMPONENTS} OPTIONAL_COMPONENTS mpi)
endif()

add_library (hicops-core SHARED ${c_sources} ${cxx_sources})

# to include the generated config.hpp
include_directories(${CMAKE_BINARY_DIR})

if (USE_TIMEMORY)
    target_link_libraries(hicops-core timemory-hicops-core)
endif()

if(USE_MPI AND MPI_FOUND)
    target_compile_options(hicops-core PUBLIC -DUSE_MPI)
endif()

set_target_properties(hicops-core
    PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)

install(TARGETS hicops-core DESTINATION lib)

INSTALL (
    DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.h*")
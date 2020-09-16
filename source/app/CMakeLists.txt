
project(hicops LANGUAGES C CXX)

# timemory headers are available
add_executable(hicops ${_EXCLUDE}
    ${CMAKE_CURRENT_LIST_DIR}/hicops.cpp)

# include core/include
target_include_directories(hicops PRIVATE ${CMAKE_CURRENT_LIST_DIR}/../core/include)

# to include the generated config.hpp
include_directories(${CMAKE_BINARY_DIR})

# link appropriate libraries
target_link_libraries(hicops hicops-core ${MPI_LIBRARIES})

# link to MPIP_LIBRARY if needed
if (USE_MPIP_LIBRARY)
    target_link_libraries(hicops ${MPIP_LIBRARY})
endif()

set_target_properties(hicops
    PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
        INSTALL_RPATH_USE_LINK_PATH ON
)

# installation
install(TARGETS hicops DESTINATION ${CMAKE_INSTALL_BINDIR})

#
# ---------------------------------------------------------------------------------------
#

if(NOT PYTHON_EXECUTABLE)
    set(PYTHON_EXECUTABLE "/usr/bin/env python")
endif()

configure_file(${CMAKE_CURRENT_LIST_DIR}/hicops_config.py
    ${CMAKE_BINARY_DIR}/scripts/wrappers/hicops_config @ONLY)


install(
    FILES
        ${CMAKE_BINARY_DIR}/scripts/wrappers/hicops_config

    DESTINATION bin
    PERMISSIONS
        OWNER_EXECUTE OWNER_READ OWNER_WRITE
        GROUP_EXECUTE GROUP_READ
        WORLD_EXECUTE WORLD_READ
    OPTIONAL)
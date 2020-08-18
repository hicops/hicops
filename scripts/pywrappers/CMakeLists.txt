
#
# ---------------------------------------------------------------------------------------
#

if(NOT PYTHON_EXECUTABLE)
    set(PYTHON_EXECUTABLE "/usr/bin/env python")
endif()

message(STATUS "Installing in.. " ${CMAKE_BINARY_DIR})

configure_file(${CMAKE_CURRENT_LIST_DIR}/hicops_comet.py
    ${CMAKE_BINARY_DIR}/wrappers/hicops_comet @ONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/psm2excel.py
    ${CMAKE_BINARY_DIR}/wrappers/psm2excel @ONLY)


install(
    FILES
        ${CMAKE_BINARY_DIR}/wrappers/hicops_comet
        ${CMAKE_BINARY_DIR}/wrappers/psm2excel

    DESTINATION bin/wrappers
    PERMISSIONS
        OWNER_EXECUTE OWNER_READ OWNER_WRITE
        GROUP_EXECUTE GROUP_READ
        WORLD_EXECUTE WORLD_READ
    OPTIONAL)
project(tools LANGUAGES C CXX)

#
# ---------------------------------------------------------------------------------------
#

if(NOT PYTHON_EXECUTABLE)
    set(PYTHON_EXECUTABLE "/usr/bin/env python")
endif()

configure_file(${CMAKE_CURRENT_LIST_DIR}/hicops_comet.py
    ${CMAKE_BINARY_DIR}/scripts/tools/hicops_comet @ONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/psm2excel.py
    ${CMAKE_BINARY_DIR}/scripts/tools/psm2excel @ONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/gen_expts.sh
    ${CMAKE_BINARY_DIR}/scripts/tools/gen_expts @ONLY)

install(
    FILES
        ${CMAKE_BINARY_DIR}/scripts/tools/hicops_comet
        ${CMAKE_BINARY_DIR}/scripts/tools/psm2excel

    DESTINATION bin/tools
    PERMISSIONS
        OWNER_EXECUTE OWNER_READ OWNER_WRITE
        GROUP_EXECUTE GROUP_READ
        WORLD_EXECUTE WORLD_READ
    OPTIONAL)

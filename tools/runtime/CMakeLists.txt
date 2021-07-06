project(runtime LANGUAGES C CXX)

#
# ---------------------------------------------------------------------------------------
#

if(NOT PYTHON_EXECUTABLE)
    set(PYTHON_EXECUTABLE "/usr/bin/env python")
endif()

if(NOT BASH_EXECUTABLE)
    set(BASH_EXECUTABLE "/bin/bash")
endif()

configure_file(${CMAKE_CURRENT_LIST_DIR}/hicops_comet.py
    ${CMAKE_BINARY_DIR}/tools/runtime/hicops_comet @ONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/hicops_expanse.py
    ${CMAKE_BINARY_DIR}/tools/runtime/hicops_expanse @ONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/psm2excel.py
    ${CMAKE_BINARY_DIR}/tools/runtime/psm2excel @ONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/psm2tsv.py
    ${CMAKE_BINARY_DIR}/tools/runtime/psm2tsv @ONLY)
    
configure_file(${CMAKE_CURRENT_LIST_DIR}/gen_expts.sh
    ${CMAKE_BINARY_DIR}/tools/runtime/gen_expts @ONLY)

install(
    FILES
        ${CMAKE_BINARY_DIR}/tools/runtime/hicops_comet
        ${CMAKE_BINARY_DIR}/tools/runtime/hicops_expanse
        ${CMAKE_BINARY_DIR}/tools/runtime/psm2excel
        ${CMAKE_BINARY_DIR}/tools/runtime/psm2tsv

    DESTINATION bin/tools
    PERMISSIONS
        OWNER_EXECUTE OWNER_READ OWNER_WRITE
        GROUP_EXECUTE GROUP_READ
        WORLD_EXECUTE WORLD_READ
    OPTIONAL)

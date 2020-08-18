
project(dbtools LANGUAGES C CXX)

#
# ---------------------------------------------------------------------------------------
#

if(NOT PYTHON_EXECUTABLE)
    set(PYTHON_EXECUTABLE "/usr/bin/env python")
endif()

if(NOT BASH_EXECUTABLE)
    set(BASH_EXECUTABLE "/bin/bash")
endif()

configure_file(${CMAKE_CURRENT_LIST_DIR}/dbprocessorLinux.py
    ${CMAKE_BINARY_DIR}/scripts/dbtools/dbprocessorLinux @ONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/dbprocessorWindows.py
    ${CMAKE_BINARY_DIR}/scripts/dbtools/dbprocessorWindows @ONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/cluster.sh
    ${CMAKE_BINARY_DIR}/scripts/dbtools/cluster @ONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/kthocc.sh
    ${CMAKE_BINARY_DIR}/scripts/dbtools/kthocc @ONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/sep_by_len.sh
    ${CMAKE_BINARY_DIR}/scripts/dbtools/sep_by_len @ONLY)

install(
    FILES
        ${CMAKE_BINARY_DIR}/scripts/dbtools/dbProcessorLinux
        ${CMAKE_BINARY_DIR}/scripts/dbtools/dbProcessorWindows
        ${CMAKE_BINARY_DIR}/scripts/dbtools/sep_by_len
        ${CMAKE_BINARY_DIR}/scripts/dbtools/cluster
        ${CMAKE_BINARY_DIR}/scripts/dbtools/kthocc

    DESTINATION dbtools
    PERMISSIONS
        OWNER_EXECUTE OWNER_READ OWNER_WRITE
        GROUP_EXECUTE GROUP_READ
        WORLD_EXECUTE WORLD_READ
    OPTIONAL)
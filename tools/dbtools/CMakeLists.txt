
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

configure_file(${CMAKE_CURRENT_LIST_DIR}/dbprep_linux.py
    ${CMAKE_BINARY_DIR}/tools/dbtools/dbprep_linux @ONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/dbprep_win.py
    ${CMAKE_BINARY_DIR}/tools/dbtools/dbprep_win @ONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/db_prep.sh
    ${CMAKE_BINARY_DIR}/tools/dbtools/db_prep @ONLY)

install(
    FILES
        ${CMAKE_BINARY_DIR}/tools/dbtools/db_prep
        ${CMAKE_BINARY_DIR}/tools/dbtools/dbprep_linux
        ${CMAKE_BINARY_DIR}/tools/dbtools/dbprep_win

    DESTINATION tools/dbtools
    PERMISSIONS
        OWNER_EXECUTE OWNER_READ OWNER_WRITE
        GROUP_EXECUTE GROUP_READ
        WORLD_EXECUTE WORLD_READ
    OPTIONAL)
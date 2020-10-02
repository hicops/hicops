project(ms2prep LANGUAGES C CXX)

#
# ---------------------------------------------------------------------------------------
#

if(NOT PYTHON_EXECUTABLE)
    set(PYTHON_EXECUTABLE "/usr/bin/env python")
endif()

if(NOT BASH_EXECUTABLE)
    set(BASH_EXECUTABLE "/bin/bash")
endif()

configure_file(${CMAKE_CURRENT_LIST_DIR}/extractms2.py
    ${CMAKE_BINARY_DIR}/tools/ms2prep/extractms2 @ONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/raw2ms2.sh
    ${CMAKE_BINARY_DIR}/tools/ms2prep/raw2ms2 @ONLY)

install(
    FILES
        ${CMAKE_BINARY_DIR}/tools/ms2prep/extractms2
        ${CMAKE_BINARY_DIR}/tools/ms2prep/raw2ms2

    DESTINATION tools/ms2prep
    PERMISSIONS
        OWNER_EXECUTE OWNER_READ OWNER_WRITE
        GROUP_EXECUTE GROUP_READ
        WORLD_EXECUTE WORLD_READ
    OPTIONAL)
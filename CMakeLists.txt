
##########################################################################################
#       CMake settings
##########################################################################################

cmake_minimum_required (VERSION 3.11 FATAL_ERROR)
project(hicops)

if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    set(MSG "")
    message(STATUS "Warning! Building from the source directory is not recommended")
    message(STATUS "If unintented, please remove 'CMakeCache.txt' and 'CMakeFiles'")
    message(STATUS "and build from a separate directory")
    message(WARNING "In-source build")
endif()

option(TIMEMORY_CCACHE_BUILD "Enable ccache build" OFF)
mark_as_advanced(TIMEMORY_CCACHE_BUILD)
if(TIMEMORY_CCACHE_BUILD)
    find_program(TIMEMORY_CCACHE_EXE ccache PATH_SUFFIXES bin)
    if(TIMEMORY_CCACHE_EXE)
        if(NOT EXISTS "${TIMEMORY_CCACHE_EXE}")
            message(WARNING
                "TIMEMORY_CCACHE_BUILD is ON but TIMEMORY_CCACHE_EXE (${TIMEMORY_CCACHE_EXE}) does not exist!")
        else()
            set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${TIMEMORY_CCACHE_EXE}")
        endif()
    endif()
endif()

cmake_policy(SET CMP0048 NEW)
cmake_policy(SET CMP0042 NEW)
if(NOT CMAKE_VERSION VERSION_LESS 3.13)
    cmake_policy(SET CMP0079 NEW)
endif()
if(NOT CMAKE_VERSION VERSION_LESS 3.14)
    cmake_policy(SET CMP0082 NEW)
endif()
# set these as the defaults
set(CMAKE_ENABLE_EXPORTS ON CACHE BOOL "Executable exports symbols for loadable modules")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON CACHE BOOL "Append directories in the linker search path")
set(CMAKE_POSITION_INDEPENDENT_CODE ON CACHE BOOL "Build position independent code")

##########################################################################################
#       CMake settings
##########################################################################################
option(USE_OMP "Enable OpenMP in hicops" ON)
option(USE_MPI "Enable MPI in hicops" ON)

find_package(OpenMP)
find_package(MPI)

if(OpenMP_FOUND AND USE_OMP)
    list(APPEND _OMP OpenMP::OpenMP_CXX)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -DUSE_OMP")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DUSE_OMP")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -DUSE_OMP")
endif()

if(USE_MPI AND MPI_FOUND)
    list(APPEND _MPI MPI::MPI_CXX)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
    target_compile_definitions(MPI::MPI_CXX INTERFACE USE_MPI)
endif()

#----------------------------------------------------------------------------------------#
#   hicops source
#----------------------------------------------------------------------------------------#

message(STATUS "Adding hicops source...")
add_subdirectory(source)

#----------------------------------------------------------------------------------------#
#  dbtool and ms2 scripts
#----------------------------------------------------------------------------------------#

message(STATUS "Adding python and shell scripts...")
add_subdirectory(scripts)
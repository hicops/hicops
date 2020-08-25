
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

# Set a default build type if none was specified
set(default_build_type "RelWithDebInfo")
 
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

##########################################################################################
#       CMake settings
##########################################################################################
option(USE_OMP "Enable OpenMP in hicops" ON)
option(USE_MPI "Enable MPI in hicops" ON)
option(USE_TIMEMORY "Enable Timemory instrumentation" OFF)

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

if (USE_TIMEMORY)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_TIMEMORY")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_TIMEMORY")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -DUSE_TIMEMORY")
    file(MAKE_DIRECTORY ${HOME}/tmp/mplconfigdir)
    set (ENV{MPLCONFIGDIR} ${HOME}/tmp/mplconfigdir)
endif()

#----------------------------------------------------------------------------------------#
#   configure config.hpp
#----------------------------------------------------------------------------------------#
message(STATUS "Configuring...")

option(TAILFIT "Use Tailfit method instead of Gumbelfit for e(x)" ON)
option(PROGRESS "Enable Progress Marks" ON)

# Allowed maximum peptide sequence length
if (NOT MAX_SEQ_LEN)
    set(MAX_SEQ_LEN   60)
endif()

# Maximum top k peaks to choose
if (NOT QALEN)
    set(QALEN   100)
endif()

# Maximum spectra batch size
if (NOT QCHUNK)
    set(QCHUNK   10000)
endif()

# Maximum hyperscore allowed
if (NOT MAX_HYPERSCORE)
    set(MAX_HYPERSCORE   100)
endif()

# configure the file
configure_file(config.hpp.in config.hpp @ONLY)

# install the configured file as well
install(FILES ${CMAKE_BINARY_DIR}/config.hpp DESTINATION include)

#----------------------------------------------------------------------------------------#
#   hicops source
#----------------------------------------------------------------------------------------#

message(STATUS "Adding hicops source...")
add_subdirectory(source)

#----------------------------------------------------------------------------------------#
#  dbtool, wrappers, ms2 scripts
#----------------------------------------------------------------------------------------#

message(STATUS "Adding tools and wrappers...")
add_subdirectory(scripts)
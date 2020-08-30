
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
option(USE_MPIP_LIBRARY "Enable MPIP instrumentation via Timemory" OFF)

##########################################################################################
#       MPI and OpenMP
##########################################################################################
find_package(OpenMP)
find_package(MPI)

if(OpenMP_FOUND AND USE_OMP)
    list(APPEND _OMP OpenMP::OpenMP_CXX)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -DUSE_OMP")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DUSE_OMP")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -DUSE_OMP")
    set (TMP ompt)
endif()

if(USE_MPI AND MPI_FOUND)
    list(APPEND _MPI MPI::MPI_CXX)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
    target_compile_definitions(MPI::MPI_CXX INTERFACE USE_MPI)
    set (TMP mpi ${TMP})
endif()

##########################################################################################
#       Configure Instrumentation
##########################################################################################

if (USE_TIMEMORY)
    set(timemory_FIND_COMPONENTS_INTERFACE timemory-hicops)
    set(COMPONENTS headers compile-options analysis-tools papi ${TMP}
        #cpu-roofline
        OPTIONAL_COMPONENTS cxx)

    # Find Timemory package
    find_package(timemory COMPONENTS ${COMPONENTS} OPTIONAL_COMPONENTS ${OPTIONAL_COMPONENTS})

    # if timemory package found then find mpip library
    if (timemory_DIR)
        file(MAKE_DIRECTORY $ENV{HOME}/mplconfigdir)
        message(STATUS "Timemory component interface enabled")
    else()
        message(WARNING "Timemory not found => Setting: USE_TIMEMORY=OFF")
        set(USE_TIMEMORY OFF CACHE BOOL "Timemory not found, disabling instrumentation" FORCE)
    endif()
else()
    message(STATUS "Timemory component interface disabled")
endif()

# safety check for MPIP library
if (USE_MPIP_LIBRARY AND NOT (USE_TIMEMORY AND USE_MPI))
    message(WARNING "MPIP instrumentation disabled. USE_MPIP_LIBRARY=ON requires USE_TIMEMORY=ON and USE_MPI=ON")
    set (USE_MPIP_LIBRARY OFF CACHE BOOL "Disabling MPIP instrumentation." FORCE)
endif()

# Find timemory-mpip.so library if required
if (USE_MPIP_LIBRARY)
    find_library(MPIP_LIBRARY
                 NAMES timemory-mpip 
                 HINTS ENV PATH
                       ENV LD_LIBRARY_PATH
                       ENV CMAKE_PREFIX_PATH
                 PATH_SUFFIXES lib lib64
                 DOC "Timemory MPIP instrumentation library"
                )
    if (MPIP_LIBRARY)
        message(STATUS "Timemory-MPIP library interface enabled")
    else ()
        message(WARNING "Timemory-MPIP library not found => Setting: USE_MPIP_LIBRARY=OFF")
        set (USE_MPIP_LIBRARY OFF CACHE BOOL "Timemory-MPIP library not found" FORCE)
    endif()
else()
    message(STATUS "Timemory-MPIP library interface disabled")
endif()

##########################################################################################
#       Configure Header files
##########################################################################################
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

##########################################################################################
#       Add sub-directories
##########################################################################################

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
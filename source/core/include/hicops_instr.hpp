/*
 * Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is licensed under the
 * Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
 * See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
 */

#pragma once

#include "config.hpp"
#include <chrono>

using time_point_t = std::chrono::system_clock::time_point;

#if defined (USE_TIMEMORY)

#    include "timemory/defines.h"
#    include "timemory/extern.hpp"

#    include "timemory/api.hpp"
#    include "timemory/components.hpp"  // 5.0
#    include "timemory/settings.hpp"    // 3.1

#    include "timemory/variadic/component_bundle.hpp"
#    include "timemory/variadic/auto_bundle.hpp"
#    include "timemory/variadic/auto_tuple.hpp"

#    include "timemory/types.hpp"                  // 3.5
#    include "timemory/macros.hpp"        // 3.2
#    include "timemory/config.hpp"

#    include "timemory/config/definition.hpp"
#    include "timemory/hash/definition.hpp"
#    include "timemory/operations/definition.hpp"
#    include "timemory/settings/definition.hpp"
#    include "timemory/storage/definition.hpp"
#    include "timemory/variadic/definition.hpp"

//--------------------------------------------------------------------------------------
//
//      Operating System
//
//--------------------------------------------------------------------------------------

#if defined(_WIN32) || defined(_WIN64) || defined(WIN32) || defined(WIN64)
#    if !defined(_WINDOWS)
#        define _WINDOWS
#    endif
#elif defined(__APPLE__) || defined(__MACH__)
#    if !defined(_MACOS)
#        define _MACOS
#    endif
#    if !defined(_UNIX)
#        define _UNIX
#    endif
#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#    if !defined(_LINUX)
#        define _LINUX
#    endif
#    if !defined(_UNIX)
#        define _UNIX
#    endif
#elif defined(__unix__) || defined(__unix) || defined(unix) || defined(_)
#    if !defined(_UNIX)
#        define _UNIX
#    endif
#endif

// shorthand
using namespace tim::component;

using papi_events_t = tim::component::papi_vector;
using wall_clock_t = tim::component::wall_clock;

struct hicops_tag
{};

// custom user bundle
using bundle_t          = user_bundle<10, hicops_tag>;

//
// Create tuple aliases
//
using time_tuple_t      = tim::auto_tuple<wall_clock, cpu_util>;
using mem_tuple_t       = tim::auto_tuple<written_bytes, read_bytes>;
using wall_tuple_t      = tim::auto_tuple<wall_clock>;

using prep_tuple_t      = tim::auto_tuple<wall_clock, read_bytes>;
using hw_counters_t     = tim::component_tuple<papi_events_t>;
using search_tuple_t    = tim::component_tuple<wall_clock, cpu_util, bundle_t>;
using comm_tuple_t      = tim::auto_tuple<wall_clock, written_bytes>;
using merge_tuple_t     = tim::auto_tuple<wall_clock, cpu_util, read_bytes>;

// 
// Mark templates as extern for faster compile times
//
extern template class tim::auto_tuple<wall_clock, cpu_util>;
extern template class tim::auto_tuple<written_bytes, read_bytes>;
extern template class tim::auto_tuple<wall_clock>;
extern template class tim::auto_tuple<wall_clock, read_bytes>;
extern template class tim::component_tuple<papi_events_t>;
extern template class tim::component_tuple<wall_clock, cpu_util, bundle_t>;
extern template class tim::auto_tuple<wall_clock, written_bytes>;
extern template class tim::auto_tuple<wall_clock, cpu_util, read_bytes>;

//
// MARKERS for time measurements
//
#define MARK_START(mark)          wall_clock_t mark; \
                                  mark.start()
#define MARK_END(mark)            mark.stop()
#define ELAPSED_SECONDS(mark)     mark.get()

#else

#define MARK_START(mark)          time_point_t mark = std::chrono::system_clock::now()
#define MARK_END(mark)            
#define ELAPSED_SECONDS(mark)    std::chrono::duration<double>(std::chrono::system_clock::now() - mark).count()

#endif // USE_TIMEMORY


//
// MACRO for printing elapsed time
//
#define PRINT_ELAPSED(es)         std::cout << "Elapsed Time: " << es << "s" << std::endl << std::endl

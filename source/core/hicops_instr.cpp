/*
 * Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is licensed under the
 * Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
 * See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
 */

#include "hicops_instr.hpp"

#if defined (USE_TIMEMORY)
#    define TIMEMORY_USE_MANAGER_EXTERN
#    include <timemory/timemory.hpp>

//
// Instantiate template classes
//
template class tim::auto_tuple<wall_clock, cpu_util>;
template class tim::auto_tuple<written_bytes, read_bytes>;
template class tim::auto_tuple<wall_clock>;
template class tim::auto_tuple<wall_clock, read_bytes>;
template class tim::component_tuple<papi_events_t>;
template class tim::component_tuple<wall_clock, cpu_util, bundle_t>;
template class tim::auto_tuple<wall_clock, written_bytes>;
template class tim::auto_tuple<wall_clock, cpu_util, read_bytes>;

#endif // USE_TIMEMORY
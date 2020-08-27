/*
 * Copyright (C) 2019  Muhammad Haseeb, Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "hicops_instr.hpp"

#if defined (USE_TIMEMORY)
#    define TIMEMORY_USE_MANAGER_EXTERN
#    include <timemory/timemory.hpp>

//
// Instantiate template classes
//
template class tim::auto_tuple<wall_clock, cpu_util, thread_cpu_util>;
template class tim::auto_tuple<written_bytes, read_bytes>;
template class tim::auto_tuple<wall_clock>;
template class tim::auto_tuple<wall_clock, thread_cpu_util, read_bytes>;
template class tim::component_tuple<papi_events_t>;
template class tim::component_tuple<wall_clock, cpu_util, bundle_t>;
template class tim::auto_tuple<wall_clock, thread_cpu_util, written_bytes>;
template class tim::auto_tuple<wall_clock, cpu_util, read_bytes>;

#endif // USE_TIMEMORY
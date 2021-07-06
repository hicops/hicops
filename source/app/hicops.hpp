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

#pragma once

#include "lbe.h"
#include "hicops_instr.hpp"
#include "argparse/argparse.hpp"

#if defined (USE_TIMEMORY)
#    include <timemory/timemory.hpp>
#endif

//--------------------------------------------------------------------------------------
//
//      timemory MPIP
//
//--------------------------------------------------------------------------------------

#if defined (USE_MPIP_LIBRARY)
extern "C"
{
    extern uint64_t timemory_start_mpip();
    extern uint64_t timemory_stop_mpip(uint64_t);
}
#endif // USE_MPIP_LIBRARY

#ifdef USE_MPI
#include <mpi.h>
#endif /* USE_MPI */
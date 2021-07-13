/*
 * Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is licensed under the
 * Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
 * See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
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
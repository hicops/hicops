/*
 * Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is licensed under the
 * Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
 * See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
 */

#pragma once

#include "common.hpp"
#include "msquery.hpp"
#include "lwqueue.h"
#include "lwbuff.h"
#include "hicops_instr.hpp"

//
// Superstep 2
//

namespace hcp
{
namespace ms2
{
// synchronize superstep 2
status_t synchronize();

// get instance of ptrs
MSQuery **& get_instance();

// initialize MS2 data index
status_t initialize(lwqueue<MSQuery *>**, int_t&, int_t&);

// delete the index
void deinitialize();

} // namespace ms2
} // namespace hcp

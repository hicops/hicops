/*
 * Copyright (C) 2021  Muhammad Haseeb, and Fahad Saeed
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
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
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

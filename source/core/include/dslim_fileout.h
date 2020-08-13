/*
 * Copyright (C) 2019 Muhammad Haseeb, Fahad Saeed
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

#ifndef INCLUDE_DSLIM_FILEOUT_H_
#define INCLUDE_DSLIM_FILEOUT_H_

#include <sstream>
#include "common.h"
#include "slm_dsts.h"
#include "slmerr.h"

/* Function Definitions */
status_t    DFile_PrintPartials(uint_t specid, Results *resPtr);
status_t    DFile_PrintScore(Index *index, uint_t specid, float_t pmass,
                           hCell* psm, double_t e_x, uint_t npsms);
status_t    DFile_InitFiles();
status_t    DFile_DeinitFiles();

#endif /* DSLIM_FILEOUT */

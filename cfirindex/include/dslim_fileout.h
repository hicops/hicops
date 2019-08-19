/*
 * Copyright (C) 2019 Muhammad Haseeb
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

/* Most beautiful includes I've ever seen */
#include <sstream>
#include "common.h"
#include "slm_dsts.h"
#include "slmerr.h"

/* Function Definitions */
STATUS  DFile_PrintPartials(UINT specid, Results *resPtr);
STATUS    DFile_PrintScore(Index *index, UINT specid, FLOAT pmass,
                           hCell* psm, DOUBLE e_x, UINT npsms);
STATUS    DFile_InitFiles();
STATUS    DFile_DeinitFiles();
ULONGLONG DFile_Factorial(ULONGLONG n);
STRING    DFile_Datetime();

#endif /* DSLIM_FILEOUT */

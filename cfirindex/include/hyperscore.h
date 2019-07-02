/*
 *  Copyright (C) 2019  Muhammad Usman Tariq, Muhammad Haseeb, Fahad Saeed
 *  Florida International University, Miami, FL
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef INCLUDE_HYPERSCORE_H_
#define INCLUDE_HYPERSCORE_H_

/* Most beautiful includes I've ever seen */
#include <sstream>
#include "common.h"
#include "slm_dsts.h"
#include "slmerr.h"

/* Function Definitions */
STATUS HYPERSCORE_Calculate(UINT specid, INT psid, FLOAT maxhv);
STATUS HYPERSCORE_WriteHeaders ();
STATUS HS_InitFile();
STATUS HS_DeinitFile();
STATUS HYPERSCORE_WriteToFile();
ULONGLONG HYPERSCORE_Factorial(ULONGLONG n);
std::string HYPERSCORE_Datetime();

#endif

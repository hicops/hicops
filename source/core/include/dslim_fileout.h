/*
 * Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is licensed under the
 * Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
 * See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
 */

#pragma once

#include <sstream>
#include "common.hpp"
#include "slm_dsts.h"
#include "slmerr.h"

/* Function Definitions */
status_t    DFile_PrintPartials(uint_t specid, Results *resPtr);
status_t    DFile_PrintScore(Index *index, uint_t specid, 
                             float_t pmass, hCell *psm, double_t e_x, uint_t npsms);
status_t    DFile_InitFiles();
status_t    DFile_DeinitFiles();
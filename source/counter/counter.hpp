/*
 * Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is licensed under the
 * Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
 * See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
 */

#pragma once 

#include <vector>
#include "common.hpp"
#include "slm_dsts.h"
#include "utils.h"

#define UNUSED_PARAM(x)              (void)(x)

status_t DBCounter(char_t *filename);

float_t UTILS_CalculatePepMass(AA *seq, uint_t len);

status_t UTILS_InitializeModInfo(SLM_vMods *vMods);

float_t UTILS_CalculateModMass(AA *seq, uint_t len, uint_t vModInfo);
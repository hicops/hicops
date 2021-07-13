/*
 * Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is licensed under the
 * Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
 * See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
 */

#pragma once

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <sstream>
#include "common.hpp"
#include "utils.h"

#define MAX_COMBS                          64

/*
 * FUNCTION: MODS_ModCounter
 *
 * DESCRIPTION: Counts the number of modifications for
 *              all peptides in SLM Peptide Index
 *
 * INPUT:
 * @threads   : Number of parallel threads
 * @conditions: Conditions of mod generation
 *
 * OUTPUT:
 * @cumulative: Number of mods
 */
ull_t  MODS_ModCounter();

/*
 * FUNCTION: MODS_GenerateMods
 *
 * DESCRIPTION: Generate modEntries for all peptide
 *              sequences in the SLM Peptide Index
 *
 * INPUT:
 * @threads   : Number of parallel threads
 * @modCount  : Number of mods that should be generated
 * @conditions: Conditions of mod generation
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t MODS_GenerateMods(Index * index);

status_t MODS_Initialize();

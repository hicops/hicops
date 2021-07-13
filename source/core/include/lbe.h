/*
 * Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is licensed under the
 * Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
 * See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
 */

#pragma once

#include <chrono>
#include <ctime>
#include <dirent.h>
#include "common.hpp"
#include "utils.h"
#include "mods.h"
#include "slmerr.h"
#include "dslim.h"
#include "msquery.hpp"

/*
 * FUNCTION: LBE_Initialize
 *
 * DESCRIPTION: Initialize internal peptides
 *              database from FASTA file
 *
 * INPUT:
 * @threads      : Number of parallel threads
 * @modconditions: String with mod conditions
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t LBE_Initialize(Index *index);

/*
 * FUNCTION: LBE_Deinitialize
 *
 * DESCRIPTION: Deallocate all memory and
 *              reset variables
 *
 * INPUT: none
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t LBE_Deinitialize(Index *index);

status_t LBE_GeneratePeps(Index *index);
/*
 * FUNCTION: LBE_Distribute
 *
 * DESCRIPTION: Apply Distribution Policy on peptides
 *
 * INPUT:
 * @threads   : Number of parallel threads
 * @policy    : Policy ID (enum)
 * @slm_chunks: Number of distribution chunks
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t LBE_Distribute(Index *index);

/*
 * FUNCTION: LBE_CountPeps
 *
 * DESCRIPTION: Count peptides in FASTA and the
 *              number of mods that will be generated
 *
 * INPUT:
 * @threads      : Number of parallel threads
 * @filename     : Path to FASTA file
 * @modconditions: Mod generation conditions
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t LBE_CountPeps(char_t *filename, Index *index, uint_t explen);

status_t LBE_CreatePartitions(Index *index);

BOOL LBE_ApplyPolicy(Index *index,  BOOL pepmod, uint_t key);

/*
 * FUNCTION: LBE_PrintHeader
 *
 * DESCRIPTION: Prints the LBE header
 *
 * INPUT : none
 * OUTPUT: none
 */
VOID   LBE_PrintHeader();
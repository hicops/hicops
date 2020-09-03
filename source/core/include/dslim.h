/*
 * Copyright (C) 2019  Muhammad Haseeb, and Fahad Saeed
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

#include "common.h"
#include "utils.h"
#include "lbe.h"
#include "slm_dsts.h"
#include "dslim_comm.h"
#include "expeRT.h"

/* Macros for SLM bitmask operations */
#define BYTE1                              8
#define BITNUM(z)                          ((z) % BYTE1)
#define IS_INIT(x,y)                       (((x)>>BITNUM(y)) & 0x1)
#define INIT(x,y)                          (x | (1 << BITNUM(y)))

#define NIBUFFS                            20

/* FUNCTION: DSLIM_Construct
 *
 * DESCRIPTION: Construct DSLIM chunks
 *
 * INPUT:
 * @threads: Number of parallel threads to launch
 * @modInfo: DSLIM mods Information
 * @peplen : Peptide Lengths to expect
 *
 * OUTPUT:
 * @status: status of execution
 */
status_t DSLIM_Construct(Index *index);

int_t DSLIM_GenerateIndex(Index *index, uint_t key);

status_t DSLIM_InitializeScorecard(Index *index, uint_t idxs);
/*
 * FUNCTION: DSLIM_AllocateMemory
 *
 * DESCRIPTION: Allocate memory for DSLIM chunks
 *
 * INPUT:
 * @chsize: Chunk size
 * @Chunks: Number of chunks
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t DSLIM_AllocateMemory(Index *index);

/*
 * FUNCTION: DSLIM_ConstructChunk
 *
 * INPUT:
 * @threads:      Number of parallel threads
 * @chunk_number: Chunk Index
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t DSLIM_ConstructChunk(uint_t threads, Index *index, uint_t chunk_number);

/*
 * FUNCTION: DSLIM_SLMTransform
 *
 * DESCRIPTION: Constructs SLIM Transform
 *
 * INPUT:
 * @threads     : Number of parallel threads
 * @chunk_number: Chunk Index
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t DSLIM_SLMTransform(uint_t threads, Index *index, uint_t chunk_number);

/*
 * FUNCTION: DSLIM_Optimize
 *
 * DESCRIPTION: Constructs SLIM Transform
 *
 * INPUT:
 * @threads     : Number of parallel threads
 * @index       : The SLM Index
 * @chunk_number: Chunk Index
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t DSLIM_Optimize(Index *index, uint_t chunk_number);

/*
 * FUNCTION: DSLIM_InitializeSC
 *
 * DESCRIPTION: Initialize Scorecard for DSLIM
 *
 * INPUT:
 * @threads: Number of parallel threads
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t DSLIM_InitializeSC(uint_t threads, Index *index);

/*
 * FUNCTION: DSLIM_Analyze
 *
 * DESCRIPTION: Analyze the DSLIM distribution
 *
 * INPUT:
 * @threads: Number of parallel threads
 * @avg    : Pointer to DSLIM mean load
 * @std    : Pointer to DSLIM load distribution
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t DSLIM_Analyze(uint_t threads, double_t &mean, double_t &std);

/*
 * FUNCTION: DSLIM_Deinitialize
 *
 * DESCRIPTION: Deallocate DSLIM memory
 *
 * INPUT: none
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t DSLIM_Deinitialize(Index *);

status_t DSLIM_DeallocateIonIndex(Index *);

status_t DSLIM_InitDistScore();

#ifdef USE_MPI
// declare this class here
class DSLIM_Comm;

status_t DSLIM_CarryForward(Index *, DSLIM_Comm *, expeRT *, hCell *, int_t);
#endif /* USE_MPI */

status_t DSLIM_DistScoreManager();

status_t DSLIM_DeallocatePepIndex(Index *);

status_t DSLIM_DeallocateSC();

status_t DSLIM_DeallocateSpecArr();

status_t DSLIM_SearchManager(Index *);

/* FUNCTION: DSLIM_QuerySpectrum
 *
 * DESCRIPTION: Query the DSLIM for all query peaks
 *              and count the number of hits per chunk
 *
 * INPUT:
 * @QA     : Query Spectra Array
 * @len    : Number of spectra in the array
 * @Matches: Array to fill in the number of hits per chunk
 * @threads: Number of parallel threads to launch
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t DSLIM_QuerySpectrum(Queries *ss, Index *index, uint_t indexchunks);

/* FUNCTION: DSLIM_WriteLIBSVM
 *
 * DESCRIPTION: Write the MS/MS spectra data in libsvm format
 *
 * INPUT:
 * @path   : Path for output data
 * @chno   : Chunk number
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t DSLIM_WriteLIBSVM(string_t path, uint_t peplen, uint_t chno);

status_t DSLIM_WriteCSV(string_t path, uint_t peplen, uint_t chno);

#ifdef DIAGNOSE
int_t DSLIM_TestBData();
#endif /* DIAGNOSE */

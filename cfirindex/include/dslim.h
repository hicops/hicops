/*
 * This file is part of PCDSFrame software
 * Copyright (C) 2019  Muhammad Haseeb, Fahad Saeed
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

#ifndef DSLIM_H_
#define DSLIM_H_

#include "common.h"
#include "utils.h"
#include "lbe.h"

/* Macros for SLM bitmask operations */
#define BYTE1                              8
#define BITNUM(z)                          ((z) % BYTE1)
#define IS_INIT(x,y)                       (((x)>>BITNUM(y)) & 0x1)
#define INIT(x,y)                          (x | (1 << BITNUM(y)))

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
STATUS DSLIM_Construct(Index *index);

INT DSLIM_GenerateIndex(Index *index, UINT key);

STATUS DSLIM_InitializeScorecard(Index *index, UINT idxs);
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
STATUS DSLIM_AllocateMemory(Index *index);

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
STATUS DSLIM_ConstructChunk(UINT threads, Index *index, UINT chunk_number);

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
STATUS DSLIM_SLMTransform(UINT threads, Index *index, UINT chunk_number);

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
STATUS DSLIM_Optimize(Index *index, UINT chunk_number);

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
STATUS DSLIM_InitializeSC(UINT threads, Index *index);

STATUS DSLIM_ModelSurvivalFunction(Results *resPtr);
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
STATUS DSLIM_Analyze(UINT threads, DOUBLE &mean, DOUBLE &std);

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
STATUS DSLIM_Deinitialize(Index *index);
STATUS DSLIM_DeallocateSC();

STATUS DSLIM_DeallocateSpecArr();

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
STATUS DSLIM_QuerySpectrum(Queries *ss, Index *index, UINT indexchunks);

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
STATUS DSLIM_WriteLIBSVM(STRING path, UINT peplen, UINT chno);



STATUS DSLIM_WriteCSV(STRING path, UINT peplen, UINT chno);
#endif /* DSLIM_H_ */

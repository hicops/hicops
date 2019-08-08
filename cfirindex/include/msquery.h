/*
 * This file is part of Load Balancing Algorithm for DSLIM
 *  Copyright (C) 2019  Muhammad Haseeb, Fahad Saeed
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

#include "common.h"
#include "utils.h"
#include "MSToolkitTypes.h"
#include "MSReader.h"
#include "MSObject.h"
#include "Spectrum.h"

/*
 * FUNCTION: MSQuery_InitializeQueryFile
 *
 * DESCRIPTION: Initialize structures using the query file
 *
 * INPUT:
 * @filename : Path to query file
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS MSQuery_InitializeQueryFile(CHAR *filename);

/*
 * FUNCTION: MSQuery_InitializeQueryFile
 *
 * DESCRIPTION: Initialize structures using only
 *              "count" spectra from file
 *
 * INPUT:
 * @start   : Start index of spectra in the query file
 * @count   : Number of spectra to use for initializing
 * @filename: Path to query file
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS MSQuery_InitializeQueryFile(UINT& start, UINT& count, CHAR *filename);

/*
 * FUNCTION: MSQuery_InitializeQueryFile
 *
 * DESCRIPTION: Extract a chunk of spectra from query file
 *
 * INPUT:
 * @QA      : Pointer to Query Array
 * @threads : Number of parallel threads
 *
 * OUTPUT:
 * @size: Size of the extracted chunk
 */
INT    MSQuery_ExtractQueryChunk(UINT *QA);

/*
 * FUNCTION: MSQuery_InitializeQueryFile
 *
 * DESCRIPTION: Extract a specific chunk of spectra from query file
 *
 * INPUT:
 * @start   : Start index of first spectrum to extract
 * @count   : Number of spectra to extract
 * @QA      : Pointer to Query Array
 * @threads : Number of parallel threads
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS MSQuery_ExtractQueryChunk(UINT start, UINT count, UINT *QA);

STATUS MSQuery_ExtractQueryChunk(UINT count, Queries &expSpecs);
/*
 * FUNCTION: MSQUERY_ProcessQuerySpectrum
 *
 * DESCRIPTION: Process a Query Spectrum and extract peaks
 *
 * INPUT:
 * @filename : Path to query file
 * @QAPtr    : Pointer to Query Array (dst)
 * @threads  : Number of parallel threads
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS MSQUERY_ProcessQuerySpectrum(CHAR *filename, UINT *QAPtr);

STATUS MSQUERY_ProcessQuerySpectrum(CHAR *filename, Queries &expSpecs, UINT offset);

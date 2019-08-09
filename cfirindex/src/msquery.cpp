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
 * GNU General Public License for more detailSpectrum.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */
#include "msquery.h"

using namespace std;
using namespace MSToolkit;

/* Global Variables */
MSReader     gReader;

extern gParams params;

/* Static global variables */
static UINT   firstScan = 0;
static UINT   currScan  = 0;
static UINT     QAcount = 0;
static STRING       MS2file;
static UINT nqchunks    = 0;
static UINT curr_chunk  = 0;
static UINT running_count   = 0;

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
STATUS MSQuery_InitializeQueryFile(CHAR *filename)
{
    UINT start, count;

    return MSQuery_InitializeQueryFile(start, count, filename);
}

/*
 * FUNCTION: MSQuery_InitializeQueryFile
 *
 * DESCRIPTION: Initialize structures using only
 *              @count spectra from file
 *
 * INPUT:
 * @start   : Start index of spectra in the query file
 * @count   : Number of spectra to use for initializing
 * @filename: Path to query file
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS MSQuery_InitializeQueryFile(UINT& start, UINT& count, CHAR *filename)
{
    STATUS status = SLM_SUCCESS;

    /* Read into local variables only */
    Spectrum Spectrum;
    MSReader tempReader;

    /* Initialize the variables */
    count = 0;
    start = 0;

    /* Add filters */
    gReader.setFilter(MS1);
    gReader.addFilter(MS2);
    gReader.addFilter(MSX);

    /* Add same filters for tempReader */
    tempReader.setFilter(MS1);
    tempReader.addFilter(MS2);
    tempReader.addFilter(MSX);

    /* Read the file */
    if (tempReader.readFile(filename, Spectrum))
    {

        while (Spectrum.getScanNumber() != 0)
        {
            /* Store the first scan number
             * into start */
            if (start == 0)
            {
                start = Spectrum.getScanNumber();
            }
            else
            {
                (VOID) Spectrum.getScanNumber();
            }
            count++;

            tempReader.readFile(NULL, Spectrum);
        }
    }
    else
    {
        status = ERR_INVLD_PARAM;
        count = 0;
        start = 0;
    }

    if (status == SLM_SUCCESS)
    {
        /* Set the Global Variables */
        firstScan = start;
        currScan  = firstScan;
        QAcount = count;
        MS2file = filename;
        curr_chunk = 0;
        running_count = 0;
        nqchunks = std::ceil(((double) count / QCHUNK));
    }

    return status;
}

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
INT MSQuery_ExtractQueryChunk(UINT *QA)
{
    STATUS status = SLM_SUCCESS;
    UINT *QAPtr = NULL;
    INT chunksize = 0;
    UINT spec = 0;
    UINT startspec = (curr_chunk * QCHUNK);
    UINT endspec = (curr_chunk + 1) * QCHUNK;

    /* Sanity checks */
    if (startspec >= QAcount)
    {
        chunksize = 0;
        status = ERR_INVLD_SIZE;
    }

    /* Query spectra available */
    if (status == SLM_SUCCESS)
    {
        if (endspec > QAcount)
        {
            endspec = QAcount;
        }

        for (spec = startspec; spec < endspec && status == SLM_SUCCESS; spec++)
        {
            /* Set the QAPtr */
            QAPtr = QA + ((spec - startspec) * QALEN);

            /* Process the Query Spectrum */
            status = MSQUERY_ProcessQuerySpectrum((CHAR *) MS2file.c_str(),
                                                  QAPtr);
        }

        if (status == SLM_SUCCESS)
        {
            /* Set the chunksize that was processed */
            chunksize = (spec - startspec);
            /* Increment the current chunk number */
            curr_chunk++;

        }
#ifdef DEBUG
        else
        {
            std::cout << "\nFATAL: There was an error processing MS2 spectra"<< endl
                      << "Chunk:\t" << curr_chunk << "\tSpec:\t" << spec << endl;
        }
#endif /* DEBUG */
    }

    return chunksize;
}

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
STATUS MSQuery_ExtractQueryChunk(UINT start, UINT count, UINT *QA)
{
    STATUS status = SLM_SUCCESS;
    UINT *QAPtr = NULL;
    UINT startspec = start;
    UINT endspec = start + count;

    if (startspec >= QAcount)
    {
        status = ERR_INVLD_SIZE;
    }

    if (status == SLM_SUCCESS)
    {
        if (endspec > QAcount)
        {
            endspec = QAcount;
        }

#ifdef _OPENMP
#pragma omp parallel for num_threads(params.threads) schedule(static)
#endif /* _OPENMP */
        for (UINT spec = startspec; spec < endspec; spec++)
        {
            QAPtr = QA + ((spec - startspec) * QALEN);

            status = MSQUERY_ProcessQuerySpectrum((CHAR *) MS2file.c_str(),
                                                 QAPtr);

        }
    }

    return status;
}

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
STATUS MSQUERY_ProcessQuerySpectrum(CHAR *filename, UINT *QAPtr)
{
    STATUS status = SLM_SUCCESS;
    Spectrum Spectrum;
    UINT threads = params.threads;

#ifndef _OPENMP
    LBE_UNUSED_PARAM(threads);
#endif

    /* TODO: FUTURE: There can be multiple spectra with different Z
     * in one MS/MS (MS2) file. How to deal/separate those
     */

    if (gReader.readFile(filename, Spectrum, currScan) && status == SLM_SUCCESS)
    {
        UINT SpectrumSize = Spectrum.size();
        UINT dIntArr[SpectrumSize];
        UINT  mzArray[SpectrumSize];

        for (UINT j = 0; j < SpectrumSize; j++)
        {
            mzArray[j] = (UINT)(Spectrum.at(j).mz * params.scale);
            dIntArr[j] = (UINT)(Spectrum.at(j).intensity * 1000);
        }

#ifdef _OPENMP
        /* Sort m/z based on Intensities */
        KeyVal_Parallel<UINT, UINT>(dIntArr, mzArray, SpectrumSize, threads);
#else
        KeyVal_Serial<UINT, UINT>(dIntArr, mzArray, SpectrumSize);
#endif /* _OPENMP */

        /* Check the size of spectrum */
        if (SpectrumSize >= QALEN)
        {
            /* Copy the last QALEN elements to QAPtr */
            std::memcpy(QAPtr, (mzArray + SpectrumSize - QALEN), (QALEN * sizeof(UINT)));
        }
        else
        {
            /* Fill in zeros which are not treated as trivial queries */
            std::memset(QAPtr, 0x0, ((QALEN - SpectrumSize) * sizeof(UINT)));

            /* Fill in rest of the spectrum */
            std::memcpy(QAPtr + (QALEN - SpectrumSize), mzArray, (SpectrumSize * sizeof(UINT)));
        }

#ifdef DEBUG
        for (UINT test = 0; test < QALEN; test++)
        {
            /* Test for integrity of QA */
            while (QAPtr[test] > params.max_mass * params.scale);
        }
#endif /* DEBUG */

        if(gReader.nextSpectrum(Spectrum))
        {
            currScan = (UINT)Spectrum.getScanNumber();
        }
    }
    else
    {
        status = ERR_INVLD_PARAM;
    }

    return status;

}

/*
 * FUNCTION: MSQuery_ExtractQueryChunk
 *
 * DESCRIPTION: Extract a specific chunk of spectra from query file
 *
 * INPUT:
 * @count   : Number of spectra to extract starting from the last spectrum read
 * @*expSpecs      : Pointer to experimental spectra struct
 * @threads : Number of parallel threads
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS MSQuery_ExtractQueryChunk(UINT count, Queries &expSpecs)
{
    STATUS status = SLM_SUCCESS;
    //UINT *QAPtr = NULL;

    /* half open interval [startspec, endspec) */
    UINT startspec = running_count;
    UINT endspec = running_count + count;

    if (startspec >= QAcount)
    {
        status = ERR_INVLD_SIZE;
    }

    if (status == SLM_SUCCESS)
    {
        if (endspec > QAcount)
        {
            endspec = QAcount;
            count = endspec - startspec;
        }
    }

    expSpecs.numSpecs = count;
    expSpecs.idx = new UINT[expSpecs.numSpecs + 1];
    expSpecs.precurse = new FLOAT[expSpecs.numSpecs];
    expSpecs.idx[0] = 0; //Set starting point to zero.

    MSReader tempReader;
    Spectrum spectrum;
    UINT l_peaks = 0;

    /* Add same filters for tempReader */
    tempReader.setFilter(MS1);
    tempReader.addFilter(MS2);
    tempReader.addFilter(MSX);

    if (tempReader.readFile(MS2file.c_str(), spectrum, currScan))
    {
        UINT l_count = 0;

        while (spectrum.getScanNumber() != 0 && l_count < count)
        {
            expSpecs.precurse[l_count] = spectrum.atZ(0).mh;

            l_peaks += spectrum.size() < QALEN? spectrum.size() : QALEN;
            l_count++;
//          expSpecs.idx[l_count] = l_peaks;

            tempReader.nextSpectrum(spectrum);
        }
    }

    if (status == SLM_SUCCESS)
    {
        expSpecs.numPeaks = l_peaks;
        expSpecs.moz = new UINT[expSpecs.numPeaks];
        expSpecs.intensity = new UINT[expSpecs.numPeaks];

        for (UINT spec = startspec; spec < endspec; spec++)
        {
            UINT index = spec - startspec; //0, 1, 2, ...., count

            status = MSQUERY_ProcessQuerySpectrum((CHAR *) MS2file.c_str(), expSpecs, index);
        }
    }

    running_count += count;

    return status;
}

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

STATUS MSQUERY_ProcessQuerySpectrum(CHAR *filename, Queries &expSpecs, UINT where)
{
    STATUS status = SLM_SUCCESS;
    Spectrum Spectrum;
    UINT threads = params.threads;

#ifndef _OPENMP
    LBE_UNUSED_PARAM(threads);
#endif
    /* TODO: FUTURE: There can be multiple spectra with different Z
     * in one MS/MS (MS2) file. How to deal/separate those
     */

    if (gReader.readFile(filename, Spectrum, currScan) && status == SLM_SUCCESS)
    {
        INT SpectrumSize = Spectrum.size();
        UINT dIntArr[SpectrumSize];
        UINT mzArray[SpectrumSize];

        for (INT j = 0; j < SpectrumSize; j++)
        {
            mzArray[j] = (UINT)(Spectrum.at(j).mz * params.scale);
            dIntArr[j] = (UINT)(Spectrum.at(j).intensity * 1000);
        }

#ifdef _OPENMP
        /* Sort m/z based on Intensities */
        KeyVal_Parallel<UINT, UINT>(dIntArr, mzArray, SpectrumSize, threads);
#else
        KeyVal_Serial<UINT, UINT>(dIntArr, mzArray, (UINT)SpectrumSize);
#endif /* _OPENMP */


        DOUBLE factor = ((DOUBLE)params.base_int/dIntArr[SpectrumSize - 1]);

        /* Set the highest peak to base intensity */
        dIntArr[SpectrumSize - 1] = params.base_int;
        UINT speclen = 1;

        /* Scale the rest of the peaks to the base peak */
        for (INT j = SpectrumSize - 2; j >= (SpectrumSize - QALEN) && j >= 0; j--)
        {
            dIntArr[j] *= factor;

            if (dIntArr[j] >= (UINT)params.min_int)
            {
                speclen++;
            }
        }

        /* Update the indices */
        UINT offset = expSpecs.idx[where];
        expSpecs.idx[where + 1] = expSpecs.idx[where] + speclen;

        /* Check the size of spectrum */
        if (speclen >= QALEN)
        {
            /* Copy the last QALEN elements to QAPtr */
            std::memcpy(&expSpecs.moz[offset], (mzArray + SpectrumSize - QALEN), (QALEN * sizeof(UINT)));
            std::memcpy(&expSpecs.intensity[offset], (dIntArr + SpectrumSize - QALEN), (QALEN * sizeof(UINT)));
        }
        else
        {
            /* Fill in zeros which are treated as trivial queries */
            std::memcpy(&expSpecs.moz[offset], (mzArray + SpectrumSize - speclen), (speclen * sizeof(UINT)));
            std::memcpy(&expSpecs.intensity[offset], (dIntArr + SpectrumSize - speclen), (speclen * sizeof(UINT)));
        }

#ifdef DEBUG
        for (UINT test = 0; test < QALEN; test++)
        {
            /* Test for integrity of QA */
            while (expSpecs->moz[test] > MAX_MASS * SCALE);
        }
#endif /* DEBUG */

        if(gReader.nextSpectrum(Spectrum))
        {
            currScan = (UINT)Spectrum.getScanNumber();
        }
    }
    else
    {
        status = ERR_INVLD_PARAM;
    }

    return status;

}

//Author: Usman
//needs fixing
UINT* UTILS_ExtractSpectra(CHAR *filename, UINT len)
{
    UINT *QA = new UINT[len];

    STATUS status = SLM_SUCCESS;

    //for (UINT qf = 0; qf < queryfiles.size(); qf++)
    {
        std::cout << std::endl << "Query File: " << filename << std::endl;

        /* Initialize Query MS/MS file */
        status = MSQuery_InitializeQueryFile(filename);

        /* DSLIM Query Algorithm */
        if (status == SLM_SUCCESS)
        {
            /* Extract a chunk of MS/MS spectra and
             * query against DSLIM Index */
            //for (; QA != NULL;)
            {
                /* Extract a chunk and return the chunksize */
                UINT ms2specs = MSQuery_ExtractQueryChunk(QA);

                /* If the chunksize is zero, all done */
                if (ms2specs <= 0)
                {
                    //break;
                }
            }
        }
    }

    return QA;
}

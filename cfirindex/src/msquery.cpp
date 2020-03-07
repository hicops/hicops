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

extern gParams params;

MSQuery::MSQuery()
{
    qfile = NULL;
    currPtr = 0;
    QAcount = 0;
    maxslen = 0;
    MS2file = new STRING;
    nqchunks = 0;
    curr_chunk = 0;
    running_count = 0;
    spectrum.intn = NULL;
    qfileIndex = 0;
    spectrum.mz = NULL;
    spectrum.SpectrumSize = 0;
    spectrum.prec_mz = 0;
    spectrum.Z = 0;
}

MSQuery::~MSQuery()
{
    currPtr = 0;
    QAcount = 0;
    nqchunks = 0;
    qfileIndex = 0;
    curr_chunk = 0;
    running_count = 0;
    maxslen = 0;

    if (MS2file != NULL)
    {
        delete MS2file;
        MS2file = NULL;
    }

    if (qfile != NULL)
    {
        delete qfile;
        qfile = NULL;
    }

    if (spectrum.intn != NULL)
    {
        delete[] spectrum.intn;
        spectrum.intn = NULL;
    }

    if (spectrum.mz != NULL)
    {
        delete[] spectrum.mz;
        spectrum.mz = NULL;
    }

    spectrum.SpectrumSize = 0;
    spectrum.prec_mz = 0;
    spectrum.Z = 0;
}

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
STATUS MSQuery::InitQueryFile(STRING *filename, INT fno)
{
    STATUS status = SLM_SUCCESS;

    /* Get a new ifstream object and open file */
    ifstream *qqfile = new ifstream(*filename);

    INT largestspec = 0;
    INT count = 0;
    INT specsize = 0;

    /* Check allocation */
    if (qqfile == NULL)
    {
        status = ERR_INVLD_PARAM;
    }

    /* Check if file opened */
    if (qqfile->is_open() && status == SLM_SUCCESS)
    {
        STRING line;

        /* While we still have lines in MS2 file */
        while (!qqfile->eof())
        {
            /* Read one line */
            getline(*qqfile, line);

            /* Empty line */
            if (line.empty())
            {
                continue;
            }
            /* Scan: (S) */
            else if (line[0] == 'S')
            {
                count++;
                largestspec = max(specsize, largestspec);
                specsize = 0;
            }
            /* Header: (H) */
            else if (line[0] == 'H' || line[0] == 'I' || line[0] == 'D' || line[0] == 'Z')
            {
                /* TODO: Decide what to do with header */
                continue;
            }
            /* MS/MS data: [m/z] [int] */
            else
            {
                specsize++;
            }
        }

        /* Check the largestspecsize */
        if (largestspec < 1)
        {
            status = ERR_INVLD_SIZE;
        }

        /* Initialize the file related params */
        if (status == SLM_SUCCESS)
        {
            currPtr  = 0;
            QAcount = count;
            *MS2file = *filename;
            curr_chunk = 0;
            running_count = 0;
            nqchunks = std::ceil(((double) QAcount / QCHUNK));
            qfileIndex = fno;
            maxslen = max(specsize, largestspec);

            /* Initialize to largest spectrum in file */
            spectrum.intn = new UINT[maxslen + 1];
            spectrum.mz = new UINT[maxslen + 1];
        }

        /* Close the file */
        qqfile->close();

        delete qqfile;
    }

    /* Return the status */
    return status;
}

STATUS MSQuery::ExtractQueryChunk(UINT count, Queries *expSpecs, INT &rem)
{
    STATUS status = SLM_SUCCESS;

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

    expSpecs->numSpecs = count;
    expSpecs->idx[0] = 0; //Set starting point to zero.

    if (qfile == NULL || qfile->is_open() == false)
    {
        /* Get a new ifstream object and open file */
        qfile = new ifstream(*MS2file);

        /* Check allocation */
        if (qfile == NULL)
        {
            status = ERR_INVLD_PARAM;
        }
    }

    /* Check if file opened */
    if (qfile->is_open() && status == SLM_SUCCESS)
    {
        for (UINT spec = startspec; spec < endspec; spec++)
        {
            ReadSpectrum();
            status = ProcessQuerySpectrum(expSpecs);
        }
    }

    if (status == SLM_SUCCESS)
    {
        /* Update the runnning count */
        running_count += count;

        /* Set the number of remaining spectra count */
        rem = QAcount - running_count;
    }

    return status;
}

VOID MSQuery::ReadSpectrum()
{
    STRING line;
    UINT speclen = 0;

    /* Check if this is the first spectrum in file */
    if (currPtr == 0)
    {
        BOOL scan = false;

        while (!qfile->eof())
        {
            /* Read one line */
            getline(*qfile, line);

            /* Empty line */
            if (line.empty() || line[0] == 'H' || line[0] == 'I' ||
                line[0] == 'D')
            {
                continue;
            }
            else if ( line[0] == 'Z')
            {
                CHAR *mh = strtok((CHAR *) line.c_str(), " \t");
                mh = strtok(NULL, " \t");
                STRING val = "1";

                if (mh != NULL)
                {
                    val = STRING(mh);
                    spectrum.Z = std::atoi(val.c_str());
                }

                val = "0.01";
                mh = strtok(NULL, " \t");

                if (mh != NULL)
                {
                    val = STRING(mh);
                    spectrum.prec_mz = (DOUBLE)std::atof(val.c_str());
                }
            }
            else if (line[0] == 'S')
            {
                if (scan == true)
                {
                    spectrum.SpectrumSize = speclen;
                    break;
                }
                else
                {
                    scan = true;
                }
            }
            /* Values */
            else
            {
                /* Split line into two DOUBLEs
                 * using space as delimiter */

                CHAR *mz1 = strtok((CHAR *) line.c_str(), " ");
                CHAR *intn1 = strtok(NULL, " ");
                STRING mz = "0.01";
                STRING intn = "0.01";

                if (mz1 != NULL)
                {
                    mz = STRING(mz1);
                }

                if (intn1 != NULL)
                {
                    intn = STRING(intn1);
                }

                spectrum.mz[speclen] = (UINT)((DOUBLE)std::atof(mz.c_str()) * params.scale);
                spectrum.intn[speclen] = (UINT)((DOUBLE)std::atof(intn.c_str()) * 1000);

                speclen++;
            }
        }
    }
    /* Not the first spectrum in file */
    else
    {
        while (!qfile->eof())
        {
            /* Read one line */
            getline(*qfile, line);

            /* Empty line */
            if (line.empty() || line[0] == 'H' || line[0] == 'I' ||
                line[0] == 'D')
            {
                continue;
            }
            else if ( line[0] == 'Z')
            {
                CHAR *mh = strtok((CHAR *) line.c_str(), " \t");
                mh = strtok(NULL, " \t");
                STRING val = "1";

                if (mh != NULL)
                {
                    val = STRING(mh);
                    spectrum.Z = std::atoi(val.c_str());
                }

                val = "0.01";
                mh = strtok(NULL, " \t");

                if (mh != NULL)
                {
                    val = STRING(mh);
                    spectrum.prec_mz = (DOUBLE)std::atof(val.c_str());
                }
            }
            else if (line[0] == 'S')
            {
                spectrum.SpectrumSize = speclen;
                break;
            }
            /* Values */
            else
            {
                /* Split line into two DOUBLEs
                 * using space as delimiter */
                CHAR *mz1 = strtok((CHAR *) line.c_str(), " ");
                CHAR *intn1 = strtok(NULL, " ");

                STRING mz = "0.0";
                STRING intn = "0.0";

                if (mz1 != NULL)
                {
                    mz = STRING(mz1);
                }
                if (intn1 != NULL)
                {
                    intn = STRING(intn1);
                }

                spectrum.mz[speclen] = (UINT)((DOUBLE)std::atof(mz.c_str()) * params.scale);
                spectrum.intn[speclen] = (UINT)((DOUBLE)std::atof(intn.c_str()) * 1000);

                speclen++;
            }
        }

        spectrum.SpectrumSize = speclen;
    }
}

STATUS MSQuery::ProcessQuerySpectrum(Queries *expSpecs)
{
    UINT *dIntArr = spectrum.intn;
    UINT *mzArray = spectrum.mz;
    INT SpectrumSize = spectrum.SpectrumSize;

    expSpecs->precurse[currPtr - running_count] = spectrum.prec_mz;

    KeyVal_Parallel<UINT, UINT>(dIntArr, mzArray, (UINT)SpectrumSize, 1);

    UINT speclen = 0;
    DOUBLE factor = 0;

    if (SpectrumSize > 0)
    {
        factor = ((DOUBLE) params.base_int / dIntArr[SpectrumSize - 1]);

        /* Set the highest peak to base intensity */
        dIntArr[SpectrumSize - 1] = params.base_int;
        speclen = 1;

        /* Scale the rest of the peaks to the base peak */
        for (INT j = SpectrumSize - 2; j >= (SpectrumSize - QALEN) && j >= 0; j--)
        {
            dIntArr[j] *= factor;

            if (dIntArr[j] >= (UINT) params.min_int)
            {
                speclen++;
            }
        }
    }

    /* Update the indices */
    UINT offset = expSpecs->idx[currPtr - running_count];
    expSpecs->idx[currPtr - running_count + 1] = expSpecs->idx[currPtr - running_count] + speclen;

    /* Check the size of spectrum */
    if (speclen >= QALEN)
    {
        /* Copy the last QALEN elements to expSpecs */
        std::memcpy(&expSpecs->moz[offset], (mzArray + SpectrumSize - QALEN), (QALEN * sizeof(UINT)));
        std::memcpy(&expSpecs->intensity[offset], (dIntArr + SpectrumSize - QALEN), (QALEN * sizeof(UINT)));
    }
    else
    {
        /* Copy the last speclen items to expSpecs */
        std::memcpy(&expSpecs->moz[offset], (mzArray + SpectrumSize - speclen), (speclen * sizeof(UINT)));
        std::memcpy(&expSpecs->intensity[offset], (dIntArr + SpectrumSize - speclen), (speclen * sizeof(UINT)));
    }

    currPtr += 1;

    return SLM_SUCCESS;
}

STATUS MSQuery::DeinitQueryFile()
{
    currPtr = 0;
    QAcount = 0;
    nqchunks = 0;
    curr_chunk = 0;
    running_count = 0;
    qfileIndex = 0;
    maxslen = 0;

    if (qfile != NULL)
    {
        qfile->close();

        delete qfile;
        qfile = NULL;
    }

    if (spectrum.intn != NULL)
    {
        delete[] spectrum.intn;
        spectrum.intn = NULL;
    }

    if (spectrum.mz != NULL)
    {
        delete[] spectrum.mz;
        spectrum.mz = NULL;
    }

    spectrum.SpectrumSize = 0;
    spectrum.prec_mz = 0;

    return SLM_SUCCESS;
}

BOOL MSQuery::isDeInit()
{
    return ((qfile == NULL) && (QAcount == 0));
}

/* Operator Overload - To copy to and from the work queue */
MSQuery& MSQuery::operator=(const MSQuery &rhs)
{
    this->MS2file = rhs.MS2file;
    this->QAcount = rhs.QAcount;
    this->currPtr = rhs.currPtr;
    this->curr_chunk = rhs.curr_chunk;
    this->maxslen = rhs.maxslen;
    this->nqchunks = rhs.nqchunks;
    this->qfile = rhs.qfile;
    this->running_count = rhs.running_count;
    this->spectrum = rhs.spectrum;
    this->qfileIndex = rhs.qfileIndex;

    return *this;
}

UINT MSQuery::getQfileIndex()
{
    return qfileIndex;
}

UINT MSQuery::getQAcount()
{
    return QAcount;
}

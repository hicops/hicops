/*
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
    MS2file = new string_t;
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
status_t MSQuery::InitQueryFile(string_t *filename, int_t fno)
{
    status_t status = SLM_SUCCESS;

    /* Get a new ifstream object and open file */
    ifstream *qqfile = new ifstream(*filename);

    int_t largestspec = 0;
    int_t count = 0;
    int_t specsize = 0;

    /* Check allocation
    if (qqfile == NULL)
    {
        status = ERR_INVLD_PARAM;
    }*/

    /* Check if file opened */
    if (qqfile->is_open() /*&& status == SLM_SUCCESS*/)
    {
        string_t line;

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

        largestspec = max(specsize, largestspec);

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
            spectrum.intn = new uint_t[maxslen + 1];
            spectrum.mz = new uint_t[maxslen + 1];
        }

        /* Close the file */
        qqfile->close();

        delete qqfile;
    }

    /* Return the status */
    return status;
}

status_t MSQuery::ExtractQueryChunk(uint_t count, Queries *expSpecs, int_t &rem)
{
    status_t status = SLM_SUCCESS;

    /* half open interval [startspec, endspec) */
    uint_t startspec = running_count;
    uint_t endspec = running_count + count;

    /*if (startspec >= QAcount)
    {
        status = ERR_INVLD_SIZE;
    }

    if (status == SLM_SUCCESS) */
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
        qfile = new ifstream;

        /* Check allocation
        if (qfile == NULL)
        {
            status = ERR_INVLD_PARAM;
        }*/

        qfile->open(MS2file->c_str());
    }

    /* Check if file opened */
    if (qfile->is_open() /*&& status == SLM_SUCCESS*/)
    {
        for (uint_t spec = startspec; spec < endspec; spec++)
        {
            ReadSpectrum();
            status = ProcessQuerySpectrum(expSpecs);
        }
    }

    //if (status == SLM_SUCCESS)
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
    string_t line;
    uint_t speclen = 0;
    char_t *saveptr;

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
            else if (line[0] == 'Z')
            {
                char_t *mh = strtok_r((char_t *) line.c_str(), " \t", &saveptr);
                mh = strtok_r(NULL, " \t", &saveptr);
                string_t val = "1";

                if (mh != NULL)
                {
                    val = string_t(mh);
                    spectrum.Z = std::atoi(val.c_str());
                }

                val = "0.01";
                mh = strtok_r(NULL, " \t", &saveptr);

                if (mh != NULL)
                {
                    val = string_t(mh);
                    spectrum.prec_mz = (double_t)std::atof(val.c_str());

                    // divide by charge to get the precursor m/z
                    spectrum.prec_mz /= MAX(1, spectrum.Z);
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

                char_t *mz1 = strtok_r((char_t *) line.c_str(), " ", &saveptr);
                char_t *intn1 = strtok_r(NULL, " ", &saveptr);
                string_t mz = "0.01";
                string_t intn = "0.01";

                if (mz1 != NULL)
                {
                    mz = string_t(mz1);
                }

                if (intn1 != NULL)
                {
                    intn = string_t(intn1);
                }

                spectrum.mz[speclen] = (uint_t)((double_t)std::atof(mz.c_str()) * params.scale);
                spectrum.intn[speclen] = (uint_t)((double_t)std::atof(intn.c_str()) * 1000);

                speclen++;
            }
        }

        spectrum.SpectrumSize = speclen;
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
                char_t *mh = strtok_r((char_t *) line.c_str(), " \t", &saveptr);
                mh = strtok_r(NULL, " \t", &saveptr);
                string_t val = "1";

                if (mh != NULL)
                {
                    val = string_t(mh);
                    spectrum.Z = std::atoi(val.c_str());
                }

                val = "0.01";
                mh = strtok_r(NULL, " \t", &saveptr);

                if (mh != NULL)
                {
                    val = string_t(mh);
                    spectrum.prec_mz = (double_t)std::atof(val.c_str());

                    // divide by charge to get the precursor m/z
                    spectrum.prec_mz /= MAX(1, spectrum.Z);
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
                char_t *mz1 = strtok_r((char_t *) line.c_str(), " ", &saveptr);
                char_t *intn1 = strtok_r(NULL, " ", &saveptr);

                string_t mz = "0.0";
                string_t intn = "0.0";

                if (mz1 != NULL)
                {
                    mz = string_t(mz1);
                }
                if (intn1 != NULL)
                {
                    intn = string_t(intn1);
                }

                spectrum.mz[speclen] = (uint_t)((double_t)std::atof(mz.c_str()) * params.scale);
                spectrum.intn[speclen] = (uint_t)((double_t)std::atof(intn.c_str()) * 1000);

                speclen++;
            }
        }

        spectrum.SpectrumSize = speclen;
    }
}

status_t MSQuery::ProcessQuerySpectrum(Queries *expSpecs)
{
    uint_t *dIntArr = spectrum.intn;
    uint_t *mzArray = spectrum.mz;
    int_t SpectrumSize = spectrum.SpectrumSize;

    expSpecs->precurse[currPtr - running_count] = spectrum.prec_mz;

    KeyVal_Parallel<uint_t, uint_t>(dIntArr, mzArray, (uint_t)SpectrumSize, 1);

    uint_t speclen = 0;
    double_t factor = 0;

    if (SpectrumSize > 0)
    {
        factor = ((double_t) params.base_int / dIntArr[SpectrumSize - 1]);

        /* Set the highest peak to base intensity */
        dIntArr[SpectrumSize - 1] = params.base_int;
        speclen = 1;

        /* Scale the rest of the peaks to the base peak */
        for (int_t j = SpectrumSize - 2; j >= (SpectrumSize - QALEN) && j >= 0; j--)
        {
            dIntArr[j] *= factor;

            if (dIntArr[j] >= (uint_t) params.min_int)
            {
                speclen++;
            }
        }
    }

    /* Update the indices */
    uint_t offset = expSpecs->idx[currPtr - running_count];
    expSpecs->idx[currPtr - running_count + 1] = expSpecs->idx[currPtr - running_count] + speclen;

    /* Check the size of spectrum */
    if (speclen >= QALEN)
    {
        /* Copy the last QALEN elements to expSpecs */
        std::memcpy(&expSpecs->moz[offset], (mzArray + SpectrumSize - QALEN), (QALEN * sizeof(uint_t)));
        std::memcpy(&expSpecs->intensity[offset], (dIntArr + SpectrumSize - QALEN), (QALEN * sizeof(uint_t)));
    }
    else
    {
        /* Copy the last speclen items to expSpecs */
        std::memcpy(&expSpecs->moz[offset], (mzArray + SpectrumSize - speclen), (speclen * sizeof(uint_t)));
        std::memcpy(&expSpecs->intensity[offset], (dIntArr + SpectrumSize - speclen), (speclen * sizeof(uint_t)));
    }

    expSpecs->numPeaks += speclen;

    currPtr += 1;

    return SLM_SUCCESS;
}

status_t MSQuery::DeinitQueryFile()
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

BOOL MSQuery::isDeInit() { return ((qfile == NULL) && (QAcount == 0)); }

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

MSQuery& MSQuery::operator=(const int_t &rhs)
{
    this->QAcount = rhs;
    this->currPtr = rhs;
    this->curr_chunk = rhs;
    this->maxslen = rhs;
    this->nqchunks = rhs;
    this->running_count = rhs;
    this->qfileIndex = rhs;

    return *this;
}

uint_t MSQuery::getQfileIndex() { return qfileIndex; }

uint_t MSQuery::getQAcount() { return QAcount; }

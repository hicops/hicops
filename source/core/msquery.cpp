/*
 * Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is licensed under the
 * Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
 * See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
 */

#include "msquery.hpp"

using namespace std;

extern gParams params;

MPI_Datatype MPI_info;
MPI_File fh;

MSQuery::MSQuery()
{
    qfile = NULL;
    currPtr = 0;
    curr_chunk = 0;
    running_count = 0;
    spectrum.intn = NULL;
    qfileIndex = 0;
    m_isinit = false;
    spectrum.mz = NULL;
    spectrum.SpectrumSize = 0;
    spectrum.prec_mz = 0;
    spectrum.Z = 0;
    spectrum.rtime = 0;

}

MSQuery::~MSQuery()
{
    currPtr = 0;
    qfileIndex = 0;
    curr_chunk = 0;
    running_count = 0;
    m_isinit = false;

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
    spectrum.rtime = 0;
}

/*
 * FUNCTION:
 *
 * DESCRIPTION: Initialize structures using the query file
 *
 * INPUT:
 * @filename : Path to query file
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t MSQuery::initialize(string_t *filename, int_t fno)
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
    if (qqfile->is_open() /*&& status =SLM_SUCCESS*/)
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
            info.QAcount = count;
            MS2file = *filename;
            curr_chunk = 0;
            running_count = 0;
            info.nqchunks = std::ceil(((double) info.QAcount / QCHUNK));
            qfileIndex = fno;
            info.maxslen = max(specsize, largestspec);

            /* Initialize to largest spectrum in file */
            spectrum.intn = new uint_t[info.maxslen + 1];
            spectrum.mz = new uint_t[info.maxslen + 1];

            m_isinit = true;
        }

        /* Close the file */
        qqfile->close();

        delete qqfile;
    }
    else
        status = ERR_FILE_ERROR;

    /* Return the status */
    return status;
}

//
// info initialized at remote process, initialize rest here
//
void MSQuery::vinitialize(string_t *filename, int_t fno)
{
    // reset these variables
    currPtr  = 0;
    curr_chunk = 0;
    running_count = 0;

    // set the file information
    MS2file = *filename;
    qfileIndex = fno;

    // init to the largest spectrum in file
    spectrum.intn = new uint_t[info.maxslen + 1];
    spectrum.mz = new uint_t[info.maxslen + 1];
    m_isinit = true;
}

status_t MSQuery::init_index()
{
    string fname = params.datapath + "/summary.io";

    // create a MPI data type
    MPI_Type_contiguous((int_t)(sizeof(info_t) / sizeof(int_t)),
                        MPI_INT,
                        &MPI_info);

    MPI_Type_commit(&MPI_info);

    // open the file
    status_t err = MPI_File_open(MPI_COMM_WORLD, fname.c_str(), (MPI_MODE_CREATE | MPI_MODE_WRONLY), MPI_INFO_NULL, &fh);

    return err;
}

status_t MSQuery::write_index() { return MPI_File_close(&fh); }

status_t MSQuery::read_index(info_t *findex, int_t size)
{
    // file name
    string fname = params.datapath + "/summary.io";
    MPI_File fh2;

    // open the file
    status_t status = MPI_File_open(MPI_COMM_WORLD, fname.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh2);

    // read the index
    status = MPI_File_read_all(fh2, findex, size, MPI_info, MPI_STATUS_IGNORE);

    // close the file
    MPI_File_close(&fh2);

    return status;
}

status_t MSQuery::archive(int_t index) { return MPI_File_write_at(fh, sizeof(info_t)*(index), &info, 1, MPI_info, MPI_STATUS_IGNORE); }

status_t MSQuery::extractbatch(uint_t count, Queries *expSpecs, int_t &rem)
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
        if (endspec > info.QAcount)
        {
            endspec = info.QAcount;
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

        qfile->open(MS2file.c_str());
    }

    /* Check if file opened */
    if (qfile->is_open() /*&& status =SLM_SUCCESS*/)
    {
        for (uint_t spec = startspec; spec < endspec; spec++)
        {
            readspectrum();
            status = pickpeaks(expSpecs);
        }
    }

    //if (status == SLM_SUCCESS)
    //{
        /* Update the runnning count */
        running_count += count;

        /* Set the number of remaining spectra count */
        rem = info.QAcount - running_count;
    //}

    return status;
}

VOID MSQuery::readspectrum()
{
    string_t line;
    uint_t speclen = 0;
    char_t *saveptr;
    char_t *Isave;

    /* Check if this is the first spectrum in file */
    if (currPtr == 0)
    {
        BOOL scan = false;

        while (!qfile->eof())
        {
            /* Read one line */
            getline(*qfile, line);

            /* Empty line */
            if (line.empty() || line[0] == 'H' || line[0] == 'D')
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
                }

                spectrum.Z = std::atoi(val.c_str());

                val = "0.01";
                mh = strtok_r(NULL, " \t", &saveptr);

                if (mh != NULL)
                {
                    val = string_t(mh);
                }

                spectrum.prec_mz = (double_t)std::atof(val.c_str());

            }
            else if (line[0] == 'I')
            {
                char_t *mh = strtok_r((char_t *) line.c_str(), " \t", &Isave);
                mh = strtok_r(NULL, " \t", &Isave);
                string_t val = "";

                if (mh != NULL)
                {
                    val = string_t(mh);
                }

                if (val.compare("RTime") == 0)
                {
                    val = "0.00";
                    mh = strtok_r(NULL, " \t", &Isave);

                    if (mh != NULL)
                    {
                        val = string_t(mh);
                    }

                    spectrum.rtime = (double_t)std::atof(val.c_str());
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
            if (line.empty() || line[0] == 'H' || line[0] == 'D')
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
                }

                spectrum.Z = std::atoi(val.c_str());

                val = "0.01";
                mh = strtok_r(NULL, " \t", &saveptr);

                if (mh != NULL)
                {
                    val = string_t(mh);
                }

                spectrum.prec_mz = (double_t)std::atof(val.c_str());

            }
            else if (line[0] == 'I')
            {
                char_t *mh = strtok_r((char_t *) line.c_str(), " \t", &Isave);
                mh = strtok_r(NULL, " \t", &Isave);
                string_t val = "";

                if (mh != NULL)
                {
                    val = string_t(mh);
                }

                if (val.compare("RTime") == 0)
                {
                    val = "0.00";
                    mh = strtok_r(NULL, " \t", &Isave);

                    if (mh != NULL)
                    {
                        val = string_t(mh);
                    }

                    spectrum.rtime = (double_t)std::atof(val.c_str());
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

status_t MSQuery::pickpeaks(Queries *expSpecs)
{
    uint_t *dIntArr = spectrum.intn;
    uint_t *mzArray = spectrum.mz;
    int_t SpectrumSize = spectrum.SpectrumSize;

    expSpecs->precurse[currPtr - running_count] = spectrum.prec_mz;
    expSpecs->charges[currPtr - running_count] = MAX(1, spectrum.Z);
    expSpecs->rtimes[currPtr - running_count] = MAX(0.0, spectrum.rtime);

    KeyVal_Parallel<uint_t, uint_t>(dIntArr, mzArray, (uint_t)SpectrumSize, 1);

    uint_t speclen = 0;
    double_t factor = 0;

    if (SpectrumSize > 0)
    {
        factor = ((double_t) params.base_int / dIntArr[SpectrumSize - 1]);

        /* Set the highest peak to base intensity */
        dIntArr[SpectrumSize - 1] = params.base_int;
        speclen = 1;

        uint_t l_min_int = params.min_int; //0.01 * dIntArr[SpectrumSize - 1];

        /* Scale the rest of the peaks to the base peak */
        for (int_t j = SpectrumSize - 2; j >= (SpectrumSize - QALEN) && j >= 0; j--)
        {
            dIntArr[j] *= factor;

            if (dIntArr[j] >= l_min_int)
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
    info.QAcount = 0;
    info.nqchunks = 0;
    curr_chunk = 0;
    running_count = 0;
    qfileIndex = 0;
    info.maxslen = 0;

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
    spectrum.Z = 0;
    spectrum.rtime = 0;

    return SLM_SUCCESS;
}

BOOL MSQuery::isDeInit() { return ((qfile == NULL) && (info.QAcount == 0)); }

/* Operator Overload - To copy to and from the work queue */
MSQuery& MSQuery::operator=(const MSQuery &rhs)
{
    this->MS2file = rhs.MS2file;
    this->info.QAcount = rhs.info.QAcount;
    this->currPtr = rhs.currPtr;
    this->curr_chunk = rhs.curr_chunk;
    this->info.maxslen = rhs.info.maxslen;
    this->info.nqchunks = rhs.info.nqchunks;
    this->qfile = rhs.qfile;
    this->running_count = rhs.running_count;
    this->spectrum = rhs.spectrum;
    this->qfileIndex = rhs.qfileIndex;

    return *this;
}

MSQuery& MSQuery::operator=(const int_t &rhs)
{
    this->info.QAcount = rhs;
    this->currPtr = rhs;
    this->curr_chunk = rhs;
    this->info.maxslen = rhs;
    this->info.nqchunks = rhs;
    this->running_count = rhs;
    this->qfileIndex = rhs;

    return *this;
}

uint_t MSQuery::getQfileIndex() { return qfileIndex; }

uint_t MSQuery::getQAcount() { return info.QAcount; }

uint_t& MSQuery::Nqchunks() { return info.nqchunks; }

uint_t& MSQuery::Curr_chunk() { return curr_chunk; }

info_t& MSQuery::Info() { return info; }

bool_t MSQuery::isinit() { return m_isinit; }

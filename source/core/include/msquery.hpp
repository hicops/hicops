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
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include <string>
#include "common.hpp"
#include "utils.h"

/* Spectrum */
struct _Spectrum
{
    uint_t *mz;
    uint_t *intn;
    uint_t SpectrumSize;
    double_t prec_mz;
    uint_t Z;
    double_t rtime;

    /* Overload the = operator - Required by MSQuery */
    _Spectrum &operator=(const _Spectrum &rhs)
    {
        this->mz = rhs.mz;
        this->intn = rhs.intn;
        this->SpectrumSize = rhs.SpectrumSize;
        this->prec_mz = rhs.prec_mz;
        this->Z = rhs.Z;
        this->rtime = rhs.rtime;

        return *this;
    }
} ;

struct _info
{
    uint_t maxslen;       // to be archived
    uint_t nqchunks;      // to be archived
    uint_t QAcount;       // to be archived

    _info() = default;
    ~_info() = default;

    /* Overload the = operator */
    _info &operator=(const _info &rhs)
    {
        this->maxslen = rhs.maxslen;
        this->nqchunks = rhs.nqchunks;
        this->QAcount = rhs.QAcount;

        return *this;
    }
};

using info_t = _info;
using spectrum_t = _Spectrum;

class MSQuery
{
private:
    /* Global Variables */
    uint_t currPtr;
    uint_t running_count;
    uint_t curr_chunk;
    info_t info;
    std::ifstream *qfile;
    uint_t qfileIndex;
    string_t MS2file;
    spectrum_t spectrum;
    bool_t m_isinit;

    VOID readspectrum();
    status_t pickpeaks(Queries *);

public:

    MSQuery();
    virtual ~MSQuery();
    uint_t getQAcount();
    status_t initialize(string_t *, int_t);
    void vinitialize(string_t *, int_t);
    static status_t init_index();
    static status_t write_index();
    static status_t read_index(info_t *, int);
    status_t archive(int_t);
    status_t extractbatch(uint_t, Queries *, int_t &);
    status_t DeinitQueryFile();
    BOOL isDeInit();
    uint_t getQfileIndex();
    MSQuery &operator=(const MSQuery &);
    MSQuery &operator=(const int_t &);

    uint_t& Curr_chunk();
    uint_t& Nqchunks();
    info_t& Info();

    bool_t isinit();

};
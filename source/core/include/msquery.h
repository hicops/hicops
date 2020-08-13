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

#ifndef MSQUERY_H_
#define MSQUERY_H_

#include <string>
#include "common.h"
#include "utils.h"

/* Spectrum */
typedef struct _Spectrum
{
    uint_t *mz;
    uint_t *intn;
    uint_t SpectrumSize;
    double_t prec_mz;
    uint_t Z;

    /* Overload the = operator - Required by MSQuery */
    _Spectrum &operator=(const _Spectrum &rhs)
    {
        this->mz = rhs.mz;
        this->intn = rhs.intn;
        this->SpectrumSize = rhs.SpectrumSize;
        this->prec_mz = rhs.prec_mz;
        this->Z = rhs.Z;

        return *this;
    }
} Spectrum;

class MSQuery
{
private:
    /* Global Variables */
    uint_t currPtr;
    uint_t running_count;
    uint_t maxslen;
    std::ifstream *qfile;
    uint_t qfileIndex;
    string_t *MS2file;
    Spectrum spectrum;

    VOID ReadSpectrum();
    status_t ProcessQuerySpectrum(Queries *);

public:
    uint_t QAcount;
    uint_t curr_chunk;
    uint_t nqchunks;

    MSQuery();
    virtual ~MSQuery();
    uint_t getQAcount();
    status_t InitQueryFile(string_t *filename, int_t fno);
    status_t ExtractQueryChunk(uint_t count, Queries *expSpecs, int_t &rem);
    status_t DeinitQueryFile();
    BOOL isDeInit();
    uint_t getQfileIndex();
    MSQuery &operator=(const MSQuery &rhs);
    MSQuery &operator=(const int_t &rhs);

};

#endif /* MSQUERY_H_ */

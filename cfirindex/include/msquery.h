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

#ifndef MSQUERY_H_
#define MSQUERY_H_

#include <string>
#include "common.h"
#include "utils.h"

/* Spectrum */
typedef struct _Spectrum
{
    UINT *mz;
    UINT *intn;
    UINT SpectrumSize;
    DOUBLE prec_mz;
} Spectrum;

class MSQuery
{
private:
    /* Global Variables */
    UINT currPtr;
    UINT QAcount;
    UINT nqchunks;
    UINT curr_chunk;
    UINT running_count;
    std::ifstream *qfile;
    STRING *MS2file;
    Spectrum spectrum;

    VOID ReadSpectrum();
    STATUS ProcessQuerySpectrum(Queries *);

public:

    MSQuery();
    virtual ~MSQuery();

    STATUS InitQueryFile(CHAR *filename);
    STATUS ExtractQueryChunk(UINT count, Queries *expSpecs, INT &rem);
    STATUS DeinitQueryFile();

};

#endif /* MSQUERY_H_ */

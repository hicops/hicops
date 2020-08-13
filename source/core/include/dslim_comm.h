/*
 * Copyright (C) 2020  Muhammad Haseeb, and Fahad Saeed
 * Florida International University (FIU), Miami, FL
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

#include "common.h"
#include "slm_dsts.h"
#include "expeRT.h"
#include "dslim.h"

#ifndef DSLIM_COMM_H_
#define DSLIM_COMM_H_

#ifdef USE_MPI

class DSLIM_Comm
{
private:
    int_t nBatches;

    int_t *sizeArray;
    int_t *fileArray;

    int_t myRXsize;

public:

    friend status_t DSLIM_CarryForward(Index *index, DSLIM_Comm *CommHandle, expeRT *ePtr, hCell *CandidatePSMS, int_t cpsmSize);
    DSLIM_Comm();
    DSLIM_Comm(int_t);
    virtual ~DSLIM_Comm();
    status_t AddBatch(int_t, int_t, int_t);
};

#endif /* USE_MPI */

#endif /* DSLIM_COMM_H_ */

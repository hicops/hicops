/*
 * Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is licensed under the
 * Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
 * See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
 */

#pragma once

#include "common.hpp"
#include "slm_dsts.h"
#include "expeRT.h"
#include "dslim.h"

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
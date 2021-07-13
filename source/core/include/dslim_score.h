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
#include <thread>
#include <unistd.h>
#include "config.hpp"
#include "slm_dsts.h"
#include "dslim.h"
#include "slmerr.h"
#include "utils.h"
#include "expeRT.h"

typedef struct _BorrowedData
{
    /* These pointers will be borrowed */
    expeRT *ePtr;
    hCell *heapArray;
    Index *index;
    int_t *sizeArray;
    int_t *fileArray;

    /* Dataset size */
    int_t cPSMsize;
    int_t nBatches;

    _BorrowedData()
    {
        ePtr = NULL;
        heapArray = NULL;
        index = NULL;
        sizeArray = NULL;
        fileArray = NULL;
        cPSMsize = 0;
        nBatches = 0;
    }

} BData;

#ifdef USE_MPI

class DSLIM_Score
{
private:

    /* Variables */
    int_t       threads;

    /* Dataset sizes */
    int_t       nSpectra;
    int_t       nBatches;
    int_t       myRXsize;

    /* These pointers will be borrowed */
    int_t      *sizeArray;
    int_t      *fileArray;
    expeRT   *ePtr;
    hCell    *heapArray;
    Index    *index;
    std::thread comm_thd;

    /* Data size that I expect to
     * receive from other processes */
    int_t      *rxSizes;
    int_t      *txSizes;

    /* key-values */
    int_t      *keys;

    fResult  *TxValues;
    fResult  *RxValues;

public:

    DSLIM_Score();
    DSLIM_Score(BData *bd);
    virtual  ~DSLIM_Score();

    status_t   CombineResults();

    status_t   ScatterScores();

    status_t   TXSizes(MPI_Request *, int_t *);
    status_t   RXSizes(MPI_Request *, int_t *);

    status_t   TXResults(MPI_Request *, int_t*);
    status_t   RXResults(MPI_Request *, int_t*);

    status_t   DisplayResults();

    status_t   Wait4RX();

    status_t   InitDataTypes();
    status_t   FreeDataTypes();
};

#endif // USE_MPI
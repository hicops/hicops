/*
 * Copyright (C) 2020 Muhammad Haseeb, and Fahad Saeed
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

#ifndef DSLIM_SCORE_H_
#define DSLIM_SCORE_H_

#include "common.h"
#include <unistd.h>
#include "config.h"
#include "slm_dsts.h"
#include "dslim.h"
#include "slmerr.h"
#include "utils.h"

class DSLIM_Score
{
private:

    /* Variables */
    INT       threads;

    /* Dataset sizes */
    INT       nSpectra;
    INT       nBatches;
    INT       myRXsize;

    /* These pointers will be borrowed */
    INT      *sizeArray;
    INT      *fileArray;
    BYICount *scPtr;
    partRes  *resPtr;
    hCell    *heapArray;
    Index    *index;
    THREAD   *comm_thd;

    /* Data size that I expect to
     * receive from other processes */
    INT      *rxSizes;
    INT      *txSizes;

    /* key-values */
    INT      *keys;

    fResult  *TxValues;
    fResult  *RxValues;

public:

    DSLIM_Score();
    DSLIM_Score(BData *bd);
    virtual  ~DSLIM_Score();

    STATUS   ComputeDistScores();

    STATUS   ScatterScores();

    STATUS   TXSizes(MPI_Request *, INT *);
    STATUS   RXSizes(MPI_Request *, INT *);

    STATUS   TXResults(MPI_Request *, INT*);
    STATUS   RXResults(MPI_Request *, INT*);

    STATUS   DisplayResults();

    STATUS   Wait4RX();

    STATUS   InitDataTypes();
    STATUS   FreeDataTypes();
};

#endif /* DSLIM_SCORE_H_ */

/*
 * Copyright (C) 2020 Muhammad Haseeb, Fahad Saeed
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

#include <semaphore.h>
#include <unistd.h>
#include "config.h"
#include "slm_dsts.h"
#include "dslim.h"
#include "slmerr.h"
#include "utils.h"

class DSLIM_Score
{
private:
    /* These pointers will be borrowed */
    Results **resPtr;
    hCell *heap;
    Index **index;
    INT *batchsizes;
    INT *filenumbers;

    /* Local variables */
    INT threads;
    BOOL isInit;
    THREAD comm_thd;

    /* Dataset size */
    INT dsetsize;
    INT num_batches;

    /* Data size that I expect to
     * receive from other processes */
    INT *rxSizes;
    INT *txSizes;

    /* key-values */
    INT *keys;
    fResult *values;

public:
    DSLIM_Score();
    virtual ~DSLIM_Score();

};

#endif /* DSLIM_SCORE_H_ */

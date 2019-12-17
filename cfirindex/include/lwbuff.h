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

#ifndef LWBUFF_H_
#define LWBUFF_H_

#include "common.h"
#include "lwqueue.h"
#include <semaphore.h>

using namespace std;

#define DEF_SIZE                    20

template <class T>
class lwbuff
{
private:

    INT cap;
    INT threshold;
    sem_t lockr;
    sem_t lockw;

    /* NOTE: The sizeof(readyQ) must be at least
     *       sizeof(workQ) + 1
     */
    lwqueue<T*> *readyQ;
    lwqueue<T*> *waitQ;

public:
    lwbuff()
    {
        cap = DEF_SIZE;

        threshold = cap/3;

        readyQ = new lwqueue<T *>(DEF_SIZE + 1, false);
        waitQ = new lwqueue<T *>(DEF_SIZE, false);

        sem_init(&lockr, 0, 1);
        sem_init(&lockw, 0, 1);
    }

    lwbuff(INT dcap, INT thr)
    {
        cap = dcap;
        threshold = thr;
        readyQ = new lwqueue<T*>(dcap + 1, false);
        waitQ = new lwqueue<T*>(dcap, false);

        sem_init(&lockr, 0, 1);
        sem_init(&lockw, 0, 1);
    }

    virtual ~lwbuff()
    {
        cap = 0;
        threshold = 0;

        vEmpty();

        delete readyQ;
        delete waitQ;

        sem_destroy(&lockr);
        sem_destroy(&lockw);

    }

    VOID Add(T *item)
    {
        waitQ->push(item);
    }

    VOID vEmpty()
    {
        /* Make sure both are empty */
        for (INT ii = 0; ii < readyQ->size(); ii ++)
        {
            readyQ->pop();
            waitQ->pop();
        }
    }

    VOID IODone(T *ptr)
    {
        readyQ->push(ptr);
    }

    VOID Replenish(T *ptr)
    {
        waitQ->push(ptr);
    }

    T *getIOPtr()
    {
        T *rtn = NULL;

        if (waitQ->size() > 0)
        {
            rtn = waitQ->front();
            waitQ->pop();
        }

        return rtn;
    }

    T *getWorkPtr()
    {
        T *rtn = NULL;

        if (readyQ->size() > 0)
        {
            rtn = readyQ->front();
            readyQ->pop();
        }

        return rtn;
    }

    T *releaseIOPtr(T *ptr)
    {
        T *rtn = NULL;

        if (ptr != NULL)
        {
            waitQ->push(ptr);
        }

        return rtn;
    }




    INT len()
    {
        return cap;
    }

    BOOL isEmptyReadyQ()
    {
        return readyQ->isEmpty();
    }

    BOOL isFullReadyQ()
    {
        return readyQ->isFull();
    }

    BOOL isEmptyWaitQ()
    {
        return waitQ->isEmpty();
    }

    BOOL isFullWaitQ()
    {
        return waitQ->isFull();
    }

    BOOL isUnderWaitQ()
    {
        return (waitQ->size() < threshold);
    }


    inline INT addOne(INT idd)
    {
        return ((idd+1) % cap);
    }

    STATUS lockr_()
    {
        return sem_wait(&lockr);
    }

    STATUS unlockr_()
    {
        return sem_post(&lockr);
    }

    STATUS lockw_()
    {
        return sem_wait(&lockw);
    }

    STATUS unlockw_()
    {
        return sem_post(&lockw);
    }

};

#endif /* LWBUFF_H_ */

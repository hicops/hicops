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
#include "lwqueue.h"
#include <semaphore.h>

using namespace std;

#define DEF_SIZE                    20

template <class T>
class lwbuff
{
private:

    int_t cap;
    int_t thr_low;
    int_t thr_high;
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

        thr_low = cap/4;
        thr_high = cap - thr_low;

        readyQ = new lwqueue<T *>(DEF_SIZE, false);
        waitQ = new lwqueue<T *>(DEF_SIZE, false);

        sem_init(&lockr, 0, 1);
        sem_init(&lockw, 0, 1);
    }

    lwbuff(int_t dcap, int_t lo, int_t hi)
    {
        cap = dcap;
        thr_low = lo;
        thr_high = hi;

        readyQ = new lwqueue<T*>(dcap, false);
        waitQ = new lwqueue<T*>(dcap, false);

        sem_init(&lockr, 0, 1);
        sem_init(&lockw, 0, 1);
    }

    virtual ~lwbuff()
    {
        cap = 0;
        thr_low = 0;
        thr_high = 0;

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
        for (int_t ii = 0; ii < readyQ->size(); ii ++)
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

    int_t len()
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

    int_t readyQStatus()
    {
        int_t sz = readyQ->size();

        if (sz < thr_low)
        {
            return -1;
        }
        else if (sz > thr_low && sz <= thr_high)
        {
            return 0;
        }

        return 1;
    }

    inline int_t addOne(int_t idd)
    {
        return ((idd+1) % cap);
    }

    status_t lockr_()
    {
        return sem_wait(&lockr);
    }

    status_t unlockr_()
    {
        return sem_post(&lockr);
    }

    status_t lockw_()
    {
        return sem_wait(&lockw);
    }

    status_t unlockw_()
    {
        return sem_post(&lockw);
    }

};
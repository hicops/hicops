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
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <unistd.h>
#include "scheduler.h"

extern VOID *DSLIM_IO_Thread_Entry(VOID *argv);
extern VOID *DSLIM_ExtraIO_Thread_Entry(VOID *argv);

using namespace std;

Scheduler::Scheduler()
{
    /* Queues to track threads */
    xtraIO = 2;
    nxtra = 0;

    dump = new lwqueue <THREAD *> (50, false);

    /* Lock for above queues */
    sem_init(&manage, 0, 1);
    sem_init(&dumpQ, 0, 1);

    /* Thresholds */
    maxpenalty = 10;
    stopXtra = 0;
    minrate = 0.45;
    waitSincelast = 0;

    Ftplus1 = 0;   /* Forecast  */
    t = 0;         /* Time instance */
    yt = 0;        /* Current observation */
    ytminus1 = 0;

    /* Intermediate equation variables */
    St = 0;
    bt = 0;
    Stminus1 = 0;
    btminus1 = 0;

    alpha = alpha1 = 0.5;
    gamma = gamma1 = 0.8;

    /* The main IO thread */
    ioThd = 1;
    pthread_create(&ioThread, NULL, &DSLIM_IO_Thread_Entry, NULL);
}

Scheduler::Scheduler(INT maxio, INT dumpsize)
{
    /* Queues to track threads */
    xtraIO = maxio;
    nxtra = 0;
    stopXtra = 0;
    minrate = 0.45;
    dump = new lwqueue <THREAD *> (dumpsize, false);

    /* Lock for above queues */
    sem_init(&manage, 0, 1);
    sem_init(&dumpQ, 0, 1);

    /* Thresholds */
    maxpenalty = 10;
    waitSincelast = 0;

    Ftplus1 = 0;   /* Forecast  */
    t = 0;         /* Time instance */
    yt = 0;        /* Current observation */
    ytminus1 = 0;

    /* Intermediate equation variables */
    St = 0;
    bt = 0;
    Stminus1 = 0;
    btminus1 = 0;

    alpha = alpha1 = 0.5;
    gamma = gamma1 = 0.8;

    /* The main IO thread */
    ioThd = 1;
    pthread_create(&ioThread, NULL, &DSLIM_IO_Thread_Entry, NULL);
}

Scheduler::~Scheduler()
{
    /* Thresholds */
    maxpenalty = 0;
    xtraIO = 0;
    nxtra = 0;
    waitSincelast = 0;

    Ftplus1 = 0;   /* Forecast  */
    t = 0;         /* Time instance */
    yt = 0;        /* Current observation */
    ytminus1 = 0;

    /* Intermediate equation variables */
    St = 0;
    bt = 0;
    Stminus1 = 0;
    btminus1 = 0;

    alpha = alpha1 = 0;
    gamma = gamma1 = 0;

    /* The two main threads */
    waitForThread(&ioThread);

    /* Wait for any active threads */
    while (getNumActivThds() != 0)
    {
        sleep(0.1);
    }

    flushDumpQueue();

    sem_destroy(&manage);
    sem_destroy(&dumpQ);

    delete dump;
}

VOID Scheduler::flushDumpQueue()
{
    sem_wait(&dumpQ);

    while (!dump->isEmpty())
    {
        THREAD * pt = dump->front();
        dump->pop();
        waitForThread(pt);
        delete pt;
    }

    sem_post(&dumpQ);
}

DOUBLE Scheduler::forecastLASP(DOUBLE yt)
{
    /* Increase the time interval value */
    t++;

    /* Check the interval number */
    if (t > 2)
    {
        /* Set the St-1 and bt-1 to St and bt */
        Stminus1 = St;
        btminus1 = bt;

        /* Compute new St and  bt */
        St = alpha * yt + (1- alpha) * (Stminus1 + btminus1);
        bt = gamma * (St - Stminus1) + (1 - gamma) * (btminus1);
    }
    else if (t == 2)
    {
        /* Set the St-1 and bt-1 to St and bt */
        Stminus1 = St;
        btminus1 = (yt - ytminus1); // Using suggestion: b1 = y2 - y1

        /* Compute new St and bt */
        St = alpha * (yt) + (1- alpha) * (Stminus1 + btminus1);
        bt = gamma * (St - Stminus1) + (1 - gamma) * (btminus1);
    }
    else /* t == 1 */
    {
        St = yt;
    }

    /* Compute forecast */
    Ftplus1 = St + bt;

    return Ftplus1;
}

STATUS Scheduler::waitForThread(THREAD *thd)
{
    VOID *ptr = NULL;
    return pthread_join(*thd, &ptr);
}

DOUBLE Scheduler::forecastLASP(DOUBLE yt, DOUBLE deltaS)
{
    /* Increase the time interval value */
    t++;

    /* Check the interval number */
    if (t > 2)
    {
        /* Set the St-1 and bt-1 to St and bt */
        Stminus1 = St;
        btminus1 = bt;

        alpha = alpha/(alpha + pow(1-alpha1, deltaS));
        gamma = gamma/(gamma + pow(1-gamma1, deltaS));

        /* Compute new St and  bt */
        St = alpha * (yt) + (1- alpha) * (Stminus1 + (deltaS * btminus1));
        bt = (gamma * (St - Stminus1) / (deltaS)) + (1 - gamma) * (btminus1);
    }
    else if (t == 2)
    {
        /* Set the St-1 and bt-1 to St and bt */
        Stminus1 = St;
        btminus1 = (yt - ytminus1); // Using suggestion: b1 = y2 - y1

        alpha = alpha/(alpha + pow(1-alpha1, deltaS));
        gamma = gamma/(gamma + pow(1-gamma1, deltaS));

        /* Compute new St and bt */
        St = alpha * (yt) + (1- alpha) * (Stminus1 + (deltaS * btminus1));
        bt = (gamma * (St - Stminus1) / (deltaS)) + (1 - gamma) * (btminus1);
    }
    else /* t == 1 */
    {
        St = yt;
    }

    /* Compute forecast */
    Ftplus1 = St + bt;

    return Ftplus1;
}

BOOL Scheduler::makeDecisions(DOUBLE yt)
{
    BOOL decision = false;

    DOUBLE rate = (yt - ytminus1)/(ytminus1);

    /* Add the next predicted rate of change */
    rate += (Ftplus1 - yt) / (yt);

    /* Divide by 2 for average */
    rate /= 2;

    waitSincelast += yt;

    if (t <= 1)
    {
        if ((waitSincelast + Ftplus1) >= maxpenalty)
        {
            decision = true;

            if (stopXtra)
            {
                stopXtra -=1;
            }
        }
    }
    else
    {
        /* Increasing very fast or too much accumulated */
        /* FIXME: maxpenalty needs to be normalized according to the resources in use */
        if (rate >= minrate && (waitSincelast + Ftplus1) >= maxpenalty)
        {
            decision = true;

            if (stopXtra)
            {
                stopXtra -=1;
            }

        }
        /* Very low rate - Stop a thread */
        else if (rate < minrate)
        {
            decision = false;

            if (nxtra && stopXtra < nxtra)
            {
                stopXtra += 1;
            }

        }
        else
        {
            /* Nominal rate and not too much accumulated,
             * Keep the state and don't stop any threads */
            decision = false;
        }
    }

    /* Set ytminus1 to yt */
    ytminus1 = yt;

    return decision;
}

STATUS Scheduler::dispatchThread()
{
    STATUS status = SLM_SUCCESS;

    if (nxtra < xtraIO)
    {
        nxtra += 1;

        /* Schedule the thread */
        THREAD *ptr = new THREAD;

        if (ptr != NULL)
        {
            /* Pass the reference to thread block as argument */
            status = pthread_create(ptr, NULL, &DSLIM_ExtraIO_Thread_Entry, (VOID *) ptr);
            ptr = NULL;
        }
        else
        {
            status = ERR_BAD_MEM_ALLOC;
            delete ptr;
        }
    }

    return status;
}

STATUS Scheduler::takeControl(VOID *argv)
{
    STATUS status = SLM_SUCCESS;

    THREAD *ptr = (THREAD *) argv;

    sem_wait(&dumpQ);

    if (dump->isFull())
    {
        flushDumpQueue();
    }

    /* Mark the thread block for deletion */
    dump->push(ptr);

    sem_wait(&manage);

    nxtra -= 1;

    if (stopXtra > 0)
    {
        stopXtra -= 1;
    }

    sem_post(&manage);

    sem_post(&dumpQ);

    return status;
}

VOID Scheduler::ioComplete()
{
    sem_wait(&manage);

    ioThd = 0;

    sem_post(&manage);

    return;
}

STATUS Scheduler::runManager(DOUBLE yt)
{
    STATUS status = SLM_SUCCESS;

    /* The new value of signal is yt + yt -1
     * since input yt is only difference from
     * the last value */
    yt = yt + ytminus1;

    /* Use double exponential smoothing forecasting
     * (LASP) to predict future */
    (VOID) this->forecastLASP(yt);

    sem_wait(&manage);

    /* Make decisions */
    if (this->makeDecisions(yt))
    {
        status = dispatchThread();
        waitSincelast = 0;
    }

    sem_post(&manage);

    return status;
}

STATUS Scheduler::runManager(DOUBLE yt, INT qchunk)
{
    STATUS status = SLM_SUCCESS;

    // FIXME: The max qchunk size can be 50000 = 2.5
    DOUBLE dS = qchunk / 25000;

    /* Use double exponential smoothing forecasting
     * (LASP) to predict future */
    (VOID) this->forecastLASP(yt, dS);

    sem_wait(&manage);

    /* Make decisions */
    if (this->makeDecisions(yt))
    {
        status = dispatchThread();
        waitSincelast = 0;
    }

    sem_post(&manage);

    return status;
}

INT Scheduler::getNumActivThds()
{
    sem_wait(&manage);

    INT val = nxtra + ioThd;

    sem_post(&manage);

    return val;
}

BOOL Scheduler::checkDecisions()
{
    sem_wait(&manage);

    BOOL ret = stopXtra;

    sem_post(&manage);

    return ret;
}

VOID Scheduler::waitForCompletion()
{
    while (true)
    {
        sem_wait(&manage);

        /* Number of extra threads running */
        if (nxtra == 0)
        {
            sem_post(&manage);
            break;
        }

        sem_post(&manage);

        /* Sleep for 0.1s */
        sleep(0.1);
    }
}

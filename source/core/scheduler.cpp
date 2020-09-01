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
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <unistd.h>
#include "scheduler.h"
#include "slm_dsts.h"

extern gParams params;
extern VOID DSLIM_IO_Threads_Entry();

using namespace std;

Scheduler::Scheduler()
{
    nIOThds = 0;
    eSignal = false;

    /* Set the total threads for preprocessing */
    maxIOThds = std::max((int_t)1, (int_t)params.maxprepthds);

    /* Lock for above queues */
    sem_init(&manage, 0, 1);

    /* Thresholds */
    maxpenalty = 2;
    stopXtra = false;
    minrate = 0.3;
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

    /* Create an IO thread */
    dispatchThread();
}

Scheduler::Scheduler(int_t maxio)
{
    /* Queues to track threads */
    maxIOThds = maxio;
    nIOThds = 0;
    eSignal = false;
    stopXtra = 0;
    minrate = 0.30;

    /* Lock for above queues */
    sem_init(&manage, 0, 1);

    /* Thresholds */
    maxpenalty = 2;
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
    dispatchThread();
}

Scheduler::~Scheduler()
{
    /* Thresholds */
    maxpenalty = 0;
    maxIOThds = 0;
    nIOThds = 0;
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

    for (auto &itr : thread_pool)
        itr.join();

    thread_pool.clear();

    sem_destroy(&manage);
}

double_t Scheduler::forecastLASP(double_t yt)
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

double_t Scheduler::forecastLASP(double_t yt, double_t deltaS)
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

BOOL Scheduler::makeDecisions(double_t yt, int_t dec)
{
    BOOL decision = false;

    /* Calculate the current rate of change.
     * Add 0.01s in denominator to attenuate the small yt changes */
    double_t rate = (yt - ytminus1)/(ytminus1 + 0.01);

    /* Add the next predicted rate of change */
    rate += (Ftplus1 - yt) / (yt + 0.01);

    /* Divide by 2 for mean */
    rate /= 2;

    waitSincelast += yt;

    /* Special case for time step 1 */
    if (t <= 1)
    {
        if ((waitSincelast + Ftplus1) >= maxpenalty)
        {
            decision = true;
            stopXtra = false;
        }
    }
    else
    {
        if (dec == -1 && nIOThds < 1)
        {
            decision = true;
            stopXtra = false;
        }
        else if (dec == -1 && nIOThds >= 1)
        {
            /* Increasing very fast or too much accumulated */
            /* FIXME: maxpenalty needs to be normalized according to the resources in use */
            if ((waitSincelast >= 6 * maxpenalty) || (rate >= minrate && (waitSincelast + Ftplus1) >= maxpenalty))
            {
                decision = true;
                stopXtra = false;
            }
        }
        else if (dec > -1)
        {
            /* Very low rate - Stop a thread */
            if (rate < minrate)
            {
                decision = false;

                /* Stop a thread if more than 1 threads running,
                 * else the one thread will stop naturally
                 */
                if (nIOThds > 1)
                {
                    stopXtra = true;
                }
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

status_t Scheduler::dispatchThread()
{
    if (nIOThds < maxIOThds)
    {
        nIOThds += 1;

        /* Pass the reference to thread block as argument */
        thread_pool.push_back(std::move(std::thread(DSLIM_IO_Threads_Entry)));
    }

    return SLM_SUCCESS;
}

status_t Scheduler::takeControl()
{
    sem_wait(&manage);

    nIOThds -= 1;

    stopXtra = false;

    sem_post(&manage);

    return SLM_SUCCESS;
}

VOID Scheduler::ioComplete()
{
    /* Raise the eSignal if this is
     * the last thread running */
    sem_wait(&manage);

    if (nIOThds == 1)
    {
        eSignal = true;
    }

    sem_post(&manage);

    return;
}

BOOL Scheduler::checkSignal()
{
    sem_wait(&manage);

    BOOL val = eSignal;

    sem_post(&manage);

    return val;
}

status_t Scheduler::runManager(double_t yt, int_t dec)
{
    status_t status = SLM_SUCCESS;

    /* The new value of signal is yt + yt -1
     * since input yt is only difference from
     * the last value */
    yt = yt + ytminus1;

    /* Use double exponential smoothing forecasting
     * (LASP) to predict future */
    (VOID) this->forecastLASP(yt);

    sem_wait(&manage);

    /* Make decisions */
    if (!eSignal && this->makeDecisions(yt, dec))
    {
        status = dispatchThread();
        waitSincelast = 0;
    }

    sem_post(&manage);

    return status;
}

int_t Scheduler::getNumActivThds()
{
    sem_wait(&manage);

    int_t val = nIOThds;

    sem_post(&manage);

    return val;
}

BOOL Scheduler::checkPreempt()
{
    sem_wait(&manage);

    BOOL ret = stopXtra;

    stopXtra = false;

    sem_post(&manage);

    return ret;
}

VOID Scheduler::waitForCompletion()
{
    while (true)
    {
        sem_wait(&manage);

        /* Number of IO threads running */
        if (nIOThds == 0)
        {
            sem_post(&manage);
            break;
        }

        sem_post(&manage);

        /* Sleep for 0.1s */
        sleep(0.1);
    }
}

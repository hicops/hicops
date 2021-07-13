/*
 * Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is licensed under the
 * Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
 * See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
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

    // Create at most IO threads
    auto ts = std::min(maxIOThds, 2);
    for (auto t = 0; t < ts; t++)
        dispatchThread();
}

Scheduler::Scheduler(int_t maxio)
{
    /* Queues to track threads */
    maxIOThds = maxio;
    nIOThds = 0;
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

    // Create at most IO threads
    auto ts = std::min(maxIOThds, 2);
    for (auto t = 0; t < ts; t++)
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
        if (dec == -1)
        {
            if (nIOThds < 1)
            {
                decision = true;
                stopXtra = false;
            }
            else if (yt >= maxpenalty/20 && waitSincelast >= 2 * maxpenalty)
            {
                decision = true;
                stopXtra = false;
            }
        }
        else if (dec == 1)
        {
            if (nIOThds > 1)
                stopXtra = true;

            decision = false;
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
    if (this->makeDecisions(yt, dec))
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

    BOOL ret = (nIOThds > 1 && stopXtra);

    if (ret)
    {
        stopXtra = false;
        nIOThds -= 1;
    }

    sem_post(&manage);

    return ret;
}

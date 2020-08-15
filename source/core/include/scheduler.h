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

#ifndef SCHEDULER_H_
#define SCHEDULER_H_

#include "common.h"
#include "lwqueue.h"

class Scheduler
{
private:

    /* Number of IO threads */
    int_t nIOThds;
    int_t maxIOThds;

    BOOL eSignal;

    /* Queues to track threads */
    lwqueue <thread_t *> *dump;

    /* Lock for above queues */
    lock_t manage;
    lock_t dumpQ;

    /* Thresholds */
    BOOL stopXtra;
    double_t maxpenalty;
    double_t minrate;
    double_t waitSincelast;

    /* Variables for forecasting penalties  */
    double_t Ftplus1;   /* Forecast  */
    double_t t;         /* Time instance */
    double_t yt;        /* Current observation */
    double_t ytminus1;  /* Last Observation */

    /* Intermediate equation variables */
    double_t St;
    double_t bt;
    double_t Stminus1;
    double_t btminus1;

    /* LASP hyper-parameters */
    double_t alpha;
    double_t alpha1;
    double_t gamma;
    double_t gamma1;


    /* Private Functions */
    status_t waitForThread(thread_t *);
    VOID   flushDumpQueue();
    double_t forecastLASP(double_t yt);
    double_t forecastLASP(double_t yt, double_t deltaS);
    status_t dispatchThread();
    BOOL   makeDecisions(double_t yt, int_t decisions);

public:
    Scheduler();
    Scheduler(int_t, int_t);
    virtual ~Scheduler();

    VOID   ioComplete();
    int_t    getNumActivThds();
    BOOL   checkPreempt();
    status_t takeControl(VOID *argv);
    status_t runManager(double_t yt, int_t dec);
    VOID   waitForCompletion();
    BOOL   checkSignal();
};

#endif /* SCHEDULER_H_ */

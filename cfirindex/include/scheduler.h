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

#ifndef SCHEDULER_H_
#define SCHEDULER_H_

#include "common.h"
#include "lwqueue.h"

class Scheduler
{
private:
    /* The two main threads */
    THREAD ioThread;
    INT ioThd;

    /* Queues to track threads */
    lwqueue <THREAD *> *dump;

    /* Lock for above queues */
    LOCK manage;
    LOCK dumpQ;

    /* Thresholds */
    INT xtraIO;
    INT nxtra;
    INT stopXtra;

    DOUBLE maxpenalty;
    DOUBLE minrate;

    DOUBLE waitSincelast;

    /* Variables for forecasting penalties  */
    DOUBLE Ftplus1;   /* Forecast  */
    DOUBLE t;         /* Time instance */
    DOUBLE yt;        /* Current observation */
    DOUBLE ytminus1;  /* Last Observation */

    /* Intermediate equation variables */
    DOUBLE St;
    DOUBLE bt;
    DOUBLE Stminus1;
    DOUBLE btminus1;

    /* LASP hyper-parameters */
    DOUBLE alpha;
    DOUBLE alpha1;
    DOUBLE gamma;
    DOUBLE gamma1;


    /* Private Functions */
    STATUS waitForThread(THREAD *);
    VOID   flushDumpQueue();
    DOUBLE forecastLASP(DOUBLE yt);
    DOUBLE forecastLASP(DOUBLE yt, DOUBLE deltaS);
    STATUS dispatchThread();
    BOOL   makeDecisions(DOUBLE yt);

public:
    Scheduler();
    Scheduler(INT, INT);
    virtual ~Scheduler();

    VOID ioComplete();
    INT    getNumActivThds();
    BOOL   checkDecisions();
    STATUS takeControl(VOID *argv);
    STATUS runManager(DOUBLE yt);
    STATUS runManager(DOUBLE yt, INT qchunk);
    VOID   waitForCompletion();
};

#endif /* SCHEDULER_H_ */

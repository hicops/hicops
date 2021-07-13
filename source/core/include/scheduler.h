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
#include <vector>
#include <thread>

class Scheduler
{
private:

    /* Number of IO threads */
    int_t nIOThds;
    int_t maxIOThds;

    std::vector<std::thread> thread_pool;

    /* Lock for above queues */
    lock_t manage;

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
    double_t forecastLASP(double_t yt);
    double_t forecastLASP(double_t yt, double_t deltaS);
    BOOL   makeDecisions(double_t yt, int_t decisions);

public:
    Scheduler();
    Scheduler(int_t);
    virtual ~Scheduler();

    status_t dispatchThread();
    int_t    getNumActivThds();
    BOOL   checkPreempt();
    status_t takeControl();
    status_t runManager(double_t yt, int_t dec);
    VOID   waitForCompletion();
};
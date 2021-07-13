/*
 * Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is licensed under the
 * Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
 * See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
 */

#pragma once

#include <vector>
#include <valarray>
#include <algorithm>
#include "lwvector.h"
#include "config.hpp"
#include "common.hpp"
#include "slm_dsts.h"

using dvector = std::vector<double_t>;
using darray  = std::valarray<double_t>;

template <class T>
BOOL isNegative(T i) { return i < 0.0; }

template <class T>
BOOL isZeroNegative(T i) { return i <= 0.0; }

template <class T>
BOOL isNotNegative(T i) { return i > 0.0; }

template <class T>
BOOL isLargerThan1(T i) { return i > 1.0; }

template <class T>
VOID printVecArray(T v)
{
    std::cout << "\nPRINT:" << std::endl;

    for (auto i : v)
    {
        std::cout << i << ", ";
    }

    std::cout << std::endl;
}

class expeRT
{
private:

    /* Size of histogram */
    const int_t SIZE = 2 + (MAX_HYPERSCORE * 10);

    /* Index number set 1 */
    int_t stt1;
    int_t end1;

    /* Index number set 2 */
    int_t stt;
    int_t ends;

    int_t hyp;

    int_t vaa;

    double_t mu_t;
    double_t beta_t;

    double_t* yy = NULL;

    /* p_x will contain the probability density function */
    lwvector<double_t> *p_x = NULL;
    lwvector<double_t> *sx = NULL;
    lwvector<double_t> *X = NULL;

    /* pdata will contain the partial logWeibull data */
    lwvector<double_t> *pdata = NULL;
    int_t pN;

    /* Constructs a log-Weibull distribution */
    VOID logWeibullResponse(double_t, double_t, int_t ,int_t);
    darray alogWeibullResponse(double_t, double_t, int_t, int_t);

    /* Learning Function - Optimizing */
    double_t logWeibullFit(lwvector<double_t> *, int_t, int_t, int_t niter=6000, double_t lr=0.12, double_t cutoff=1e-3);

    inline double_t MeanSqError(const darray &);

    template <class T>
    VOID LinearFit(T& x, T& y, int_t n, double_t &a, double_t &b);

    template <class T>
    inline int_t argmax(T &data, int_t i1, int_t i2, double_t value);

    template <class T>
    inline int_t rargmax(T &data, int_t i1, int_t i2, double_t value);

    template <class T>
    inline int_t largmax(T &data, int_t i1, int_t i2, double_t value);

    dvector vrange(int_t, int_t);
    darray  arange(int_t, int_t);

public:

    /* Constructor */
    expeRT();

    /* Destructor */
    virtual ~expeRT();

    /* Function to reset the data */
    VOID ResetPartialVectors();

    status_t Reconstruct(ebuffer *ebs, int_t specno, partRes *fR);

    /* Add distibution data */
    status_t AddlogWeibull(int_t, double_t, double_t, int_t, int_t);

    status_t StoreIResults(Results *, int_t, ebuffer *);

    /* Model using log-Weibull in DISTMEM */
    status_t ModelSurvivalFunction(double_t &, const int_t);

    /* Model using log-Weibull in SHM */
    status_t ModelSurvivalFunction(Results *);

    /* Model using log-Weibull in DISTMEM */
    status_t ModelTailFit(double_t &, const int_t);

    /* Model using log-Weibull in SHM */
    status_t ModelTailFit(Results *);


    /* Model the partial distribution using logWeibull */
    status_t Model_logWeibull(Results *);

#if 0
    /* In case of distributed memory, we will call this */
    status_t ModelSurvivalFunction(double_t &, int_t);
#endif
};
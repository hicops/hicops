/*
 * Copyright (C) 2020  Muhammad Haseeb, and Fahad Saeed
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

#ifndef EXPERT_H_
#define EXPERT_H_

#include <vector>
#include <valarray>
#include <algorithm>
#include "lwvector.h"
#include "config.h"
#include "common.h"
#include "slm_dsts.h"

using dvector = std::vector<DOUBLE>;
using darray  = std::valarray<DOUBLE>;

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
    const INT SIZE = 2 + (MAX_HYPERSCORE * 10);

    /* Index number set 1 */
    INT stt1;
    INT end1;

    /* Index number set 2 */
    INT stt;
    INT ends;

    INT hyp;

    INT vaa;

    DOUBLE mu_t;
    DOUBLE beta_t;

    DOUBLE* yy = NULL;

    /* p_x will contain the probability density function */
    lwvector<DOUBLE> *p_x = NULL;
    lwvector<DOUBLE> *sx = NULL;
    lwvector<DOUBLE> *X = NULL;

    /* pdata will contain the partial logWeibull data */
    lwvector<DOUBLE> *pdata = NULL;
    INT pN;

    /* Constructs a log-Weibull distribution */
    VOID logWeibullResponse(DOUBLE, DOUBLE, INT ,INT);
    darray alogWeibullResponse(DOUBLE, DOUBLE, INT, INT);

    /* Learning Function - Optimizing */
    DOUBLE logWeibullFit(lwvector<DOUBLE> *, INT, INT, INT niter=6000, DOUBLE lr=0.12, DOUBLE cutoff=1e-3);

    inline DOUBLE MeanSqError(const darray &);

    template <class T>
    VOID LinearFit(T& x, T& y, INT n, DOUBLE &a, DOUBLE &b);

    template <class T>
    inline INT argmax(T &data, INT i1, INT i2, DOUBLE value);

    template <class T>
    inline INT rargmax(T &data, INT i1, INT i2, DOUBLE value);

    template <class T>
    inline INT largmax(T &data, INT i1, INT i2, DOUBLE value);

    dvector vrange(INT, INT);
    darray  arange(INT, INT);

public:

    /* Constructor */
    expeRT();

    /* Destructor */
    virtual ~expeRT();

    /* Function to reset the data */
    VOID ResetPartialVectors();

    STATUS Reconstruct(ebuffer *ebs, INT specno, partRes *fR);

    /* Add distibution data */
    STATUS AddlogWeibull(INT, DOUBLE, DOUBLE, INT, INT);

    STATUS StoreIResults(Results *, INT, ebuffer *);

    /* Model using log-Weibull in DISTMEM */
    STATUS ModelSurvivalFunction(DOUBLE &, const INT);

    /* Model using log-Weibull in SHM */
    STATUS ModelSurvivalFunction(Results *);

    /* Model using log-Weibull in DISTMEM */
    STATUS ModelTailFit(DOUBLE &, const INT);

    /* Model using log-Weibull in SHM */
    STATUS ModelTailFit(Results *);


    /* Model the partial distribution using logWeibull */
    STATUS Model_logWeibull(Results *);

#if 0
    /* In case of distributed memory, we will call this */
    STATUS ModelSurvivalFunction(DOUBLE &, INT);
#endif
};

#endif /* EXPERT_H_ */

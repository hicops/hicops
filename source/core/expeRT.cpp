/*
 * Copyright (C) 2020  Muhammad Haseeb, and Fahad Saeed
 * Florida International University (FIU), Miami, FL
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
#include <vector>
#include <fstream>
#include <numeric>
#include "expeRT.h"
#include "sgsmooth.h"

using namespace std;

expeRT::expeRT()
{
    stt1 = stt = 0;
    end1 = ends = SIZE - 1;

    yy  = NULL;

    p_x = new lwvector<double_t>(SIZE);
    sx = new lwvector<double_t>(SIZE);
    X = new lwvector<double_t>(SIZE);

    /* Remember the constructor for valarrays */
    pdata = new lwvector<double_t>(SIZE, 0.0);
    pN = 0;

    mu_t = 0.0;
    beta_t = 4.0;

    hyp = vaa = 0;
}

expeRT::~expeRT()
{
    /* Delete all vectors */
    yy  = NULL;
    stt1= stt = 0;
    end1 = ends = SIZE - 1;

    /* Clear all vectors */
    if (p_x != NULL)
    {
        delete p_x;
        p_x = NULL;
    }

    if (pdata != NULL)
    {
        delete pdata;
        pdata = NULL;
    }

    if (sx != NULL)
    {
        delete sx;
        sx = NULL;
    }

    if (X != NULL)
    {
        delete X;
        X = NULL;
    }

    pN = 0;

    mu_t = beta_t = 0.0;

    hyp = vaa = 0;
}

inline dvector expeRT::vrange(int_t stt, int_t end)
{
    dvector xx (end-stt+1);

    int_t k = stt;

    for (auto &r : xx)
    {
        r = k;
        k++;
    }

    return xx;
}

inline darray expeRT::arange(int_t stt, int_t end)
{
    darray xx (end-stt+1);

    int_t k = stt;

    for (auto &r : xx)
    {
        r = k;
        k++;
    }

    return xx;
}

status_t expeRT::ModelSurvivalFunction(Results *rPtr)
{
    status_t status = SLM_SUCCESS;

    /* Assign to internal variables */
    yy = rPtr->survival;
    hyp = rPtr->maxhypscore;

    /*if (yy == NULL)
    {
        status = ERR_INVLD_PARAM;
        mu_t = 0;
        beta_t = 100;
    }*/

    if (status == SLM_SUCCESS)
    {
        /* Find the curve region */
        end1 = rargmax<double_t *>(yy, 0, hyp-1, 1.0);
        stt1 = argmax<double_t *>(yy, 0, end1, 1.0);

        /* To handle special cases */
        if (stt1 == end1)
        {
            stt1 = end1;
            end1 += 1;
        }

        /* Slice off yyt between stt1 and end1 */
        p_x->Assign(yy + stt1, yy + end1 + 1);

        lwvector<double_t> *yyt = p_x;

        lwvector<double_t> *yhat = NULL;

        /* Database size
         * vaa = accumulate(yy, yy + hyp + 1, 0); */
        vaa = rPtr->cpsms;

        /* Check if no distribution data except for hyp */
        if (vaa < 1)
        {
            mu_t = 0;
            beta_t = 100;
            stt = stt1;
            ends = end1;
            status = ERR_NOT_ENOUGH_DATA;
        }

        /* Enough distribution data - Proceed */
        if (status == SLM_SUCCESS)
        {
            /* Smoothen the curve using Savitzky-Golay filter */
            int_t svgl = std::min(7, end1 - stt1);
            int_t pln = 5;

            /* mu estimation markers */
            auto l = std::max_element(yy + stt1, (yy + end1 + 1)) -
                    (yy + stt1);
            auto k = l;

            /* Window length = 7 or max odd number */
            if (svgl % 2 == 0)
            {
                svgl -= 1;
            }

            /* If window size > 1 */
            if (svgl > 1)
            {
                /* Polynomial = 5 or max */
                pln = std::min(5, svgl - 1);

                if ((int_t) yyt->Size() >= (svgl-1+2))
                {
                    yhat = new lwvector<double_t>(yyt->Size(), 0.0);

                    /* Smoothen the curve using SavGol: window length=7, provide win=(window-1)/2, polynomial=5 */
                    sg_smooth(yyt, yhat, std::max(1, (svgl-1)/2), pln);

                    /* Adjust for negatives */
                    std::replace_if(yhat->begin(), yhat->end(), isNegative<double_t>, 0);

                    /* Normalize yhat */
                    yhat->divide((double_t)std::max(std::accumulate(yhat->begin(), yhat->end(), 1), vaa));

                    /* Normalize yy */
                    yyt->divide((double_t)vaa);

                    /* Update the marker */
                    k = std::max_element(yhat->begin(),yhat->end()) -
                                         yhat->begin();

                    /* Mix yhat (35%) + yyt (65%) */
                    for (auto id = 0; id < (int) yyt->Size(); id++)
                    {
                        (*yyt)[id] = (*yhat)[id] * 0.35 + (*yyt)[id] * 0.65;
                    }

                    delete yhat;
                    yhat = NULL;
                }
                else
                {
                    /* Update the marker */
                    k = l;

                    /* Normalize yy */
                    yyt->divide((double_t)vaa);
                }
            }
            else
            {
                /* Normalize yhat */
                const int_t vaa2 = std::max(accumulate(yyt->begin(), yyt->end(), 1), vaa);
                yyt->divide((double_t)vaa2);
            }

            /* Training parameters */
            mu_t = (stt1 + (k+l)/2.0);

            /* Train the logWeibull model */
            (VOID) logWeibullFit(yyt, stt1, end1);

            /* Modeled response * vaa (included) */
            logWeibullResponse(mu_t, beta_t, 0, hyp);

            /* Disabled this code since we don't need to tail fit again */
#if 0
            /* Filter p_x again */
            ends = rargmax((*p_x), 0, hyp - 1, 0.99);
            stt = argmax((*p_x), 0, ends, 0.99);

            p_x->clip(stt, ends);

            /* Compute survival function s(x) */
            sx->Assign(p_x->begin(), p_x->end());

            /* cumulative_sum(sx) */
            std::partial_sum(p_x->begin(), p_x->end(), sx->begin());

            /* Survival function s(x) */
            sx->divide((double_t)vaa);
            sx->add((double_t)-1);
            sx->negative();

            /* Adjust for > 1 */
            std::replace_if(sx->begin(), sx->end(), isLargerThan1<double_t>, 0.999);

            /* Adjust for negatives */
            int_t replacement = rargmax(*sx, 0, sx->Size()-1, 1e-4);
            std::replace_if(sx->begin(), sx->end(), isZeroNegative<double_t>, (*sx)[replacement]);

            /* log10(s(x)) */
            std::transform(sx->begin(), sx->end(), sx->begin(), [](double_t& c) {   return log10(c);});

            /* Offset markers */
            auto mark = 0;
            auto mark2 = 0;
            auto hgt = (*sx)[sx->Size() - 1] - (*sx)[0];

            /* If length > 4, then find thresholds */
            if (sx->Size() > 3)
            {
                mark = largmax<lwvector<double_t>>(*sx, 0, sx->Size()-1, (*sx)[0] + hgt * 0.22) - 1;
                mark2 = rargmax<lwvector<double_t>>(*sx, 0, sx->Size()-1, (*sx)[0] + hgt*0.87);

                if (mark2 == (int_t)sx->Size())
                {
                    mark2 -= 1;
                }

                /* To handle special cases */
                if (mark >= mark2)
                {
                    mark = mark2 - 1;
                }
            }
            /* If length < 4 business as usual */
            else if (sx->Size() == 3)
            {
                /* Mark the start of the regression point */
                mark = largmax(*sx, 0, sx->Size()-1, (*sx)[0] + hgt * 0.22) - 1;
                mark2 = sx->Size() - 1;

                /* To handle special cases */
                if (mark >= mark2)
                {
                    mark = mark2 - 1;
                }
            }
            else
            {
                mark = 0;
                mark2 = sx->Size() - 1;
            }

            /* Make the x-axis */
            X->AddRange(stt+mark, stt+mark2);

            /* Make the y-axis */
            sx->clip(mark, mark2);

            LinearFit<lwvector<double_t>>(*X, *sx, sx->Size(), mu_t, beta_t);

            sx->Erase();
            X->Erase();

            //cout << "y = " << mu_t << "x + " << beta_t << endl;
            //cout << "eValue: " << pow(10, hyp * mu_t + beta_t) * vaa << endl;
#endif /* 0 */
        }
    }

    /* Disabled this code since we don't need to tail fit again */
#if 0
    rPtr->mu = mu_t * 1e6;
    rPtr->beta = beta_t *1e6;
#else
    /* Assign variables back to rPtr */
    rPtr->mu = (*p_x)[hyp] * 1e6;
#endif

    rPtr->minhypscore = stt1;
    rPtr->nexthypscore = end1;

    /* Clear arrays, vectors and variables */
    p_x->Erase();
    mu_t = 0.0;
    beta_t = 4.0;
    stt1 = stt = 0;
    end1 = ends = SIZE - 1;
    yy  = NULL;
    hyp = vaa = 0;

    return status;
}

status_t expeRT::ModelSurvivalFunction(double_t &eValue, const int_t max1)
{
    status_t status = SLM_SUCCESS;

    /* Get hyp from the packet */
    hyp = max1;

    /* Find the curve region */
    end1 = rargmax((*pdata), 0, hyp - 1, 1.0);
    stt1 = argmax((*pdata), 0, end1, 1.0);

    /* Slice off yyt between stt1 and end1 */
    p_x->Assign(pdata->begin() + stt1, pdata->begin() + end1 + 1);

    lwvector<double_t> *yyt = p_x;
    lwvector<double_t> *yhat = NULL;

    /* Database size
     * vaa = accumulate(yy, yy + hyp + 1, 0); */
    vaa = pN;

    /* Check if no distribution data except for hyp */
    if (vaa < 1)
    {
        mu_t = 0;
        beta_t = 100;
        stt = stt1;
        ends = end1;
        status = ERR_NOT_ENOUGH_DATA;
    }

    /* Enough distribution data - Proceed */
    if (status == SLM_SUCCESS)
    {
        /* Smoothen the curve using Savitzky-Golay filter */
        int_t svgl = std::min(7, end1 - stt1);
        int_t pln = 5;

        /* mu estimation markers */
        auto l = std::max_element(yyt->begin(), yyt->end()) - yyt->begin();
        auto k = l;

        /* Window length = 7 or max odd number */
        if (svgl % 2 == 0)
        {
            svgl -= 1;
        }

        /* If window size > 1 */
        if (svgl > 1)
        {
            /* Polynomial = 5 or max */
            pln = std::min(5, svgl - 1);

            if ((int_t) yyt->Size() >= (svgl - 1 + 2))
            {
                yhat = new lwvector<double_t>(yyt->Size(), 0.0);

                /* Smoothen the curve using SavGol: window length=7, provide win=(window-1)/2, polynomial=5 */
                sg_smooth(yyt, yhat, std::max(1, (svgl - 1) / 2), pln);

                /* Aor negatives */
                std::replace_if(yhat->begin(), yhat->end(), isNegative<double_t>, 0);

                /* Normalize yhat */
                yhat->divide((double_t) std::max(accumulate(yhat->begin(), yhat->end(), 1), vaa));

                /* Normalize yy */
                yyt->divide((double_t) vaa);

                /* Update the marker */
                k = std::max_element(yhat->begin(), yhat->end()) - yhat->begin();

                /* Mix yhat (35%) + yyt (65%) */
                for (auto id = 0; id < (int) yyt->Size(); id++)
                {
                    (*yyt)[id] = (*yhat)[id] * 0.35 + (*yyt)[id] * 0.65;
                }

                delete yhat;
                yhat = NULL;
            }
            else
            {
                /* Update the marker */
                k = l;

                /* Normalize yy */
                yyt->divide((double_t) vaa);
            }
        }
        else
        {
            /* Normalize yhat */
            const int_t vaa2 = std::max(accumulate(yyt->begin(), yyt->end(), 1), vaa);
            yyt->divide((double_t) vaa2);
        }

        /* Training parameters */
        mu_t = (stt1 + (k + l) / 2.0);

        /* Train the logWeibull model */
        (VOID) logWeibullFit(yyt, stt1, end1);

        /* Modeled response * vaa (included) */
        logWeibullResponse(mu_t, beta_t, 0, hyp);

        /* Disabled this code since we don't need to tail fit again */
#if 0
        /* Filter p_x again */
        ends = rargmax((*p_x), 0, hyp - 1, 0.99);
        stt = argmax((*p_x), 0, ends, 0.99);

        p_x->clip(stt, ends);

        /* Compute survival function s(x) */
        sx->Assign(p_x->begin(), p_x->end());

        /* cumulative_sum(sx) */
        std::partial_sum(p_x->begin(), p_x->end(), sx->begin());

        /* Survival function s(x) */
        sx->divide((double_t) vaa);
        sx->add((double_t) -1);
        sx->negative();

        /* Adjust for > 1 */
        std::replace_if(sx->begin(), sx->end(), isLargerThan1<double_t>, 0.999);

        /* Adjust for negatives */
        int_t replacement = rargmax(*sx, 0, sx->Size() - 1, 1e-4);
        std::replace_if(sx->begin(), sx->end(), isZeroNegative<double_t>, (*sx)[replacement]);

        /* log10(s(x)) */
        std::transform(sx->begin(), sx->end(), sx->begin(), [](double_t& c)
        {   return log10(c);});

        /* Offset markers */
        auto mark = 0;
        auto mark2 = 0;
        auto hgt = (*sx)[sx->Size() - 1] - (*sx)[0];

        /* If length > 4, then find thresholds */
        if (sx->Size() > 3)
        {
            mark = largmax<lwvector<double_t>>(*sx, 0, sx->Size() - 1, (*sx)[0] + hgt * 0.22) - 1;
            mark2 = rargmax<lwvector<double_t>>(*sx, 0, sx->Size() - 1, (*sx)[0] + hgt * 0.87);

            if (mark2 == (int_t) sx->Size())
            {
                mark2 -= 1;
            }

            /* To handle special cases */
            if (mark >= mark2)
            {
                mark = mark2 - 1;
            }
        }
        /* If length < 4 business as usual */
        else if (sx->Size() == 3)
        {
            /* Mark the start of the regression point */
            mark = largmax(*sx, 0, sx->Size() - 1, (*sx)[0] + hgt * 0.22) - 1;
            mark2 = sx->Size() - 1;

            /* To handle special cases */
            if (mark >= mark2)
            {
                mark = mark2 - 1;
            }
        }
        else
        {
            mark = 0;
            mark2 = sx->Size() - 1;
        }

        /* Make the x-axis */
        X->AddRange(stt + mark, stt + mark2);

        /* Make the y-axis */
        sx->clip(mark, mark2);

        LinearFit<lwvector<double_t>>(*X, *sx, sx->Size(), mu_t, beta_t);

        sx->Erase();
        X->Erase();

        //cout << "y = " << mu_t << "x + " << beta_t << endl;
        //cout << "eValue: " << pow(10, hyp * mu_t + beta_t) * vaa << endl;
#endif /* 0 */
    }

    /* Disabled this code since we don't need to tail fit again */
#if 0
    /* Compute the eValue */
    eValue = pow(10, hyp * mu_t + beta_t) * vaa;
#else
    /* Compute the eValue */
    eValue = (*p_x)[hyp];
#endif /* 0 */

    /* Clear arrays, vectors and variables */
    pdata->setmem(0);
    pN = 0;
    p_x->Erase();
    mu_t = 0.0;
    beta_t = 4.0;
    stt1 = stt = 0;
    end1 = ends = SIZE - 1;
    yy  = NULL;
    hyp = vaa = 0;

    return status;
}

status_t expeRT::ModelTailFit(Results *rPtr)
{
    status_t status = SLM_SUCCESS;

    /* Assign to internal variables */
    yy = rPtr->survival;
    hyp = rPtr->maxhypscore;

    /*if (yy == NULL)
    {
        status = ERR_INVLD_PARAM;
        mu_t = 0;
        beta_t = 100;
    }*/

    if (status == SLM_SUCCESS)
    {
        /* Find the curve region */
        end1 = rargmax<double_t *>(yy, 0, hyp-1, 1.0);
        stt1 = argmax<double_t *>(yy, 0, end1, 1.0);

        /* To handle special cases */
        if (stt1 == end1)
        {
            stt1 = end1;
            end1 += 1;
        }

        /* Slice off yyt between stt1 and end1 */
        p_x->Assign(yy + stt1, yy + end1 + 1);

        /* Database size
         * vaa = accumulate(yy, yy + hyp + 1, 0); */
        vaa = rPtr->cpsms;

        /* Check if no distribution data except for hyp */
        if (vaa < 1)
        {
            mu_t = 0;
            beta_t = 100;
            stt = stt1;
            ends = end1;
            status = ERR_NOT_ENOUGH_DATA;
        }

        /* Enough distribution data - Proceed */
        if (status == SLM_SUCCESS)
        {
            /* Filter p_x again */
            ends = end1;
            stt = stt1;

            /* Compute survival function s(x) */
            sx->Assign(p_x->begin(), p_x->end());

            /* cumulative_sum(sx) */
            std::partial_sum(p_x->begin(), p_x->end(), sx->begin());

            /* Survival function s(x) */
            sx->divide((double_t)vaa);
            sx->add((double_t)-1);
            sx->negative();

            /* Adjust for > 1 */
            std::replace_if(sx->begin(), sx->end(), isLargerThan1<double_t>, 0.999);

            /* Adjust for negatives */
            int_t replacement = rargmax(*sx, 0, sx->Size()-1, 1e-4);
            std::replace_if(sx->begin(), sx->end(), isZeroNegative<double_t>, (*sx)[replacement]);

            /* log10(s(x)) */
            std::transform(sx->begin(), sx->end(), sx->begin(), [](double_t& c) {   return log10(c);});

            /* Offset markers */
            auto mark = 0;
            auto mark2 = 0;
            auto hgt = (*sx)[sx->Size() - 1] - (*sx)[0];

            /* If length > 4, then find thresholds */
            if (sx->Size() > 3)
            {
                mark = largmax<lwvector<double_t>>(*sx, 0, sx->Size()-1, (*sx)[0] + hgt * 0.22) - 1;
                mark2 = rargmax<lwvector<double_t>>(*sx, 0, sx->Size()-1, (*sx)[0] + hgt*0.87);

                if (mark2 == (int_t)sx->Size())
                {
                    mark2 -= 1;
                }

                /* To handle special cases */
                if (mark >= mark2)
                {
                    mark = mark2 - 1;
                }
            }
            /* If length < 4 business as usual */
            else if (sx->Size() == 3)
            {
                /* Mark the start of the regression point */
                mark = largmax(*sx, 0, sx->Size()-1, (*sx)[0] + hgt * 0.22) - 1;
                mark2 = sx->Size() - 1;

                /* To handle special cases */
                if (mark >= mark2)
                {
                    mark = mark2 - 1;
                }
            }
            else
            {
                mark = 0;
                mark2 = sx->Size() - 1;
            }

            /* Make the x-axis */
            X->AddRange(stt+mark, stt+mark2);

            /* Make the y-axis */
            sx->clip(mark, mark2);

            LinearFit<lwvector<double_t>>(*X, *sx, sx->Size(), mu_t, beta_t);

            sx->Erase();
            X->Erase();

            //cout << "y = " << mu_t << "x + " << beta_t << endl;
            //cout << "eValue: " << pow(10, hyp * mu_t + beta_t) * vaa << endl;
        }
    }

    /* Assign variables back to rPtr */
    rPtr->mu = mu_t * 1e6;
    rPtr->beta = beta_t *1e6;
    rPtr->minhypscore = stt1;
    rPtr->nexthypscore = end1;

    /* Clear arrays, vectors and variables */
    p_x->Erase();
    mu_t = 0.0;
    beta_t = 4.0;
    stt1 = stt = 0;
    end1 = ends = SIZE - 1;
    yy  = NULL;
    hyp = vaa = 0;

    return status;
}

status_t expeRT::ModelTailFit(double_t &eValue, const int_t max1)
{
    status_t status = SLM_SUCCESS;

    /* Get hyp from the packet */
    hyp = max1;

    /* Find the curve region */
    end1 = rargmax((*pdata), 0, hyp - 1, 1.0);
    stt1 = argmax((*pdata), 0, end1, 1.0);

    /* Slice off yyt between stt1 and end1 */
    p_x->Assign(pdata->begin() + stt1, pdata->begin() + end1 + 1);

    /* Database size
     * vaa = accumulate(yy, yy + hyp + 1, 0); */
    vaa = pN;

    /* Check if no distribution data except for hyp */
    if (vaa < 1)
    {
        mu_t = 0;
        beta_t = 100;
        stt = stt1;
        ends = end1;
        status = ERR_NOT_ENOUGH_DATA;
    }

    /* Enough distribution data - Proceed */
    if (status == SLM_SUCCESS)
    {
        /* Normalize yhat */
        const int_t vaa2 = std::max(accumulate(p_x->begin(), p_x->end(), 1), vaa);

        p_x->divide((double_t) vaa2);

        /* Filter p_x again */
        ends = end1;
        stt = stt1;

        /* Compute survival function s(x) */
        sx->Assign(p_x->begin(), p_x->end());

        /* cumulative_sum(sx) */
        std::partial_sum(p_x->begin(), p_x->end(), sx->begin());

        /* Survival function s(x) */
        sx->add((double_t) -1);
        sx->negative();

        /* Adjust for > 1 */
        std::replace_if(sx->begin(), sx->end(), isLargerThan1<double_t>, 0.999);

        /* Adjust for negatives */
        int_t replacement = rargmax(*sx, 0, sx->Size() - 1, 1e-4);
        std::replace_if(sx->begin(), sx->end(), isZeroNegative<double_t>, (*sx)[replacement]);

        /* log10(s(x)) */
        std::transform(sx->begin(), sx->end(), sx->begin(), [](double_t& c)
        {   return log10(c);});

        /* Offset markers */
        auto mark = 0;
        auto mark2 = 0;
        auto hgt = (*sx)[sx->Size() - 1] - (*sx)[0];

        /* If length > 4, then find thresholds */
        if (sx->Size() > 3)
        {
            mark = largmax<lwvector<double_t>>(*sx, 0, sx->Size() - 1, (*sx)[0] + hgt * 0.22) - 1;
            mark2 = rargmax<lwvector<double_t>>(*sx, 0, sx->Size() - 1, (*sx)[0] + hgt * 0.87);

            if (mark2 == (int_t) sx->Size())
            {
                mark2 -= 1;
            }

            /* To handle special cases */
            if (mark >= mark2)
            {
                mark = mark2 - 1;
            }
        }
        /* If length < 4 business as usual */
        else if (sx->Size() == 3)
        {
            /* Mark the start of the regression point */
            mark = largmax(*sx, 0, sx->Size() - 1, (*sx)[0] + hgt * 0.22) - 1;
            mark2 = sx->Size() - 1;

            /* To handle special cases */
            if (mark >= mark2)
            {
                mark = mark2 - 1;
            }
        }
        else
        {
            mark = 0;
            mark2 = sx->Size() - 1;
        }

        /* Make the x-axis */
        X->AddRange(stt + mark, stt + mark2);

        /* Make the y-axis */
        sx->clip(mark, mark2);

        LinearFit<lwvector<double_t>>(*X, *sx, sx->Size(), mu_t, beta_t);

        sx->Erase();
        X->Erase();

        //cout << "y = " << mu_t << "x + " << beta_t << endl;
        //cout << "eValue: " << pow(10, hyp * mu_t + beta_t) * vaa << endl;
    }

    /* Compute the eValue */
    eValue = pow(10, hyp * mu_t + beta_t) * vaa;

    /* Clear arrays, vectors and variables */
    pdata->setmem(0);
    pN = 0;
    p_x->Erase();
    mu_t = 0.0;
    beta_t = 4.0;
    stt1 = stt = 0;
    end1 = ends = SIZE - 1;
    yy  = NULL;
    hyp = vaa = 0;

    return status;
}

VOID expeRT::ResetPartialVectors()
{
    /* Clear arrays, vectors, variables */
    pdata->setmem(0);
    mu_t = 0.0;
    beta_t = 4.0;
    stt = 0;
    ends = SIZE - 1;
    yy  = NULL;
    hyp = vaa = pN = 0;

}

status_t expeRT::StoreIResults(Results *rPtr, int_t spec, ebuffer *ofs)
{
    status_t status = 0;

    int_t curptr = spec * 128 * sizeof(ushort_t);
    yy = rPtr->survival;

    if (yy == NULL)
    {
        status = ERR_INVLD_PARAM;
    }

    if (status == SLM_SUCCESS)
    {
        /* Find the curve region */
        ends = rargmax<double_t *>(yy, 0, SIZE - 1, 0.99);
        stt = argmax<double_t *>(yy, 0, ends, 0.99);

        //rPtr->mu = curptr;

        for (auto ii = stt; ii <= ends; ii++)
        {
            ushort_t k = (yy[ii]);

            /* Encode into 65500 levels */
            if (rPtr->cpsms > 65500)
            {
                k = (ushort_t)(((double_t)(k * 65500))/rPtr->cpsms);
            }

            memcpy(ofs->ibuff + curptr, (const VOID *) &k, sizeof(k));
            curptr += sizeof(k);
        }

        rPtr->minhypscore = stt;
        rPtr->nexthypscore = ends;
        //rPtr->beta = rPtr->mu + 128 * sizeof(ushort_t);
    }

    yy = NULL;
    stt = 0;
    ends = 0;

    return status;
}

status_t expeRT::Reconstruct(ebuffer *ebs, int_t specno, partRes *fR)
{
    status_t status = SLM_SUCCESS;

    auto min  = fR->min;
    auto max2 = fR->max2;

    pN += fR->N;

    char_t *buffer = ebs->ibuff + (specno * 256);

    for (auto jj = min; jj <= max2; jj++)
    {
        ushort_t *val = (ushort_t*) (buffer + (jj - min) * 2);

        double_t val1 = (*val);

        /* Decode from 65500 levels */
        if (fR->N > 65500)
        {
            val1 = (val1/65500) * fR->N;
        }

        (*pdata)[jj] = (*pdata)[jj] + val1;
    }

    return status;
}

status_t expeRT::Model_logWeibull(Results *rPtr)
{
    status_t status = SLM_SUCCESS;

    /* */
    yy = rPtr->survival;
    hyp = rPtr->maxhypscore;

    /*if (yy == NULL)
    {
        status = ERR_INVLD_PARAM;
        mu_t = 0;
        beta_t = 100;
    }*/

    /*if (status == SLM_SUCCESS)*/
    {
        /* Find the curve region */
        ends = rargmax<double_t *>(yy, 0, SIZE-1, 1.0);
        stt = argmax<double_t *>(yy, 0, ends, 1.0);

        /* Slice off yyt between stt1 and end1 */
        p_x->Assign(yy + stt, yy + ends + 1);
        lwvector<double_t> *yyt = p_x;
        lwvector<double_t> *yhat = NULL;

        /* Database size
         * vaa = accumulate(yy, yy + hyp + 1, 0); */
        vaa = rPtr->cpsms;

        /* Sanity checking */
        if (stt == ends || vaa < 1)
        {
            mu_t = stt;
            beta_t = 0;
            status = ERR_NOT_ENOUGH_DATA;
        }

        /* mu estimation markers */
        auto l = std::max_element(yy + stt, (yy + ends + 1)) -
                                 (yy + stt);
        auto k = l;

        /* Enough distribution data - Proceed */
        if (status == 0)
        {
            /* Smoothen the curve using Savitzky-Golay filter */
            int_t svgl = std::min(7, ends - stt);
            int_t pln = 3;

            /* Window length = 7 or max odd number */
            if (svgl % 2 == 0)
            {
                svgl -= 1;
            }

            /* If window size > 1 */
            if (svgl > 1)
            {
                /* Polynomial = 3 or max */
                pln = std::min(3, svgl - 1);

                if ((int_t) yyt->Size() >= (svgl-1+2))
                {
                    yhat = new lwvector<double_t>(yyt->Size(), 0.0);

                    /* Smoothen the curve using SavGol: window length=7, provide win=(window-1)/2, polynomial=5 */
                    sg_smooth(yyt, yhat, std::max(1, (svgl-1)/2), pln);

                    /* Normalize yhat */
                    yhat->divide((double_t)vaa);

                    /* Normalize yyt */
                    yyt->divide((double_t)vaa);

                    /* Update the marker */
                    k = std::max_element(yhat->begin(),yhat->end()) -
                                         yhat->begin();

                    /* Mix yhat (35%) + yyt (65%) */
                    for (auto id = 0; id < (int) yyt->Size(); id++)
                    {
                        (*yyt)[id] = (*yhat)[id] * 0.4 + (*yyt)[id] * 0.6;
                    }

                    delete yhat;
                    yhat = NULL;
                }
                else
                {
                    /* Update the marker */
                    k = l;

                    /* Normalize yyt */
                    yyt->divide((double_t)vaa);
                }
            }
            else
            {
                k = l;

                /* Normalize yyt */
                yyt->divide((double_t)vaa);
            }

            /* Training parameters */
            mu_t = (stt + (k+l)/2.0);

            /* Train the logWeibull model */
            (VOID) logWeibullFit(yyt, stt, ends);
        }
    }

    /* Assign variables back to rPtr */
    rPtr->mu = mu_t;
    rPtr->beta = beta_t;
    rPtr->minhypscore = stt;
    rPtr->nexthypscore = ends;
    rPtr->maxhypscore = hyp;

    /* Reset variables */
    p_x->Erase();
    stt1 = stt = 0;
    end1 = ends = SIZE - 1;
    mu_t = 0;
    vaa = 0;
    beta_t = 4.0;
    yy = NULL;

    return status;
}

status_t expeRT::AddlogWeibull(int_t N, double_t mu, double_t beta, int_t Min, int_t Max)
{
    status_t status = 0;

    if (beta > 0)
    {
        if (pdata->Size() < SIZE)
        {
            /* Code from logWeibullResponse */
            // x = arange(st, en)
            p_x->MakeRange(0, SIZE - 1);

            //z = (xx - mu)/beta;
            p_x->add(-mu);
            p_x->divide(beta);

            for (auto ii = 0; ii < p_x->Size(); ii++)
            {
                (*p_x)[ii] = (N) * (1 / beta) * exp(-((*p_x)[ii] + exp(-(*p_x)[ii])));
            }

            /* Add to the current response */
            pdata->add(*p_x);
            pN += N;
        }
        else
        {
            Min = std::max(Min-20, 0);
            Max = std::min(Max+20, SIZE-1);

            /* Code from logWeibullResponse */
            // x = arange(st, en)
            p_x->MakeRange(Min, Max);

            //z = (xx - mu)/beta;
            p_x->add(-mu);
            p_x->divide(beta);

            for (auto ii = 0; ii < p_x->Size(); ii++)
            {
                (*p_x)[ii] = (N) * (1 / beta) * exp(-((*p_x)[ii] + exp(-(*p_x)[ii])));
                (*pdata)[ii+Min] = (*pdata)[ii+Min] + (*p_x)[ii];
            }

            pN += N;
        }

        /* Maybe erase p_x for safety purposes */
        //p_x->Erase();
    }
    else if (beta < 1e-5)
    {
        (*pdata)[(int_t)roundf(mu)] += (double_t)N;
        pN += N;
    }
    else
    {
        ;
    }

    return status;
}

inline double_t expeRT::MeanSqError(const darray &y)
{
    return (y * y).sum();
}

double_t expeRT::logWeibullFit(lwvector<double_t> *yy, int_t s, int_t e, int_t niter, double_t lr, double_t cutoff)
{
    double_t curerr = INFINITY;

    beta_t = 4.0;

    const darray y(yy->data(), yy->Size());
    const darray X1(arange(s, e));

    for (auto i = 0; i < niter; i++)
    {
        /* Gumbel distribution response */
        darray h_x(alogWeibullResponse(mu_t, beta_t, s, e));

        /* Difference */
        darray diff = y - h_x;

        /* Current error */
        curerr = MeanSqError(diff);

        /* Check for break condition */
        if (curerr < cutoff)
        {
            break;
        }

        /* Compute the partial derivatives */
        darray b = -h_x/(beta_t);
        darray c = (mu_t - X1)/beta_t - (mu_t-X1)/beta_t * exp((mu_t - X1)/beta_t);

        b = b + b * c;

        auto d = (diff * b).sum();

        darray e = h_x/beta_t;
        e = (e - e * exp((mu_t - X1)/beta_t));

        auto ee = (diff * e).sum();

        /* Update the mu and beta */
        mu_t   += lr * ee;
        beta_t += lr * d;
    }

    return curerr;
}

inline VOID expeRT::logWeibullResponse(double_t mu, double_t beta, int_t st, int_t en)
{
    // x = arange(st, en)
    p_x->MakeRange(st, en);

    //z = (xx - mu)/beta;
    p_x->add(-mu);
    p_x->divide(beta);

    for (auto ii = st; ii <= en; ii++)
    {
        (*p_x)[ii-st] =  (vaa) * (1/beta) * exp(-((*p_x)[ii-st] + exp(-(*p_x)[ii-st])));
    }
}

inline darray expeRT::alogWeibullResponse(double_t mu, double_t beta, int_t st, int_t en)
{
    darray xx = arange(st, en);

    darray z = (xx - mu)/beta;

    xx =  (1/beta) * exp(-(z + exp(-z)));

    return xx;
}

template <class T>
inline int_t expeRT::rargmax(T &data, int_t i1, int_t i2, double_t value)
{
    int_t rv = i2;

    for (auto p = i2; p >= i1; p--)
    {
        if (data[p]>= value)
        {
            rv = p;
            break;
        }
    }

    return rv;
}

template <class T>
inline int_t expeRT::argmax(T &data, int_t i1, int_t i2, double_t value)
{
    int_t rv = i1;

    for (auto p = i1; p <= i2; p++)
    {
        if (data[p] >= value)
        {
            rv = p;
            break;
        }
    }

    return rv;
}

template <class T>
inline int_t expeRT::largmax(T &data, int_t i1, int_t i2, double_t value)
{
    int_t rv = i1;

    for (auto p = i1; p <= i2; p++)
    {
        if (data[p] <= value)
        {
            rv = p;
            break;
        }
    }

    return rv;
}

/****************************************************************************
 *
 *  Source: https://people.sc.fsu.edu/~jburkardt/cpp_src/llsq/llsq.html
 *
 *
 *  Purpose:
 *
 *    LLSQ solves a linear least squares problem matching a line to data.
 *
 *  Discussion:
 *
 *    A formula for a line of the form Y = A * X + B is sought, which
 *    will minimize the root-mean-square error to N data points ( X[I], Y[I] );
 *
 *  Licensing:
 *
 *    This code is distributed under the GNU LGPL license.
 *
 *  Modified:
 *
 *    17 July 2011
 *
 *  Author:
 *
 *    John Burkardt
 *
 *  Parameters:
 *
 *    Input, int N, the number of data values.
 *
 *    Input, double X[N], Y[N], the coordinates of the data points.
 *
 *    Output, double &A, &B, the slope and Y-intercept of the least-squares
 *    approximate to the data.
 *
 ****************************************************************************
 *    Last Modified: April 12, 2020 Muhammad Haseeb (mhaseeb@fiu.edu)
 ****************************************************************************/
template <class T>
VOID expeRT::LinearFit(T& x, T& y, int_t n, double_t &a, double_t &b)
{
    int_t i;
    double_t bot;
    double_t top;
    double_t xbar;
    double_t ybar;
//
//  Special case.
//
    if (n == 1)
    {
        a = 0.0;
        b = y[0];
        return;
    }
//
//  Average X and Y.
//
    xbar = 0.0;
    ybar = 0.0;

    for (i = 0; i < n; i++)
    {
        xbar = xbar + x[i];
        ybar = ybar + y[i];
    }

    xbar = xbar / (double_t) n;
    ybar = ybar / (double_t) n;
//
//  Compute Beta.
//
    top = 0.0;
    bot = 0.0;

    for (i = 0; i < n; i++)
    {
        top = top + (x[i] - xbar) * (y[i] - ybar);
        bot = bot + (x[i] - xbar) * (x[i] - xbar);
    }

    a = top / bot;

    b = ybar - a * xbar;

    return;
}

#if 0
status_t expeRT::ModelSurvivalFunction(double_t &eValue, int_t max)
{
    status_t status = SLM_SUCCESS;

    /* Get hyp from the packet */
    hyp = max;

    /* Find the curve region */
    ends = rargmax(*pdata, 0, hyp - 1, 0.99);
    stt = argmax(*pdata, 0, ends, 0.99);

    /* To handle special cases */
    if (stt == ends)
    {
        stt = ends;
        ends += 1;
    }

    /* Candidate PSMs */
    vaa = pN;

    /* Check if no distribution data except for hyp */
    if (vaa < 1)
    {
        mu_t = 0;
        beta_t = 100;
        status = -2;
    }

    /*  */
    if (status == SLM_SUCCESS)
    {
        /* Normalize the data */
        pdata->divide(vaa);

        /* Compute survival function s(x) */
        sx->Assign(pdata->begin() + stt, pdata->begin() + ends + 1);

        //sx->print();

        /* Cumulative Sum */
        std::partial_sum(begin(*pdata) + stt, begin(*pdata) + ends + 1, sx->begin());

        /* Survival function s(x) */
        sx->add(-1.0);
        sx->negative();

        /* Adjust for > 1 */
        std::replace_if(sx->begin(), sx->end(), isLargerThan1<double_t>, 0.999);

        /* Adjust for negatives */
        int_t replacement = rargmax(*sx, 0, sx->Size()-1, 1e-4);
        std::replace_if(sx->begin(), sx->end(), isZeroNegative<double_t>, (*sx)[replacement]);

        /* log10(s(x)) */
        std::transform(sx->begin(), sx->end(), sx->begin(), [](double_t& c) { return log10(c); });

        /* Offset markers */
        auto mark = 0;
        auto mark2 = 0;
        auto hgt = (*sx)[sx->Size() - 1] - (*sx)[0];

        /* If length > 4, then find thresholds */
        if (sx->Size() > 3)
        {
            mark = largmax<lwvector<double_t>>(*sx, 0, sx->Size() - 1, (*sx)[0] + hgt * 0.15) - 1;
            mark2 = rargmax<lwvector<double_t>>(*sx, 0, sx->Size() - 1,
                    (*sx)[0] + hgt * 0.85);

            if (mark2 == (int_t) sx->Size())
            {
                mark2 -= 1;
            }

            /* To handle special cases */
            if (mark >= mark2)
            {
                mark = mark2 - 1;
            }
        }
        /* If length < 4 business as usual */
        else if (sx->Size() == 3)
        {
            /* Mark the start of the regression point */
            mark = largmax<lwvector<double_t>>(*sx, 0, sx->Size() - 1, (*sx)[0] + hgt * 0.22) - 1;
            mark2 = sx->Size() - 1;

            /* To handle special cases */
            if (mark >= mark2)
            {
                mark = mark2 - 1;
            }
        }
        else
        {
            mark = 0;
            mark2 = sx->Size() - 1;
        }

        X->AddRange(stt+mark, stt+mark2);

        /* Make the y-axis */
        sx->clip(mark, mark2);

        LinearFit<lwvector<double_t>>(*X, *sx, sx->Size(), mu_t, beta_t);

        //cout << "y = " << mu_t << "x + " << beta_t << endl;
        //cout << "eValue: " << pow(10, hyp * mu_t + beta_t) * vaa << endl;

        X->Erase();
        sx->Erase();
    }

    /* Compute the eValue */
    eValue = pow(10, hyp * mu_t + beta_t) * vaa;

    /* Clear arrays, vectors, variables */
    pdata->setmem(0);
    mu_t = 0.0;
    beta_t = 4.0;
    stt = 0;
    ends = SIZE - 1;
    yy  = NULL;
    hyp = vaa = pN = 0;

    return status;
}
#endif /* 0 */

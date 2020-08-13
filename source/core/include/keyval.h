/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 Eduard L�pez
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * MODIFIED: Muhammad Haseeb, 2020
 *
*/

#ifndef KEYVAL_H_
#define KEYVAL_H_

#include <omp.h>
#include "common.h"

#ifdef USE_OMP
/*
 * FUNCTION: KeyVal_Parallel_Internal
 *
 * DESCRIPTION: Internal Parallel Key Value Sort
 *
 * INPUT:
 * @key   : Key array
 * @val   : Value array
 * @left  : Start index of arrays
 * @right : End index of arrays
 * @cutoff: Sort in serial if right-left < cutoff
 *
 * OUTPUT: none
 */
template<class K, class V>
VOID KeyVal_Parallel_Internal(K *key, V *val, int_t left, int_t right, int_t cutoff)
{
    int_t i = left, j = right;
    K tmp;
    V tmp2;

    K pivotkey = key[(left + right) / 2];
    {
        /* PARTITION PART */
        while (i <= j)
        {
            while (key[i] < pivotkey)
                i++;
            while (key[j] > pivotkey)
                j--;
            if (i <= j)
            {
                tmp = key[i];
                tmp2 = val[i];

                key[i] = key[j];
                val[i] = val[j];

                key[j] = tmp;
                val[j] = tmp2;
                i++;
                j--;
            }
        }

    }

    if (((right - left) < cutoff))
    {
        if (left < j)
        {
            KeyVal_Parallel_Internal<K, V>(key, val, left, j, cutoff);
        }
        if (i < right)
        {
            KeyVal_Parallel_Internal<K, V>(key, val, i, right, cutoff);
        }

    }
    else
    {
#pragma omp task
        {
            KeyVal_Parallel_Internal<K, V>(key, val, left, j, cutoff);
        }
#pragma omp task
        {
            KeyVal_Parallel_Internal<K, V>(key, val, i, right, cutoff);
        }
    }

}

/*
 * FUNCTION: KeyVal_Parallel
 *
 * DESCRIPTION: Internal Parallel Key Value Sort
 *
 * INPUT:
 * @key       : Key array
 * @val       : Value array
 * @lenArray  : Length of arrays
 * @numThreads: Number of parallel threads
 *
 * OUTPUT: none
 */
template<class K, class V>
VOID KeyVal_Parallel(K* key, V *val, uint_t lenArray, uint_t numThreads)
{
    int_t cutoff = 1000;

#pragma omp parallel num_threads(numThreads)
    {
#pragma omp single nowait
        {
            KeyVal_Parallel_Internal<K, V>(key, val, 0, lenArray - 1, cutoff);
        }
    }
}

#endif /* USE_OMP */

/*
 * FUNCTION: KeyVal_Serial_Internal
 *
 * DESCRIPTION: Internal Serial Key Value Sort
 *
 * INPUT:
 * @key  : Key array
 * @val  : Value array
 * @left : Start index of arrays
 * @right: End index of arrays
 *
 * OUTPUT: none
 */
template<class K, class V>
VOID KeyVal_Serial_Internal(K *key, V *val, int_t low, int_t high)
{
    int_t i = low;
    int_t j = high;
    K pivotkey = key[(i + j) / 2];
    K tmpkey;
    V tmpval;

    while (i <= j)
    {
        while (key[i] < pivotkey)
            i++;
        while (key[j] > pivotkey)
            j--;

        if (i <= j)
        {
            tmpkey = key[i];
            tmpval = val[i];

            key[i] = key[j];
            val[i] = val[j];

            key[j] = tmpkey;
            val[j] = tmpval;

            i++;
            j--;
        }
    }

    if (j > low)
        KeyVal_Serial_Internal(key, val, low, j);
    if (i < high)
        KeyVal_Serial_Internal(key, val, i, high);
}

/*
 * FUNCTION: KeyVal_Serial
 *
 * DESCRIPTION: Internal Parallel Key Value Sort
 *
 * INPUT:
 * @key     : Key array
 * @val     : Value array
 * @lenArray: Length of arrays
 *
 * OUTPUT: none
 */
template<class K, class V>
VOID KeyVal_Serial(K *key, V *val, uint_t lenArray)
{
    KeyVal_Serial_Internal<K, V>(key, val, 0, lenArray-1);
}


#endif /* KEYVAL_H_ */

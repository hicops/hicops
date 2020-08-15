/*
 * This file is part of SLM-Transform
 *  Copyright (C) 2019  Muhammad Haseeb, Fahad Saeed
 *  Florida International University, Miami, FL
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef INCLUDE_UTILS_H_
#define INCLUDE_UTILS_H_

#include <omp.h>
#include <parallel/algorithm>
#include <bits/stdc++.h>
#include <functional>
#include <cstdlib>
#include "common.h"
#include "common.h"
#include "slm_dsts.h"

#define BITSET(x, bit)                      ((x) |  ((ull_t)((ull_t)1 << (bit)))
#define BITCLR(x, bit)                      ((x) & ~((ull_t)((ull_t)1 << (bit)))
#define ISBITSET(x,bit)                     ((x) &  ((ull_t)((ull_t)1 << (bit))))


/*
 * FUNCTION: UTILS_Sort
 *
 * DESCRIPTION: Sort an array
 *
 * INPUT:
 * @data      : Array to sort
 * @size      : Size of array
 * @descending: Sort in descending order?
 *
 * OUTPUT: none
 */
template <class Num>
void   UTILS_Sort(Num *data, uint_t size, bool descending = false)
{
    if (descending == false)
    {
        std::sort(data, data + size);
    }
    else
    {
        std::sort(data, data+size, std::greater<Num>());
    }
}

/*
 * FUNCTION: UTILS_ParallelSort
 *
 * DESCRIPTION: Parallel Sort an array
 *
 * INPUT:
 * @data      : Array to sort
 * @size      : Size of array
 * @descending: Sort in descending order?
 *
 * OUTPUT: none
 */
template <class Num>
void   UTILS_ParallelSort(Num *data, uint_t size, bool descending = false)
{
    if (descending == false)
    {
        __gnu_parallel::sort(data, data + size);
    }
    else
    {
        __gnu_parallel::sort(data, data+size, std::greater<Num>());
    }
}

/*
 * FUNCTION: UTILS_GetNumProcs
 *
 * DESCRIPTION: Get the number of processors
 *              from the environment variable
 *
 * INPUT: none
 *
 * OUTPUT:
 * @procs: Number of processors
 */
uint_t   UTILS_GetNumProcs();

/*
 * FUNCTION: UTILS_Shuffle
 *
 * DESCRIPTION: Shuffles an array

 *
 * INPUT:
 * @arr : Array to shuffle
 * @N   : Size of array
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t UTILS_Shuffle(uint_t *old, uint_t N);

/*
 * FUNCTION: UTILS_ShuffleI
 *
 * DESCRIPTION: Shuffles an array using specific seed
 *
 * INPUT:
 * @arr : Array to shuffle
 * @N   : Size of array
 * @seed: The seed to use
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t UTILS_ShuffleI(uint_t *arr, uint_t N, ull_t seed);

/*
 * FUNCTION: UTILS_GenerateSpectrum
 *
 * DESCRIPTION: Generates theoretical spectrum of a peptide
 *
 * INPUT:
 * @seq     : Peptide sequence
 * @len     : Length of peptide
 * @Spectrum: Pointer to the theoretical spectrum
 *
 * OUTPUT:
 * @mass: Precursor mass of peptide
 */
float_t  UTILS_GenerateSpectrum(char_t *seq, uint_t len, uint_t *Spectrum);

/*
 * FUNCTION: UTILS_CalculatePepMass
 *
 * DESCRIPTION: Calculate precursor mass of a peptide
 *
 * INPUT:
 * @seq: Peptide sequence
 * @len: Length of peptide
 *
 * OUTPUT:
 * @mass: Precursor mass of peptide
 */
float_t UTILS_CalculatePepMass(AA *seq, uint_t len);

#ifdef VMODS

/*
 * FUNCTION: UTILS_InitializeModInfo
 *
 * DESCRIPTION: Initialize the gModInfo structure
 *
 * INPUT:
 * @modconditions: The initialization information
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t UTILS_InitializeModInfo(string_t modconditions);

/*
 * FUNCTION: UTILS_InitializeModInfo
 *
 * DESCRIPTION: Initialize the gModInfo structure
 *
 * INPUT:
 * @vMods: gModInfo information
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t UTILS_InitializeModInfo(SLM_vMods *vMods);

/*
 * FUNCTION: UTILS_CalculateModMass
 *
 * DESCRIPTION: Calculate precursor mass of a mod
 *
 * INPUT:
 * @seq     : Modified Peptide sequence
 * @len     : Length of Modified Peptide
 * @vmodInfo: Information about the modifications
 *
 * OUTPUT:
 * @mass: Precursor mass of modified peptide
 */
float_t UTILS_CalculateModMass(AA *seq, uint_t len, uint_t vModInfo);

/*
 * FUNCTION: UTILS_GenerateModSpectrum
 *
 * DESCRIPTION: Generates theoretical spectrum for a mod
 *
 * INPUT:
 * @seq     : Modified peptide sequence
 * @len     : Length of modified peptide
 * @Spectrum: Pointer to the theoretical spectrum
 * @modInfo : Modified peptide information
 *
 * OUTPUT:
 * @mass: Precursor mass of modified peptide
 */
float_t UTILS_GenerateModSpectrum(char_t *seq, uint_t len, uint_t *Spectrum, modAA modInfo);

#endif /* VMODS */

#endif /* INCLUDE_UTILS_H_ */

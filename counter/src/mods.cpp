/*
 * This file is part of SLM-Transform
 *  Copyright (C) 2019  Fatima Afzali, Muhammad Haseeb, Fahad Saeed
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
#include "mods.h"
#include "lbe.h"

using namespace std;

/* Global Variables */
map<AA, INT> condLookup;
UINT lclcntr;
UINT partcntr;
STRING *mods;
vector<STRING> tokens;
UINT limit = 0;
ULONGLONG Comb[MAX_COMBS][MAX_COMBS];

/* External Variables */
extern gParams params;
/* Peptide Sequences */
extern vector<STRING>     Seqs;

/* Static Functions */
static LONGLONG count(STRING s, STRING conditions);
static VOID MODS_GenCombinations();

/*
 * FUNCTION: MODS_GenCombinations
 *
 * DESCRIPTION: Generates all nCk required
 *
 * INPUT: none
 *
 * OUTPUT: none
 */
static VOID MODS_GenCombinations()
{
    //run this at start of main to fill Comb with the proper values
    Comb[0][0] = 1;
    Comb[1][0] = 1;
    Comb[1][1] = 1;
    for (int i = 1; i < 64; i++)
    {
        Comb[i][0] = 1;
        for (int j = 1; j <= (i + 1) / 2; j++)
        {
            Comb[i][j] = Comb[i - 1][j - 1] + Comb[i - 1][j];
        }
        for (int j = i / 2; j < i; j++)
        {
            Comb[i][j] = Comb[i][i - j];
        }
        Comb[i][i] = 1;
    }
}

/*
 * FUNCTION: combine
 *
 * DESCRIPTION: Generates nCk
 *
 * INPUT:
 * @n: n
 * @k: k
 *
 * OUTPUT:
 * @Comb: nCk
 */
ULONGLONG combine(int n, int k) {
    //return factorial(n) / factorial(k) / factorial(n-k);
    return Comb[n][k];
}

/*
 * FUNCTION: partition2
 *
 * DESCRIPTION:
 *
 * INPUT:
 * @a: Count of AAs can be modified in seq
 * @b: Count of mods that can be made of condition in A.
 *
 * OUTPUT:
 * @sum: Number of ways to pick upto b elements
 *       from multiset A.
 */
LONGLONG partition2(vector<INT> a, INT b)
{
    INT sum = 0;
    vector<INT> a2;

    /* Set up basic conditions */
    for (INT t : a)
    {
        sum += t;
    }

    if (b > sum)
    {
        b = sum;
    }

    if (b <= 0)
    {
        return (LONGLONG) (b == 0);
    }

    if (a.size() == 0)
    {
        return 1;
    }

    sum = 0;

    for (UINT i = 1; i < a.size(); i++)
    {
        a2.push_back(a[i]);
    }

    for (INT i = 0; i < a[0] + 1; i++)
    {
        sum += combine(a[0], i) * partition2(a2, b - i);
    }

    return sum;
}

/*
 * FUNCTION: partition2
 *
 * DESCRIPTION:
 *
 * INPUT:
 * @A    : List of AA counts for each mod condition
 * @B    : Count of mods per modcondition
 * @limit: Max mods per peptide sequence
 *
 * OUTPUT:
 * @sum:
 */
LONGLONG partition3(vector<vector<INT> > A, vector<INT> B, INT limit)
{
    if (A.size() == 1)
    {
        if (B[0] <= limit)
        {
            return partition2(A[0], B[0]);
        }
        else
        {
            return partition2(A[0], limit);
        }
    }

    LONGLONG sum = 0;
    vector<vector<INT> > A2;
    vector<INT> B2;

    for (UINT i = 1; i < A.size(); i++)
    {
        A2.push_back(A[i]);
        B2.push_back(B[i]);
    }

    for (INT i = 0; i < B[0] + 1; i++)
    {
        sum += (partition2(A[0], i) - partition2(A[0], i - 1)) * partition3(A2, B2, limit - i);
    }

    return sum;
}

/*
 * FUNCTION: count
 *
 * DESCRIPTION: Main partition method
 *
 * INPUT:
 * @s         : Peptide sequence
 * @conditions: Mod conditions
 *
 * OUTPUT:
 * @nmods: Number of mods generated for @s
 */
static LONGLONG count(STRING s, STRING conditions)
{
    map<CHAR, INT> AAcounts;
    vector<vector<INT>> A;
    vector<INT> B;
    vector<INT> temp;

    for (auto c : s)
    {
        AAcounts[c] += 1;
    }

    for (UINT i = 0; i < (tokens.size() - 1) / 2; i++)
    {
        for (UINT j = 0; j < tokens[2 * i + 1].length(); j++)
        {
            temp.push_back(AAcounts[tokens[2 * i + 1][j]]);
        }

        A.push_back(temp);
        temp.clear();
        B.push_back(stoi(tokens[2 * i + 2]));
    }

    return partition3(A, B, limit);
}

/*
 * FUNCTION: MODS_ModCounter
 *
 * DESCRIPTION: Counts the number of modifications for
 *              all peptides in the index
 *
 * INPUT:
 * @threads   : Number of parallel threads
 * @conditions: Conditions of mod generation
 *
 * OUTPUT:
 * @cumulative: Number of mods
 */
ULONGLONG MODS_ModCounter()
{
    ULONGLONG cumulative = 0;

    UINT threads = params.threads;
    STRING conditions = params.modconditions;

#ifdef VMODS
    STRING token;
    stringstream ss(conditions);

    while (ss >> token)
    {
        tokens.push_back(token);
    }

    limit = stoi(tokens[0]);

    /* Return if no mods to generate */
    if (limit > 0)
    {
        /* Generate all possible combinations pre-handed */
        (VOID) MODS_GenCombinations();

        /* Parallel modcounter */
#ifdef _OPENMP
            /* The parallel for loop */
#pragma omp parallel for num_threads (threads) schedule(static) reduction(+: cumulative)
            for (UINT i = 0; i < Seqs.size(); i++)
            {
                cumulative += count(Seqs.at(i), conditions) - 1;
            }

#else
        LBE_UNUSED_PARAM(threads);

       for (UINT i = 0; i < Seqs.size(); i++)
       {
            cumulative += count(Seqs.at(i), conditions) - 1;
       }

#endif /* _OPENMP */

#endif /* VMODS */

#ifndef VMODS
        LBE_UNUSED_PARAM(conditions);
#endif /* VMODS */
    }

    lclcntr = 0;
    partcntr = 0;
    tokens.clear();
    Seqs.clear();
    limit = 0;

    return cumulative;

}

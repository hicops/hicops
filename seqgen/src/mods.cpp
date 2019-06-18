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

using namespace std;

/* Global Variables */
map<AA, INT> condLookup;
UINT lclcntr;
STRING *mods;
vector<STRING> tokens;
UINT limit = 0;
ULONGLONG Comb[MAX_COMBS][MAX_COMBS];

/* External Variables */
#ifdef VMODS
extern vector<STRING>     Seqs; /* Peptide Sequences */
extern varEntry    *modEntries; /* SLM Mods Index    */

/* Static Functions */
static LONGLONG count(STRING s, STRING conditions);
static VOID MODS_ModList(STRING peptide, vector<INT> conditions,
                         INT total, varEntry container, INT letter,
                         bool novel, INT modsSeen, UINT refid);
static VOID MODS_GenCombinations(VOID);
#endif /* VMODS */


INT cmpvarEntries(const VOID* lhs, const void *rhs)
{
    varEntry *a = (varEntry *) lhs;
    varEntry *b = (varEntry *) rhs;

    /* The signs are reversed on purpose
     * since larger 1st bit number means
     * smaller the modEntry item.
     */
    if (*a < *b)
    {
        return 1;
    }

    /* b is larger */
    if (*a > *b)
    {
        return -1;
    }

    return 0;
}

/*
 * FUNCTION: MODS_GenCombinations
 *
 * DESCRIPTION: Generates all nCk required
 *
 * INPUT: none
 *
 * OUTPUT: none
 */
static VOID MODS_GenCombinations(VOID)
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

#ifdef VMODS

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
 * FUNCTION: MODS_ModList
 *
 * DESCRIPTION: Populates varModEntries for given peptide sequence
 *
 * INPUT:
 * @peptide   : Peptide sequence
 * @conditions: Mod conditions
 * @total     :
 * @container :
 * @letter    :
 * @novel     :
 * @modsSeen  :
 *
 * OUTPUT: none
 */
static VOID MODS_ModList(STRING peptide, vector<INT> conditions,
                         INT total, varEntry container, INT letter,
                         bool novel, INT modsSeen, UINT refid)
{
    if (novel && letter != 0)
    {
        modEntries[lclcntr] = container;
        lclcntr++;
#ifdef DEBUG
        mods[lclcntr] = peptide;
#endif /* DEBUG */
    }

    if (total == 0 || letter >= (INT) peptide.length())
    {
        return;
    }

    if (condLookup[peptide[letter]] != -1 and conditions[condLookup[peptide[letter]]] > 0)
    {
        varEntry dupContainer = container;
        vector<INT> dupConditions;

        dupContainer.sites.sites |= ((ULONGLONG)1 << (letter));
        dupContainer.sites.modNum += (1 << (4 * modsSeen)) * (1 + condLookup[peptide[letter]]);
        dupContainer.seqID = refid;

        for (INT i : conditions)
        {
            dupConditions.push_back(i);
        }

        dupConditions[condLookup[peptide[letter]]] -= 1;
        STRING dupPeptide = peptide;
        dupPeptide[letter] += 32;

        MODS_ModList(dupPeptide, dupConditions, total - 1, dupContainer, letter + 1, true, modsSeen + 1, refid);
    }

    if (letter < (INT) peptide.size())
    {
        MODS_ModList(peptide, conditions, total, container, letter + 1, false, modsSeen, refid);
    }

    return;
}
#endif /* VMODS */

/*
 * FUNCTION: MODS_ModCounter
 *
 * DESCRIPTION: Counts the number of modifications for
 *              all peptides in SLM Peptide Index
 *
 * INPUT:
 * @threads   : Number of parallel threads
 * @conditions: Conditions of mod generation
 *
 * OUTPUT:
 * @cumulative: Number of mods
 */
LONGLONG MODS_ModCounter(UINT threads, STRING conditions)
{
    LONGLONG cumulative = 0;

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
#pragma omp parallel for schedule(static) reduction(+: cumulative)
            for (UINT i = 0; i < Seqs.size(); i++)
            {
                cumulative += count(Seqs.at(i), conditions) - 1;
            }


#else

        LBE_UNUSED_PARAM(threads);

        for (UINT i = 0; i < Seqs.size(); i++)
        {
            UINT mods =  count(Seqs.at(i), conditions) - 1;
            cumulative += mods;
        }
#endif /* _OPENMP */

#endif /* VMODS */

#ifndef VMODS
        LBE_UNUSED_PARAM(conditions);
#endif /* VMODS */
    }

    return cumulative;

}

/*
 * FUNCTION: MODS_GenerateMods
 *
 * DESCRIPTION: Generate modEntries for all peptide
 *              sequences in the SLM Peptide Index
 *
 * INPUT:
 * @threads   : Number of parallel threads
 * @modCount  : Number of mods that should be generated
 * @conditions: Conditions of mod generation
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS MODS_GenerateMods(UINT threads, UINT modCount, STRING conditions)
{
    STATUS status = SLM_SUCCESS;

    LBE_UNUSED_PARAM(threads);

#ifdef VMODS
    STRING allLetters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    varEntry container;
    vector<INT> condList;

    /* Reset the local counter */
    lclcntr = 0;

    /* Reset conditions for all letters */
    for (INT i = 0; i < 26; i++)
    {
        condLookup[allLetters[i]] = -1;
    }

    /* Set up conditions lookup table */
    for (UINT i = 0; i < (tokens.size() - 1) / 2; i++)
    {
        for (UINT j = 0; j < tokens[(2 * i) + 1].size(); j++)
        {
            /* Fill up the condLookup */
            condLookup[tokens[(2 * i) + 1][j]] = i;
        }

        /* Push the allowed modification numbers into condList */
        condList.push_back(stoi(tokens[(2 * i) + 2]));
    }

#ifdef DEBUG
    /* Print the number of modifications */
    cout << numMods << '\t' << modCount << endl;
    STRING *mods = new STRING[modCount];
    mods = new STRING[modCount];
#endif /* DEBUG */

    for (UINT i = 0; i < Seqs.size(); i++)
    {
#ifdef DEBUG
        cout << Seqs.at(i).length() << '\t' << lclcntr << '\t' << Seqs.at(i) << '\t' << modCount << endl;
#endif /* DEBUG */
        UINT stt = lclcntr;
        MODS_ModList(Seqs[i], condList, limit, container, 0, false, 0, i);
        UINT ssz = lclcntr - stt;

        std::qsort((void *)(modEntries+stt), ssz, sizeof(varEntry), cmpvarEntries);

    }

    if (lclcntr != modCount)
    {
        status = ERR_INVLD_SIZE;
    }

#ifdef DEBUG
    /* Print for DEBUG information */
    for(UINT i = 0; i < modCount; i++)
    {
        cout << i << ": " << mods[i] << '\t';
        printf("0x%08x\t", entries[i].sites.modNum);
        cout << bitset<MAX_SEQ_LEN>(entries[i].sites.sites) << endl;
    }

    delete[] mods;

#endif /* DEBUG */

#endif /* VMODS */

    lclcntr = 0;
    tokens.clear();
    limit = 0;

    return status;
}

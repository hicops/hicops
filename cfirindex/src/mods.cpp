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
STRING *mods;
Index *lclindex;
/* SLM Mods Index   */
pepEntry    *modEntries;

/* Initialized only once global */
map<AA, INT> condLookup;
vector<STRING> tokens;
UINT limit = 0;
ULONGLONG Comb[MAX_COMBS][MAX_COMBS];
vector<INT> condList;

/* External Variables */
extern gParams params;
extern UINT *varCount;

#ifdef VMODS

/* Peptide Sequences */
extern vector<STRING> Seqs;

/* Static Functions */
static ULONGLONG count(STRING s, STRING conditions);

static VOID MODS_ModList(STRING peptide, vector<INT> conditions,
                         INT total, pepEntry container, INT letter,
                         bool novel, INT modsSeen, UINT refid,
                         UINT &global, UINT &local);

static VOID MODS_GenCombinations();

#endif /* VMODS */

/*
 * FUNCTION: MODS_Initialize
 *
 * DESCRIPTION: Initialize global variables and conditions
 *
 * INPUT: none
 *
 * OUTPUT: status of execution
 */
STATUS MODS_Initialize()
{
#ifdef VMODS
    STRING conditions = params.modconditions;
    STRING allLetters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    pepEntry container;
    STRING token;
    stringstream ss(conditions);

    while (ss >> token)
    {
        tokens.push_back(token);
    }

    limit = stoi(tokens[0]);

    /* Generate all possible combinations pre-handed */
    (VOID) MODS_GenCombinations();

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

#endif /* VMODS */

    /* Return SLM_SUCCESS */
    return SLM_SUCCESS;
}

INT cmpvarEntries(const VOID* lhs, const void *rhs)
{
    pepEntry *a = (pepEntry *) lhs;
    pepEntry *b = (pepEntry *) rhs;

    /* The signs are reversed on purpose
     * since larger 1st bit number means
     * smaller the modEntry item.
     */
    if (*a << *b)
    {
        return 1;
    }

    /* b is larger */
    if (*a >> *b)
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
static VOID MODS_GenCombinations()
{
    //run this at start of main to fill Comb with the proper values
    Comb[0][0] = 1;
    Comb[1][0] = 1;
    Comb[1][1] = 1;

    for (int i = 1; i < MAX_COMBS; i++)
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
ULONGLONG combine(int n, int k)
{
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
    LONGLONG sum = 0;
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
 * @sequence  : Peptide sequence
 * @conditions: Mod conditions
 *
 * OUTPUT:
 * @nmods: Number of mods generated for @seq
 */
static ULONGLONG count(STRING s, STRING conditions)
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
                         INT total, pepEntry container, INT letter,
                         bool novel, INT modsSeen, UINT refid, UINT &global, UINT &local)
{
    if (novel && letter != 0)
    {
        if (LBE_ApplyPolicy(lclindex, true, global) == true)
        {
            modEntries[local] = container;
            modEntries[local].Mass = UTILS_CalculateModMass((AA *)Seqs.at(modEntries[local].seqID).c_str(),
                                                                Seqs.at(0).length(),
                                                                modEntries[local].sites.modNum);
            local++;
        }

        global++;

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
        pepEntry dupContainer = container;
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

        MODS_ModList(dupPeptide, dupConditions, total - 1, dupContainer, letter + 1, true, modsSeen + 1, refid, global, local);
    }

    if (letter < (INT) peptide.size())
    {
        MODS_ModList(peptide, conditions, total, container, letter + 1, false, modsSeen, refid, global, local);
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
ULONGLONG MODS_ModCounter(Index *index)
{
    ULONGLONG cumulative = 0;

    UINT threads = params.threads;
    STRING conditions = params.modconditions;

#ifdef VMODS

    /* Return if no mods to generate */
    if (limit > 0)
    {
        /* Parallel modcounter */
#ifdef _OPENMP
            /* The parallel for loop */
#pragma omp parallel for num_threads (threads) schedule(static) reduction(+: cumulative)
            for (UINT i = 0; i < Seqs.size(); i++)
            {
                varCount[i] = count(Seqs.at(i), conditions) - 1;
                cumulative += varCount[i];
            }

#else
        LBE_UNUSED_PARAM(threads);

       for (UINT i = 0; i < Seqs.size(); i++)
       {
            varCount[i] = count(Seqs.at(i), conditions) - 1;
            cumulative += varCount[i];
       }

       if (varCount[Seqs.size()]  != index->modCount)
       {
           cumulative = (UINT)(-1);
       }

#endif /* _OPENMP */

#endif /* VMODS */

#ifndef VMODS
        LBE_UNUSED_PARAM(conditions);
#endif /* VMODS */
    }

    UINT count = varCount[0];
    varCount[0] = 0;

    for (UINT ii = 1; ii <= Seqs.size(); ii++)
    {
        UINT tmpcount = varCount[ii];
        varCount[ii] = varCount[ii - 1] + count;
        count = tmpcount;
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
STATUS MODS_GenerateMods(Index *index)
{
    STATUS status = SLM_SUCCESS;
    pepEntry* idx = index->pepEntries + index->lclpepCnt;
    STRING conditions = params.modconditions;

    modEntries = idx;

#ifdef VMODS
    STRING allLetters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    pepEntry container;

    /* Copy the global condList into the local one */
    vector<INT> lclcondList = condList;

    lclindex = index;

#ifdef DEBUG
    /* Print the number of modifications */
    cout << numMods << '\t' << modCount << endl;
    STRING *mods = new STRING[modCount];
    mods = new STRING[modCount];
#endif /* DEBUG */

#ifdef _OPENMP
#pragma omp parallel for num_threads(params.threads) schedule (static)
#endif /* _OPENMP */
    for (UINT i = 0; i < Seqs.size(); i++)
    {
#ifdef DEBUG
        cout << Seqs.at(i).length() << '\t' << lclcntr << '\t' << Seqs.at(i) << '\t' << modCount << endl;
#endif /* DEBUG */

        /* Make global and local index */
        UINT globalidx = varCount[i];
        UINT localidx  = UINT(varCount[i] / params.nodes);

        UINT stt = localidx;
        MODS_ModList(Seqs[i], lclcondList, limit, container, 0, false, 0, i, globalidx, localidx);
        UINT ssz = localidx - stt;

        std::qsort((void *)(modEntries + stt), ssz, sizeof(pepEntry), cmpvarEntries);
    }

/*
    if (partcntr != index->lclmodCnt)
    {
        status = ERR_INVLD_SIZE;
    }
*/

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

    /* Return the status */
    return status;
}

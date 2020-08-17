/*
 * Copyright (C) 2019  Fatima Afzali, Muhammad Haseeb, Fahad Saeed
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
#include "mods.h"
#include "lbe.h"

using namespace std;

/* Global Variables */
string_t *mods;
Index *lclindex;
/* SLM Mods Index   */
pepEntry    *modEntries;

/* Initialized only once global */
map<AA, int_t> condLookup;
vector<string_t> tokens;
uint_t limit = 0;
ull_t Comb[MAX_COMBS][MAX_COMBS];
vector<int_t> condList;

/* External Variables */
extern gParams params;
extern uint_t *varCount;

#ifdef VMODS

/* Peptide Sequences */
extern vector<string_t> Seqs;

/* Static Functions */
static ull_t count(string_t s);

static VOID MODS_ModList(string_t peptide, vector<int_t> conditions,
                         int_t total, pepEntry container, int_t letter,
                         bool novel, int_t modsSeen, uint_t refid,
                         uint_t &global, uint_t &local);

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
status_t MODS_Initialize()
{
#ifdef VMODS
    string_t conditions = params.modconditions;
    string_t allLetters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    pepEntry container;
    string_t token;
    stringstream ss(conditions);

    while (ss >> token)
    {
        tokens.push_back(token);
    }

    limit = stoi(tokens[0]);

    /* Generate all possible combinations pre-handed */
    (VOID) MODS_GenCombinations();

    /* Reset conditions for all letters */
    for (int_t i = 0; i < 26; i++)
    {
        condLookup[allLetters[i]] = -1;
    }

    /* Set up conditions lookup table */
    for (uint_t i = 0; i < (tokens.size() - 1) / 2; i++)
    {
        for (uint_t j = 0; j < tokens[(2 * i) + 1].size(); j++)
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

int_t cmpvarEntries(const VOID* lhs, const void *rhs)
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
ull_t combine(int n, int k)
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
longlong_t partition2(vector<int_t> a, int_t b)
{
    longlong_t sum = 0;
    vector<int_t> a2;

    /* Set up basic conditions */
    for (int_t t : a)
    {
        sum += t;
    }

    if (b > sum)
    {
        b = sum;
    }

    if (b <= 0)
    {
        return (longlong_t) (b == 0);
    }

    if (a.size() == 0)
    {
        return 1;
    }

    sum = 0;

    for (uint_t i = 1; i < a.size(); i++)
    {
        a2.push_back(a[i]);
    }

    for (int_t i = 0; i < a[0] + 1; i++)
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
longlong_t partition3(vector<vector<int_t> > A, vector<int_t> B, int_t limit)
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

    longlong_t sum = 0;
    vector<vector<int_t> > A2;
    vector<int_t> B2;

    for (uint_t i = 1; i < A.size(); i++)
    {
        A2.push_back(A[i]);
        B2.push_back(B[i]);
    }

    for (int_t i = 0; i < B[0] + 1; i++)
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
static ull_t count(string_t s)
{
    map<char_t, int_t> AAcounts;
    vector<vector<int_t>> A;
    vector<int_t> B;
    vector<int_t> temp;

    for (auto c : s)
    {
        AAcounts[c] += 1;
    }

    for (uint_t i = 0; i < (tokens.size() - 1) / 2; i++)
    {
        for (uint_t j = 0; j < tokens[2 * i + 1].length(); j++)
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
static VOID MODS_ModList(string_t peptide, vector<int_t> conditions,
                         int_t total, pepEntry container, int_t letter,
                         bool novel, int_t modsSeen, uint_t refid, uint_t &global, uint_t &local)
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

    if (total == 0 || letter >= (int_t) peptide.length())
    {
        return;
    }

    if (condLookup[peptide[letter]] != -1 and conditions[condLookup[peptide[letter]]] > 0)
    {
        pepEntry dupContainer = container;
        vector<int_t> dupConditions;

        dupContainer.sites.sites |= ((ull_t)1 << (letter));
        dupContainer.sites.modNum += (1 << (4 * modsSeen)) * (1 + condLookup[peptide[letter]]);
        dupContainer.seqID = refid;

        for (int_t i : conditions)
        {
            dupConditions.push_back(i);
        }

        dupConditions[condLookup[peptide[letter]]] -= 1;
        string_t dupPeptide = peptide;
        dupPeptide[letter] += 32;

        MODS_ModList(dupPeptide, dupConditions, total - 1, dupContainer, letter + 1, true, modsSeen + 1, refid, global, local);
    }

    if (letter < (int_t) peptide.size())
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
ull_t MODS_ModCounter()
{
    ull_t cumulative = 0;

    uint_t threads = params.threads;
    string_t conditions = params.modconditions;

#ifdef VMODS

    /* Return if no mods to generate */
    if (limit > 0)
    {
        /* Parallel modcounter */
#ifdef USE_OMP
            /* The parallel for loop */
#pragma omp parallel for num_threads (threads) schedule(static) reduction(+: cumulative)
            for (uint_t i = 0; i < Seqs.size(); i++)
            {
                varCount[i] = count(Seqs.at(i)) - 1;
                cumulative += varCount[i];
            }

#else
        LBE_UNUSED_PARAM(threads);

       for (uint_t i = 0; i < Seqs.size(); i++)
       {
            varCount[i] = count(Seqs.at(i)) - 1;
            cumulative += varCount[i];
       }

#endif /* USE_OMP */

#endif /* VMODS */

#ifndef VMODS
        LBE_UNUSED_PARAM(conditions);
#endif /* VMODS */


        uint_t count = varCount[0];
        varCount[0] = 0;

        for (uint_t ii = 1; ii <= Seqs.size(); ii++)
        {
            uint_t tmpcount = varCount[ii];
            varCount[ii] = varCount[ii - 1] + count;
            count = tmpcount;
        }

        if (varCount[Seqs.size()] != cumulative)
        {
            cumulative = (uint_t) (-1);
        }
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
status_t MODS_GenerateMods(Index *index)
{
    status_t status = SLM_SUCCESS;
    pepEntry* idx = index->pepEntries + index->lclpepCnt;
    string_t conditions = params.modconditions;

    modEntries = idx;

#ifdef VMODS
    const string_t allLetters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    pepEntry container;

    /* Copy the global condList into the local one */
    vector<int_t> lclcondList = condList;

    lclindex = index;

#ifdef DEBUG
    /* Print the number of modifications */
    cout << numMods << '\t' << modCount << endl;
    string_t *mods = new string_t[modCount];
    mods = new string_t[modCount];
#endif /* DEBUG */

#ifdef USE_OMP
#pragma omp parallel for num_threads(params.threads) schedule (static)
#endif /* USE_OMP */
    for (uint_t i = 0; i < Seqs.size(); i++)
    {
#ifdef DEBUG
        cout << Seqs.at(i).length() << '\t' << lclcntr << '\t' << Seqs.at(i) << '\t' << modCount << endl;
#endif /* DEBUG */

        /* Make global and local index */
        uint_t globalidx = varCount[i];
        uint_t localidx  = uint_t(varCount[i] / params.nodes);

        uint_t stt = localidx;
        MODS_ModList(Seqs[i], lclcondList, limit, container, 0, false, 0, i, globalidx, localidx);
        uint_t ssz = localidx - stt;

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
    for(uint_t i = 0; i < modCount; i++)
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

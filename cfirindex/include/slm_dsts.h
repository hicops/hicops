/*
 * This file is part of PCDSFrame software
 * Copyright (C) 2019  Muhammad Haseeb, Fahad Saeed
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

#ifndef INCLUDE_SLM_DSTS_H_
#define INCLUDE_SLM_DSTS_H_

#include "common.h"
#include "minheap.h"
#include <cstring>

/* Types of modifications allowed by SLM_Mods     */
#define MAX_MOD_TYPES                        15
#define MAX_SEQ_LEN                          60

/************************* Common DSTs ************************/

/* Add distribution policies */
typedef enum _DistPolicy
{
    _chunk,
    _cyclic,
    _zigzag,

} DistPolicy;

typedef struct _SLM_varAA
{
    AA     residues[4]   ; /* Modified AA residues in this modification - Upto 4 */
    UINT  modMass        ; /* Scaled mass of the modification                    */
    USHORT aa_per_peptide; /* Allowed modified residues per peptide sequence     */

    _SLM_varAA()
    {
        residues[0]= '\0';
        modMass = 0;
        aa_per_peptide = 0;
    }

    _SLM_varAA& operator=(const _SLM_varAA& rhs)
    {
        /* Check for self assignment */
        if (this != &rhs)
        {
            this->aa_per_peptide = rhs.aa_per_peptide;
            this->modMass = rhs.modMass;

            this->residues[0] = rhs.residues[0];
            this->residues[1] = rhs.residues[1];
            this->residues[2] = rhs.residues[2];
            this->residues[3] = rhs.residues[3];
        }

        return *this;
    }

} SLM_varAA;

typedef struct _SLM_Mods
{
    USHORT       vmods_per_pep    ; /* Total allowed modified residues per sequence */
    USHORT            num_vars    ; /* Number of types of modifications added to index - Max: 7 */
    SLM_varAA vmods[MAX_MOD_TYPES]; /* Information for each modification added to index */

    _SLM_Mods()
    {
        vmods_per_pep = 0;
        num_vars = 0;
    }

    /* Overload = operator */
    _SLM_Mods& operator=(const _SLM_Mods& rhs)
    {
        /* Check for self assignment */
        if (this != &rhs)
        {
            this->vmods_per_pep = rhs.vmods_per_pep;
            this->num_vars = rhs.num_vars;

            for (UINT i = 0; i < MAX_MOD_TYPES; i++)
            {
                this->vmods[i] = rhs.vmods[i];
            }
        }

        return *this;
    }

} SLM_vMods;

typedef struct _pepSeq
{
    AA        *seqs; /* Stores peptide sequence, could store as strings as well */
    USHORT   peplen; /* Stores sequence length */
    UINT        AAs; /* Total number of characters */

    _pepSeq()
    {
        seqs = NULL;
        peplen = 0;
        AAs = 0;
    }

} PepSeqs;

typedef struct _modAA
{
    ULONGLONG  sites; /* maxlen(pep) = 60AA + 2 bits (termini mods)      */
    UINT  modNum    ; /* 4 bits per mods num, Max 8 mods allowed per pep */

    _modAA()
    {
        modNum = 0x0;
        sites = 0x0;
    }

    /* Overload = operator */
    _modAA& operator=(const _modAA& rhs)
    {
        /* Check for self assignment */
        if (this != &rhs)
        {
            this->sites = rhs.sites;
            this->modNum = rhs.modNum;
        }

        return *this;
    }

} modAA;

typedef struct _pepEntry
{
    FLOAT  Mass; /* Mass of Peptide            */
    IDX   seqID; /* Normal Peptide Sequence ID */
    modAA sites; /* Modified AA information    */

    /* Overload = operator */
    _pepEntry& operator=(const _pepEntry& rhs)
    {
        /* Check for self assignment */
        if (this != &rhs)
        {
            this->Mass = rhs.Mass;
            this->seqID = rhs.seqID;
            this->sites.modNum = rhs.sites.modNum;
            this->sites.sites = rhs.sites.sites;
        }
        return *this;
    }

    BOOL operator >(const _pepEntry& rhs)
    {
        return this->Mass > rhs.Mass;
    }

    BOOL operator >=(const _pepEntry& rhs)
    {
        return this->Mass >= rhs.Mass;
    }

    BOOL operator <(const _pepEntry& rhs)
    {
        return this->Mass < rhs.Mass;
    }

    BOOL operator <=(const _pepEntry& rhs)
    {
        return this->Mass <= rhs.Mass;
    }

    BOOL operator ==(const _pepEntry& rhs)
    {
        return this->Mass == rhs.Mass;
    }

    BOOL operator>>(const _pepEntry& rhs)
    {
        INT s1 = 0;
        INT s2 = 0;

        INT s3= MAX_SEQ_LEN;
        INT s4= MAX_SEQ_LEN;

        /* compute distance from n-term */
        while ((this->sites.sites >> s1 & 0x1) != 0x1 && s1++ < MAX_SEQ_LEN);
        while ((rhs.sites.sites   >> s2 & 0x1) != 0x1 && s2++ < MAX_SEQ_LEN);

        /* compute distance from c-term */
        while ((this->sites.sites >> s3 & 0x1) != 0x1 && s3-- > s1);
        while ((rhs.sites.sites   >> s4 & 0x1) != 0x1 && s4-- > s2);

       /* return which has larger distance */
        return (s1 + MAX_SEQ_LEN - s3) > (s2 + MAX_SEQ_LEN - s4);

    }

    BOOL operator>>=(const _pepEntry& rhs)
    {
        INT s1 = 0;
        INT s2 = 0;

        INT s3= MAX_SEQ_LEN;
        INT s4= MAX_SEQ_LEN;

        /* compute distance from n-term */
        while ((this->sites.sites >> s1 & 0x1) != 0x1 && s1++ < MAX_SEQ_LEN);
        while ((rhs.sites.sites   >> s2 & 0x1) != 0x1 && s2++ < MAX_SEQ_LEN);

        /* compute distance from c-term */
        while ((this->sites.sites >> s3 & 0x1) != 0x1 && s3-- > s1);
        while ((rhs.sites.sites   >> s4 & 0x1) != 0x1 && s4-- > s2);

       /* return which has larger distance */
        return (s1 + MAX_SEQ_LEN - s3) >= (s2 + MAX_SEQ_LEN - s4);

    }

    BOOL operator<<(const _pepEntry& rhs)
    {
        INT s1 = 0;
        INT s2 = 0;

        INT s3= MAX_SEQ_LEN;
        INT s4= MAX_SEQ_LEN;

        /* compute distance from n-term */
        while ((this->sites.sites >> s1 & 0x1) != 0x1 && s1++ < MAX_SEQ_LEN);
        while ((rhs.sites.sites   >> s2 & 0x1) != 0x1 && s2++ < MAX_SEQ_LEN);

        /* compute distance from c-term */
        while ((this->sites.sites >> s3 & 0x1) != 0x1 && s3-- > s1);
        while ((rhs.sites.sites   >> s4 & 0x1) != 0x1 && s4-- > s2);

       /* return which has larger distance */
        return (s1 + MAX_SEQ_LEN - s3) < (s2 + MAX_SEQ_LEN - s4);
    }

    BOOL operator<<=(const _pepEntry& rhs)
    {
        INT s1 = 0;
        INT s2 = 0;

        INT s3= MAX_SEQ_LEN;
        INT s4= MAX_SEQ_LEN;

        /* compute distance from n-term */
        while ((this->sites.sites >> s1 & 0x1) != 0x1 && s1++ < MAX_SEQ_LEN);
        while ((rhs.sites.sites   >> s2 & 0x1) != 0x1 && s2++ < MAX_SEQ_LEN);

        /* compute distance from c-term */
        while ((this->sites.sites >> s3 & 0x1) != 0x1 && s3-- > s1);
        while ((rhs.sites.sites   >> s4 & 0x1) != 0x1 && s4-- > s2);

       /* return which has larger distance */
        return (s1 + MAX_SEQ_LEN - s3) <= (s2 + MAX_SEQ_LEN - s4);
    }

    /* Default constructor */
    _pepEntry()
    {
        Mass = 0;
        seqID = 0;
        sites.modNum = 0;
        sites.sites = 0;
    }

} pepEntry;

/************************* SLM Index DSTs ************************/
/*
 * Structure to store the sum of matched b and y ions and
 * their summed intensities for a given experimental spectrum
 * against all the candidate peptides.
 */
typedef struct _BYC
{
    UCHAR   bc;        /* b ion count */
    UCHAR   yc;        /* y ion count */

    _BYC()
    {
        bc = 0;
        yc = 0;
    }

} BYC;

typedef struct _iBYC
{
    UINT ibc;
    UINT iyc;

    _iBYC()
    {
        ibc = 0;
        iyc = 0;
    }
} iBYC;

typedef struct _SLMchunk
{
    UINT    *iA; /* Ions Array (iA)   */
    UINT    *bA; /* Bucket Array (bA) */

    _SLMchunk()
    {
        iA = NULL;
        bA = NULL;
    }

#ifdef FUTURE
    UCHAR *bits = NULL; /* Scorecard bits    */
#endif /* FUTURE */
} SLMchunk;

/* Structure for each pep file */
typedef struct _Index
{
    UINT pepCount     ;
    UINT modCount     ;
    UINT totalCount   ;

    UINT lclpepCnt    ;
    UINT lclmodCnt    ;
    UINT lcltotCnt    ;
    UINT nChunks      ;
    UINT chunksize    ;
    UINT lastchunksize;

    PepSeqs          pepIndex;
    pepEntry      *pepEntries;
    SLMchunk        *ionIndex;

    _Index()
    {
        pepCount = 0;
        modCount = 0;
        totalCount = 0;

        lclpepCnt = 0;
        lclmodCnt = 0;
        lcltotCnt = 0;
        nChunks = 0;
        chunksize = 0;
        lastchunksize = 0;

        pepEntries = NULL;
        ionIndex = NULL;
    }
} Index;

/* Structure for global Parameters */
typedef struct _globalParams
{
    UINT threads;
    UINT min_len;
    UINT max_len;
    UINT maxz;
    UINT topmatches;
    UINT scale;
    UINT min_shp;
    UINT min_cpsm;
    UINT nodes;
    UINT myid;
    UINT spadmem;

    UINT min_mass;
    UINT max_mass;
    UINT dF;

    INT  base_int;
    INT  min_int;

    DOUBLE dM;
    DOUBLE res;
    DOUBLE expect_max;

    DOUBLE *perf;

    STRING dbpath;
    STRING datapath;
    STRING workspace;
    STRING modconditions;

    DistPolicy policy;

    SLM_vMods vModInfo;

    _globalParams()
    {
        threads = 1;
        min_len = 6;
        max_len = 40;
        maxz = 3;
        topmatches = 10;
        scale = 100;
        expect_max = 20;
        min_shp = 4;
        min_cpsm = 4;
        base_int = 100000;
        min_int = 0.01 * base_int;
        nodes = 1;
        myid = 0;
        spadmem = 2048;
        min_mass = 500;
        max_mass = 5000;
        dF = (UINT) (0.02 * scale);
        dM = 500.0;
        res = 0.01;
        perf = NULL;
        policy = _cyclic;
	}
}gParams;


/* Experimental MS/MS spectra data */
typedef struct _queries
{
    UINT        *moz; /* Stores the m/z values of the spectra */
    UINT  *intensity; /* Stores the intensity values of the experimental spectra */
    UINT        *idx; /* Row ptr. Starting index of each row */
    FLOAT  *precurse; /* Stores the precursor mass of each spectrum. */
    UINT    numPeaks;
    UINT    numSpecs; /* Number of theoretical spectra */

    _queries()
    {
        moz = NULL;
        intensity = NULL;
        idx = NULL;
        precurse = NULL;
        numPeaks = 0;
        numSpecs = 0;
    }

} Queries;

/* Score entry that goes into the heap */
typedef struct _heapEntry
{
    /* The index * + offset */
    USHORT   idxoffset;

    /* Number of shared ions in the spectra */
    USHORT sharedions;

    /* Total ions in spectrum */
    USHORT totalions;

    /* Parent spectrum ID in the respective chunk of index */
    UINT  psid;

    /* Computed hyperscore */
    FLOAT hyperscore;

    /* Constructor */
    _heapEntry()
    {
        idxoffset = 0;
        psid = 0;
        hyperscore = 0;
        sharedions = 0;
        totalions = 0;
    }

    /* Overload = operator */
    _heapEntry& operator=(const _heapEntry& rhs)
    {
        /* Check for self assignment */
        if (this != &rhs)
        {
            this->idxoffset = rhs.idxoffset;
            this->psid = rhs.psid;
            this->hyperscore = rhs.hyperscore;
            this->sharedions = rhs.sharedions;
            this->totalions = rhs.totalions;
        }

        return *this;
    }

    BOOL operator <=(const _heapEntry& rhs)
    {
        return this->hyperscore <= rhs.hyperscore;
    }

    BOOL operator <(const _heapEntry& rhs)
    {
        return this->hyperscore < rhs.hyperscore;
    }

    BOOL operator >=(const _heapEntry& rhs)
    {
        return this->hyperscore >= rhs.hyperscore;
    }

    BOOL operator >(const _heapEntry& rhs)
    {
        return this->hyperscore > rhs.hyperscore;
    }

    BOOL operator ==(const _heapEntry& rhs)
    {
        return this->hyperscore == rhs.hyperscore;
    }

} hCell;

/* This structure will be made per thread */
typedef struct _Results
{
    /* Number of candidate PSMs (n) */
    UINT cpsms;

    /* Min heap to keep track of top
     */
    minHeap<hCell> topK;

    /************************
     * Required variables per
     * query for expect score
     ************************/

    /* The y = mx + b form
     * for linear regression
     */
    FLOAT weight;
    FLOAT bias;

    INT minhypscore;
    INT maxhypscore;
    INT nexthypscore;

    /* Survival function s(x) vs log(score) */
    DOUBLE *survival;
    DOUBLE *xaxis;

    /* Constructor */
    _Results()
    {
        cpsms = 0;
        weight = 0;
        bias = 0;
        minhypscore = 0;
        maxhypscore = 0;
        nexthypscore = 0;
        survival = NULL;
        xaxis = NULL;
    }

    void reset()
    {
        cpsms = 0;
        weight = 0;
        bias = 0;
        minhypscore = 0;
        maxhypscore = 0;
        nexthypscore = 0;

        std::memset(survival, 0x0, sizeof(UINT) * (2 + MAX_HYPERSCORE * 10));

        topK.heap_reset();
    }

} Results;

typedef struct _BYICount
{
    BYC     *byc;       /* Both counts */
    iBYC   *ibyc;       /* Sum of b/y ion intensities */
    Results        res;

    _BYICount()
    {
        byc = NULL;
        ibyc = NULL;
    }

} BYICount;

#define BYISIZE                 (sizeof(UCHAR) * 2 + sizeof(UINT) * 2)

#endif /* INCLUDE_SLM_DSTS_H_ */

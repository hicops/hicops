/*
 * Copyright (C) 2020  Muhammad Haseeb, Fahad Saeed
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
    uint_t  modMass        ; /* Scaled mass of the modification                    */
    ushort_t aa_per_peptide; /* Allowed modified residues per peptide sequence     */

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
    ushort_t       vmods_per_pep    ; /* Total allowed modified residues per sequence */
    ushort_t            num_vars    ; /* Number of types of modifications added to index - Max: 7 */
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

            for (uint_t i = 0; i < MAX_MOD_TYPES; i++)
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
    ushort_t   peplen; /* Stores sequence length */
    uint_t        AAs; /* Total number of characters */

    _pepSeq()
    {
        seqs = NULL;
        peplen = 0;
        AAs = 0;
    }


} PepSeqs;
typedef struct _modAA
{
    ull_t  sites; /* maxlen(pep) = 60AA + 2 bits (termini mods)      */
    uint_t  modNum    ; /* 4 bits per mods num, Max 8 mods allowed per pep */

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
    float_t  Mass; /* Mass of Peptide            */
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
        int_t s1 = 0;
        int_t s2 = 0;

        int_t s3= MAX_SEQ_LEN;
        int_t s4= MAX_SEQ_LEN;

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
        int_t s1 = 0;
        int_t s2 = 0;

        int_t s3= MAX_SEQ_LEN;
        int_t s4= MAX_SEQ_LEN;

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
        int_t s1 = 0;
        int_t s2 = 0;

        int_t s3= MAX_SEQ_LEN;
        int_t s4= MAX_SEQ_LEN;

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
        int_t s1 = 0;
        int_t s2 = 0;

        int_t s3= MAX_SEQ_LEN;
        int_t s4= MAX_SEQ_LEN;

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
    ushort_t   bc;        /* b ion count */
    ushort_t   yc;        /* y ion count */
    uint_t    ibc;
    uint_t    iyc;

    _BYC()
    {
        bc = 0;
        yc = 0;
        ibc = 0;
        iyc = 0;
    }
} BYC;

typedef struct _SLMchunk
{
    uint_t    *iA; /* Ions Array (iA)   */
    uint_t    *bA; /* Bucket Array (bA) */

    _SLMchunk()
    {
        iA = NULL;
        bA = NULL;
    }

#ifdef FUTURE
    uchar_t *bits = NULL; /* Scorecard bits    */
#endif /* FUTURE */
} SLMchunk;

/* Structure for each pep file */
typedef struct _Index
{
    uint_t pepCount     ;
    uint_t modCount     ;
    uint_t totalCount   ;

    uint_t lclpepCnt    ;
    uint_t lclmodCnt    ;
    uint_t lcltotCnt    ;
    uint_t nChunks      ;
    uint_t chunksize    ;
    uint_t lastchunksize;

    PepSeqs     pepIndex;
    pepEntry *pepEntries;
    SLMchunk   *ionIndex;

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
    uint_t threads;
    uint_t min_len;
    uint_t max_len;
    uint_t maxz;
    uint_t topmatches;
    uint_t scale;
    uint_t min_shp;
    uint_t min_cpsm;
    uint_t nodes;
    uint_t myid;
    uint_t spadmem;

    uint_t min_mass;
    uint_t max_mass;
    uint_t dF;

    int_t  base_int;
    int_t  min_int;

    double_t dM;
    double_t res;
    double_t expect_max;

    string_t dbpath;
    string_t datapath;
    string_t workspace;
    string_t modconditions;

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
        dF = 0;
        dM = 500.0;
        res = 0.01;
        policy = DistPolicy::_cyclic;
	}
}gParams;


/* Experimental MS/MS spectra data */
typedef struct _queries
{
    uint_t        *moz; /* Stores the m/z values of the spectra */
    uint_t  *intensity; /* Stores the intensity values of the experimental spectra */
    uint_t        *idx; /* Row ptr. Starting index of each row */
    float_t  *precurse; /* Stores the precursor mass of each spectrum. */
    int_t     numPeaks;
    int_t     numSpecs; /* Number of theoretical spectra */
    int_t     batchNum;

    void reset()
    {
        numPeaks = 0;
        numSpecs = 0;
        batchNum = 0;
    }

    _queries()
    {
        this->idx       = NULL;
        this->precurse  = NULL;
        this->moz       = NULL;
        this->intensity = NULL;
        numPeaks        = 0;
        numSpecs        = 0;
        batchNum        = 0;
    }

    VOID init()
    {
        this->idx       = new uint_t[QCHUNK + 1];
        this->precurse  = new float_t[QCHUNK];
        this->moz       = new uint_t[QCHUNK * QALEN];
        this->intensity = new uint_t[QCHUNK * QALEN];
        numPeaks        = 0;
        numSpecs        = 0;
        batchNum = 0;
    }

    VOID deinit()
    {
        numPeaks = -1;
        numSpecs = -1;
        batchNum = -1;

        /* Deallocate the memory */
        if (this->moz != NULL)
        {
            delete[] this->moz;
            this->moz = NULL;
        }

        if (this->intensity != NULL)
        {
            delete[] this->intensity;
            this->intensity = NULL;
        }

        if (this->precurse != NULL)
        {
            delete[] this->precurse;
            this->precurse = NULL;
        }

        if (this->idx != NULL)
        {
            delete[] this->idx;
            this->idx = NULL;
        }
    }

    ~_queries()
    {
        numPeaks = 0;
        numSpecs = 0;
        batchNum = 0;

        /* Deallocate the memory */
        if (this->moz != NULL)
        {
            delete[] this->moz;
            this->moz = NULL;
        }

        if (this->intensity != NULL)
        {
            delete[] this->intensity;
            this->intensity = NULL;
        }

        if (this->precurse != NULL)
        {
            delete[] this->precurse;
            this->precurse = NULL;
        }

        if (this->idx != NULL)
        {
            delete[] this->idx;
            this->idx = NULL;
        }
    }

} Queries;

/* Score entry that goes into the heap */
typedef struct _heapEntry
{
    /* Number of shared ions in the spectra */
    uchar_t   sharedions;

    /* The index * + offset */
    uchar_t    idxoffset;

    /* Total ions in spectrum */
    ushort_t totalions;

    /* Parent spectrum ID in the respective chunk of index */
    int_t        psid;

    float_t     pmass;

    /* Computed hyperscore */
    float_t hyperscore;

    /* Constructor */
    _heapEntry()
    {
        idxoffset  = 0;
        psid       = 0;
        hyperscore = 0;
        sharedions = 0;
        totalions  = 0;
        pmass      = 0;
    }

    /* Copy constructor */
    _heapEntry(const _heapEntry &obj)
    {
        memcpy(this, &obj, sizeof(_heapEntry));

        /*idxoffset  = obj.idxoffset;
        psid       = obj.psid;
        hyperscore = obj.hyperscore;
        sharedions = obj.sharedions;
        totalions  = obj.totalions;
        pmass      = obj.pmass;*/
    }

    /* Overload = operator */
    _heapEntry& operator=(const _heapEntry& rhs)
    {
        /* Check for self assignment */
        if (this != &rhs) 
        {
            memcpy(this, &rhs, sizeof(_heapEntry));
        }

        return *this;
    }

    /* Overload = operator */
    _heapEntry& operator=(const int_t& rhs)
    {
        this->idxoffset = rhs;
        this->psid = rhs;
        this->hyperscore = rhs;
        this->sharedions = rhs;
        this->totalions = rhs;
        this->pmass = rhs;

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

/* Structure to contain a communication request */
typedef struct _commRqst
{
    uint_t btag;
    uint_t bsize;
    uint_t buff;

    VOID _comRqst()
    {
        btag = 0;
        bsize = 0;
        buff = -1;
    }

    /* Overload = operator */
    _commRqst& operator=(const _commRqst& rhs)
    {
        /* Check for self assignment */
        if (this != &rhs)
        {
            memcpy(this, &rhs, sizeof(_commRqst));
        }

        return *this;
    }

} commRqst;


/* This structure will be made per thread */
typedef struct _Results
{
    /* Number of candidate PSMs (n) */
    uint_t cpsms;

    /* Min heap to keep track of top
     */
    minHeap<hCell> topK;

    /* The y ~ logWeibull(X, mu, beta)
     * for data fit
     */
    int_t mu;
    int_t beta;

    int_t minhypscore;
    int_t maxhypscore;
    int_t nexthypscore;

    /* Survival function s(x) vs log(score) */
    double_t *survival;

    /* Constructor */
    _Results()
    {
        cpsms = 0;
        mu = 0;
        beta = 0;
        minhypscore = 0;
        maxhypscore = 0;
        nexthypscore = 0;
        survival = NULL;;
    }

    void reset()
    {
        cpsms = 0;
        mu = 0;
        beta = 0;
        minhypscore = 0;
        maxhypscore = 0;
        nexthypscore = 0;

        std::memset(survival, 0x0, sizeof(double_t) * (2 + MAX_HYPERSCORE * 10));

        topK.heap_reset();
    }

    void reset2()
    {
        cpsms = 0;
        mu = 0;
        beta = 0;
        minhypscore = 0;
        maxhypscore = 0;
        nexthypscore = 0;

        std::memset(survival, 0x0, sizeof(double_t) * (2 + MAX_HYPERSCORE * 10));
    }

} Results;

/* Data structure for partial result Tx/Rx */
typedef struct _partResult
{
    /* Convert: x10 + 0.5 */
    ushort_t min;
    ushort_t max2;
    int_t max;

    /* Total number of samples scored */
    int_t N;
    int_t qID;

    /* Default contructor */
    _partResult()
    {
        min = 0;
        max2 = 0;
        N  = 0;
        max = 0;
        qID = 0;
    }

    _partResult(int_t def)
    {
        min = def;
        max2 = def;
        N  = def;
        max = def;
        qID = 0;
    }

    /* Destructor */
    ~_partResult()
    {
        min = 0;
        max = 0;
        N  = 0;
        max2 = 0;
        qID = 0;
    }

    _partResult& operator=(const int_t& rhs)
    {
        /* Check for self assignment */
            min = rhs;
            max2 = rhs;
            N = rhs;
            max = rhs;
            qID = rhs;

        return *this;
    }

    _partResult& operator=(const _partResult& rhs)
    {
        /* Check for self assignment */
        if (this != &rhs)
        {
            min = rhs.min;
            max2 = rhs.max2;
            N = rhs.N;
            max = rhs.max;
            qID = rhs.qID;
        }

        return *this;
    }

    BOOL operator==(const _partResult& rhs)
    {
        BOOL val = false;

        /* Check for self assignment */
        if (this != &rhs)
        {
            val = (min == rhs.min &&
                  max == rhs.max &&
                    N == rhs.N &&
                 max2 == rhs.max2 &&
                  qID == rhs.qID);
        }
        else
        {
            val= true;
        }

        return val;
    }

    BOOL operator==(const int_t rhs)
    {
        BOOL val = false;

        val = (min == rhs && max == rhs && N == rhs && max2 == rhs);

        return val;
    }

} partRes;

typedef struct _BYICount
{
    BYC     *byc;       /* Both counts */
    Results  res;

    _BYICount()
    {
        byc = NULL;
    }

} BYICount;

#define BYISIZE                 (sizeof(ushort_t) * 2 + sizeof(uint_t) * 2)

typedef struct _fResult
{
    int_t eValue;
    int_t specID;
    int_t npsms;

    _fResult()
    {
        eValue = 0;
        specID = 0;
        npsms = 0;
    }

    ~_fResult()
    {
        this->eValue = 0;
        this->specID = 0;
        this->npsms = 0;
    }

    _fResult &operator=(const _fResult& rhs)
    {
        /* Check for self assignment */
        if (this != &rhs)
        {
            memcpy(this, &rhs, sizeof(_fResult));
        }

        return *this;
    }

    _fResult &operator=(const int_t& rhs)
    {
        /* Check for self assignment */
        this->eValue = rhs;
        this->specID = rhs;
        this->npsms = rhs;

        return *this;
    }

} fResult;

#define psize          256

typedef struct _ebuffer
{
    char_t *ibuff;
    partRes *packs;
    int_t currptr;
    int_t batchNum;
    BOOL isDone;

    _ebuffer()
    {
        packs = new partRes[QCHUNK];
        ibuff = new char_t[(psize * sizeof (ushort_t)) * QCHUNK];
        currptr = 0;
        batchNum = -1;
        isDone = true;
    }

    ~_ebuffer()
    {
        if (packs != NULL)
        {
            delete[] packs;
            packs = NULL;
        }

        if (ibuff != NULL)
        {
            delete[] ibuff;
            ibuff = NULL;
        }

        currptr = 0;
        batchNum = -1;
        isDone = true;
    }

} ebuffer;

#endif /* INCLUDE_SLM_DSTS_H_ */

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

#ifndef INCLUDE_SLM_DSTS_H_
#define INCLUDE_SLM_DSTS_H_

#include "common.h"

/* Types of modifications allowed by SLM_Mods     */
#define MAX_MOD_TYPES                        15

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
    AA     residues[4]   = ""; /* Modified AA residues in this modification - Upto 4 */
    uint_t  modMass         = 0; /* Scaled mass of the modification                    */
    ushort_t aa_per_peptide = 0; /* Allowed modified residues per peptide sequence     */

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
    ushort_t       vmods_per_pep = 0; /* Total allowed modified residues per sequence */
    ushort_t            num_vars = 0; /* Number of types of modifications added to index - Max: 7 */
    SLM_varAA vmods[MAX_MOD_TYPES]; /* Information for each modification added to index */

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
    AA        *seqs = NULL; /* Stores peptide sequence, could store as strings as well */
    ushort_t   peplen    = 0; /* Stores sequence length */
    uint_t        AAs    = 0; /* Total number of characters */
} PepSeqs;

typedef struct _modAA
{
    ull_t  sites = 0x0; /* maxlen(pep) = 60AA + 2 bits (termini mods)      */
    uint_t  modNum =     0x0; /* 4 bits per mods num, Max 8 mods allowed per pep */

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
    uchar_t   bc  = 0;        /* b ion count */
    uchar_t   yc  = 0;        /* y ion count */
} BYC;

typedef struct _iBYC
{
    uint_t ibc   = 0;
    uint_t iyc   = 0;
} iBYC;

typedef struct _BYICount
{
    BYC     *byc = NULL;       /* Both counts */
    iBYC   *ibyc = NULL;       /* Sum of b/y ion intensities */
} BYICount;

#define BYISIZE                 (sizeof(uchar_t) * 2 + sizeof(uint_t) * 2)

typedef struct _SLMchunk
{
    uint_t    *iA = NULL; /* Ions Array (iA)   */
    uint_t    *bA = NULL; /* Bucket Array (bA) */

#ifdef FUTURE
    uchar_t *bits = NULL; /* Scorecard bits    */
#endif /* FUTURE */
} SLMchunk;

/* Structure for each pep file */
typedef struct _Index
{
    uint_t pepCount         = 0;
    uint_t modCount         = 0;
    uint_t totalCount       = 0;

    uint_t lclpepCnt        = 0;
    uint_t lclmodCnt        = 0;
    uint_t lcltotCnt        = 0;
    uint_t nChunks          = 0;
    uint_t chunksize        = 0;
    uint_t lastchunksize    = 0;

    PepSeqs          pepIndex;
    pepEntry      *pepEntries;
    SLMchunk        *ionIndex;
} index_t;

/* Structure for global Parameters */
typedef struct _globalParams
{
    uint_t threads = 1;
    uint_t min_len = 4;
    uint_t max_len = 40;
    uint_t maxz     = 3;
    uint_t topmatches = 10;
    uint_t scale      = 100;
    uint_t min_shp    = 4;
    uint_t nodes      = 1;
    uint_t myid       = 0;
    uint_t spadmem    = 2048;

    uint_t min_mass = 500.0;
    uint_t max_mass = 5000.0;
    uint_t dF       = (uint_t)(0.02 * scale);

    double_t dM       = 500.0;
    double_t res      = 0.01;

    double_t *perf    = NULL;

    string_t dbpath;
    string_t datapath;
    string_t modconditions;

    DistPolicy policy = _cyclic;

    SLM_vMods vModInfo;

}gParams;


/* Same as specSeqs below but has intensity values for experimental spectra */
typedef struct _eSpecSeqs
{
    uint_t                *moz;       /* Stores the m/z values of the spectra */
    uint_t                *intensity; /* Stores the intensity values of the experimental spectra */
    //BOOL              *iType;     /* Stores the ion type of the coresponding peak in miz */
    uint_t                *idx;       /* Row ptr. Starting index of each row */
    float_t               *precurse;  /* Stores the precursor mass of each spectrum. */
    uint_t                 numPeaks;  /* Total length of moz array i.e. total number of peaks */
    uint_t                 numSpecs;  /* Number of theoretical spectra */
} ESpecSeqs;

/************************* SLM Query DSTs ************************/
typedef struct _Query
{
    /* Raw chosen fragments/peaks from the input MS/MS spectrum   */
    PEAK Peaks[QALEN] = {0};

/* HM: Enable if intensity information is also required in future */
#if (defined(REQUIRE_INTENSITY))
    intensity_t Intensities[QUERYPK] = {0};
#endif /* (defined(REQUIRE_INTENSITY)) */
} Query;

#endif /* INCLUDE_SLM_DSTS_H_ */

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
    UINT  modMass         = 0; /* Scaled mass of the modification                    */
    USHORT aa_per_peptide = 0; /* Allowed modified residues per peptide sequence     */

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
    USHORT       vmods_per_pep = 0; /* Total allowed modified residues per sequence */
    USHORT            num_vars = 0; /* Number of types of modifications added to index - Max: 7 */
    SLM_varAA vmods[MAX_MOD_TYPES]; /* Information for each modification added to index */

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
    AA        *seqs = NULL; /* Stores peptide sequence, could store as strings as well */
    USHORT   peplen    = 0; /* Stores sequence length */
    UINT        AAs    = 0; /* Total number of characters */
} PepSeqs;

typedef struct _modAA
{
    ULONGLONG  sites = 0x0; /* maxlen(pep) = 60AA + 2 bits (termini mods)      */
    UINT  modNum =     0x0; /* 4 bits per mods num, Max 8 mods allowed per pep */

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
typedef struct _BYICount
{
    UCHAR   *bc;        /* b ion count */
    UINT    *ibc;       /* Sum of b ion intensities */
    UCHAR   *yc;        /* y ion count */
    UINT    *iyc;       /* Sum of y ion intensities */
    UINT     size;      /* Length of score card (above declared arrays) */
    INT     especid;   /* Experimental spectrum id */
} BYICount;

typedef struct _SLMchunk
{
    UINT    *iA = NULL; /* Ions Array (iA)   */
    UINT    *bA = NULL; /* Bucket Array (bA) */
//    BYICount sc;
#ifdef FUTURE
    UCHAR *bits = NULL; /* Scorecard bits    */
#endif /* FUTURE */
} SLMchunk;

/* Structure for each pep file */
typedef struct _Index
{
    UINT pepCount         = 0;
    UINT modCount         = 0;
    UINT totalCount       = 0;

    UINT lclpepCnt        = 0;
    UINT lclmodCnt        = 0;
    UINT lcltotCnt        = 0;
    UINT nChunks          = 0;
    UINT chunksize        = 0;
    UINT lastchunksize    = 0;

    PepSeqs          pepIndex;
    pepEntry      *pepEntries;
    SLMchunk        *ionIndex;
} Index;

/* Structure for global Parameters */
typedef struct _globalParams
{
    UINT threads = 1;
    UINT min_len = 4;
    UINT max_len = 40;
    UINT maxz     = 3;
    UINT topmatches = 10;
    UINT scale      = 100;
    UINT min_shp    = 4;
    UINT nodes      = 1;
    UINT myid       = 0;
    UINT spadmem    = 2048;

    UINT min_mass = 500.0;
    UINT max_mass = 5000.0;
    UINT dF       = (UINT)(0.02 * scale);

    DOUBLE dM       = 500.0;
    DOUBLE res      = 0.01;

    DOUBLE *perf    = NULL;

    STRING dbpath;
    STRING datapath;
    STRING modconditions;

    DistPolicy policy = _cyclic;

    SLM_vMods vModInfo;

}gParams;


/* Same as specSeqs below but has intensity values for experimental spectra */
typedef struct _eSpecSeqs
{
    UINT                *moz;       /* Stores the m/z values of the spectra */
    UINT                *intensity; /* Stores the intensity values of the experimental spectra */
    //BOOL              *iType;     /* Stores the ion type of the coresponding peak in miz */
    UINT                *idx;       /* Row ptr. Starting index of each row */
    FLOAT               *precurse;  /* Stores the precursor mass of each spectrum. */
    UINT                 numPeaks;  /* Total length of moz array i.e. total number of peaks */
    UINT                 numSpecs;  /* Number of theoretical spectra */
} ESpecSeqs;

/************************* SLM Query DSTs ************************/
typedef struct _Query
{
    /* Raw chosen fragments/peaks from the input MS/MS spectrum   */
    PEAK Peaks[QALEN] = {0};

/* HM: Enable if intensity information is also required in future */
#if (defined(REQUIRE_INTENSITY))
    INTENSITY Intensities[QUERYPK] = {0};
#endif /* (defined(REQUIRE_INTENSITY)) */
} Query;

#endif /* INCLUDE_SLM_DSTS_H_ */

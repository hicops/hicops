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

#include <thread>
#include "utils.h"
#include "slm_dsts.h"


/* Calculates index in Amino Acid mass array */
#define AAidx(x)                 (x - 'A')

/* Mass of water - Added to each peptide mass */
#define H2O                      18.015f
#define PROTON                   1.00727647

/* Not an Amino Acid (NAA) mass */
#define NAA                      -20000

#ifdef VMODS
/* Global Mods Info  */
SLM_vMods      gModInfo;
#endif /* VMODS */

extern gParams params;

/* Amino Acids Masses */
FLOAT AAMass[26] = {
                    71.03712,   // A
                    NAA,        // B
                    103.00919,  // C
                    115.030,    // D
                    129.0426,   // E
                    147.068,    // F
                    57.02146,   // G
                    137.060,    // H
                    113.084,    // I
                    NAA,        // J
                    128.094,    // K
                    113.084,    // L
                    131.0405,   // M
                    114.043,    // N
                    NAA,        // O
                    97.0527,    // P
                    128.05858,  // Q
                    156.1012,   // R
                    87.032,     // S
                    101.0476,   // T
                    NAA,        // U
                    99.06841,   // V
                    186.0793,   // W
                    NAA,        // X
                    163.0633,   // Y
                    NAA,        // Z
                    };

/* Static Mods for Amino Acids */
FLOAT StatMods[26] = {
                      0,        // A
                      0,        // B
                      57.021464,// C + 57.02
                      0,        // D
                      0,        // E
                      0,        // F
                      0,        // G
                      0,        // H
                      0,        // I
                      0,        // J
                      0,        // K
                      0,        // L
                      0,        // M
                      0,        // N
                      0,        // O
                      0,        // P
                      0,        // Q
                      0,        // R
                      0,        // S
                      0,        // T
                      0,        // U
                      0,        // V
                      0,        // W
                      0,        // X
                      0,        // Y
                      0,        // Z
                      };

/* Macros to extract AA masses */
#define GETAA(x,z)                 ((AAMass[AAidx(x)]) + (StatMods[AAidx(x)]) + ((PROTON) * (z)))


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
UINT UTILS_GetNumProcs()
{
    UINT procs = 1;

#if (_WIN32 == 1)
    std::cout <<"HostOS: Windows\n";
    char var[] = "NUMBER_OF_PROCESSORS";
    char *val = std::getenv(var);
    procs = std::atoi(val);
#else
    std::cout <<"HostOS: Linux\n";
    procs = std::thread::hardware_concurrency();
#endif /* _WIN32 */

    return procs;
}

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
STATUS UTILS_Shuffle(UINT *arr, UINT N)
{
    STATUS status = SLM_SUCCESS;

    /* Check if not already initialized */
    if (arr != NULL)
    {
        /* Create the default seed */
        ULONGLONG defseed = 0xdefceed;

        /* Provide the default seed to UINT_Shuffle */
        return UTILS_ShuffleI(arr, N, defseed);
    }
    else
    {
        status = ERR_INVLD_PTR;
    }

    return status;
}

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
STATUS UTILS_ShuffleI(UINT *arr, UINT N, ULONGLONG seed)
{
    STATUS status = SLM_SUCCESS;
    UINT *indices = arr;

    /* Check if not already initialized */
    if (indices != NULL)
    {
        /* Check if default seed,
         * then create a random seed */
        if (seed == 0xdefceed)
        {
            seed = time(0);
        }

        /* Shuffle the sequence using the provided seed or use time(); */
        std::shuffle(indices, (indices + N), std::default_random_engine(seed));
    }
    else
    {
        status = ERR_INVLD_PTR;
    }

    return status;
}

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
FLOAT UTILS_GenerateSpectrum(CHAR *seq, UINT len, UINT *Spectrum)
{
    /* Calculate Peptide sequences Mass */
    FLOAT mass = UTILS_CalculatePepMass(seq, len);

    UINT maxz = params.maxz;
    UINT scale = params.scale;

    /* If there is a non-AA char, the mass will be -ve */
    /* FIXME: No stupid characters should be allowed in
     *        peptide sequence */
    if (mass > 0)
    {
        /* Set the array to zeros */
        std::memset(Spectrum, 0x0, (sizeof(UINT) * (iSERIES * maxz * (len-1))));

        /* Generate Spectrum */
        for (UINT z = 0; z < maxz; z++)
        {
            /* Indices for b and y series start */
            UINT bstart = z * (len - 1);
            UINT ystart = z * (len - 1) + maxz * (len - 1);

            /* Mass of fragment = [M + (z-1)H]/z */

            /* First b-ion */
            Spectrum[bstart] = (UINT)((GETAA(seq[0], z+1) * scale)/(z+1));
            /* First y-ion */
            Spectrum[ystart] = (UINT)(((GETAA(seq[len-1], z+1) + H2O) * scale)/(z+1));

            /* Loop until length - 1 only */
            for (UINT l = 1; l < len - 1; l++)
            {
                /* Extract b-ions */
                Spectrum[bstart + l] = Spectrum[bstart + (l-1)] +
                                       (UINT)((GETAA(seq[l], 0) * scale)/(z+1));

                /* Extract y-ions */
                Spectrum[ystart + l] = Spectrum[ystart + (l-1)] +
                                       (UINT)(((GETAA(seq[len-1-l], 0)) * scale)/(z+1));
            }
        }
    }

    return mass;
}

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
FLOAT UTILS_CalculatePepMass(AA *seq, UINT len)
{
    /* Initialize mass to H2O */
    FLOAT mass = H2O;

    /* Calculate peptide mass */
    for (UINT l = 0; l < len; l++)
    {
        mass += AAMass[AAidx(seq[l])];
    }

    return mass;
}

#ifdef VMODS

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
STATUS UTILS_InitializeModInfo(SLM_vMods *vMods)
{
    STATUS status = SLM_SUCCESS;

    gModInfo = *vMods;

    return status;
}
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
FLOAT UTILS_CalculateModMass(AA *seq, UINT len, UINT vModInfo)
{
    /* Initialize mass to H2O */
    FLOAT mass = H2O;

    /* Calculate peptide mass */
    for (UINT l = 0; l < len; l++)
    {
        mass += AAMass[AAidx(seq[l])];
    }

    /* Add the mass of modifications present in the peptide */
    UINT start = 0x0F;
    UINT modNum = vModInfo & start;

    while (modNum != 0)
    {
        mass += ((FLOAT)(gModInfo.vmods[modNum - 1].modMass)/params.scale);
        start = (start << 4);
        modNum = ((vModInfo & start) >> start);
    }


    return mass;
}

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
FLOAT UTILS_GenerateModSpectrum(CHAR *seq, UINT len, UINT *Spectrum, modAA modInfo)
{
    /* Calculate Mod Peptide Mass */
    STATUS status = SLM_SUCCESS;
    FLOAT mass = 0;

    DOUBLE minmass = params.min_mass;
    DOUBLE maxmass = params.max_mass;
    UINT maxz = params.maxz;
    UINT scale = params.scale;


    CHAR modPos[MAX_SEQ_LEN] = {};
    INT modNums[MAX_MOD_TYPES] = {};
    INT modSeen = 0;

    /* Check if valid modInfo */
    if (modInfo.sites == 0 || modInfo.modNum == 0)
    {
        status = ERR_INVLD_MOD;
        mass = NAA;
    }

    /* Compute Mod Mass */
    if (status == SLM_SUCCESS)
    {
        mass = UTILS_CalculateModMass(seq, len, modInfo.modNum);
    }

    /* Check if a valid precursor mass */
    if (mass > minmass && mass < maxmass)
    {
        for (UINT i = 0; i < MAX_MOD_TYPES; i++)
        {
            modNums[i] = ((modInfo.modNum & (0x0F << (4 * i))) >> (4 * i)) - 1; // -1 to store index instead
            if (modNums[i] != -1)
            {
                modSeen++;
            }
        }

        for (UINT i = 0; i < MAX_SEQ_LEN; i++)
        {
            modPos[i] = ISBITSET(modInfo.sites,i) ? 1 : 0;
        }

        if (mass > 0)
        {
            /* Set the array to zeros */
            //std::memset(Spectrum, 0x0, (sizeof(UINT) * (iSERIES * MAXz * (len-1))));

            /* Generate Normal Spectrum */
            for (UINT z = 0; z < maxz; z++)
            {
                /* Indices for b and y series start */
                UINT bstart = z * (len - 1);
                UINT ystart = z * (len -1)  + maxz * (len - 1);

                /* Mass of fragment = [M + (z-1)H]/z */

                /* First b-ion */
                Spectrum[bstart] = (UINT)((GETAA(seq[0], z+1) * scale));
                /* First y-ion */
                Spectrum[ystart] = (UINT)(((GETAA(seq[len-1], z+1) + H2O) * scale));

                /* Loop until length - 1 only */
                for (UINT l = 1; l < len - 1; l++)
                {
                    /* Extract b-ions */
                    Spectrum[bstart + l] = Spectrum[bstart + (l - 1)] +
                                           (UINT)((GETAA(seq[l], 0) * scale));

                    /* Extract y-ions */
                    Spectrum[ystart + l] = Spectrum[ystart + (l - 1)] +
                                           (UINT)(((GETAA(seq[len-1-l], 0)) * scale));
                }
            }

            /* Adjust b-ions with additional masses */
            for (UINT z = 0; z < maxz; z++)
            {
                /* Indices for b series start */
                UINT bstart = z * (len-1);
                UINT counter = 0;

                /* Loop until length - 1 only */
                for (UINT l = 0; l < len - 1; l++)
                {
                    counter += modPos[l];

                    for (UINT k = 0; k < counter; k++)
                    {
                        Spectrum[bstart + l] += gModInfo.vmods[modNums[k]].modMass;
                    }

                    /* Divide by charge here */
                    Spectrum[bstart + l] /= (z + 1);
                }
            }

            /*  Adjust y-ions with additional masses */
            for (UINT z = 0; z < maxz; z++)
            {
                /* Indices for y series start */
                UINT ystart = z * (len -1) + maxz * (len -1);
                UINT counter = 0;

                /* Loop until length - 1 only */
                for (INT l = (len - 1); l > 0; l--)
                {
                    counter += modPos[l];

                    for (UINT k = 0; k < counter; k++)
                    {
                        Spectrum[ystart + (len - 1) - l] += gModInfo.vmods[modNums[modSeen - 1 - k]].modMass;
                    }

                    /* Divide by charge here */
                    Spectrum[ystart + (len - 1) - l] /= (z + 1);
                }
            }
        }
    }

    return mass;
}

#endif /* VMODS */

//****************************************************************************
//
//  Source: https://people.sc.fsu.edu/~jburkardt/cpp_src/llsq/llsq.html
//
//
//  Purpose:
//
//    LLSQ solves a linear least squares problem matching a line to data.
//
//  Discussion:
//
//    A formula for a line of the form Y = A * X + B is sought, which
//    will minimize the root-mean-square error to N data points ( X[I], Y[I] );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, double X[N], Y[N], the coordinates of the data points.
//
//    Output, double &A, &B, the slope and Y-intercept of the least-squares
//    approximate to the data.
//

//****************************************************************************

void UTILS_LinearRegression(INT n, DOUBLE x[], DOUBLE y[], DOUBLE &a, DOUBLE &b)
{
    INT i;
    DOUBLE bot;
    DOUBLE top;
    DOUBLE xbar;
    DOUBLE ybar;
//
//  Special case.
//
    if (n == 1)
    {
        a = 0.0;
        b = y[0];
        return;
    }
//
//  Average X and Y.
//
    xbar = 0.0;
    ybar = 0.0;

    for (i = 0; i < n; i++)
    {
        xbar = xbar + x[i];
        ybar = ybar + y[i];
    }

    xbar = xbar / (DOUBLE) n;
    ybar = ybar / (DOUBLE) n;
//
//  Compute Beta.
//
    top = 0.0;
    bot = 0.0;
    for (i = 0; i < n; i++)
    {
        top = top + (x[i] - xbar) * (y[i] - ybar);
        bot = bot + (x[i] - xbar) * (x[i] - xbar);
    }
    a = top / bot;

    b = ybar - a * xbar;

    return;
}

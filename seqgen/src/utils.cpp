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
UINT UTILS_GetNumProcs(VOID)
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
        mass += (gModInfo.vmods[modNum - 1].modMass/SCALE);
        start = (start << 4);
        modNum = ((vModInfo & start) >> start);
    }


    return mass;
}

#endif /* VMODS */

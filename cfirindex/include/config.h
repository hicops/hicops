/*
 *  This file is part of SLM-Transform
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

#ifndef INCLUDE_CONFIG_H_
#define INCLUDE_CONFIG_H_

/************** Developer Mode options: Danger Area *****************/


/* Contains MACROS that configure SLM Index build */

/* Manually switch _OPENMP */
#undef _OPENMP

/* Enable to show SLMIndex Benchmark information  */
#undef BENCHMARK

#define MAX_SEQ_LEN              60

/* Number of fragments (peaks) per query spectrum */
#define QALEN                    100

/* How many query spectra per chunk */
#define QCHUNK                   50000

/* How many top matches required at max */
#define MAXMATCHES               99999

/* Enable DSLIM Mods construction                 */
#define VMODS

/* Chunk size in number ions per chunk (def: 2 bn)   */
#define MAX_IONS                0x80000000

/* Enable debug checkpoints                       */
#undef DEBUG

/* Recycle and reuse the already allocated memory
 * instead of new allocations when possible       */
#undef RECYCLE

/* Future: Define if the peak intensities
 * are required for the query spectra             */
#undef REQUIRE_INTENSITY

/* TODO:  Only b- and -y-ions used in current
 *               implementation                   */

/* Define the Ion Series and add it to iSERIES    */
#define aIONS                     0
#define bIONS                     1
#define cIONS                     0
#define xIONS                     0
#define yIONS                     1
#define zIONS                     0
#define NLOSS                     0
#define IMMONIUM                  0

/* Number of ion series to generate per spectrum  */
#define iSERIES                   (aIONS + \
                                   bIONS + \
                                   cIONS + \
                                   xIONS + \
                                   yIONS + \
                                   zIONS + \
                                   NLOSS + \
                                   IMMONIUM)

/*******************************************************************/

#endif /* INCLUDE_CONFIG_H_ */

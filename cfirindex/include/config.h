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

#ifndef INCLUDE_CONFIG_H_
#define INCLUDE_CONFIG_H_

/************** Developer Mode options: Danger Area *****************/

/* Contains MACROS that configure SLM Index build */

/* Manually switch _OPENMP */
//#undef _OPENMP

/* Manually switch MPI */
//#define DISTMEM

/* Switch the benchmarking */
//#define BENCHMARK

#define DIAGNOSE

#define MAX_SEQ_LEN              60

/* Number of fragments (peaks) per query spectrum */
#define QALEN                    100

/* How many query spectra per batch */
#define QCHUNK                   10000

/* Precision to inter-convert between FLOAT and INT
 * for Tx/Rx purposes
 */
#define PRECISION                1000000

/* Max RX buffer size = 512MB */
#define  RXBUFFERSIZE              (512 * 1024 * 1024)

/* Maximum hyperscore that is possible */
#define MAX_HYPERSCORE           100

/* Max number of RX instances before buffer flush */
#define MAX_RX_INST              500

/* Enable DSLIM Mods construction                 */
#define VMODS

/* Chunk size in number ions per chunk (def: 2 bn)   */
#define MAX_IONS                0x80000000

/* Enable debug checkpoints                       */
#undef DEBUG

/* Define to analyze the distributed partial results
 * instead of dumping them into file
 */
#undef ANALYSIS

/* FUTURE:  Only b- and y-ions used in current
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

/* Number of ion series generated per spectrum  */
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

/*
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

#pragma once

// include the CMake configured macros
#include <config.hpp>

//
// Universal Constants: Do not change them
//

// -------------------------------------------------------------------------------- //

/* Max RX buffer size = 512MB */
#define  RXBUFFERSIZE              (512 * 1024 * 1024)

/* Chunk size in number ions per chunk (def: 2 bn)   */
#define MAX_IONS                0x80000000

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

// -------------------------------------------------------------------------------- //
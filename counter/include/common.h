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

#ifndef INCLUDE_COMMON_H_
#define INCLUDE_COMMON_H_

#include <string>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iomanip>

#include "config.h"
#include "slmerr.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* Set/Get the bit for modification site in a peptide */
#define MODSITE(x)                       (1 << (x))

/* Handle unusued parameters to avoid compiler warnings */
#define LBE_UNUSED_PARAM(x)              ((void *)(&x))

#define MIN(x,y)                         ((x < y)? x: y)
#define MAX(x,y)                         ((x > y)? x: y)

/* typedef the data types */
typedef int STATUS;
typedef int INT;
typedef long long LONGLONG;
typedef short SHORT;
typedef long LONG;
typedef char  CHAR;
typedef void VOID;
typedef bool BOOL;
typedef double DOUBLE;
typedef std::string STRING;

/* Unsigned data types */
typedef unsigned int UINT;
typedef unsigned int IDX;
typedef unsigned int ID;
typedef unsigned long ULONG;
typedef unsigned long long ULONGLONG;
typedef unsigned short USHORT;
typedef unsigned char UCHAR;
typedef float FLOAT;

typedef unsigned int CHG;
typedef unsigned short* SLM_SC;

typedef char  AA;         // Amino Acid
typedef double PEAK;      // m/z of a fragement (peak)
typedef double INTENSITY; // Intensity of a fragment (peak)
typedef double* SPECTRUM; // Theoretical Spectrum

#endif /* INCLUDE_COMMON_H_ */

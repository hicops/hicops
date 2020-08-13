/*
 * Copyright (C) 2019  Muhammad Haseeb, and Fahad Saeed
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
#include <pthread.h>
#include <semaphore.h>

#include "config.h"
#include "slmerr.h"

#ifdef USE_OMP
#include <omp.h>
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif /* USE_MPI */

/* Set/Get the bit for modification site in a peptide */
#define MODSITE(x)                       (1 << (x))

/* Handle unusued parameters to avoid compiler warnings */
#define LBE_UNUSED_PARAM(x)              ((void *)(&x))

#define MIN(x,y)                         ((x < y)? x: y)
#define MAX(x,y)                         ((x > y)? x: y)
#define True                             true
#define False                            false

/* alias the data types */
using status_t = int;
using int_t = int;
using longlong_t = long long;
using short_t = short;
using LONG = long;
using char_t = char;
using VOID = void;
using BOOL = bool;
using double_t = double;
using string_t = std::string;

/* Unsigned data types */
using uint_t = unsigned int;
using IDX  = unsigned int;
using ID =  unsigned int;
using ulong_t  = unsigned long;
using ull_t = unsigned long long;
using ushort_t = unsigned short;
using uchar_t = unsigned char;
using float_t  = float;

using CHG = unsigned int;
using SLM_SC = unsigned short*;

using AA = char;          // Amino Acid
using PEAK = double;      // m/z of a fragement (peak)
using intensity_t = double; // Intensity of a fragment (peak)
using SPECTRUM = double*; // Theoretical Spectrum

using  thread_t = pthread_t;
using  lock_t = sem_t;

#endif /* INCLUDE_COMMON_H_ */

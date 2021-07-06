/*
 * From : https://raw.githubusercontent.com/thatchristoph/vmd-cvs-github/master/plugins/signalproc/src/sgsmooth.C
 *
 * Sliding window signal processing (and linear algebra toolkit).
 *
 * supported operations:
 * <ul>
 * <li> Savitzky-Golay smoothing.
 * <li> computing a numerical derivative based of Savitzky-Golay smoothing.
 * <li> required linear algebra support for SG smoothing using STL based
 *      vector/matrix classes
 * </ul>
 *
 * \brief Linear Algebra "Toolkit".
 *
 * modified by Rob Patro, 2016
 * Also modified by Muhammad Haseeb, 2020 (https://GitHub.com/mhaseeb123)
 *
 */

#pragma once

#include <vector>
#include "common.hpp"
#include "lwvector.h"

// savitzky golay smoothing.
void sg_smooth(lwvector<double> *v, lwvector<double> *ptr, const int w, const int deg);

//! numerical derivative based on savitzky golay smoothing.
std::vector<double> sg_derivative(const std::vector<double> &v, const int w,
                                const int deg, const double h=1.0);

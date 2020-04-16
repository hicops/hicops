#ifndef __SGSMOOTH_HPP__
#define __SGSMOOTH_HPP__

#include <vector>
#include "lwvector.h"

// savitzky golay smoothing.
void sg_smooth(lwvector<double> *v, lwvector<double> *ptr, const int w, const int deg);

//! numerical derivative based on savitzky golay smoothing.
std::vector<double> sg_derivative(const std::vector<double> &v, const int w,
                                const int deg, const double h=1.0);

#endif // __SGSMOOTH_HPP__

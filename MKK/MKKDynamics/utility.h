// Some utility functions

#ifndef  UTILITY_H
#define  UTILITY_H

#include <chrono>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdint.h>

#include "dim.h"

// MACROS
#define VOLUMESPHERE pow(M_PI,((double)(DIM))/2.)/exp(lgamma(1+((double)(DIM))/2.)) // volume prefactor for sphere
#define DBL_EPSILON  2.2204460492503131e-016 // smallest # such that 1.0+DBL_EPSILON!=1.0
#define DBL_LARGE 1.0e8

using clock_type = std::chrono::time_point<std::chrono::steady_clock>;
using duration_type = std::chrono::duration<double>;

// 1-cubic; 2-Dn
#define BOX_TYPE 2

// random number
extern std::mt19937_64 rg;
extern std::uniform_real_distribution <> rrand;
// to set seed: rg.seed(long seed);
// to get random double: rrand(rg);

void setminimg(double x[DIM]);

// get diameter interval (1-x) - (1+x) from polydispersity K
double diameterRange(double polydispersity);
double averagesprvolume(double rx);
double gerRfrompfrac(double vspheres, double rx, int N);
// minimum image under Dn periodic boundary condition
void DnLatticePoint(double p[DIM], int32_t mirror[DIM]);
#endif

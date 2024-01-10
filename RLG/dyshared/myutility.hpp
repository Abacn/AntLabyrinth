//
//  myutility.hpp
//  RLGDynamics
//
//  Created by Yi Hu on 6/13/19.
//

#ifndef myutility_hpp
#define myutility_hpp

#include <random>
#include "commondef.hpp"

typedef void(*PBCFunc)(Vector &);

/** Some common utility functions */
namespace MyUtility{
// global random number generator
extern std::mt19937 rg;
// sgn function
template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}
// initialize the random number generator
void init(unsigned int seed);
// Volume of unit sphere at the dimention d (d is defined by DIM macro)
double getVunitsphere();
// uniform distributed random number in [0, 1)
double rand();
// normal distribution random number
double randn();
// phi squared distribution random number
double randc();
double norm_squared(const point &pi, const point &pj);
double norm_squared(const Vector &vec);
// square equation root solver, return the smaller root (when A>0)
double rootsolverA(double A, double B, double C);
// square equation root solver, return the larger root (when A>0)
double rootsolverB(double A, double B, double C);
void pbc(Vector &vec);
};
#endif /* myutility_hpp */

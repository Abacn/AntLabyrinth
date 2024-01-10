//
//  myutility.cpp
//  RLGDynamics
//
//  Created by Yi Hu on 6/13/19.
//

#include <cmath>
#include <random>
#include "myutility.hpp"

namespace{
  std::uniform_real_distribution <> rrand; // uniform distribution in [0, 1)
  std::normal_distribution<> rnorm;        // normal distribution device
  std::chi_squared_distribution<> rchisqr((double)DIM); // chi-squared distribution (for v^2)
  double Vunitshpere;
}

namespace MyUtility{
  std::mt19937 rg;  // random generator
  void init(unsigned int seed)
  {
    // initialize random seed
    rg.seed(seed);
    // volume of unit sphere
    Vunitshpere = pow(M_PI, DIM*0.5)/tgamma(1.0+DIM*0.5);
  }
  
  double getVunitsphere()
  {
    return Vunitshpere;
  }
  
  double rand()
  {
    return rrand(rg);
  }
  
  double randn()
  {
    return rnorm(rg);
  }

  double randc()
  {
    return rchisqr(rg);
  }
  
  // L2 norm
  double norm_squared(const Vector &pi, const Vector &pj)
  {
    double dtmp, result = 0.;
    for(int i=0; i<DIM; ++i)
    {
      dtmp = pi[i] - pj[i];
      result += dtmp*dtmp;
    }
    return result;
  }
  
  double norm_squared(const Vector &vec)
  {
    return vec.norm_squared();
  }
  
  // square equation root solver, return the smaller root (when A>0)
  double rootsolverA(double A, double B, double C)
  {
    double dt = B*B - 4.*A*C;
    if(dt < 0.) return NAN;
    else return -0.5*(B+sqrt(dt))/A;
  }
  
  // square equation root solver, return the larger root (when A>0)
  double rootsolverB(double A, double B, double C)
  {
    double dt = B*B - 4.*A*C;
    if(dt < 0.) return NAN;
    else return 0.5*(sqrt(dt)-B)/A;
  }
};

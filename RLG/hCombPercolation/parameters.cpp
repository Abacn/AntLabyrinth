//
//  parameters.cpp
//  VoidPercolation
//
// Global parameters that have to be used cross the classes

#include <cstdint>
#include <cmath>

#include "parameters.hpp"

namespace Parameter{
  const double BoxLength = 1.0;
  const double HalfBoxL = 0.5*BoxLength;
  const double PhiCutD[]
  = {0, 0, 1.128,    // 0-2D
    3.51, 6.249, 9.17,   // 3-5D
    12.22, 15.4, 18.6,  // 6-8D
    22., 25.2, 30     // 9-11D
  };
  const double slopes[]
  = {
    0., 0., -3., // 0-2D
    -8., -15.1, -20., // 3-5D
    -26.9, -24.7, -15., // 6-8D
    -16., -15., -15.  // 9-11D
  };
  
  const double corrlengths[] // nu
  = { 1., 1.5, 1.33,
    0.8774, 0.6852, 0.5723
  };
  const unsigned nPhiCutD = 11;  // check bound
  
  double getPhiCut(uint64_t N)
  {
    double targetphi, targetslp, nu;
    if(DIM<nPhiCutD)
    {
      targetphi = PhiCutD[DIM];
      targetslp = slopes[DIM];
    }
    else
    {
      targetphi = 2.7*DIM;
      targetslp = -3.2*DIM;
    }
    if(DIM<6)
    {
      nu = corrlengths[DIM];
    }
    else
    {
      nu = 0.5;
    }
    return targetphi*0.99 + targetslp*pow(N, -1.0/(DIM*nu)); // 99.5% percentile
  }
};

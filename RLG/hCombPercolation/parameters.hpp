//
//  parameters.hpp
//  VoidPercolation
//
//  Created by Yi Hu on 10/6/18.

#ifndef parameters_hpp
#define parameters_hpp

#include <cstdint>
#include "dim.h"

namespace Parameter{
  extern const double BoxLength;
  extern const double HalfBoxL;
  extern const double PhiCutD[];
  extern const unsigned nPhiCutD;
  double getPhiCut(uint64_t N);
};

#endif /* parameters_hpp */

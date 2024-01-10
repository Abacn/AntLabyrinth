//
//  hcombpbc.hpp
//  RLGDynamics
//
//  Created by Yi Hu on 11/4/19.
//

#ifndef hcombpbc_hpp
#define hcombpbc_hpp


#include <cstdint>

#include "commondef.hpp"

#if DIM>=8 && DIM%2==0
// Test shows that Dn+ quantizer costs twice of the time for Dn.
// Consider that the system size can be reduced by factor of 2, which makes no difference
// Conclusion: not necessary to implement Dn+ instead of Dn in d=8,10,12,...
// #define DNPLUSFLAG
#endif

#ifndef DNPLUSFLAG
#if DIM==24
#define HCOMBPBC_IS_LEECH
#endif
#endif

namespace HCombPBC{
  // periodic boundary condition
  // D_n symmetry
#ifdef DNPLUSFLAG
  int latticepoint(double p[DIM], int32_t mirror[DIM]);
  void getpbcvec(double pbcvec[DIM], int32_t mirror[DIM], int hfflag);
#else
  // Take a dD coords p, put the nearest lattice point in mirror,
  // also subtract mirror from p.
  void latticepoint(double p[DIM], int32_t mirror[DIM]);
  void getpbcvec(double pbcshift[DIM], int32_t mirror[DIM]);
#endif
  extern const int Vcell;
  void pbc(Vector &vec);
  double norm_squared_PBC(const Vector &va, const Vector &vb);
};
#endif /* hcombpbc_hpp */

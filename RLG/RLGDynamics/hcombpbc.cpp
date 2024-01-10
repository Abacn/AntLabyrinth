//
//  hcombpbc.cpp
//  RLGDynamics
//
//  Created by Yi Hu on 11/4/19.
//

#include "hcombpbc.hpp"

#ifdef HCOMBPBC_IS_LEECH
extern "C" void quantizer_L24(const double *t, int32_t *lpoint);

constexpr double SCALE_FACTOR = 2 * M_SQRT2;
#endif

/**
 * Utility for Dn and related periodic conditions.
 *
 * If macro DNPLUSFLAG is not defined, Dn PBC is compiled.
 * If macro DNPLUSFLAG is defined, Dn+/Dn0+ PBC (D8+=E8, D90+=L9, dense packing
 * in d=8,9) is compiled; in d=24, the Leech lattice PBC is compiled.
 */
namespace HCombPBC{
  // periodic boundary condition
  // D_n symmetry
#ifdef DNPLUSFLAG
  void latticepoint_inner(double p[DIM], int32_t mirror[DIM])
#else
  void latticepoint(double p[DIM], int32_t mirror[DIM])
#endif
  {
#ifndef HCOMBPBC_IS_LEECH
    int rp, gind = 0, summr = 0;
    double maxdelta = 0, dtmp;
    for(rp=0; rp<DIM; ++rp)
    {
      // fast round eqv mirror[rp] = round(p[rp]);
      dtmp = p[rp] + 6755399441055744.0;
      mirror[rp] = reinterpret_cast<int32_t&>(dtmp);
      summr += mirror[rp];
      p[rp] -= mirror[rp];
    }
    if(summr%2 != 0)
    {
      for(rp=0; rp<DIM; ++rp)
      {
        dtmp = fabs(p[rp]);
        if(maxdelta < dtmp)
        {
          maxdelta = dtmp;
          gind = rp;
        }
      }
      // use next nearest neighbor image
      if(p[gind] < 0)
      {
        --mirror[gind];
        ++p[gind];
      }
      else
      {
        ++mirror[gind];
        --p[gind];
      }
    }
#else
    // leech lattice
    int rp;
    double in[DIM];
    // scale by a factor of 4 to get integer lattice points
    for(rp=0; rp<DIM; ++rp)
    {
      in[rp] = p[rp] * SCALE_FACTOR;
    }
    quantizer_L24(in, mirror);
    for(rp=0; rp<DIM; ++rp)
    {
      p[rp] = (in[rp] - mirror[rp]) / SCALE_FACTOR;
    }
#endif
  }

#ifdef DNPLUSFLAG
  int latticepoint(double p[DIM], int32_t mirror[DIM])
  {
    int rp, rst = 0;
    double norm1=0., norm2 = 0.;
    double pb[DIM];
    int32_t mirrorb[DIM];
    latticepoint_inner(p, mirror);
    for(rp=0; rp<DIM; ++rp)
    {
      pb[rp] = p[rp] - 0.5;
    }
    latticepoint_inner(pb, mirrorb);
    for(rp=0; rp<DIM; ++rp)
    {
      norm1 += p[rp]*p[rp];
      norm2 += pb[rp]*pb[rp];
    }
    if(norm2<norm1)
    {
      // copy pb to p
      for(rp=0; rp<DIM; ++rp)
      {
        p[rp] = pb[rp];
        mirror[rp] = mirrorb[rp];
      }
      rst = 2; // shift (1/2)
    }
    return rst;
  }
  
  void getpbcvec(double pbcshift[DIM], int32_t mirror[DIM], int hfflag)
  {
    if(2 == hfflag)
    {
      for(int rp=0; rp<DIM; ++rp)
      {
        pbcshift[rp] = mirror[rp] + 0.5;
      }
    }
    else
    {
      for(int rp=0; rp<DIM; ++rp)
      {
        pbcshift[rp] = mirror[rp];
      }
    }
  }
#else

  void getpbcvec(double pbcshift[DIM], int32_t mirror[DIM])
  {
    for(int rp=0; rp<DIM; ++rp)
    {
#ifdef HCOMBPBC_IS_LEECH
      pbcshift[rp] = mirror[rp] / SCALE_FACTOR;
#else
      pbcshift[rp] = mirror[rp];
#endif // HCOMBPBC_IS_LEECH
    }
  }
#endif // DNPLUSFLAG

#if defined(DNPLUSFLAG) || defined(HCOMBPBC_IS_LEECH)
  const int Vcell = 1; // Dn+ or leech lattice
#else
  const int Vcell = 2; // Dn lattice
#endif

  void pbc(Vector &vec)
  {
    int32_t mirror[DIM];
    latticepoint(vec.x, mirror);
  }

  // norm under periodic boundary condition
  double norm_squared_PBC(const Vector &va, const Vector &vb)
  {
    Vector vec_ = va - vb;
    pbc(vec_);
    return vec_.norm_squared();
  }
};

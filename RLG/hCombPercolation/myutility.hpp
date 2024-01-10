//
//  myutility.hpp
//  QhullPercolation
//
//  Created by Yi Hu on 3/21/19.

#ifndef myutility_hpp
#define myutility_hpp

#include <list>
#include "commondef.hpp"

namespace MyUtility
{
  // get the lattice point of p, put result into mirror
  // return type: void
  typedef int mirrorvec_type;
#ifdef DNPLUSFLAG
  // Dn+
#define DNVCELL 1.0
  template <class T>
  void latticepoint_inner(T &p, mirrorvec_type mirror[DIM])
#else
  // Dn
#define DNVCELL 2.0
  template <class T>
  void latticepoint(T &p, mirrorvec_type mirror[DIM])
#endif
  {
    int rp, gind = 0, summr = 0;
    double maxdelta = 0, dtmp;
    for(rp=0; rp<DIM; ++rp)
    {
      // fast round eqv mirror[rp] = round(p[rp]);
      dtmp = p[rp] + 6755399441055744.0;
      mirror[rp] = reinterpret_cast<int32_t&>(dtmp);
      p[rp] -= mirror[rp];
      dtmp = fabs(p[rp]);
      if(maxdelta < dtmp)
      {
        maxdelta = dtmp;
        gind = rp;
      }
      summr += mirror[rp];
    }
    if(summr%2 != 0)
    {
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
  }
  
#ifdef DNPLUSFLAG
  // return value: 0 - is nearest mirror;
  // 2 - nearest mirror differs by (1/2)+mirror
  template <class T>
  int latticepoint(T &p, mirrorvec_type mirror[DIM])
  {
    int rp, rst=0;
    double norm1=0., norm2 = 0.;
    double pb[DIM];
    mirrorvec_type mirrorb[DIM];
#ifndef DNLAMFLAG
    // even dimension
    for(rp=0; rp<DIM; ++rp)
#else
    // odd dimension
    pb[DIM-1] = p[DIM-1];
    for(rp=0; rp<DIM-1; ++rp)
#endif
    {
      pb[rp] = p[rp] - 0.5;
    }
    latticepoint_inner<T>(p, mirror);
    latticepoint_inner<double[DIM]>(pb, mirrorb);
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
#endif

// will change the value in p but the value useless. Parsing the regerence for performance only
#ifndef DNPLUSFLAG
  template <class T>
  double norm_squared_distance_PBC(T &p)
  {
    mirrorvec_type mirror[DIM];
    int rp;
    double norm=0.;
    latticepoint<T>(p, mirror);
    for(rp=0; rp<DIM; ++rp)
    {
      norm += p[rp]*p[rp];
    }
    return norm;
  }
#else
  template <class T>
  double norm_squared_distance_PBC(T &p)
  {
    mirrorvec_type mirror[DIM];
    int rp;
    double norm1=0., norm2 = 0.;
    double pb[DIM];
    latticepoint_inner<T>(p, mirror);
#ifndef DNLAMFLAG
    // even dimension
    for(rp=0; rp<DIM; ++rp)
#else
    // odd dimension
    pb[DIM-1] = p[DIM-1];
    for(rp=0; rp<DIM-1; ++rp)
#endif
    {
      pb[rp] = p[rp] - 0.5;
    }
    latticepoint_inner<double[DIM]>(pb, mirror);
    for(rp=0; rp<DIM; ++rp)
    {
      norm1 += p[rp]*p[rp];
      norm2 += pb[rp]*pb[rp];
    }
    if(norm2<norm1)
    {
      return norm2;
    }
    else
    {
      return norm1;
    }
  }
#endif

  double norm_squared_distance(const point &pa);
  double norm_squared_distance(const Vector &pa);
  double norm_squared_distance(const point &pa, const point &pb);
  double norm_squared_distance(const Vector &pa, const Vector &pb);
  double norm_squared_distance_PBC(const point &pa, const point &pb);
  
  /** get all nearest neighbor copies within the shell of size sqrcut
   */
  std::list<std::pair<point, double> > get_point_images(point &pa, double sqrcut);
  
  void pbc(point &pa);
  int pbc_is_incell(point &pa);
  point &point_add(point &pa, const point &pb);
  Vector dif_points(const point &point1, const point &point2);
  int latticepoint(Vector &p, int mirror[DIM]);
  Vector dif_vec_PBC(Vector &vec);
  Vector dif_vec_PBC(const Vector &vec);
  Vector dif_points_PBC(const point &point1, const point &point2);
  Vector point2vec(const point &p);
  template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
  }
  void setseed(unsigned int seed);
  double rand();
  double randn();
  void gauss_eliminate(double *a, double *b, double *x, const int n);
}

#endif /* myutility_hpp */

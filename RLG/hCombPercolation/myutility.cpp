//
//  myutility.cpp
//  QhullPercolation
//
//  Created by Yi Hu on 3/21/19.

#include <cmath>
#include <random>
#include "dim.h"
#include "myutility.hpp"
#include "parameters.hpp"

// module internal variables and methods
namespace{
  std::mt19937 rg;  // random generator
  std::uniform_real_distribution <> rrand; // uniform distribution in (0, 1]
  std::normal_distribution<> rnorm;

  /** return all nearest neghbor copies (shifted by (±1,±1,0,...)) of pa within shell sqrcut
   */
  std::list<std::pair<point, double> > get_point_images_dn(const point &pa, const std::vector<int> &idx, double nnnorm, double sqrcut)
  {
    double incp, incq, n2norm, n3norm;
    std::list<std::pair<point, double> > results;
    int rp, rq;
    // (1, 1, 0, 0, ...)
    for(rp=0; rp<DIM-1; ++rp)
    {
      incp = 1.0 - 2*fabs(pa[idx[rp]]);
      n2norm = nnnorm+incp;
      for(rq=rp+1; rq<DIM; ++rq)
      {
        incq = 1.0 - 2*fabs(pa[idx[rq]]);
        n3norm = n2norm+incq;
        if(n3norm <= sqrcut)
        {
          point pb = pa;
          pb[idx[rp]] -= copysign(1.0, pb[idx[rp]]);
          pb[idx[rq]] -= copysign(1.0, pb[idx[rq]]);
          results.push_back({pb, n3norm});
        }
        else
        {
          break;
        }
      }
      if(rq==rp+1)
      {
        break;
      }
    }
    if(sqrcut>=1.0)
    {
      int rr, rs;
      double incr, incs, n4norm, n5norm;
      // (1, 0, 0, ..., -1)
      for(rp=DIM-1; rp>0; --rp)
      {
        incp = 1.0 + 2*fabs(pa[idx[rp]]);
        n2norm = nnnorm+incp;
        if(n2norm>sqrcut) break;
        for(rq=0; rq<rp; ++rq)
        {
          incq = 1.0 - 2*fabs(pa[idx[rq]]);
          n3norm = n2norm+incq;
          if(n3norm>sqrcut) break;
          point pb = pa;
          pb[idx[rp]] += copysign(1.0, pb[idx[rp]]);
          pb[idx[rq]] -= copysign(1.0, pb[idx[rq]]);
          results.push_back({pb, n3norm});
        }
        if(rq==0) break;
      }
      // (2, 0, 0, ...)
      for(rp=0; rp<DIM; ++rp)
      {
        incp = 4.0 - 4*fabs(pa[idx[rp]]);
        n2norm = nnnorm+incp;
        if(n2norm>sqrcut) break;
        point pb = pa;
        pb[idx[rp]] -= copysign(2.0, pb[idx[rp]]);
        results.push_back({pb, n2norm});
      }
      // (1, 1, 1, 1, 0, ...)
      for(rp=0; rp<DIM-3; ++rp)
      {
        incp = 1.0 - 2*fabs(pa[idx[rp]]);
        n2norm = nnnorm+incp;
        if(n2norm>sqrcut) break;
        for(rq=rp+1; rq<DIM-2; ++rq)
        {
          incq = 1.0 - 2*fabs(pa[idx[rq]]);
          n3norm = n2norm+incq;
          if(n3norm>sqrcut) break;
          for(rr=rq+1; rr<DIM-1; ++rr)
          {
            incr = 1.0 - 2*fabs(pa[idx[rr]]);
            n4norm = n3norm+incr;
            if(n4norm>sqrcut) break;
            for(rs=rr+1; rs<DIM; ++rs)
            {
              incs = 1.0 - 2*fabs(pa[idx[rs]]);
              n5norm = n4norm+incs;
              if(n5norm>sqrcut) break;
              // add
              point pb = pa;
              pb[idx[rp]] -= copysign(1.0, pb[idx[rp]]);
              pb[idx[rq]] -= copysign(1.0, pb[idx[rq]]);
              pb[idx[rr]] -= copysign(1.0, pb[idx[rr]]);
              pb[idx[rs]] -= copysign(1.0, pb[idx[rs]]);
              results.push_back({pb, n5norm});
            }
            if(rs==rr+1) break;
          }
          if(rr==rq+1) break;
        }
        if(rq==rp+1) break;
      }
    }
    return results;
  }
}

namespace MyUtility
{
  // distance to the origin
  double norm_squared_distance(const point &pa)
  {
    double result = 0.0, dtmp;
    for(int rp=0; rp<DIM; ++rp)
    {
      dtmp = pa[rp];
      result += dtmp*dtmp;
    }
    return result;
  }
  
  double norm_squared_distance(const Vector &pa)
  {
    double result = 0.0, dtmp;
    for(int rp=0; rp<DIM; ++rp)
    {
      dtmp = pa[rp];
      result += dtmp*dtmp;
    }
    return result;
  }
  
  // distance of two points
  double norm_squared_distance(const point &pa, const point &pb)
  {
    double result = 0.0, dtmp;
    int rp;
    for(rp=0; rp<DIM; ++rp)
    {
      dtmp = pa[rp] - pb[rp];
      result += dtmp*dtmp;
    }
    return result;
  }
  
  // vector form
  double norm_squared_distance(const Vector &pa, const Vector &pb)
  {
    double result = 0.0, dtmp;
    int rp;
    for(rp=0; rp<DIM; ++rp)
    {
      dtmp = pa[rp] - pb[rp];
      result += dtmp*dtmp;
    }
    return result;
  }
  
  // apply pbc condition
  void pbc(point &pa)
  {
    mirrorvec_type mirror[DIM];
    latticepoint<point>(pa, mirror);
  }
  
  // return value: 0 - in cell
  // 1 - not in cell, diff by integer
  // 2 - not in cell, diff by integer + 1/2
  int pbc_is_incell(point &pa)
  {
    mirrorvec_type mirror[DIM];
    int retval = 0;
#ifdef DNPLUSFLAG
    retval = latticepoint<point>(pa, mirror);
#else
    latticepoint<point>(pa, mirror);
#endif
    if(0 == retval)
    {
      for(int rp=0; rp<DIM; ++rp)
      {
        if(0 != mirror[rp])
        {
          retval = 1;
          break;
        }
      }
    }
    return retval;
  }
  
  point &point_add(point &pa, const point &pb)
  {
    int i;
    for(i=0; i<DIM; ++i)
    {
      pa[i] += pb[i];
    }
    return pa;
  }
  
  Vector dif_points(const point &point1, const point &point2)
  {
    Vector a(DIM);
    int i;
    for(i=0; i<DIM; ++i)
    {
      a[i] = point1[i] - point2[i];
    }
    return a;
  }

  Vector dif_vec_PBC(Vector &vec)
  {
    mirrorvec_type mirror[DIM];
    latticepoint<Vector>(vec, mirror);
    return vec;
  }
  
  Vector dif_vec_PBC(const Vector &vec_)
  {
    mirrorvec_type mirror[DIM];
    Vector vec = vec_;
    latticepoint<Vector>(vec, mirror);
    return vec;
  }
  
  Vector dif_points_PBC(const point &point1, const point &point2)
  {
    // point1 - point2
    Vector vec = dif_points(point1, point2);
    return dif_vec_PBC(vec);
  }
  
  double norm_squared_distance_PBC(const point &pa, const point &pb)
  {
    double dis[DIM];
    int rp;
    for(rp=0; rp<DIM; ++rp) dis[rp] = pa[rp] - pb[rp];
    return norm_squared_distance_PBC<double[DIM]>(dis);
  }
  
  /** get all nearest and next nearest neighbor copies within the shell of size sqrcut
   */
  std::list<std::pair<point, double> > get_point_images(point &pa, double sqrcut)
  {
    mirrorvec_type mirror[DIM];
    std::vector<int> idx(DIM);
    std::list<std::pair<point, double> > results;
    // first check Dn images
#ifdef DNPLUSFLAG
    latticepoint_inner(pa, mirror);
#else
    latticepoint(pa, mirror);
#endif
    double sqdis = norm_squared_distance(pa);
    if(sqdis <= sqrcut)
    {
      results.push_back({pa, sqdis});
      if(sqrcut>=0.5)
      {
        iota(idx.begin(), idx.end(), 0);
        // argsort coordinate with descend abs value
        std::stable_sort(idx.begin(), idx.end(), [&pa](int i1, int i2) {return fabs(pa[i1]) > fabs(pa[i2]);});
        results.splice(results.end(), get_point_images_dn(pa, idx, sqdis, sqrcut));
      }
    }
#ifdef DNPLUSFLAG
    // shifts 1/2
#ifndef DNLAMFLAG
    // even dimension
    for(int rp=0; rp<DIM; ++rp)
#else
    for(int rp=0; rp<DIM-1; ++rp)
#endif
    {
      pa[rp] += 0.5;
    }
    latticepoint_inner(pa, mirror);
    sqdis = norm_squared_distance(pa);
    if(sqdis <= sqrcut)
    {
      results.push_back({pa, sqdis});
      if(sqrcut>=0.5)
      {
        iota(idx.begin(), idx.end(), 0);
        // argsort coordinate with descend abs value
        std::stable_sort(idx.begin(), idx.end(), [&pa](int i1, int i2) {return fabs(pa[i1]) > fabs(pa[i2]);});
        results.splice(results.end(), get_point_images_dn(pa, idx, sqdis, sqrcut));
      }
    }
#endif
    return results;
  }
  
  Vector point2vec(const point &p)
  {
    Vector vec(0., DIM);
    int rp;
    for(rp=0; rp<DIM; ++rp) vec[rp] = p[rp];
    return vec;
  }
  void setseed(unsigned int seed)
  {
    rg.seed(seed);
  }
  
  double rand()
  {
    return rrand(rg);
  }
  
  double randn()
  {
    return rnorm(rg);
  }

  // adapt gaussian elimination on size DIM linear equation
  void gauss_eliminate(double *a, double *b, double *x, const int n)
  {
    int j, col, row, max_row, dia;
    double max, tmp;
    
    for (dia = 0; dia < n; ++dia) {
      double *vrow = a + dia*n;
      max_row = dia; max = fabs(vrow[dia]);
      row = dia;
      while (++row < n)
      {
        tmp = fabs(*(a+row*n+dia));
        if (max < tmp)
        {
          max_row = row;
          max = std::move(tmp);
        }
      }
      // incomplete swap row. Only swap those columns i>=dia
      if(dia != max_row)
      {
        double *row1 = a+dia*n, *row2 = a+max_row*n;
        for (int i = dia; i < n; ++i) {
          tmp = row1[i];
          row1[i] = row2[i];
          row2[i] = tmp;
        }
        tmp = b[dia]; b[dia] = b[max_row]; b[max_row] = tmp;
      }
      row = dia;
      while (++row < n) {
        double *vrow2 = a + row*n;
        tmp = vrow2[dia] / vrow[dia];
        col = dia;
        while (++col < n)
          vrow2[col] -= tmp * vrow[col];
        b[row] -= tmp * b[dia];
      }
    }
    row = n;
    while (--row >= 0) {
      double *vrow2 = a + row*n;
      tmp = b[row];
      j = n;
      while (--j > row)
        tmp -= x[j] * vrow2[j];
      x[row] = tmp / vrow2[row];
    }
  }
}

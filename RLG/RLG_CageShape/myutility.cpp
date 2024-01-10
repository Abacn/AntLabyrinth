//
//  myutility.cpp
//  VoidPercolation
//

#include "myutility.hpp"

// Utility class includes some static functions

// calculated the squared distance of two points
namespace MyUtility
{
  
  // distance to the origin
  double norm_squared_distance(const Point &pa)
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
  double norm_squared_distance(const Point &pa, const Point &pb)
  {
    double result = 0.0, dtmp;
    for(int rp=0; rp<DIM; ++rp)
    {
      dtmp = pa[rp] - pb[rp];
      result += dtmp*dtmp;
    }
    return result;
  }
  
  Vector point2vec(const Point &p)
  {
    Vector vec({0.});
    for(int rp=0; rp<DIM; ++rp) vec[rp] = p[rp];
    return vec;
  }
  
  // Calculate the squared circumcircle radius
  double facet_weight(const Eigen::Matrix<double, DIM+1, DIM+1> matA, int covertex)
  {
    using namespace Eigen;
    Matrix<double, DIM+1, DIM+1> matB;  // augmented matrix
    int rp, rq, t_rp, t_rq;
    double weight;
    for(rp=0, t_rp=1; rp<=DIM; ++rp)
    {
      matB(rp, rp) = 0.0;
      if(rp == covertex) continue;
      matB(t_rp, 0) = matB(0, t_rp) = 1.0;
      for(rq=rp+1, t_rq=t_rp+1; rq<=DIM; ++rq)
      {
        if(rq == covertex) continue;
        matB(t_rp, t_rq) = matB(t_rq, t_rp) = matA(rp, rq);
        ++t_rq;
      }
      ++t_rp;
    }
    auto matinv = matB.inverse();
    weight = -0.5*matinv(0, 0);
    return weight;
  }
  
}


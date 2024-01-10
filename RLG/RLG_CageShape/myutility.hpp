//
//  myutility.hpp
//  VoidPercolation
//
//  Created by Yi Hu on 8/7/18.
//

#ifndef myutility_h
#define myutility_h

#include <ctime>
#include "graph_construct_cage.hpp"

namespace MyUtility
{
  double norm_squared_distance(const Point &pa);
  double norm_squared_distance(const Point &pa, const Point &pb);
  double facet_weight(const Eigen::Matrix<double, DIM+1, DIM+1> matA, int covertex);
  Vector point2vec(const Point &p);
  static class{
  private:
    clock_t t0;
  public:
    void set(){t0 = clock();}
    float get(){return (float)(clock()-t0)/CLOCKS_PER_SEC;}
  } timer;
}

#endif /* myutility_h */

//
//  neilist_hcomb.hpp
// Neighbor list, specified and optimized for honeyComb periodic box
//  RLGDynamics
//
//  Created by Yi Hu on 11/4/19.
//

#ifndef neilist_hcomb_h
#define neilist_hcomb_h

#include <algorithm>
#include <vector>

#include "myutility.hpp"
#include "hcombpbc.hpp"

/** neighbor list of one particle in Dn and related simulation box. */
struct NeighborList_HComb{
  // members
  int status; // 0: uninitialized; 1-initialized
  std::vector<std::pair<Point_handle, Vector> > neighbors;
  // default constructor
  NeighborList_HComb(): status(0){}
  NeighborList_HComb(Point_handle idx, double sqrcut, const points &all_points)
  {
    construct(idx, sqrcut, all_points);
  }
  NeighborList_HComb(const point &tracer, double sqrcut, const points &all_points)
  {
    construct(tracer, sqrcut, all_points);
  }
  
  Point_handle construct(const point &tracer, double sqrcut, const points &all_points)
  {
    Point_handle i;
    int32_t mirror[DIM];
    double norm;
    Vector pbcshift, diffvec;
    if(!neighbors.empty()) neighbors.clear();
    for(i=0; i<all_points.size(); ++i)
    {
      // calc diffvec
      diffvec.assub(tracer, all_points[i]);
#ifdef DNPLUSFLAG
      hfflag = HCombPBC::latticepoint(diffvec.x, mirror);
#else
      HCombPBC::latticepoint(diffvec.x, mirror);
#endif
      norm = diffvec.norm_squared();
      if(norm < sqrcut)
      {
        // tracer - obstacle = (-vec) + pbcshift
        // vec = obstacle - tracer +  pbcshift
#ifdef DNPLUSFLAG
        HCombPBC::getpbcvec(pbcshift.x, mirror, hfflag);
#else
        HCombPBC::getpbcvec(pbcshift.x, mirror);
#endif
        neighbors.push_back(std::make_pair(i, pbcshift));
      }
    }
    status = 1;
    return (Point_handle)neighbors.size();
  }
  
  void deconstruct()
  {
    if(!neighbors.empty())
    {
      neighbors.clear();
    }
    status = 0;
  }
};

#endif /* neilist_hcomb_h */

//
//  neilist.hpp
//  RLGDynamics neighbor list
//
//  Created by Yi Hu on 6/13/19.
//

#ifndef neilist_hpp
#define neilist_hpp

#include <vector>
#include "commondef.hpp"

typedef double(*NormFunc)(const Vector &, const Vector &);

/* Neighbor list of one particle. */
template<NormFunc norm_func>
struct NeighborList{
  // members
  int status; // 0: uninitialized; 1-initialized
  std::vector<Point_handle> neighbors;
  // default constructor
  NeighborList(): status(0){}
  NeighborList(Point_handle idx, double sqrcut, const points &all_points)
  {
    construct(idx, sqrcut, all_points);
  }
  NeighborList(const point &tracer, double sqrcut, const points &all_points)
  {
    construct(tracer, sqrcut, all_points);
  }
  
  Point_handle construct(Point_handle idx, double sqrcut, const points &all_points)
  {
    Point_handle i;
    double norm;
    Vector vec;
    if(!neighbors.empty()) neighbors.clear();
    for(i=0; i<all_points.size(); ++i)
    {
      if(i==idx) continue;
      norm = norm_func(all_points[idx], all_points[i]);
      if(norm < sqrcut)
      {
        neighbors.push_back(i);
      }
    }
    status = 1;
    return (Point_handle)neighbors.size();
  }
  
  Point_handle construct(const point &tracer, double sqrcut, const points &all_points)
  {
    Point_handle i;
    double norm;
    Vector vec;
    if(!neighbors.empty()) neighbors.clear();
    for(i=0; i<all_points.size(); ++i)
    {
      norm = norm_func(tracer, all_points[i]);
      if(norm < sqrcut)
      {
        neighbors.push_back(i);
      }
    }
    status = 1;
    return (Point_handle)neighbors.size();
  }
  
  Point_handle construct(Point_handle idx, double sqrcut, const points &all_points, const std::vector<Point_handle> &possible_neighbors)
  {
    double norm;
    Vector vec;
    if(!neighbors.empty()) neighbors.clear();
    for(auto i: possible_neighbors)
    {
      if(i==idx) continue;
      norm = norm_func(all_points[idx], all_points[i]);
      if(norm < sqrcut)
      {
        neighbors.push_back(i);
      }
    }
    status = 1;
    return (Point_handle)neighbors.size();
  }
  
  Point_handle construct(const point &tracer, double sqrcut, const points &all_points, const std::vector<Point_handle> &possible_neighbors)
  {
    double norm;
    Vector vec;
    if(!neighbors.empty()) neighbors.clear();
    for(auto i: possible_neighbors)
    {
      norm = norm_func(tracer, all_points[i]);
      if(norm < sqrcut)
      {
        neighbors.push_back(i);
      }
    }
    status = 1;
    return (Point_handle)neighbors.size();
  }
  
  void deconstruct()
  {
    if(!neighbors.empty()) neighbors.clear();
    status = 0;
  }
};
#endif /* neilist_hpp */

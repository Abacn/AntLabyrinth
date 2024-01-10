//
//  box.hpp
//  RLGDynamics
//
//  Created by Yi Hu on 6/13/19.
//

#ifndef box_hpp
#define box_hpp

#include <fstream>

#include "commondef.hpp"
#include "read_input.hpp"
#include "neilist.hpp"
#include "msdrecorder.hpp"
#include "basebox.hpp"

#include "gridlist.hpp"
// use cell list only in small dimension
#if DIM<=4
#define GRIDLIST_ON
#endif

/** Tracer dynamics simulation on (hyper)cubic periodic box. */
class Box: public BaseBox{
public:
  // static method
  // periodic boundary condition
  static void pbc(Vector &vec)
  {
    for(int rp=0; rp<DIM; ++rp)
    {
      if(vec[rp]<-0.5) vec[rp] += 1.;
      else if(0.5<=vec[rp]) vec[rp] -= 1.;
    }
  }
  // norm under periodic boundary condition
  static double norm_squared_PBC(const Vector &va, const Vector &vb)
  {
    Vector vec_ = va - vb;
    pbc(vec_);
    return vec_.norm_squared();
  }
private:
  // member
  ReadInputDy input;
  points obstacles;
  Tracer tracer;
  Brownian br;
  uint32_t count;
  double scale, sqscale, quadscale, invscale, invsqscale;
  double rcut, sqrcut, rtol, rconstC;   // true rcut; rtol=(rcut-scale)^2, rconstC=r^2-(rcut-r)^2
  int nsampletime;        // number of time stamps
  double *timearray, *r2array, *r4array, *chidisarray, *rcovarray;
  double *logr2array, *logr4array;
  double maxt, nowt, maxsinglet;
  int *samplecount;
  int nextrecord;         // next record index
  MSDRecorder *msdrecorder;
  NeighborList<norm_squared_PBC> *nlists;
  Vector nlisttracer;     // store the tracer position when type 2 neighbor list is built
  // additional info
  long ncollision, nescape; // number of collision, escape from neighbor shell
  GridList grids;
  std::ofstream ofd;
  // method
  void genrandomPoints();
  void genrandomVelocity();
  Point_handle getnextcollide(Point_handle lastcollide, double &t_elapse);
  void record(std::vector<double> const &vals);
  double noncollide(Point_handle lastcollide);
  void updatetracer(Point_handle thiscollide, double t_elapse);
  bool samplefinished();
  std::ofstream oftraj; // trajectory, debug usage
public:
  Box(ReadInputDy &input);
  ~Box();
  int process();
  void dump();
};

#endif /* box_hpp */


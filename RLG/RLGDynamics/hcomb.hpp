//
//  hcomb.hpp
//  RLGDynamics
//
//  Created by Yi Hu on 8/5/19.
//

#ifndef hcomb_hpp
#define hcomb_hpp

#include <cstdint>
#include <fstream>

#include "commondef.hpp"
#include "read_input.hpp"
#include "neilist_hcomb.hpp"
#include "msdrecorder.hpp"
#include "basebox.hpp"
#include "hcombpbc.hpp"

/** Tracer dynamics simulation on D_n periodic box. */
class HComb: public BaseBox{
public:
  // static method
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
  NeighborList_HComb nlist;
  Vector nlisttracer;     // store the tracer position when type 2 neighbor list is built
  // additional info
  long ncollision, nescape, nlistsize; // number of collision, escape from neighbor shell, neighbor list size
  std::ofstream ofd;
  // method
  void genrandomPoints();
  void genrandomVelocity();
  Point_handle getnextcollide(Point_handle lastcollide, double &t_elapse);
  void record(std::vector<double> const &vals);
  double noncollide(Point_handle lastcollide);
  void updatetracer(Point_handle thiscollide, double t_elapse);
  void updatercut(double newrcut);
  void updatenlist();
  bool samplefinished();
  std::ofstream oftraj; // trajectory, debug usage
public:
  HComb(ReadInputDy &input);
  ~HComb();
  int process();
  void dump();
};


#endif /* hcomb_hpp */

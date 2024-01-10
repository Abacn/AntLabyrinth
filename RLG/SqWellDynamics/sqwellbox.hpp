//
//  sqwellbox.hpp
//  SqWellDynamics
//
//  Created by Yi Hu on 5/1/20.

#ifndef sqwellbox_hpp
#define sqwellbox_hpp

#include <fstream>

#include "commondef.hpp"
#include "basereadinput.hpp"
#include "neilist.hpp"
#include "msdrecorder.hpp"
#include "basebox.hpp"

#include "gridlist.hpp"
// use cell list only in small dimension
#if DIM<=4
#define GRIDLIST_ON
#endif

class ReadInputDySq: public ReadInput {
public:
  // member method
  int read(const char *inputf);
  // member variable
  System_type systype;
  int repeatrun;
  uint32_t N;
  double phi, rmax;
  int sleft, texp, sinterval, tscaleflag;
  unsigned int seed;
  std::string datafile;
  enum Output_mode outmode;
  int rcuttype;
  double rcut;
  double xi, lambda;
};

class SqwellBox: public BaseBox{
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
  ReadInputDySq input;
  points obstacles;
  Tracer tracer;
  uint32_t count;
  double sconstA;
  double scale, sqscale, quadscale, invscale, invsqscale;
  // true rcut; rtol=(rcut-scale)^2, rconsA = ((1+lambda)^d-1)*scale^d, rconstC=r^2-(rcut-r)^2
  double rcut, sqrcut, rtol, scalepdim, rconstA, rconstC;
  int nsampletime;        // number of time stamps
  double *timearray, *r2array, *r4array, *chidisarray, *rcovarray;
  double maxt, nowt, maxsinglet;
  int *samplecount;
  int nextrecord;         // next record index
  MSDRecorder *msdrecorder;
  NeighborList<norm_squared_PBC> *nlists;
  Vector nlisttracer;     // store the tracer position when type 2 neighbor list is built
  // additional info
  long ncollision, nescape, nlistsize; // number of collision, escape from neighbor shell, neighbor list size
  GridList grids;
  std::ofstream ofd;
  // method
  void genrandomPoints();
  void genrandomVelocity();
  Point_handle getnextcollide(Point_handle lastcollide, double &t_elapse);
  void record(double t_elapse);
  double noncollide(Point_handle lastcollide);
  void updatetracer(Point_handle thiscollide, double t_elapse);
  void updatercut(double newrcut);
  void updatenlist();
  bool samplefinished();
  std::ofstream oftraj; // trajectory, debug usage
public:
  SqwellBox(ReadInputDySq &input);
  ~SqwellBox();
  int process();
  void dump();
};


#endif /* sqwellbox_hpp */

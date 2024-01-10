//
//  sphere.hpp
//  RLGDynamics
//
//

#ifndef sphere_hpp
#define sphere_hpp

#include <fstream>
#include <random>

#include "commondef.hpp"
#include "read_input.hpp"
#include "myutility.hpp"
#include "neilist.hpp"
#include "msdrecorder.hpp"
#include "basebox.hpp"

/** Tracer dynamics simulation in a sphere shell. */
class Sphere: public BaseBox{
private:
  // member
  ReadInputDy input;
  points obstacles;
  Tracer tracer;
  uint32_t count;
  double expectN;
  int N; // now obstacles
  int statuscode;
  // poisson distribution
  std::poisson_distribution<> pois;
  double rcut, sqrcut, rtol, rconstA, rconstC;   // true rcut; rtol=(rcut-scale)^2, rconstC=r^2-(rcut-r)^2
  double shelltol;
  int nsampletime;        // number of time stamps
  double *timearray, *r2array, *r4array, *chidisarray, *rcovarray;
  double *logr2array, *logr4array;
  double maxt, nowt;
  int *samplecount;
  int nextrecord;         // next record index
  bool nlistflag;         // use neighbor list or not
  MSDRecorder *msdrecorder;
  NeighborList<MyUtility::norm_squared> *nlists;
  Vector nlisttracer;     // store the tracer position when type 2 neighbor list is built
  // additional info
  long ncollision, nescape; // number of collision, escape from neighbor shell
  std::ofstream ofd;
  // methods
  void genrandomPoints();
  void genrandomVelocity();
  Point_handle getnextcollide(Point_handle lastcollide, double &t_elapse);
  void record(std::vector<double> const &vals);
  double noncollide(Point_handle lastcollide);
  void updatetracer(Point_handle thiscollide, double t_elapse);
  void updatercut(double newrcut);
  bool samplefinished();
  std::ofstream oftraj; // trajectory, debug usage
public:
  Sphere(ReadInputDy &input);
  ~Sphere();
  int process();
  void dump();
};

#endif /* sphere_hpp */

/*
 * Modified from
 * Packing of hard spheres via molecular dynamics
 * Developed by Monica Skoge, 2006, Princeton University
 * Contact: Aleksandar Donev (adonev@math.princeton.edu) with questions
 * This code may be used, modified and distributed freely.
 * Please cite:
 *
 * "Packing Hyperspheres in High-Dimensional Euclidean Spaces"
 * 	M. Skoge, A. Donev, F. H. Stillinger and S. Torquato, 2006
 *
 * if you use these codes.
*/

//-----------------------------------------------------------------------------
// Box maker
//---------------------------------------------------------------------------

#ifndef  BOX_H
#define  BOX_H


#include <vector>
#include <math.h>
#include <string>
#include <chrono>

#include "dim.h"

#include "read_input.h"
#include "vector.h"
#include "event.h"
#include "sphere.h"
#include "heap.h"
#include "nlist.h"
#include "utility.h"
#include "displacements.h"

#define M 1.0

class Box;
class Displacements;

//---------------------------------------------------------------------------
// Class neighbor
//---------------------------------------------------------------------------
class Neighbor
{
public:
  int i;

  Neighbor(int i_i);

public:
  virtual void Operation(int j, const vector<> &shift) = 0;
};

//---------------------------------------------------------------------------
// Predicts collisions, inherits neighbor operation
//---------------------------------------------------------------------------
class collision : public Neighbor
{
public:

  Box *b;
  double ctime;
  int cpartner;
  vector<> cpartnerpboffset;

public:
  collision(int i_i, Box *b);

  virtual void Operation(int j, const vector<> &shift);
  virtual ~collision() {};

};


class Box {

public:

  // constructor and destructor
  Box(int N_i, double r_i, double rmin_i, double rmax_i, double growthrate_i, double maxpf_i, double shiftscale, int seed);
  Box(read_input const& input);
  ~Box();

  // Creating configurations
  void CreateSpheres(double temp);
  void RecreateSpheres(const std::string fileprefix, double temp);

  // Callibrate Vecolities
  double Velocity(double temp);
  void VelocityGiver(double temp);

  // find collisions from all neighbor
  // neighbor list version
  void ForAllNeighbors(int, Neighbor&);
  void PredictCollision(int i, int j, const vector<> &shift, collision* ccollision);
// Processing an event
  void Process(int n, int option=0, double nextt=0.);
  void Statistics(double ctime);
  void Synchronize(bool rescale);
  void StartMeasure(double nextsampletime, double sampletimedelt);
  void Reset();

// Debugging
  void OutputEvents();

// Statistics
  double Energy();
  void CallibVelocity();
  std::vector<double> Thermodynamics();
  double PackingFraction();
  void PrintStatistics(int mode=0);
  void RunTime();
  // append configuration to wconfigfile
  void WriteConfiguration(const std::string wconfigfile);
  // write configuration (include velocity to wconfigfile)
  void WriteLastConfiguration(const std::string wconfigfile);
  void WriteVelocities(const std::string wconfigfile);
  void WriteShifts(const std::string shiftfile) const;
  void CopyPositions(const Box &b, double temp);

  // MKK shifting
  // reset is min img vector: result = minimg(i-j)
  void getshift(int i, int j, vector<>& result);
  // i-j = result + shift
  void getshift(int i, int j, vector<>& result, vector<>& shift);
  // dx shifted
  vector<> getshift(const int i, const int j);
  // input index, output particle coordinates
  void vectoridx(vector<>& x, const int i);
  static double rootsolverA(double A, double B, double C);

  // variables used outside the class
  const int N;                   // number of spheres
  const double boxvolume;
  const double maxpf;
  double rscale, rmeanfin;             // intermediate radius, final mean radius
  //array
  EventHeap h;                         // event heap

  double gtime;                  // this is global clock
  double rtime;                  // reset time, total time = rtime + gtime

  // arrays
  Sphere *s;                      // array of spheres

  // neighbor list
  nlist* neighs;
  vector<> **allshifts;             // shift table

  // snapshot position
  vector<> *x0;

private:
  void CreateSphere(int Ncurrent, double givenR = 0.);      // Create single sphere
  void CreateShift();   // create MKK shift
  void ReadPositions(const std::string filename, const std::string shiftname = "");
  int ReadVelocities(const std::string filename);
  void SetInitialEvents();

  // Predicting next event
  Event FindNextEvent(int i);
  void CollisionChecker(Event c);
  Event FindNextCollision(int i);
  Event FindNextTransfer(int i);
  int ProcessEvent();
  void Collision(Event e);

  // status variables
  double rx;                     // radius range is (1-rx)r to (1+rx)r
  double rinitial, rfinal;    // max min average radius
  double growthrate;             // growthrate of the spheres
  double grsq;                   // grsq = (rinitial*growthrate)^2;
  double shiftscale;
  double plasttime;              // last time stamp in computing pressure
  double xmomentum;              // exchanged momentum
  double pf;                     // packing fraction
  double nextsampletime;
  double sampletimedelt;
  // statistics
  double energy;                 // kinetic energy
  double energychange;
  double collisionrate;          // average rate of collision between spheres
  double* sums_pmsd;             // particle MSD sums
  double sum_pmqd;               // particle MQD sum
  uint64_t* dist_sd;             // distribution of displacement
  uint64_t ncollisions;               // number of collisions
  uint64_t ntransfers;                // number of transfers
  uint64_t nchecks;                   // number of checks
  uint64_t toteventscheck, toteventscoll, toteventstransfer; // total number of events
  uint64_t ncycles;                   // counts # cycles for output
  uint64_t nlastbuildnlist;           // count number of nlist build
  clock_type tstart, tlastresize;     // run time of program
  duration_type tnlist;
  double visaccu[DIM][DIM];           // viscosity accumulant

  // internal variables
  const int calcpmsd;                 // flag calculate particle MSD (1) or not (0), set to (2) also output distribution
  static const int sd_left = -19, sd_delta = 4, sd_right = 7; // mark the 2^(-19) to 2^(7), 2^(1/4) space, total 106
  static const int sd_nbin = (sd_right - sd_left) * sd_delta + 2;
  uint64_t create_counter;            // counter number of attempt in create sphere
  uint64_t pmsd_counter;              // counter number of particle MSD sample
  static const uint64_t COUNT_LIMIT = 10000000ULL;  // average limit number of create sphere attempt per sphere
  Displacements* disp_stat;
};

#endif

#ifndef DISPLACEMENTS_H
#define DISPLACEMENTS_H

#include <cstdint>

#include "box.h"

class Box;

// Displacement distribution
class Displacements
{
public:
  Displacements(Box* box, double t_start, int t_n, const char* fprefix);
  ~Displacements();
  void initial_snapshot();
  int record(double ctime);
  // dump remaining
  void dump();
private:
  // method
  // dump one
  void dump_t(int now_tn);
  // members
  Box* box;
  double t_start; // smallest time interval
  const int t_n; // number of time applied. t_start, 2 * t_start, 4 * t_start, ..., 2^(t_n - 1) * t_start;
  int* counts;
  const double bin_left; // left most displacement
  const int n_bin;       // number of bins. The bins are bin_left, bin_left * 10^(1/10), bin_left * 10^(2/10)
  const int max_count;   // max count of snapshot per distribution, has to be 2^N
  const char* fprefix;
  int64_t** bins; // distribution. Size is t_n * (2 + n_bin)

  int N; // number of particles
  vector<> **snapshots; // start snapshots
  double *lasttimes; // start times

  // members initialized during initial_snapshot
  double t_scale; // \hat t = internal_t / t_scale
  double Dt_scale; // \hat \Delta = internal_x * internal_x / Dt_scale

  // on-the-fly vars
  double next_dump_t;
  int next_dump_tn;
};

#endif /* DISPLACEMENTS_H */

#ifndef PRESSURETENSOR_H
#define PRESSURETENSOR_H

#include <cstdint>
#include <array>

#include "box.h"

class Box;

typedef std::array<std::array<double, DIM>, DIM> PresureTensorType;

// Pressure tendor related statistics
class PressureTensor
{
public:
  PressureTensor(Box* box, double t_interval, double t_max, const char* fname);
  ~PressureTensor();
  int record();
  int record_remain();
  void reset();
  // dump remaining
  void dump() const;

  void accumulate(int i, int j, double val) // inline
  {
    visaccu[i][j] += val;
  }

  double next_recordtime() const
  {
    return t_nextrecord;
  }
private:
  double viscosity_snapshot(const PresureTensorType& a, const PresureTensorType& b) const;
  Box* box;
  const double t_dumpinterval;
  double t_nextrecord, t_nextexprecord;
  const int t_n, n_exp;
  int next_dump_tn, i_exp, i_bin;
  int* counts;
  double* vals;
  const double t_bin_interval;
  const int max_count;
  const char* fname;
  PresureTensorType visaccu;           // viscosity accumulant
  PresureTensorType* tsnapshot;        // tensor snapshot
  long long debugcount;
};

#endif /* PRESSURETENSOR_H */

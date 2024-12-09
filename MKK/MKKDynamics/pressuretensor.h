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
  PressureTensor(Box* box, double t_interval, int t_max, const char* fname);
  ~PressureTensor();
  int record();
  void reset();
  // dump remaining
  void dump() const;

  void accumulate(int i, int j, double val) // inline
  {
    visaccu[i][j] += val;
  }

private:
  Box* box;
  const double t_dumpinterval;
  double t_nextrecord;
  const int t_n;
  int next_dump_tn;
  int* counts;
  double* vals;
  const char* fname;
  const int max_count;
  PresureTensorType visaccu;           // viscosity accumulant
  PresureTensorType* tsnapshot;        // tensor snapshot
  long long debugcount;
};

#endif /* PRESSURETENSOR_H */

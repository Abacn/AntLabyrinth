#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <string>

#include "pressuretensor.h"
#include "box.h"

PressureTensor::PressureTensor(Box* box, double t_interval, int t_n, const char* fname):
box(box),
t_dumpinterval(t_interval),
t_n(t_n),
next_dump_tn(0),
max_count(1024),
fname(fname)
{
  counts = new int[t_n];
  vals = new double[t_n];
  tsnapshot = new PresureTensorType[t_n];

  for (int i = 0; i < DIM; ++i)
    for (int j = 0; j < DIM; ++j)
      visaccu[i][j] = 0.0;

  for (int i = 0; i < t_n; ++i)
  {
    counts[i] = 0;
    vals[i] = 0.0;
    for (int j = 0; j < DIM; ++j)
      for (int k = 0; k < DIM; ++k)
        tsnapshot[i][j][k] = 0.0;
    this->t_nextrecord = t_dumpinterval;
  }
}

PressureTensor::~PressureTensor()
{
  delete[] counts;
  delete[] vals;
  delete[] tsnapshot;
}


int PressureTensor::record()
{
  if (t_nextrecord <= 0.0) return -1;
  double ttime = box->rtime;
  bool will_dump = false;
  bool will_print = false; // debug;
  while (t_nextrecord <= ttime && next_dump_tn < t_n)
  {
    for (int i = next_dump_tn; i < t_n; ++i)
    {
      double t_nextrecordn = t_dumpinterval * (1LL << i) * (counts[i] + 1);
      if (t_nextrecordn > ttime) break;

      // resolve viscosity
      double accu = 0.;
      double diag = 0.;
      for (int j = 0; j < DIM; ++j)
      {
        diag += visaccu[j][j] - tsnapshot[i][j][j];
      }
      diag /= DIM;
      for (int j = 0; j < DIM; ++j)
      {
        for (int k = 0; k < DIM; ++k)
        {
          double entry = visaccu[j][k] - tsnapshot[i][j][k];
          if (j == k)
            accu += (entry - diag) * (entry - diag);
          else
            accu += entry * entry;
        }
      }
      vals[i] += accu;

      // write back to snapshot
      for (int j = 0; j < DIM; ++j)
        for (int k = 0; k < DIM; ++k)
          tsnapshot[i][j][k] = visaccu[j][k];

      will_print = true;
      // increment counts
      if (++counts[i] == max_count)
      {
        ++next_dump_tn;
        if (next_dump_tn & 1) will_dump = true; // output results on-the-fly on every other time interval filled
      }
    }
    if (next_dump_tn < t_n)
      t_nextrecord = t_dumpinterval * (1LL << next_dump_tn) * (counts[next_dump_tn] + 1);
  }
  if (will_dump)
  {
    dump();
  }
  return 0;
}


void PressureTensor::dump() const
{
  std::ofstream ofs(fname);
  double sigma = box->rmeanfin * 2.0;
  for (int i = 0; i < t_n; ++i)
  {
    if (counts[i] == 0) break;
    double t = t_dumpinterval * (1LL << i);
    // box volume set unit sphere diameter
    double viscosity = vals[i] / (2.0 * (DIM - 1) * (DIM + 2) * box->boxvolume * t * counts[i]) * pow(sigma, DIM - 1);
    ofs << t / sigma << "\t" << viscosity << "\n";
  }
}

void PressureTensor::reset()
{
  for (int i = 0; i < t_n; ++i)
  {
    counts[i] = 0;
    vals[i] = 0;
  }
  for (int i = 0; i < DIM; ++i)
    for (int j = 0; j < DIM; ++j)
      visaccu[i][j] = 0.0;
}

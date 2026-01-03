#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <string>

#include "pressuretensor.h"
#include "box.h"

PressureTensor::PressureTensor(Box* box, double t_min, double t_max, const char* fname):
box(box),
t_dumpinterval(t_min),
t_n(ceil(log2(t_max / t_min) + 1.0)),
n_exp(t_n - 1 - 10),
next_dump_tn(0),
t_bin_interval(t_min* pow(2.0, n_exp)),
max_count(1024), // needs to be 2^x
fname(fname)
{
  // t_min, t_max;
  double sigma = box->rmeanfin * 2.0;
  // std::cout << " DBG " << t_min / sigma << " " << t_max / sigma << std::endl;
  // t_n, n_exp, t_dumpinterval, t_bin_interval
  // std::cout << " DBG " << t_n << " " << n_exp << " " << t_dumpinterval / sigma << " " << t_bin_interval / sigma << std::endl;
  counts = new int[t_n];
  vals = new double[t_n];
  tsnapshot = new PresureTensorType[max_count + 1];

  for (int i = 0; i < DIM; ++i)
    for (int j = 0; j < DIM; ++j)
      visaccu[i][j] = 0.0;

  for (int i = 0; i < t_n; ++i)
  {
    counts[i] = 0;
    vals[i] = 0.0;
  }

  for (int i = 0; i < max_count; ++i)
  {
    for (int j = 0; j < DIM; ++j)
      for (int k = 0; k < DIM; ++k)
        tsnapshot[i][j][k] = 0.0;
  }
  t_nextexprecord = t_dumpinterval;
  t_nextrecord = t_dumpinterval;
  i_exp = i_bin = 0;
}

PressureTensor::~PressureTensor()
{
  delete[] counts;
  delete[] vals;
  delete[] tsnapshot;
}

double PressureTensor::viscosity_snapshot(const PresureTensorType& a, const PresureTensorType& b) const
{
  // resolve viscosity
  double accu = 0.;
  double diag = 0.;
  for (int j = 0; j < DIM; ++j)
  {
    diag += b[j][j] - a[j][j];
  }
  diag /= DIM;
  for (int j = 0; j < DIM; ++j)
  {
    for (int k = 0; k < DIM; ++k)
    {
      double entry = b[j][k] - a[j][k];
      if (j == k)
        accu += (entry - diag) * (entry - diag);
      else
        accu += entry * entry;
    }
  }
  return accu;
}

bool dump_snapshot()
{
  static std::chrono::time_point<std::chrono::steady_clock> last = std::chrono::steady_clock::now();
  std::chrono::time_point<std::chrono::steady_clock> tp = std::chrono::steady_clock::now();
  if (std::chrono::duration<double>(tp - last).count() > 900.0) // dump at most every 15 minutes
  {
    last = tp;
    return true;
  }
  return false;
}

int PressureTensor::record()
{
  if (t_nextrecord <= 0.0) return -1;
  double ttime = box->rtime;
  bool will_dump = false;
  bool will_print = false; // debug;
  while (t_nextrecord <= ttime && i_bin < max_count)
  {
    if (i_exp < n_exp)
    {
      vals[i_exp] += viscosity_snapshot(tsnapshot[i_bin], visaccu);
      ++counts[i_exp++];
      t_nextrecord += t_nextexprecord;
      t_nextexprecord *= 2;
    }
    else
    {
      i_exp = 0;
      tsnapshot[++i_bin] = visaccu;
      // std::cout << " DBG " << t_nextrecord / (box->rmeanfin * 2.0) << " " << ttime / (box->rmeanfin * 2.0) << std::endl;
      t_nextexprecord = t_dumpinterval;
      t_nextrecord = i_bin * t_bin_interval + t_nextexprecord;
      if (dump_snapshot())
      {
        dump();
      }
    }
  }
  return 0;
}

int PressureTensor::record_remain()
{
  record();
  // record remaining time intervals
  for (int i = 0; i < i_bin; ++i)
  {
    for (int j = 1, j_bin = 0; i + j <= i_bin; j <<= 1, j_bin++)
    {
      vals[n_exp + j_bin] += viscosity_snapshot(tsnapshot[i], tsnapshot[i + j]);
      ++counts[n_exp + j_bin];
    }
  }
  dump();
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
    ofs << t / sigma << "\t" << viscosity << "\n"; // "\t" << counts[i] << "\n";
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

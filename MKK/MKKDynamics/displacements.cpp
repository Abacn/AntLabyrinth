#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <string>

#include "displacements.h"
#include "box.h"

Displacements::Displacements(Box* box, double t_start, int t_n, const char* fprefix):
box(box),
t_start(t_start),
t_n(t_n),
bin_left(-3),
n_bin(61),
max_count(1024),
fprefix(fprefix)
{
  counts = new int[t_n];
  bins = new int64_t*[t_n];
  snapshots = new vector<>*[t_n];
  lasttimes = new double[t_n];
  for (int i = 0; i < t_n; ++i)
  {
    counts[i] = 0;
    bins[i] = new int64_t[n_bin + 2];
    std::memset(bins[i], 0, (n_bin + 2) * sizeof(int64_t));
    snapshots[i] = nullptr;
    lasttimes[i] = 0.;
  }
}

void Displacements::initial_snapshot()
{
  if (counts[0] != 0 || snapshots[0] != nullptr)
  {
    std::cerr << "Initialize snapshot has to be called on clean (never recorded) displacement stat!" << std::endl;
    exit(1);
  }
  snapshots[0] = new vector<>[box->N];
  for (int i = 0; i < box->N; ++i)
  {
    snapshots[0][i] = box->s[i].x;
  }
  t_scale = box->rmeanfin*2.0 / DIM;
  Dt_scale = box->rmeanfin*box->rmeanfin*4.0 / DIM;
  next_dump_t = t_start * t_scale;
  next_dump_tn = 0;
}

Displacements::~Displacements()
{
  for (int i = 0; i < t_n; ++i)
  {
    delete[] bins[i];
    if (snapshots[i] != nullptr)
    {
      if (i + 1 < t_n && snapshots[i] != snapshots[i + 1])
      {
        delete[] snapshots[i];
      }
      snapshots[i] = nullptr;
    }
  }
  delete[] lasttimes;
  delete[] snapshots;
  delete[] bins;
  delete[] counts;
}


int Displacements::record(double ctime)
{
  double ttime = box->rtime + ctime;
  if (next_dump_t == 0.0) return -1;
  int ret_val = 0;
  while ( next_dump_t < ttime && next_dump_tn < t_n)
  {
    vector<>* endx = new vector<>[box->N];
    for (int i = 0; i < box->N; ++i)
    {
      endx[i] = box->s[i].x + box->s[i].v * (next_dump_t - (box->rtime + box->s[i].lutime));
    }
    double target_t_interval = next_dump_t - lasttimes[next_dump_tn];
    double next_t_interval = target_t_interval;
    for (int i = next_dump_tn; i < t_n; ++i)
    {
      if (fabs((next_dump_t - lasttimes[i]) / target_t_interval - 1.0) < 1e-10)
      {
        // std::cout << "record " << i << " at " << next_dump_t << " " << counts[i] << std::endl;
        // record i
        for (int j = 0; j < box->N; ++j)
        {
          double disp = vector<>::norm_squared(endx[j], snapshots[i][j]) / Dt_scale;
          int nind = floor((log10(disp) - bin_left) * 10) + 1;
          if (nind < 0) nind = 0;
          else if (nind > n_bin + 1) nind = n_bin + 1;
          ++bins[i][nind];
        }
        ++counts[i];
        // lazy additin to next time-interval bins
        if (1 == counts[i] && i + 1 < t_n)
        {
          snapshots[i + 1] = snapshots[i];
          lasttimes[i + 1] = lasttimes[i];
        }
        // current first time-interval bins at capacity
        if (max_count == counts[i])
        {
          delete[] snapshots[i];
          snapshots[i] = nullptr;
          ++next_dump_tn;
          next_t_interval *= 2;
          dump_t(i);
        }
        else
        {
          if (i + 1 < t_n && snapshots[i] != snapshots[i + 1])
          {
            delete[] snapshots[i];
          }
          snapshots[i] = endx;
        }
        lasttimes[i] = next_dump_t;
      }
      else
      {
        break;
      }
      target_t_interval *= 2;
    }
    next_dump_t += next_t_interval;
    ret_val++;
  }
  return ret_val;
}


void Displacements::dump()
{
  for (int i = next_dump_tn; i < t_n; ++i)
  {
    if (counts[i] == 0)
    {
      break;
    }
    dump_t(i);
  }
}

void Displacements::dump_t(int now_tn)
{
  std::stringstream sst;
  sst << "t = " << t_start << " * 2 ^ " << now_tn;

  if (counts[now_tn] == 0)
  {
    std::cerr << "count is zero for " << sst.str() << "-th , won't dump displacements." << std::endl;
    return;
  }
   std::string fname = std::string(fprefix) + "_" + std::to_string(now_tn) + ".dat";
   std::ofstream ofs(fname);
   ofs << "#" << sst.str() << "; Dt0 = 10^" << bin_left << "; count = " << counts[now_tn] << " * " << box->N << std::endl;
   for (int i = 0; i < n_bin + 2; ++i) ofs << (double) bins[now_tn][i] / counts[now_tn] << std::endl;
}

//
//  brownian.cpp
//  RLGDynamics
//
//  Created by Yi Hu on 8/27/22.
//

#include "basebox.hpp"

const int Brownian::PER_COL = 10;

void Brownian::setup(double dtmin, double dtmax, int n_mult/*=10*/)
{
  t_gridmin = dtmin; t_deltmax = dtmax; n_repeat = PER_COL * (1ULL << n_mult);
  reset();
}

void Brownian::reset()
{
  t_delt = t_gridmin / PER_COL; // PER_COL collisions per sample
  t_next = t_delt;
  n_thisrepeat = 1;
}

int Brownian::reached(double t)
{
  return t > t_next;
}

double Brownian::inc()
{
  double t = t_next;
  // calculate t_next
  if(n_thisrepeat < n_repeat)
  {
    t_next += t_delt;
    ++n_thisrepeat;
  }
  else if(t_delt < t_deltmax)
  {
    t_delt *= 2;
    if(t_delt > t_deltmax)
    {
      t_delt = t_deltmax;
    }
    else
    {
      n_thisrepeat = n_repeat/2 + 1;
    }
    t_next += t_delt;
  }
  else
  {
    // t_delt == t_deltmax
    t_next += t_delt;
  }
  return t;
}

double Brownian::get_delt() const
{
  return t_delt;
}

#include "box.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "vector.h"

//==============================================================
//==============================================================
//  Class Sphere:
//==============================================================
//==============================================================


//==============================================================
// Constructor
//==============================================================
Sphere::Sphere()
{
}


//==============================================================
// Constructor
//==============================================================
Sphere::Sphere(const Sphere& s)
{
  i = s.i;
  x = s.x;
  v = s.v;
  lutime = s.lutime;
  nextevent = s.nextevent;
  nextcollision = s.nextcollision;
  r = s.r;
  //gr = s.gr;
  m = s.m;
  species=s.species;
}

//==============================================================
// Constructor
//==============================================================
Sphere::Sphere(int i_i, const vector<DIM> &x_i, double lutime_i, double r_i, double gr_i, double m_i, int species_i):
  i(i_i),
  x(x_i),
  lutime(lutime_i),
  r(r_i),
  //gr(gr_i),
  m(m_i),
  species(species_i)
{
}

//==============================================================
// Destructor
//==============================================================
Sphere::~Sphere()
{

}



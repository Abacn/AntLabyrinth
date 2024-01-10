#ifndef  SPHERE_H
#define  SPHERE_H

#include <stdint.h>
#include "vector.h"


class Sphere {

 public:

  // constructor and destructor

  Sphere();
  Sphere(const Sphere& s);
  Sphere(int i_i, const vector<DIM> &x, double lutime_i, double r_i, double gr_i, double m_i, int species_i);
  ~Sphere();

 //variables

  int i;                          // sphere ID

  // impending event
  Event nextevent;                // next event...can be collision or transfer
  Event nextcollision;            // next collision if next event is transfer
  // maybe nextnext event

  // past information
  vector<DIM> x;          // position
  vector<DIM> v;          // velocity
  double lutime;          // last update time
  double r;               // sphere radius
  //double gr;            // sphere growth rate
  double m;               // sphere mass
  int species;            // species number (not used during the MD)
 // make sure efficent in memory

};

#endif

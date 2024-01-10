//
//  basebox.hpp
//  RLGDynamics
//  Base class for system
//
//  Created by Yi Hu on 7/16/19.
//

#ifndef basebox_hpp
#define basebox_hpp

#include "commondef.hpp"
#include "myutility.hpp"

/** Tracer status. */
struct Tracer{
  point coord;   // current coordinates
  Vector pbcvec; // counter of periodic box travelled
  Vector v;      // velocity
  Tracer(): coord(0.), pbcvec(0.), v(0.){}
  void setzero(){coord.setzero(); pbcvec.setzero(); v.setzero();}
  Vector getX() const {return coord+pbcvec;} // get displacement
};

/** Brownian dynamics handler. */
class Brownian {
public:
  void setup(double dtmin, double dtmax, int n_mult=10);
  void reset();          // reset the handler in preparation of next run
  int reached(double t); // if time interval bump reached
  double inc();          // return and increment next velocity reset time (t++)
  double get_delt() const; // get interval
protected:
  const static int PER_COL;
private:
  double t_gridmin, t_deltmax, t_next, t_delt;
  uint64_t n_repeat, n_thisrepeat;
};

// null collision obstacle index
#define NULL_COLLISION ((uint32_t)-1)
// brownian collision obstacle index
#define BROWNIAN_COLLISION ((uint32_t)-2)

/** Base class for simulation box and geometry. */
class BaseBox{
public:
  virtual int process() = 0; // run one simulation
  virtual void dump() = 0;   // dump results
  virtual ~BaseBox(){}
};

#endif /* basebox_hpp */

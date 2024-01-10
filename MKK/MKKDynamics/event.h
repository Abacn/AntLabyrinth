#ifndef  EVENT_H
#define  EVENT_H

#include "vector.h"
#define INF    100000000

class Event {

 public:

  // constructor and destructor
  Event(double time_i, int i_i, int j_i, vector<> shift_i);
  Event(double time_i, int i_i, int j_i);
  Event(const Event& e);
  Event();

  ~Event();

  bool operator<(const Event&) const;
  bool operator>(const Event&) const;
  void erase();

 //variables


  double time;             // time of next collision
  int i;        // collision partner with lower number
  int j;        // collision partner with higher number
  vector<> shift;        // shift

  /* 0<=j<=N                binary collision between i and j
     j=N+DIM+1+x            transfer where x=-(k+1) for left wall
                            and x=k+1 for right wall
     j=INF                 both check after event that did not altered motion of                           i and check after event that altered motion of i, i.e                           rescaling of velocities. I currently don't see need t                           o separate the two

     j=-1                  check after collision
  */

};

#endif

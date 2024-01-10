#include "event.h"
#include "utility.h"

//==============================================================
//==============================================================
//  Class Event
//==============================================================
//==============================================================


//==============================================================
// Constructor
//==============================================================
Event::Event(double time_i, int i_i, int j_i, vector<> shift_i):
  time(time_i),
  i(i_i),
  j(j_i),
  shift(shift_i)
{
}

Event::Event(double time_i, int i_i, int j_i):
  time(time_i),
  i(i_i),
  j(j_i)
{
}

Event::Event(const Event& e)
{
  time = e.time;
  i = e.i;
  j = e.j;
  shift = e.shift;
}

Event::Event()
{
}


//==============================================================
// Destructor
//==============================================================
Event::~Event()
{
}

void Event::erase()
{
  time = DBL_LARGE;
  i = 0;
  j = 0;
}

bool Event::operator<(const Event& e) const
{
  return e.time < time;
}

bool Event::operator>(const Event& e) const
{
  return e.time > time;
}



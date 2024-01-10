#ifndef  NLIST_H
#define  NLIST_H

#include <functional>
#include <utility>

#include "dim.h"
#include "vector.h"

// neighbor square distance function
typedef  double(*NormFunc)(const int, const int, vector<>&);
// input index, output particle coordinates
typedef  void(*VectorIdx)(vector<>&, const int);

// neighbor list class
class nlist
{
public:
  // constructor
  nlist(int N, int arrsize = 10);
  // destructor
  ~nlist();
  // set rv
  void setrv(double newev, double newr);
  // try shrink or bump rv. 1-bump; 2-shrink
  bool tryadjustrv(int option);
  // construct all
  int construct(VectorIdx vec, NormFunc normfunc);

  int checkstatus();
  int checkstatus(const vector<> &veci, int i);
  double getRFactor() const;
  double getsqrhfsh() const;
  std::pair<int, const int*> getnlist(int i);
  const vector<> &getstartpos(int i);
  const vector<> &getshift(int i, int idx);
private:
  // expand the arr size
  void expand(int i);
  const int N; // total number of particles
  int status; // status: 0-not updated; 1-updated
  int arrsize;
  int *arrhead;
  int **arr;
  vector<> *startpos; // starting position
  vector<> **shifts;   // shift applied to neighbor
  // rv - shell radius; half shell width=rv/2-r
  double rv, sqrrv, halfshell, sqrhs, rcontact;
};

#endif

#include "nlist.h"

#include <cassert>

#include "utility.h"

nlist::nlist(int N, int arrsize):
N(N), status(0), arrsize(arrsize)
{
  arr = new int* [N];
  arrhead = new int[N];
  startpos = new vector<>[N];
  shifts = new vector<>*[N];
  for (int rp = 0; rp < N; ++rp)
  {
    arrhead[rp] = arrsize;
    arr[rp] = new int[arrsize+1];
    // the head records number of particles in the neighbor list
    arr[rp][0] = 0;
    shifts[rp] = new vector<>[arrsize];
  }
}

nlist::~nlist()
{
  for (int rp = 0; rp < N; ++rp)
  {
    delete[]  arr[rp];
    delete[] shifts[rp];
  }
  delete[] arr;
  delete[] arrhead;
  delete[] startpos;
  delete[] shifts;
}

// set rv
void nlist::setrv(double newrv, double newr)
{
  rv = newrv;
  sqrrv = rv * rv;
  rcontact = newr;
  halfshell = rv / 2 - newr;
  assert(halfshell > 0.0);
  sqrhs = halfshell * halfshell;
  status = 0;
}

bool nlist::tryadjustrv(int option)
{
  double newrv = (1 == option ? rv * pow(2.0, 1.0 / DIM) : rv / pow(2.0, 1.0 / DIM));
  if (newrv > 0.7107 || newrv <= 2 * rcontact) return false;
  else setrv(newrv, rcontact);
  return true;
}

// get current rfactor
double nlist::getRFactor() const
{
  return rv / (2 * rcontact);
}


// get halfshell*halfshell
double nlist::getsqrhfsh() const
{
  return sqrhs;
}

// get neighbor list of particle i
std::pair<int, const int*> nlist::getnlist(int i)
{
  return std::make_pair(arr[i][0], arr[i] + 1);
}

// get start position
const vector<>& nlist::getstartpos(int i)
{
  return startpos[i];
}

// get shift vec
const vector<>& nlist::getshift(int i, int idx)
{
  return shifts[i][idx];
}

// construct all
int nlist::construct(VectorIdx vec, NormFunc normfunc)
{
  // total number of pairs
  int ttn = 0;
  for (int i = 0; i < N; ++i)
  {
    // copy vector
    vec(startpos[i], i);
    // construct neighbor list for particle i
    arr[i][0] = 0;
  }
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < i; ++j)
    {
      vector<> shiftij;
      double dsqr = normfunc(i, j, shiftij);
      if (dsqr < sqrrv)
      {
        if (arr[i][0] == arrhead[i])
        {
          expand(i);
        }
        shifts[i][arr[i][0]] = shiftij;
        arr[i][++arr[i][0]] = j;
        if (arr[j][0] == arrhead[j])
        {
          expand(j);
        }
        shifts[j][arr[j][0]] = -shiftij;
        arr[j][++arr[j][0]] = i;
        ++ttn;
      }
    }
  }
  status = 1;
  return ttn;
}

// check status
int nlist::checkstatus()
{
  return status;
}

int nlist::checkstatus(const vector<> &veci, int i)
{
  if (0 == status) return 0; // the whole neighbor list is not ready
  else
  {
    vector<> dx = startpos[i] - veci;
    setminimg(dx.x);
    if (dx.norm_squared() > sqrhs)
    {
      status = 0;
    }
    return status;
  }
}

// expand the arr size of particle i
void nlist::expand(int i)
{
  int newsize = 2 * arrhead[i];
  int *newarr = new int[newsize + 1];
  vector<> *newshifts = new vector<>[newsize];
  for (int rp = 0; rp < arrhead[i]; ++rp)
  {
    newarr[rp] = arr[i][rp];
    newshifts[rp] = shifts[i][rp];
  }
  newarr[arrhead[i]] = arr[i][arrhead[i]];
  delete[]  arr[i];
  delete[] shifts[i];
  arr[i] = newarr;
  shifts[i] = newshifts;
  arrhead[i] = newsize;
}

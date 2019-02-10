#include <string.h>
#include <math.h>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "leathsitenode.hpp"
#include "read_input.h"
#include "limited_queue.hpp"
#include "ant.hpp"

/* Ant labyrinth on a cluster which is grew by Leath algorithm
 * map version
 * used for p < p_c
 */

using namespace std;

#if APPLY_RAND_TYPE==1
AntWalker::AntWalker(const read_input_ant &inp) :
input(inp), posmin{1, -1}, randrange(*new uniform_int_distribution<int>(0, 2*DIM-1))
#else
AntWalker::AntWalker(const read_input_ant &inp) :
input(inp), posmin{1, -1}
#endif
{
  bool convergedflag = false;
  double powv, base = (DIM-1.0)/DIM;
  count=0;
  msdarray = new double[input.arraylen];
  tmparray = new double[input.arraylen];  // multiple start point
  tmparrcount = new long[input.arraylen];
  arraylen = input.arraylen;
  memset(msdarray, 0, sizeof(double)*input.arraylen);
  // srandom(input.seed);
#if APPLY_RAND_TYPE==1
  rg.seed(input.seed);
#else
  *((int*)randstate) = input.seed;  // Avoid all zero in randstate array
  initstate(input.seed, randstate, 256);
#endif
  size2_dlist = new double[input.arraylen];
  for(int rp=0; rp<arraylen; ++rp)
  {
    if(!convergedflag)
    {
      powv = pow(base, 1<<rp);
      if(powv == 0.0) convergedflag = true;
        size2_dlist[rp] = 0.5*(1-powv);
        }
    else
    {
      size2_dlist[rp] = 0.5;
    }
  }
}

AntWalker::~AntWalker()
{
  delete [] size2_dlist;
  delete [] msdarray;
  delete [] tmparray;
  delete [] tmparrcount;
#if APPLY_RAND_TYPE==1
  delete &randrange;
#endif
}

int AntWalker::antrun()
{
  long nowtime, nextrectime = 1;
  long stepview=10000000, nextview=0, csize, mapsize = 0, timeoffset, timeoffset2;
  int movidx, movdir, movtmp, recidx=0, recidx2;
  unsigned char isoccupy;
  int multdtfactor = 0;
  long multdt = 1L << multdtfactor; // 2^0 = 1 in the beginning
  int rp;
  const int max_ctfactor = 16;
  const long max_count = 1L << max_ctfactor;  // 2^16 = at most 65536 samples
  const int dqfactor = arraylen-1-max_ctfactor;       // from 2^dqfactor ~ 2^(arraylen-1), improve sampling
  // need a queue of max size 2^max_ctfactor-1
  const long dq = max_count >> 1;
  limited_queue<LeathSiteNode> tmpqueue(dq);  // limited size queue
  LeathSiteNode multants[arraylen]; // initialized as 0 by default
  map<LeathSiteNode, unsigned char> siterec = { {zerovc, 1} };
  map<LeathSiteNode, unsigned char>::iterator itbuff;
  LeathSiteNode antsite, tmpsite;  // 0 at original
  csize = 1; // a site at origin
  memset(tmparray, 0, sizeof(double)*arraylen);
  memset(tmparrcount, 0, sizeof(long)*arraylen);
  tmpqueue.push(antsite);
  for(nowtime = 1; recidx<arraylen; ++nowtime)
  {
    if(nowtime == nextrectime)
    {
      // increment the recidx
      ++recidx;
      nextrectime <<= 1;
    }
    // trial move
#if APPLY_RAND_TYPE==1
    movtmp = randrange(rg);
#else
    movtmp = random()%(2*DIM);
#endif
    movidx = movtmp >> 1;
    movdir = movtmp & 0x1;
    antsite.x[movidx] += posmin[movdir];
    itbuff = siterec.find(antsite);
    if(itbuff == siterec.end())
    {
      // need run new
#if APPLY_RAND_TYPE==1
      if(rrand(rg) < input.pf)
#else
        if((double)random()/(double)RAND_MAX < input.pf)
#endif
        {
          isoccupy = 1;
          ++csize;
        }
        else
        {
          isoccupy = 0;
        }
      siterec[antsite] = isoccupy;
      mapsize = siterec.size();
      if(mapsize == count_limit)  // reach the count limit
      {
        if(arraylen > recidx)
          arraylen = recidx; // force to stop the loop
      }
    }
    else
    {
      isoccupy = itbuff->second;
    }
    if(!isoccupy)
    {
      antsite.x[movidx] -= posmin[movdir]; // return to the last position
    }
    if((nowtime & (multdt-1)) == 0) // record
    {
      timeoffset = nowtime >> multdtfactor;
      recidx2 = multdtfactor;
      // record, small time regime
      while(recidx2 < dqfactor)
      {
        tmparray[recidx2] += (antsite - multants[recidx2]).norm_squared();
        tmparrcount[recidx2]++;
        multants[recidx2] = antsite;
        if (timeoffset & 1) break;
        timeoffset >>= 1;
        recidx2++;
      }
      // record, long time regime
      if(recidx2 == dqfactor)
      {
        for(timeoffset2 = 1; timeoffset2 <= tmpqueue.size(); timeoffset2 <<= 1)
        {
          tmparray[recidx2] += (antsite - tmpqueue[tmpqueue.size() - timeoffset2]).norm_squared();
          tmparrcount[recidx2]++;
          ++recidx2;
        }
        tmpqueue.push(antsite);
      }
      if(tmparrcount[multdtfactor] == max_count)
      {
        if(csize<=2)
        {
          // special cases, analytical solution available
          if(csize==1)
          {
            if(mapsize == 1 + 2*DIM)
            {
              // size 1 cluster
              // optimization have to use goto
              goto end_of_antrun;
            }
          }
          else
          {
            // size 2 cluster
            if(mapsize == 4*DIM)
            {
              for(rp=0; rp<arraylen; ++rp)
              {
                msdarray[rp] += size2_dlist[rp];
              }
              goto end_of_antrun;
            }
          }
        }
        ++multdtfactor;
        multdt = 1L << multdtfactor;
      }
    }
#ifdef DEBUG
    if(nowtime >= nextview)
    {
      nextview += stepview;
      cout << "\t" << nowtime << "\t" << csize << endl;
    }
#endif
  }
  // The last one
  tmparray[arraylen-1] = antsite.norm_squared();
  tmparrcount[arraylen-1] = 1;
  for(rp=0; rp<arraylen; ++rp)
  {
    msdarray[rp] += tmparray[rp] / tmparrcount[rp];
  }
end_of_antrun:
  count++;
  return 0;
}



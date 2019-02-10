#include <string.h>
#include <math.h>
#include <string>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "leathsitenode.hpp"
#include "read_input.h"
#include "limited_queue.hpp"
#include "ant_hbd.hpp"

/* Ant labyrinth on a cluster which is grew by Leath algorithm
 * hybrid version
 * Use transfer matrices for small clusters and simulation for large
 */

using namespace std;

HybridAntWalker::HybridAntWalker(const read_input_ant &inp) : AntWalker(inp)
{
  mat_limit = inp.L;
  count_stat = count_dyn = count_partdyn = 0;
  for(int i=0; i<DIM; ++i)
  {
    unitvc[i][i] = 1;
  }
}

HybridAntWalker::~HybridAntWalker()
{
}

// Master function of run one step
int HybridAntWalker::antrun()
{
  int rp, rstatus;
  bool dyn_flag;
  mat_limit_type csize;
  MapType siterec = { {zerovc, 1} };
  vector<LeathSiteNode> sitelist;    // Store the first mat_limit occupied sites
  memset(tmparray, 0, sizeof(double)*arraylen);
  csize = run_leath(siterec, sitelist);
  // cout << csize << endl;
  dyn_flag = false;
  if(csize <= mat_limit)
  {
    rstatus = run_static(siterec, sitelist);
    if(rstatus < 0) // matrix solve fail, just run dyn
    { 
      rstatus = simu_len = arraylen;
      dyn_flag = true;
      ++count_dyn;
    }
    else if(rstatus > 0)
    {
      dyn_flag = true;
      if(rstatus >= arraylen) // useless transfer matrix run happens, parameters does not optimized
      {
        rstatus = simu_len = arraylen;
        if(!optnotice_flag)
        {
          cout << "[Notice] mat_limit too large compare to arraylen" << endl;
          optnotice_flag = true;
        }
        ++count_dyn;
      }
      else
      {
        simu_len = rstatus + 10;
        if(simu_len > arraylen)
        {
          // Use tmat result for long time but also need run full simulation
          simu_len = arraylen;
          ++count_dyn;
        }
        else
        {
          ++count_partdyn;
        }
      }
    }
    else
    {
      // all points covered by static run
      ++count_stat;
    }
    for(rp=rstatus; rp<arraylen; ++rp)
    {
      msdarray[rp] += tmparray[rp];
    }
  }
  else
  {
    rstatus = simu_len = arraylen;
    dyn_flag = true;
    ++count_dyn;
  }
  if(dyn_flag)
  {
    sitelist.clear();  // By now run_dyn doesn't need sitelist. Save some memory
    run_dyn(siterec, csize);
  }
  for(rp=0; rp<rstatus; ++rp)
  {
    msdarray[rp] += tmparray[rp];
  }
  count++;
  return 0;
}

// Leath algorithm generating size 
mat_limit_type HybridAntWalker::run_leath(MapType& siterec, vector<LeathSiteNode> &sitelist)
{
  LeathSiteNode posvc2;
  mat_limit_type csize = 1;
  int rp;
  double weight;
  int qfront = 0;  // current front of the vector
  // Start with cluster size = 2 , unitvc[0]
  MapType::iterator itbuff;
  for(sitelist.push_back(zerovc);csize <= mat_limit && qfront < sitelist.size(); ++qfront)
  {
    const LeathSiteNode popNode = sitelist[qfront];
    // scan nearest neighbors
    for(rp=0; rp<DIM; ++rp)
    {
      // +1
      posvc2 = popNode + unitvc[rp];
      itbuff = siterec.find(posvc2);
      if(itbuff == siterec.end())
      {
        // Does not have the key
#if APPLY_RAND_TYPE==1
        weight = rrand(rg);
#else
        weight = (double)random()/(double)RAND_MAX;
#endif
        if(weight < input.pf)
        {
          // occupied
          sitelist.push_back(posvc2);
          siterec.insert(itbuff, make_pair(posvc2, ++csize) );
        }
        else
        {
          siterec.insert(itbuff, make_pair(posvc2, (mat_limit_type)0) );
        }
      }
      // -1
      posvc2 = popNode - unitvc[rp];
      itbuff = siterec.find(posvc2);
      if(itbuff == siterec.end())
      {
        // Does not have the key
#if APPLY_RAND_TYPE==1
        weight = rrand(rg);
#else
        weight = (double)random()/(double)RAND_MAX;
#endif
        if(weight < input.pf)
        {
          // occupied
          sitelist.push_back(posvc2);
          siterec.insert(itbuff, make_pair(posvc2, ++csize) );
        }
        else
        {
          siterec.insert(itbuff, make_pair(posvc2, (mat_limit_type)0) );
        }
      }
    }
  }
  return csize;
}

int HybridAntWalker::run_static(const MapType& siterec, vector<LeathSiteNode> &sitelist)
{
  int rp;
  mat_limit_type csize = sitelist.size();
  //size_t mapsize = siterec.size();
  if(csize<=2)
  {
    // special cases, analytical solution available
    if(csize==1)
    {
      // size 1 cluster
      // if(mapsize == 1 + 2*DIM) // always true
      return 0;
    }
    else
    {
      // size 2 cluster
      // if(mapsize == 4*DIM) // always true
      for(rp=0; rp<arraylen; ++rp)
      {
        tmparray[rp] = size2_dlist[rp];
      }
      return 0;
    }
  }
  else
  {
    return run_matrix(siterec, sitelist);
  }
}

int HybridAntWalker::run_dyn(MapType& siterec, long csize)
{
  long nowtime, nextrectime = 1;
  long stepview=10000000, nextview=stepview, mapsize = 0, timeoffset, timeoffset2;
  int movidx, movdir, movtmp, recidx=0, recidx2;
  mat_limit_type isoccupy;
  int multdtfactor = 0;
  long multdt = 1L << multdtfactor; // 2^0 = 1 in the beginning
  int rp;
  const int max_ctfactor = 16;
  const long max_count = 1L << max_ctfactor;  // 2^16 = at most 65536 samples
  int dqfactor = simu_len-1-max_ctfactor;       // from 2^dqfactor ~ 2^(simu_len-1), improve sampling
  if(dqfactor<0) dqfactor = 0;
  // need a queue of max size 2^max_ctfactor-1
  const long dq = max_count >> 1;
  limited_queue<LeathSiteNode> tmpqueue(dq);  // limited size queue
  LeathSiteNode multants[simu_len]; // initialized as 0 by default
  MapType::iterator itbuff;
  LeathSiteNode antsite, tmpsite;  // 0 at original
  csize = 1; // a site at origin
  memset(tmparrcount, 0, sizeof(long)*simu_len);
  tmpqueue.push(antsite);
  for(nowtime = 1; recidx<simu_len; ++nowtime)
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
      siterec.insert(itbuff, make_pair(antsite, (mat_limit_type)isoccupy));
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
  tmparray[simu_len-1] = antsite.norm_squared();
  tmparrcount[simu_len-1] = 1;
  for(rp=0; rp<simu_len; ++rp)
  {
    tmparray[rp] /= tmparrcount[rp];
  }
  return 0;
}

int HybridAntWalker::antdump()
{
  if(count == input.maxcluster)
  {
    double t0 = msdarray[0]/count;
    msdarray[0] = input.pf*count;
    cout << "Simulated D^2 at t=1: " << t0 << ", deviation: " << fabs(t0 - input.pf)/input.pf << endl;
  }
  return AntWalker::antdump();
}

long HybridAntWalker::getCount(int idx)
{
  long ret_val = -1;
  switch(idx)
  {
    case 0: ret_val = count; break;
    case 1: ret_val = count_stat; break;
    case 2: ret_val = count_dyn; break;
    case 3: ret_val = count_partdyn; break;
  }
  return ret_val;
}

string HybridAntWalker::getCountStr()
{
  stringstream sstr;
  string strcout;
  sstr << "(" << getCount(0) << "/" << getCount(1) << "/" << getCount(2) << "/" << getCount(3) << ")";
  sstr >> strcout;
  return strcout;
}

#include <string.h>
#include <math.h>
#include <set>
#include <queue>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "leathsitenode.hpp"
#include "invade.hpp"

/* Ant labyrinth on a cluster which is grew by Leath algorithm
 * (ordered_)map version
 * used for p > p_c
 */

using namespace std;

InvadePclt::InvadePclt(const read_input_ant &inp) :
input(inp), posmin{1, -1}
{
  count=0;
  pcarray = new double[input.arraylen];
  sqpcarray = new double[input.arraylen];
  tmparrcount = new long[input.arraylen];
  for(int i=0; i<DIM; ++i)
  {
    unitvc[i][i] = 1;
  }
  for(int i=0; i<input.arraylen; ++i)
  {
    tmparrcount[i] = (long)(pow(2, i/4.0)*100L);
  }
  maxtime = tmparrcount[input.arraylen-1];
#ifdef DEBUG
  cout << "Maxcluster size: " << maxtime << endl;
#endif
  arraylen = input.arraylen;
  memset(pcarray, 0, sizeof(double)*input.arraylen);
  memset(sqpcarray, 0, sizeof(double)*input.arraylen);
  // srandom(input.seed);
#if APPLY_RAND_TYPE==1
  rg.seed(input.seed);
#else
  initstate(input.seed, randstate, 256);
#endif
}

InvadePclt::~InvadePclt()
{
  delete [] pcarray;
  delete [] sqpcarray;
  delete [] tmparrcount;
}

int InvadePclt::ivdrun()
{
  long nextrectime = 0;
  long csize=0;
  double weight, dtmp;
  int rp;
  LeathSiteNode posvc;
  set<LeathSiteNode> siterec = { {zerovc} };
  pair<set<LeathSiteNode>::iterator,bool> retval;
  priority_queue<PQueueNode> mypque;
  mypque.push(PQueueNode(zerovc, 0.0));
  
  for(csize = 0; csize<=maxtime;)
  {
    const LeathSiteNode &popNode = mypque.top().getSiteNode();
    // pop the front element
    mypque.pop();
    // scan nearest neighbors
    for(rp=0; rp<DIM; ++rp)
    {
      // +1
      posvc = popNode + unitvc[rp];
      retval = siterec.insert(posvc);
      if(retval.second)
      {
        // Does not have the key
#if APPLY_RAND_TYPE==1
        weight = rrand(rg);
#else
        weight = (double)random()/(double)RAND_MAX;
#endif
        mypque.push(PQueueNode(posvc, weight));
      }
      // -1
      posvc = popNode - unitvc[rp];
      retval = siterec.insert(posvc);
      if(retval.second)
      {
        // Does not have the key
#if APPLY_RAND_TYPE==1
        weight = rrand(rg);
#else
        weight = (double)random()/(double)RAND_MAX;
#endif
        mypque.push(PQueueNode(posvc, weight));
      }
    }
    ++csize;
    // cout << csize << "\t" << siterec.size() << "\t" << double(csize)/siterec.size() << endl;
    if(csize == tmparrcount[nextrectime])
    {
      dtmp = (double)csize / siterec.size();
      pcarray[nextrectime] += dtmp;
      sqpcarray[nextrectime] += dtmp*dtmp;
      ++nextrectime;
    }
  }

  count++;
  return 0;
}

int InvadePclt::ivddump()
{
  ofstream output(input.datafile, ios::out);
  long rp;
  output << "count: " << count << "\n\nsize\tpc\tVar\n";
  for(int i=0; i<arraylen; ++i)
  {
    rp = tmparrcount[i];
    output << rp << "\t" << setprecision(10) << pcarray[i]/count << "\t" << (sqpcarray[i]-pcarray[i]*pcarray[i]/count)/count << "\n";
  }
  return 0;
}

// Main function
int main(int argc, const char * argv[])
{
  const char *inputf, *dftinputf = "inputfile_ivd.dat";
  long stepview, nextview, nowclu;
#ifdef DEBUG
  clock_t t1, t2;
  t1 = clock();
#endif
  if (argc < 2)
  {
    inputf = dftinputf;
  }
  else
  {
    inputf = argv[1];
  }
  read_input_ant input;
  int error = input.read(inputf);
  if (error) return error;
  InvadePclt ivdwalker(input);
  stepview = input.maxcluster / 50;
  nextview = stepview;
  cout << "Start calculation..." << endl;
  for(nowclu=0; nowclu<input.maxcluster; ++nowclu)
  {
    ivdwalker.ivdrun();
    if(nowclu >= nextview)
    {
      cout << "(" << nowclu << ") " << endl;
      ivdwalker.ivddump();
      nextview += stepview;
    }
  }
  cout << "(" << nowclu << ") \nEnd of calculation" << endl;
  ivdwalker.ivddump();
#ifdef DEBUG
  t2 = clock();
  cout << (float)(t2-t1)/CLOCKS_PER_SEC << endl;
#endif
}

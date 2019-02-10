//
//  leath_counter.cpp
//  MSDLimit
//
//  Created by Yi Hu on 5/3/18.
//

#include <string.h>
#include <math.h>
#include <queue>
#include <functional>
#include <set>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "leathsitenode.hpp"
#include "leath_counter.hpp"

using namespace std;

LeathCounter::LeathCounter(const read_input_leath &inp): Ndstb(1001), input(inp)
{
  count=countsites=percocount=0;
  stat = statsq = 0.0;
  distribution = new long[Ndstb];
  for(int i=0; i<DIM; ++i)
  {
    unitvc[i][i] = 1;
  }
  memset(distribution, 0, sizeof(long)*Ndstb);
  rg.seed(input.seed);
  recFile.open("msddump.dat", ios::out);
  recFile << DIM << "D\np:\t" << input.pf << "\nSite limit:\t" << count_limit << "\nsize\tperimeter\tmsd\n";
}


LeathCounter::~LeathCounter()
{
  delete [] distribution;
  recFile.close();
}

int LeathCounter::clugen()
{
  LeathSiteNode posvc2;
  Coordinate<DIM, long> sr2, r2s; // <r>^2 and <r^2>

  long csize = 0;
  int rp;
  double weight, result;
  queue<LeathSiteNode> myque;
  // Start with cluster size = 2 , unitvc[0]
  // unordered_set<LeathSiteNode> siterec = {zerovc};
  set<LeathSiteNode> siterec = {zerovc};  // for p close to p_p, set is better than unordered_set
  pair<set<LeathSiteNode>::iterator,bool> retval;
  myque.push(zerovc);

  while(siterec.size() < count_limit && !myque.empty())
  {
    const LeathSiteNode &popNode = myque.front();
    // scan nearest neighbors
    for(rp=0; rp<DIM; ++rp)
    {
      // +1
      posvc2 = popNode + unitvc[rp];
      retval = siterec.insert(posvc2);
      if(retval.second)
      {
        // Does not have the key
        weight = rrand(rg);
        if(weight < input.pf)
        {
          // occupied
          myque.push(posvc2);
        }
      }
      // -1
      posvc2 = popNode - unitvc[rp];
      retval = siterec.insert(posvc2);
      if(retval.second)
      {
        // Does not have the key
        weight = rrand(rg);
        if(weight < input.pf)
        {
          // occupied
          myque.push(posvc2);
        }
      }
    }
    // count
    sr2 += popNode.Long();
    r2s += popNode.Long().squared();
    ++csize;
    // pop the front element
    myque.pop();
  }
  if(siterec.size() >= count_limit)
  {
    percocount++;
    if(1==percocount)
      cout << "Size limit exceed at step " << count << endl;
  }
  result = (double)(r2s.sum()) - sr2.norm_squared()/(double)csize;  // s(<r^2> - <r>^2)
  result = result/csize;
  stat += result;
  statsq += result*result;
  if(csize < Ndstb) distribution[csize]++;
  if(count < 2000000) // Only output first 2x10^6 clusters, otherwise disk space stuck
    recFile << csize << "\t" << siterec.size()-csize << "\t" << setprecision(10) << result*2.0 << "\n";
  // result did not *=2 in early versions, but this does not change the local slope
  countsites += csize;
  count++;
  return 0;
}


int LeathCounter::cludump()
{
  double msdresult, msddevi;
  if(count<=0)
    return -1; // counterError
  ofstream output(input.datafile, ios::out);
  // First line is MDS_\infty
  msdresult = stat*2/count;
  msddevi = sqrt(statsq*4/count - msdresult*msdresult);
  output << msdresult << "\t" << msddevi;
  output << endl;
  output << "Number of cluster/percolated: " << count << "\t" << percocount << "\n";
  // Third line outputs cluster distribution information
  output << "Cluster distributions: " << "\n";
  for(int i=1; i<Ndstb; ++i)
  {
    output << i << "\t" << (double)distribution[i]/(double)count << "\n";
  }
  output << endl;
  return 0;
}


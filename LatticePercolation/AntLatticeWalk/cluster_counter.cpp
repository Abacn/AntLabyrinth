//
//  cluster_counter.cpp
//  MSDLimit
//
//  Created by Yi Hu on 5/3/18.
//

#include <string.h>
#include <math.h>
#include <memory>
#include <queue>
#include "cluster_counter.hpp"

#ifdef PERCOLATION_DETECTION
#include <map>
#endif

using namespace std;

ClusterCounter::ClusterCounter(const read_input_lattice &inp, grid_field<DIM, char> &sys): input(inp), system(sys), Ndstb(1001)
{
  constructor_type = 0;
  N=pow(input.L,DIM);
  system.set_size(input.L); //set grid side
  system.set_seed(input.seed);
  system.initialize();      // initialize cells to 0
  stat=0.0;
  statsq = 0.0;
  count=percocount=0;
  distribution = new long[Ndstb];
  memset(distribution, 0, sizeof(long)*Ndstb);
}

ClusterCounter::ClusterCounter(const read_input_lattice &inp): input(inp), system(*new grid_field<DIM, char>), Ndstb(1001)
{
  constructor_type = 1;
  N=pow(input.L,DIM);
  system.set_size(input.L); //set grid side
  system.set_seed(input.seed);
  system.initialize();      // initialize cells to 0
  stat=0.0;
  statsq = 0.0;
  count=percocount=0;
  distribution = new long[Ndstb];
  memset(distribution, 0, sizeof(long)*Ndstb);
}

ClusterCounter::~ClusterCounter()
{
  if(1==constructor_type)
    delete &system;
  delete distribution;
}

int ClusterCounter::clugen()
{
  int rp;
  long i, tmpi;
  int err = system.initialize(input);
#ifdef PERCOLATION_DETECTION
  map<long, long> siterec;
  long tmpoff;
  int ispercolated = 0;
#endif
  if(err) return err; // percolated (todo)
  
  // convenient constant vectors
  Coordinate<DIM, int> zerovc, posvc, offstvc, posvc2, offstvc2;
  Coordinate<DIM, int> unitvc[DIM];
  const Coordinate<DIM, int> grid_size = system.get_size();
  for(int i=0; i<DIM; ++i)
  {
    unitvc[i][i] = 1;
  }
  double result = 0.0, tmpd;
  
  queue<shared_ptr<SiteNode> > myque;
  
  // scan clusters
  for(i=0; i<N; ++i)
  {
    if(1==system.get(i)) // not empty
    {
      Coordinate<DIM, long> sr2, r2s; // <r>^2 and <r^2>
      int csize = 0;
      Coordinate<DIM, int> coord = system.int2coord(i);
      shared_ptr<SiteNode> st(new SiteNode(coord, zerovc));
      system.set(i, 2); // the start point
      myque.push(st);
#ifdef PERCOLATION_DETECTION
      if(!ispercolated)
      {
        siterec.clear();
        siterec[i] = 0;
      }
#endif
      while(!myque.empty())
      {
        const shared_ptr<SiteNode> &popNode = myque.front();
        posvc = popNode->coord;
        offstvc = popNode->offset;
        myque.pop();
        // find neighbors
        for(rp=0; rp<DIM; ++rp)
        {
          // +1
          posvc2 = posvc + unitvc[rp];
          system.applypbc(posvc2, rp);
          tmpi = system.coord2int(posvc2);
          if(system.get(tmpi) == 1)
          {
            shared_ptr<SiteNode> st2(new SiteNode(posvc2, offstvc+unitvc[rp]));
            system.set(tmpi, 2);
            myque.push(st2);
#ifndef PERCOLATION_DETECTION
          }
#else
            // Add the site to record
            if(!ispercolated)
            {
              siterec[tmpi] = system.coord2int(st2->offset);
            }
          }
          else if(!ispercolated && system.get(tmpi) == 2)
          {
            // Check if the recorded site has the same offset
            if(siterec[tmpi] != system.coord2int(offstvc+unitvc[rp]))
            {
              // Percolated
              ispercolated = 1;
              percocount++;
#ifdef DEBUG
              cout << "Percolated: Replica " << count+1 << endl;
#endif
            }
          }

#endif
          // -1
          posvc2 = posvc - unitvc[rp];
          system.applypbc(posvc2, rp);
          tmpi = system.coord2int(posvc2);
          if(system.get(tmpi) == 1)
          {
            shared_ptr<SiteNode> st2(new SiteNode(posvc2, offstvc-unitvc[rp]));
            system.set(tmpi, 2);
            myque.push(st2);
#ifndef PERCOLATION_DETECTION
          }
#else
            // Add the site to record
            if(!ispercolated)
            {
              siterec[tmpi] = system.coord2int(st2->offset);
            }
          }
          else if(!ispercolated && system.get(tmpi) == 2)
          {
            // Check if the recorded site has the same offset
            if(siterec[tmpi] != system.coord2int(offstvc-unitvc[rp]))
            {
              // Percolated
              ispercolated = 1;
              percocount++;
#ifdef DEBUG
              cout << "Percolated: Replica " << count+1 << endl;
#endif
            }
          }
#endif
        }
        // count
        sr2 += offstvc.Long();
        r2s += offstvc.Long().squared();
        ++csize;
      }
      result += (double)(r2s.sum()) - (double)(sr2.squared().sum())/csize;
      if(csize<Ndstb)
      {
        ++distribution[csize];
      }
      ++distribution[0];
    }
  }
  result *= 2.0/system.getOccupied();
  stat += result;
  statsq += result*result;
  count++;
  return 0;
}

int ClusterCounter::cludump()
{
  double msdresult, msddevi;
  if(count<=0)
    return -1; // counterError
  ofstream output(input.datafile, ios::out);
  // First line is MDS_\infty
  msdresult = stat/count;
  output << msdresult;
  if(count > 1)
  {
    // standard deviation = sqrt((<x^2>-<x>^2)/(count-1))
    msddevi = sqrt(statsq/count - msdresult*msdresult);
    output << "\t" << msddevi;
  }
  output << std::endl;
#ifdef PERCOLATION_DETECTION
  output << "Fraction of percolated replica: " << (double)percocount/count;
#endif
  output << std::endl;
  // Third line outputs cluster distribution information
  double allparticles = count*system.getOccupied();
  output << "Cluster distributions: " << std::endl;
  for(int i=1; i<Ndstb; ++i)
  {
    output << i << "\t" << distribution[i]/allparticles*i << "\n";
  }
  output << std::endl;
  
  return 0;
}


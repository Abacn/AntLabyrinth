//
//  leath_counter.hpp
//  AntLatticeWalk
//
//  Created by Yi Hu on 6/27/18.
//  Copyright © 2018 Yi Hu. All rights reserved.
//

#ifndef leath_counter_h
#define leath_counter_h

#include "read_input.h"
#include "leathsitenode.hpp"
#include <random>

// Cluster counter, leath algorithm
class LeathCounter
{
private:
  long count, countsites, percocount;
  const long count_limit = 40000000000ll / (DIM*4+50); // 10^10 (40G) / (DIM*4+50 bit)
  long N;
  const read_input_leath &input;
  int constructor_type;
  double stat, statsq;
  long *distribution;  //  The first element store the count of all clusters
  const int Ndstb; // the size of distribution matrix.
  LeathSiteNode zerovc, unitvc[DIM]; // help variables
  std::ofstream recFile;
protected:
  std::mt19937_64 rg;
  std::uniform_real_distribution <> rrand;
public:
  LeathCounter(const read_input_leath &inp);
  ~LeathCounter();
  int clugen();
  int cludump();
  long getCount(){return count;}
  long getCountSites(){return countsites;}
  double getMSD(){return stat*2/count;}
};



#endif /* leath_counter_h */

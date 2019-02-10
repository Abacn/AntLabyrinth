//
//  ant.h
//  AntLatticeWalk
//
//  Created by Yi Hu on 5/2/18.
//  Copyright © 2018 Yi Hu. All rights reserved.
//

#ifndef my_ant_h
#define my_ant_h

#include "leathsitenode.hpp"
#include "read_input.h"
#include "ant_global.hpp"

// rand_type:
// 1 : mt19937_64, slower but longer period: 2^19937
// 2 : use random() in RAND_MAX; random() has a period of 16*(2^31-1)
// See https://www.gnu.org/software/libc/manual/html_node/BSD-Random.html#BSD-Random

#define APPLY_RAND_TYPE 1

#if APPLY_RAND_TYPE==1
#include <random>
#elif APPLY_RAND_TYPE==2
#include <stdlib.h>
#else
#error APPLY_RAND_TYPE must be either 1 or 2
#endif

// Cluster counter, leath algorithm
class AntWalker
{
protected:
  long count;
  const long count_limit = 2000000000l / (DIM*4+50); // 10^7 (2G) / (DIM*4+50 bit)
  const read_input_ant &input;
  double *msdarray;  //  The first element store the count of all clusters
  double *tmparray;  // multiple start point
  double *size2_dlist;  // list of D^2 for size of 2
  long *tmparrcount;
  int arraylen;
  LeathSiteNode zerovc; // help variables
#if APPLY_RAND_TYPE==1
  std::mt19937 rg;
  std::uniform_real_distribution <> rrand;
  std::uniform_int_distribution<int> &randrange;
#endif
  const int posmin[2];    // utility array
public:
  AntWalker(const read_input_ant &inp);
  ~AntWalker();
  int antrun();
  int antdump();
  long getCount(){return count;}
};


#endif /* my_ant_h */

//
//  invade.hpp
//  AntLatticeWalk
//
//  Created by Yi Hu on 7/31/18.
//  Copyright © 2018 Yi Hu. All rights reserved.
//

#ifndef invade_hpp
#define invade_hpp


#include "leathsitenode.hpp"
#include "read_input.h"

// rand_type:
// 1 : mt19937_64, slower but longer period: 2^19937
// 2 : use random() in RAND_MAX; random() has a period of 16*(2^31-1)
// See https://www.gnu.org/software/libc/manual/html_node/BSD-Random.html#BSD-Random

#define APPLY_RAND_TYPE 1

#if APPLY_RAND_TYPE==1
#include <random>
#elif APPLY_RAND_TYPE==2
#include <stdlib.h>
char randstate[256];
#else
#error APPLY_RAND_TYPE must be either 1 or 2
#endif

// Cluster counter, leath algorithm
class InvadePclt
{
private:
  long count;
  const long count_limit = 2000000000l / (DIM*4+50); // 10^7 (2G) / (DIM*4+50 bit)
  const read_input_ant &input;
  double *pcarray;  //  pc
  double *sqpcarray;  // perimeter size
  long *tmparrcount;
  long maxtime;
  int arraylen;
  LeathSiteNode zerovc, unitvc[DIM]; // help variables
protected:
#if APPLY_RAND_TYPE==1
  std::mt19937_64 rg;
  std::uniform_real_distribution <> rrand;
#endif
  const int posmin[2];    // utility array
public:
  InvadePclt(const read_input_ant &inp);
  ~InvadePclt();
  int ivdrun();
  int ivddump();
  long getCount(){return count;}
};

// Priority queue node
class PQueueNode
{
public:
  inline PQueueNode(const LeathSiteNode sitenode, double weight): sitenode(sitenode), weight(weight){}
  inline PQueueNode(const PQueueNode& pq):sitenode(pq.sitenode), weight(pq.weight){}
  inline bool operator<(const PQueueNode &b) const{return(weight < b.weight);}
  inline bool operator==(const PQueueNode &b) const{return(weight == b.weight);}
  inline const LeathSiteNode getSiteNode() const{return sitenode;}
  inline operator double() const{return weight;}
private:
  LeathSiteNode sitenode;
  double weight;
};
#endif /* invade_h */

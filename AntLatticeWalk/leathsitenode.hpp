//
//  leathsitenode.hpp
//  AntLatticeWalk
//
//  Created by Yi Hu on 7/6/18.
//

#ifndef leathsitenode_h
#define leathsitenode_h

#include "dim.h"
#include "coordinate.h"

#ifndef TYPE_LEATHSITENODE
#define TYPE_LEATHSITENODE
// Site node data structure for Leath algo
typedef Coordinate<DIM, int> LeathSiteNode;
#endif

namespace std {
  // costom hash function
  template <>
  struct hash<LeathSiteNode>
  {
    std::size_t operator()(const LeathSiteNode& k) const
    {
      long hv = k.x[0];
      const long primelist[] = {31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199};
      for(int rp=1; rp<DIM; ++rp)
      {
//      hv += 47*k.x[rp]; // multiply by a prime number
#if DIM<=36
        hv = hv*primelist[rp] + k.x[rp];
#else
        hv = hv*47 +  k.x[rp];
#endif
      }
      return std::hash<long>()(hv);
//      return hv;
    }
  };
  
}

#endif /* leathsitenode_h */

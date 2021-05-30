//
//  invade.hpp
//  AntLatticeWalk
//
//  Created by Yi Hu on 7/31/18.
//  Copyright © 2018 Yi Hu. All rights reserved.
//

#ifndef invade_hpp
#define invade_hpp

#include <vector>

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

// Priority queue node
template <typename S>
class PQueueNode_T
{
public:
  inline PQueueNode_T(typename std::set<S>::const_iterator sitenode, double weight): sitenode(sitenode), weight(weight){}
  inline PQueueNode_T(const PQueueNode_T& pq):sitenode(pq.sitenode), weight(pq.weight){}
  inline bool operator<(const PQueueNode_T &b) const{return(weight < b.weight);}
  inline bool operator==(const PQueueNode_T &b) const{return(weight == b.weight);}
  inline const S& getSiteNode() const{return *sitenode;}
  inline operator double() const{return weight;}
private:
  typename std::set<S>::const_iterator sitenode;
  double weight;
};

using PQueueNode = PQueueNode_T< LeathSiteNode >;
using BQueueNode = PQueueNode_T< LeathBondNode >;

// Cluster counter, leath algorithm
class InvadePclt
{
protected:
  long count;
  const long count_limit = 2000000000l / (DIM*4+50); // 10^7 (2G) / (DIM*4+50 bit)
  const read_input_ant &input;
  double *pcarray;  //  pc
  double *sqpcarray;  // perimeter size
  long *tmparrcount;
  long maxtime;
  int arraylen;
  LeathSiteNode zerovc, unitvc[DIM]; // help variables
#if APPLY_RAND_TYPE==1
  std::mt19937_64 rg;
  std::uniform_real_distribution <> rrand;
#endif
  const int posmin[2];    // utility array
public:
  InvadePclt(const read_input_ant &inp);
  virtual ~InvadePclt();
  virtual int ivdrun() = 0;
  int ivddump();
  long getCount(){return count;}
};

// hypercubic lattice
class InvadePclt_Zn: public InvadePclt
{
public:
  InvadePclt_Zn(const read_input_ant &inp);
  virtual ~InvadePclt_Zn();
  virtual int ivdrun() override;
};

// hypercubic lattice bond percolation
class InvadePclt_Znb: public InvadePclt
{
public:
  InvadePclt_Znb(const read_input_ant &inp);
  virtual ~InvadePclt_Znb();
  virtual int ivdrun() override;
private:
  std::vector<LeathBondNode> getNeighbors(LeathSiteNode const &node);
};

// Dn (fcc) lattice, site percolation
class InvadePclt_Dn: public InvadePclt
{
public:
  InvadePclt_Dn(const read_input_ant &inp);
  virtual ~InvadePclt_Dn();
  virtual int ivdrun() override;
private:
  std::vector<LeathBondNode> getNeighbors(LeathSiteNode const &node);
};

// Dn (fcc) lattice, bond percolation
class InvadePclt_Dnb: public InvadePclt
{
public:
  InvadePclt_Dnb(const read_input_ant &inp);
  virtual ~InvadePclt_Dnb();
  virtual int ivdrun() override;
  // size of Dn base vectors
  static constexpr int size_dnbase = DIM*(DIM-1);
private:
  std::vector<LeathBondNode> getNeighbors(LeathSiteNode const &node);
  // Dn base vectors
  Coordinate<DIM, IdxType> dnbase[size_dnbase];
};

// dense packing lattice, site percolation
class InvadePclt_DS: public InvadePclt
{
public:
  InvadePclt_DS(const read_input_ant &inp);
  virtual ~InvadePclt_DS();
  virtual int ivdrun() override;
#if DIM>=6 && DIM<8
  static constexpr int ndim = 8;
#else
  static constexpr int ndim = DIM;
#endif
  using SiteNode = Coordinate<ndim, IdxType>;
  using PQueueNode = PQueueNode_T<SiteNode>;
};

// dense packing lattice, bond percolation
class InvadePclt_DSb: public InvadePclt
{
public:
  InvadePclt_DSb(const read_input_ant &inp);
  virtual ~InvadePclt_DSb();
  virtual int ivdrun() override;
#if DIM>=6 && DIM<8
  static constexpr int ndim = 8;
#else
  static constexpr int ndim = DIM;
#endif
  using SiteNode = Coordinate<ndim, IdxType>;
  using BondNode = std::pair< Coordinate<ndim, IdxType>, short>;
  using BQueueNode = PQueueNode_T<BondNode>;
private:
  std::vector<BondNode> getNeighbors(SiteNode const &node);
};

#endif /* invade_h */

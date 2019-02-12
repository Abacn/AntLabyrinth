//
//  cluster_counter.hpp
//  AntLatticeWalk
//
//  Created by Yi Hu on 5/3/18.
//

#ifndef cluster_counter_h
#define cluster_counter_h

#include "grid_field.h"

// The switch of percolation detection. Comment out if off
#define PERCOLATION_DETECTION 1

class ClusterCounter
{
private:
  int count, percocount;
  long N;
  grid_field<DIM, char> &system;
  const read_input_lattice &input;
  int constructor_type;
  double stat, statsq;
  long *distribution;  //  The first element store the count of all clusters
  const int Ndstb; // the size of distribution matrix.
public:
  ClusterCounter(const read_input_lattice &inp);
  ClusterCounter(const read_input_lattice &inp, grid_field<DIM, char> &sys);
  ~ClusterCounter();
  int clugen();
  int cludump();
  int getCount(){return count;}
  double getMSD(){return stat/count;}
};

// Node of sites, used for queue
class SiteNode
{
public:
  Coordinate<DIM, int> coord;
  Coordinate<DIM, int> offset;
  SiteNode(Coordinate<DIM, int> coord, Coordinate<DIM, int> offset)
  {
    this->coord = coord;
    this->offset = offset;
  }
};

#endif /* cluster_counter_h */

//
//  gridlist.hpp
//  QhullPercolation
//
//  Created by Yi Hu on 4/9/19.
//

#ifndef gridlist_hpp
#define gridlist_hpp

#include <vector>
#include "dim.h"
#include "commondef.hpp"

/** Grid list (cell list) of obstacles under cubic PBC. */
class GridList{
public:
  GridList();
  bool is_valid() const{return size>3;}
  int set_points(const points &pos);
  int update_rcut(double rcut);
  void get_neighbor(const point &coord, std::vector<Point_handle> &neighbors) const;
  std::size_t get_size() const{return size;}
  bool is_empty() const{return grids.empty();}
private:
  // members
  double l;      // cell length
  // number of cells per direction, number of total cells, number of neighbor cells 
  std::size_t n, size, nnbr;
  std::vector<std::vector<Point_handle> > grids;
  // functions
  Point_handle grid2idx(const long grid[DIM]) const;
  void idx2grid(Point_handle idx, long grid[DIM]) const;
  Point_handle pos2idx(const point &coord) const;
};


#endif /* gridlist_hpp */

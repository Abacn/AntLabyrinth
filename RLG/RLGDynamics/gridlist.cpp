//
//  gridlist.cpp
//  QhullPercolation
//
//  Created by Yi Hu on 4/9/19.
//

#include <vector>
#include "commondef.hpp"
#include "dim.h"
#include "gridlist.hpp"

GridList::GridList(): n(1), size(1), nnbr(1)
{}

Point_handle GridList::grid2idx(const long grid[DIM]) const
{
  Point_handle result = 0;
  long ltmp;
  for(int rp=0; rp<DIM; ++rp)
  {
    result *= n;
    ltmp = grid[rp];
    if(ltmp<0) ltmp+=n;
    else if(ltmp>=n) ltmp -= n;
    result += ltmp;
  }
  return result;
}

void GridList::idx2grid(Point_handle idx, long grid[DIM]) const
{
  for(int rp=DIM-1; rp>=0; --rp)
  {
    grid[rp] = idx%n;
    idx /= n;
  }
}

Point_handle GridList::pos2idx(const point &coord) const
{
  Point_handle result = 0;
  for(int rp=0; rp<DIM; ++rp)
  {
    result *= n;
    result += (Point_handle)((coord[rp]+0.5)*n);
  }
  return result;
}

int GridList::set_points(const points &pos)
{
  if(n<4) return -1;
  Point_handle p = 0;
  Point_handle idx;
  grids = std::vector<std::vector<Point_handle> >(size);
  for(auto it=pos.begin(); it != pos.end(); ++it, ++p)
  {
    idx = pos2idx(*it);
    grids[idx].push_back(p);
  }
  return 0;
}

/* Update rcut for grid
 * return value: 0-grid does not change; 1-grid needs change; 2-grid not needed because rcut is large
 */
int GridList::update_rcut(double rcut_)
{
  std::size_t newn = (Point_handle)(1.0/rcut_);
  int retval;
  if(newn != n)
  {
    n = newn;
    if(n<4)
    {
      grids.clear();
      retval = 2;
    }
    else
    {
      size = 1;
      nnbr = 1;
      for(int rp=0;rp<DIM;++rp){size *= n; nnbr*=3;}
      l = 1./n;
      grids.clear();
      retval = 1;
    }
  }
  else
  {
    retval = 0;
  }
  return retval;
}

void GridList::get_neighbor(const point &coord, std::vector<Point_handle> &neighbors) const
{
  int rq;
  long grid[DIM], grid_new[DIM];
  Point_handle idx = pos2idx(coord), rp, idxb;
  idx2grid(idx, grid);
  neighbors.clear();
  for(rq=0; rq<DIM; ++rq) grid_new[rq] = grid[rq]-1;
  for(rp=0; rp<nnbr; ++rp, ++grid_new[0])
  {
    for(rq = 0; grid_new[rq]-grid[rq]==2; ++grid_new[++rq])
      grid_new[rq] = grid[rq] - 1;
    idxb = grid2idx(grid_new);
    if(grids[idxb].size() > 0)
      neighbors.insert(neighbors.end(), grids[idxb].begin(), grids[idxb].end());
  }
}

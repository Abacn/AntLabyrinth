//
//  qhullinterface.cpp
//  QhullPercolation
//
//  Created by Yi Hu on 3/22/19.
// all methods that involve quickhull.hpp

#include <iostream>
#include <cstdint>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_set>
#include "dim.h"

#include "graph_construct.hpp"
#include "myutility.hpp"
// only this file includes qhull lib
#include "quickhull.hpp"

using point_iterator_type = typename points::const_iterator;
using quick_hull_type = quick_hull< DIM, point_iterator_type >;

void create_qhull(quick_hull_type& quick_hull_t, const points &points_inv, const int32_t inv_cur_size)
{
  quick_hull_t.add_points(points_inv.cbegin(), points_inv.cbegin()+inv_cur_size);
  auto const initial_simplex_ = quick_hull_t.get_affine_basis();
  quick_hull_t.create_initial_simplex(initial_simplex_.cbegin(), std::prev(initial_simplex_.cend()));
  quick_hull_t.create_convex_hull();
}

// run quickhull
void Graph_Construct::exec_qhull(int32_t currentidx, Triangulation_Data &TData, const points &points_rel_, const points &points_inv_, const int32_t cur_size, const std::vector<Point_handle> &inv_idx_, void *qhullhandle)
{
  int32_t atmp[DIM+1], btmp[DIM+1], rtmp[DIM];
  unsigned long rp, rq, nneighb;
  double dtmp;
  Cell_handle cellidx;
  double const eps_ = 1e-14; // 100*machine epsilon
  int dealloc_tag = 0;
  quick_hull_type *qtp; // quick hull type pointer
  point circumcenter;
  if(nullptr == qhullhandle)
  {
#ifdef DEBUG
    clock_t t1, t2;
    t1 = clock();
#endif
    qtp = new quick_hull_type(eps_); // dimension, eps
    create_qhull(*qtp, points_inv_, cur_size);
    nneighb=0;
    qtp->facet_make_reduce(points_inv_.begin()+currentidx, &nneighb);
    count_neighbor += nneighb;
    count_neighbor2 += nneighb*nneighb;
    ev_point++;
    dealloc_tag = 1;
#ifdef DEBUG
    t2 = clock();
    qhulltime += (double)(t2-t1)/CLOCKS_PER_SEC;
#endif
  }
  else
  {
    qtp = static_cast<quick_hull_type*>(qhullhandle);
  }
  quick_hull_type &quick_hull_t = *qtp;
  auto idx0vert = points_inv_.cbegin();
  atmp[DIM] = currentidx;
  btmp[DIM] = inv_idx_[currentidx];
  // store indices of cells
  std::vector<Cell_handle> cellidxlist(quick_hull_t.facets_reduced_.size());
  // std::cout << quick_hull_t.facets_.size() << std::endl;
  rq = 0;
  for (auto const & facet_ : quick_hull_t.facets_reduced_) // traverse facets
  {
    // get facet point indices
    auto const & vertices_ = facet_.vertices_;
    rp = 0;
    for (auto const & vertex_ : vertices_)
    {
      int32_t idx = (int32_t)(vertex_ - idx0vert); // idx in inv
      if(idx == currentidx)
      {
        std::cerr << "Error: Particle " << currentidx << " has self neighbor\n";
      }
      atmp[rp++] = idx;
    }
    // convert relative index to absolute
    for(int i=0; i<DIM; ++i)
      btmp[i] = inv_idx_[atmp[i]];
    // create cell if not exist
    Base_Cell cell1(btmp, 0);
    // the cell may be linked at this time only if the index of central point is smallest or second-smallest
    // in this version, already removed in quickhull.facet_make_reduce
    // if(cell1.head() != inv_idx_[currentidx] && cell1.get_idx(1) != inv_idx_[currentidx])
    TriMapType::iterator it1 =  TData.cellmap.lower_bound(cell1);
    if(it1== TData.cellmap.end() || it1->first != cell1)
    {
      // not exist yet
      if(optflag)
      {
        // in partial run, a cell where current idx is second largest may not be added
        if(cell1.head() == inv_idx_[currentidx] || (!singletdataflag && cell1.get_idx(1) == inv_idx_[currentidx] ) )
        {
          dtmp = cell_weight(atmp, points_rel_, circumcenter);
          if(dtmp > wcut)
          {
            MyUtility::pbc(MyUtility::point_add(circumcenter, points_original[btmp[DIM]]));
            // add cell
            cellidx = TData.add_cell(it1, cell1, circumcenter);
          }
          else
          {
            // cell too small, do not add it
            if(cell1.head() == inv_idx_[currentidx])
            {
              // but count it
              cellidx = TData.add_cell();
            }
            else
            {
              cellidx = Null_cell_handle;
            }
          }
        }
        else
        {
          // Because cell is added during the smallest point scan.
          // If it does not exist in cellmap, it means this cell has skipped (too small)
          cellidx = Null_cell_handle;
        }
      }
      else
      {
        dtmp = cell_weight(atmp, points_rel_, circumcenter);
        MyUtility::pbc(MyUtility::point_add(circumcenter, points_original[btmp[DIM]]));
        // add cell
        cellidx = TData.add_cell(it1, cell1, circumcenter);
      }
    }
    else
    {
      cellidx = it1->second;
    }
    
    cellidxlist[rq] = cellidx;
    if(Null_cell_handle == cellidx)
    {
      ++rq;
      continue; // too small cell, just skip
    }
    // scan neighbors
    int rs = -1;
    for(unsigned long rr : facet_.neighbours_)
    {
      ++rs;
      if(rr < rq)
      {
        // connect cell1 and cell2 iff the common facet is smallest
        Cell_handle hcell2 = cellidxlist[rr];
        if(hcell2 == Null_cell_handle)
        {
          // cell 2 is too small, just skip
          TData.add_facet();
          continue;
        }
        // construct neighbor facet
        for(int rt=0; rt<DIM; ++rt)
        {
          if(rt == rs)
          {
            rtmp[rt] = currentidx;
          }
          else
          {
            rtmp[rt] = (int32_t)(facet_.vertices_[rt] - idx0vert);
          }
        }
        Base_Cell::sort_facet_indices(rtmp);
        if(rtmp[0] == currentidx)             // absolute current index
        {
          // connect cell
          dtmp = facet_weight(rtmp, points_rel_);
          if(optflag && dtmp<wcut)
          {
            TData.add_facet();
          }
          else
          {
            TData.add_facet(dtmp, cellidxlist[rq], cellidxlist[rr]);
          }
        }
      }
    }
    ++rq;
  }
  if(dealloc_tag)
  {
    delete qtp;
  }
}

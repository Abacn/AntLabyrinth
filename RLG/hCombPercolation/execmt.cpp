//
//  execmt.cpp
//  D4Percolation
//
//  Created by Yi Hu on 6/5/19.
//

#ifdef _MULTITHREAD

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdint>
#include <iomanip>
#include <iterator>
#include <deque>
#include <map>
#include <sstream>
#include <string>
#include <unordered_set>
#include <future>
#include "dim.h"

#include "graph_construct.hpp"
#include "myutility.hpp"
// only this file includes qhull lib
#include "quickhull.hpp"

// check frequency in millisecond
#if DIM<=6
#define ASYNC_CHK_FREQ 1
#elif DIM==7
#define ASYNC_CHK_FREQ 5
#elif DIM==8
#define ASYNC_CHK_FREQ 20
#elif DIM==9
#define ASYNC_CHK_FREQ 50
#else
#define ASYNC_CHK_FREQ 200
#endif


using point_iterator_type = typename points::const_iterator;
using quick_hull_type = quick_hull< DIM, point_iterator_type >;

int Graph_Construct::async_fun(int32_t idx, int rp)
{
  double const eps_ = 1e-14; // 100*machine epsilon
  void create_qhull(quick_hull_type& quick_hull_t, const points &points_inv, const int32_t inv_cur_size);
  prs[rp] = calcinv(idx, point_rel_omp[rp], point_inv_omp[rp], cur_sizes[rp], inv_idxs[rp]);
  if(prs[rp]+DIM<=cur_sizes[rp])
  {
    quick_hull_type qhull(eps_);
    create_qhull(qhull, point_inv_omp[rp], cur_sizes[rp]);
    unsigned long nneighb = 0;
    qhull.facet_make_reduce(point_inv_omp[rp].begin()+prs[rp], &nneighb); // save memory
    // incrementing global variable in thread
    count_neighbor += nneighb;
    count_neighbor2 += nneighb*nneighb;
    ev_point++;
    exec_qhull(prs[rp], tdatas[rp], point_rel_omp[rp], point_inv_omp[rp], cur_sizes[rp], inv_idxs[rp], &qhull);
  }
  return rp;
}

void Graph_Construct::buildTriangulation(Triangulation_Data &TData, int32_t startidx, int32_t stepidx)
{
  int32_t endidx=input.N-DIM;
  int num_run_points = ceil((endidx+1.0-startidx)/stepidx);
  int n_thread=std::min(input.n_thread, num_run_points);
  int32_t pa, pb, rp;
  pa = startidx;
  std::deque< std::future<int> > async_queues;
  tdatas = std::vector<Triangulation_Data>(n_thread);
  if(n_thread==1 && input.outmode != OUTPUT_MODE_SAVEGRAPH)
  {
    singletdataflag = true;
  }
  else
  {
    singletdataflag = false;
  }
  // allocate inverse coordinates
  for(rp=0; rp<n_thread; ++rp)
  {
    async_queues.push_back(std::async(std::launch::async, &Graph_Construct::async_fun, this, pa, rp));
    pa+=stepidx;
  }
#if GRAPH_CONSTRUCT_QHULL_SKIPNUM>10
  int32_t plastadd=0;
  std::vector<uint64_t> nlastfacets(n_thread, 0); // trace the facet size
#endif
#if 7 <= DIM
  std::string statusfile = std::string("status_")+std::to_string(startidx)+".dat";
  auto fst = std::ofstream(statusfile); // create an empty file
  fst.close();
#endif
  bool continuing_flag = true;
  auto iterfirst = async_queues.begin();
  // Now pa is the next async job
  pb = 0;
  while(!async_queues.empty())
  {
    while(1)
    {
      // repeatedly check thread ends
      if(iterfirst == async_queues.end()) iterfirst = async_queues.begin();
      auto waitresult = iterfirst->wait_for(std::chrono::milliseconds(ASYNC_CHK_FREQ));
      if(waitresult == std::future_status::ready) break;
      ++iterfirst;
    }
    rp = iterfirst->get();
#if GRAPH_CONSTRUCT_QHULL_SKIPNUM>10
    uint64_t fsize = tdatas[rp].facetlist.size();
    if(fsize>nlastfacets[rp]) {plastadd=pb;nlastfacets[rp]=fsize;}
    else if(pb-plastadd>GRAPH_CONSTRUCT_QHULL_SKIPNUM) {continuing_flag = false;} // very unlikely there will be new facets
#endif
    iterfirst = async_queues.erase(iterfirst);
    if(continuing_flag && pa<=endidx)
    {
      async_queues.push_back(std::async(std::launch::async, &Graph_Construct::async_fun, this, pa, rp));
      iterfirst = async_queues.begin();
      pa+=stepidx;
    }
#if 7 <= DIM
    // display status when dimension is high
    if(pb % commondef_h_STATUS_INTERVAL == 0)
    {
      fst.open(statusfile);
      fst << "now/total\tmaxnbrdis\tNneighbor\n";
      fst << pb << '/' << num_run_points << '\t' << maxnbrdis << '\t' << (double)count_neighbor/ev_point << '\n';
      fst.close();
    }
    ++pb;
  }
  remove(statusfile.c_str());
#else
    ++pb;
  }
#endif
  // deallocate inverse coordinates
  for(int rp=0; rp<input.n_thread; ++rp) inv_idxs[rp].clear();
  for(int rp=0; rp<input.n_thread; ++rp) point_rel_omp[rp].clear();
  for(int rp=0; rp<input.n_thread; ++rp) point_inv_omp[rp].clear();
  TData = std::move(tdatas[0]);
  if(!singletdataflag)
  {
    uint64_t cellsize=TData.cellsize, facetsize=TData.facetsize, dupcells=0;
    for(rp=1; rp<n_thread; ++rp)
    {
      // add cell
      Cell_handle cellidx, scidx; // cell idx in TData and in tdatas[rp]
      std::map<Cell_handle, Cell_handle> chmap;
      std::vector<Base_Cell> bcells(tdatas[rp].celllist.size());
      for(const auto &val : tdatas[rp].cellmap)
      {
        bcells[val.second] = val.first;
      }
      for(scidx=0; scidx<tdatas[rp].celllist.size(); ++scidx)
      {
        auto &fc = bcells[scidx];
        TriMapType::iterator it1 = TData.cellmap.lower_bound(fc);
        if(it1== TData.cellmap.end() || it1->first != fc)
        {
          // not exist yet
          cellidx = TData.add_cell(it1, fc, tdatas[rp].celllist[scidx]);
        }
        else
        {
          // already exist
          cellidx = it1->second;
          ++dupcells;
        }
        chmap.insert({scidx, cellidx});
      }
      // add facets
      Cell_handle cell1, cell2;
      for(auto const &fct: tdatas[rp].facetlist)
      {
        cell1 = chmap[fct.cell1];
        cell2 = chmap[fct.cell2];
        TData.add_facet(fct.weight, cell1, cell2);
      }
      cellsize += tdatas[rp].cellsize;
      facetsize += tdatas[rp].facetsize;
      tdatas[rp].clear();
    }
    TData.cellsize = cellsize-dupcells;
    TData.facetsize = facetsize;
  }
  tdatas.clear();
  // sort facet by circumradius (descend)
  std::sort(TData.facetlist.rbegin(), TData.facetlist.rend());
}
#endif

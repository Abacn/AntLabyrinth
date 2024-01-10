//
//  graph_common.cpp
//  QhullPercolation
//
//  Created by Yi Hu on 3/26/19.
//
// This file contains the methods related to Graph_Construct that is the same for different targets.
// For convenience of copy-paste between different targets
// Because the class definition with the same name in different targets can be different

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
//#ifdef _MULTITHREAD
//#include <omp.h>
//#endif
#include "graph_construct.hpp"
#include "cellset.hpp"
#include "myutility.hpp"
#include "dim.h"

// Facet weight is the square of circumradius
double Graph_Construct::facet_weight(const int32_t point_arr[DIM], const points& points_rel_)
{
  int rp, rq;
  double weight, dtmp;
  double matA[DIM][DIM], vecb[DIM], vecx[DIM];
//#pragma omp parallel for private(rp,rq)
  for(rp=0; rp<DIM; ++rp)
  {
    vecb[rp] = 1.;
  }
  // build augmented matrix
  for(rp=0; rp<DIM; ++rp)
  {
    matA[rp][rp] = 0.;
    for(rq=rp+1; rq<DIM; ++rq)
    {
      dtmp = MyUtility::norm_squared_distance(points_rel_[point_arr[rp]], points_rel_[point_arr[rq]]);
      matA[rp][rq] = matA[rq][rp] = dtmp;
    }
  }
  MyUtility::gauss_eliminate((double*)matA, vecb, vecx, DIM);
  dtmp = 0.;
  for(rp=0; rp<DIM; ++rp) dtmp+=vecx[rp];
  weight = 0.5/dtmp;
  return weight;
}

// get cell weight and circumcenter. The weight is the square of circumradius
// points_rel_[poin_arr[DIM]] is zero
double Graph_Construct::cell_weight(const int32_t point_arr[DIM+1], const points& points_rel_, point &circumcenter)
{
  int rp, rq;
  double weight, dtmp;
  double matA[DIM+1][DIM+1], vecb[DIM+1], vecx[DIM+1];
//#pragma omp parallel for private(rp,rq)
  for(rp=0; rp<=DIM; ++rp)
  {
    vecb[rp]=1.;
  }
  // build augmented matrix
  for(rp=0; rp<=DIM; ++rp)
  {
    matA[rp][rp] = 0.;
    for(rq=rp+1; rq<=DIM; ++rq)
    {
      dtmp = MyUtility::norm_squared_distance(points_rel_[point_arr[rp]], points_rel_[point_arr[rq]]);
      // record the largest neighbor distance
      if(maxnbrdis<dtmp)
      {
        maxnbrdis = dtmp;
        adapt_rcut(false);
      }
      // add debug info - neighbor distribution
      matA[rp][rq] = matA[rq][rp] = dtmp;
    }
  }
  MyUtility::gauss_eliminate((double*)matA, vecb, vecx, DIM+1);
  dtmp = 0.;
  for(rp=0; rp<=DIM; ++rp) dtmp+=vecx[rp];
  dtmp = 1.0/dtmp;
  weight = 0.5*dtmp;
  for(rp=0; rp<DIM; ++rp)
  {
    circumcenter[rp] = 0.;
    for(rq=0; rq<DIM; ++rq)
    {
      circumcenter[rp] += vecx[rq]*points_rel_[point_arr[rq]][rp];
    }
    circumcenter[rp] *= dtmp;
  }
  
  if(weight>maxcellr) maxcellr=weight;
  return weight;
}

void Graph_Construct::counterReset()
{
  maxnbrdis = 0.;
  count_neighbor = count_neighbor2 = ev_point = 0;
  maxcellr = 0.;
}

// generate random distribution
void Graph_Construct::genrandomPoints()
{
  long i=0, j=0;
  double dtmp;
  for(i=0; i<input.N; ++i)
  {
    for(j=0; j<DIM; ++j)
    {
      dtmp = MyUtility::rand()*2.0 - 1.0;
      points_original[i][j] = dtmp;
    }
    MyUtility::pbc(points_original[i]);
  }
}

#ifndef _MULTITHREAD
void Graph_Construct::buildTriangulation(Triangulation_Data &TData, int32_t startidx, int32_t stepidx)
{
  int32_t pa, pr;
#if GRAPH_CONSTRUCT_QHULL_SKIPNUM>10
  int32_t plastadd=0;
  uint64_t nlastfacet=0; // trace the facet size
#endif
  // -DIM because for the last N points, all cells must have been scanned
  int32_t endidx=input.N-DIM;
#if 7 <= DIM
  int num_run_points = ceil((endidx+1.0-startidx)/stepidx);
  std::string statusfile = std::string("status_")+std::to_string(startidx)+".dat";
  std::ofstream fst(statusfile); // create an empty file
  fst.close();
#endif
  if(input.outmode == OUTPUT_MODE_SAVEGRAPH)
  {
    singletdataflag = false;
  }
  else
  {
    singletdataflag = true;
  }
  for(pa=startidx; pa<=endidx; pa+=stepidx)
  {
#ifdef DEBUG
    clock_t t1, t2, t3;
    t1 = clock();
#endif
    pr = calcinv(pa, points_rel, points_inv, inv_cur_size, inv_idx);
#ifdef DEBUG
    t2 = clock();
#endif
    // new facet can be added only when center point is at most DIM largest
    if(pr+DIM>inv_cur_size) continue;
    exec_qhull(pr, TData, points_rel, points_inv, inv_cur_size, inv_idx);
#ifdef DEBUG
    t3 = clock();
    calcinvtime += (double)(t2-t1)/CLOCKS_PER_SEC;
    execqhulltime += (double)(t3-t2)/CLOCKS_PER_SEC;
#endif
#if GRAPH_CONSTRUCT_QHULL_SKIPNUM>10
    uint64_t fsize = TData.facetlist.size();
    if(fsize>nlastfacet) {plastadd=pa;nlastfacet=fsize;}
    else if((pa-plastadd)/stepidx>GRAPH_CONSTRUCT_QHULL_SKIPNUM) {break;} // very unlikely there will be new facets
#endif
#if 7 <= DIM
    // display status when dimension is high
    if((pa-startidx)/stepidx % commondef_h_STATUS_INTERVAL == 0)
    {
      fst.open(statusfile);
      fst << "now/total\tmaxnbrdis\tNneighbor\n";
      fst << (pa-startidx)/stepidx << '/' << num_run_points << '\t' << maxnbrdis << '\t' << (double)count_neighbor/ev_point << '\n';
      fst.close();
    }
  }
  remove(statusfile.c_str());
#else
  }
#endif
  std::sort(TData.facetlist.rbegin(), TData.facetlist.rend());
}
#endif

// evaluate percolation
int Graph_Construct::eval_percolation(Triangulation_Data &TData)
{
  const int n_cri = 5;
  double weights[n_cri] = {0.0};
  bool pcflag[n_cri] = {false};
  double dtmp, dtmp2;
  Cell_disjoint_set cdset(TData.celllist);
  std::size_t countfacet = 0;
  std::uint32_t percoflag;
  std::vector<Facet_decorate>::iterator fit;
  
  for(fit = TData.facetlist.begin(); !pcflag[n_cri-1] && fit != TData.facetlist.end(); ++fit)
  {
    percoflag = cdset.merge(*fit);
    ++countfacet;
    
    for(int rp=0; rp<n_cri; ++rp)
    {
#define SET_PERCO(condition) if(condition){weights[rp] = fit->weight;pcflag[rp] = true;}
      if(pcflag[rp]) continue;
      switch(rp)
      {
        case 0:   // 0: any cluster percolates in any d;
          SET_PERCO(percoflag)
          break;
        case 1: // 1: any cluster percolates in first d;
          SET_PERCO(percoflag & 1)
          break;
        case 2: // 2: any cluster percolates in second d;
          SET_PERCO(percoflag & 2)
          break;
        case 3: // 3: the system is percolated in all d;
          SET_PERCO(cdset.isPerco_all())
          break;
        case 4: // 4: there exists a cluster percolates in all d
          SET_PERCO(percoflag == Perco_flag::perco_all)
          break;
      }
#undef SET_PERCO
    }
  }
  // tracer radius is 1
  if(countfacet < TData.facetlist.size())
  {
    threshold = weights[n_cri-1];
  }
  else
  {
    
      std::cout << "  Sample " << count+1 << " not percolated for all considered facet (num=" << countfacet <<
#if RUN_FAIL_OPTION==1
      ")...\n  Loose wcut" << std::endl;
#else
      ")...\n  Disabling optimization instead" << std::endl;
#endif
      return 1;
  }
  ++count;
  int retval = graphsum(TData, countfacet, cdset.get_cell_size());
  if(retval)
  {
    std::cout << "Error: Cannot write to output file" << std::endl;
    exit(-1);
  }
  if(input.outmode >= OUTPUT_MODE_PRODUCT)
  {
    std::stringstream dbgoutstr;
    // actual facetlist size and celllist size, number of evaluated points
    dbgoutstr << TData.facetlist.size() << "\t" << TData.celllist.size() << "\t" << ev_point;
    for(int rp=0; rp<n_cri; ++rp)
    {
      // percolation threshold in different critera
      dbgoutstr << "\t" << weight2Phi(weights[rp]);
    }
    // percolation cluster size, number of clusters
    dbgoutstr << "\t" << cdset.getmax_clu_size() << "\t" << cdset.getclustercount();
    // average number of neighbors <n>, and its variance <n^2>-<n>^2
    dtmp = (double)count_neighbor/ev_point;
    dtmp2 = (double)count_neighbor2/ev_point;
    dbgoutstr << "\t" << std::setprecision(8) << dtmp << "\t" << dtmp2-dtmp*dtmp;
    // maximum cell square circumradius
    dbgoutstr << "\t" << maxcellr;
    // maximum square neighbor distance
    dbgoutstr << "\t" << maxnbrdis;
    graphsumextend(TData, dbgoutstr.str().c_str());
  }
  if(count==1 && defaultoptflag)
  {
    adapt_rcut(true);
    std::cout << "Adjusted rcut to " << rcut << std::endl;
  }
  return 0;
}

// calculate inverse vectors of point idx
int32_t Graph_Construct::calcinv(int32_t idx, points &points_rel_, points &points_inv_, int32_t &inv_cur_size_, std::vector<Point_handle> &inv_idx_) const
{
  int32_t i=0, j=0, rp=0;
  int32_t rela_idx=0;
  if(0 == inv_idx_.size()) // cleared
  {
    inv_idx_.resize(input.N);
    points_rel_.resize(input.N);
    points_inv_.resize(input.N);
  }
  for(i=0; i<input.N; ++i)
  {
    if(i==idx)
    {
      for(j=0; j<DIM; ++j) points_rel_[rp][j] = points_inv_[rp][j] = 0.; // as origin
      rela_idx = rp;
      inv_idx_[rp++] = i;
    }
    else
    {
      for(j=0; j<DIM; ++j)
      {
        points_rel_[rp][j] = points_original[i][j] - points_original[idx][j];
      }
      auto nbs = MyUtility::get_point_images(points_rel_[rp], rcut);
      while(rp+nbs.size() >= inv_idx_.size())
      {
        inv_idx_.resize(2*inv_idx_.size());
        points_rel_.resize(2*points_rel_.size());
        points_inv_.resize(2*points_inv_.size());
      }
      for(auto const &nb: nbs) // within the cutoff (square) radius
      {
        points_rel_[rp] = nb.first;
        for(j=0; j<DIM; ++j)
        {
          points_inv_[rp][j] = points_rel_[rp][j] / nb.second;
        }
        inv_idx_[rp++] = i;
      }
    }
  }
  inv_cur_size_ = rp;
  return rela_idx;
}

void Graph_Construct::adapt_rcut(bool is_decr)
{
  //double target_rcut = facet_max_w*3.0*pow(2., 2.0/DIM); // 3.0 is empirical
  // approximate the appropriate radius of outer shell
  double target_rcut;
  if(count>0)
    target_rcut = maxnbrdis*pow(1.05, 2.0/DIM); //  5% more
  else if(input.N >= 1000)
    target_rcut = maxnbrdis*pow(1.15, 2.0/DIM); // 15% more
  else
    target_rcut = maxnbrdis*pow(1.25, 2.0/DIM); // 25% more
  if(target_rcut>1.5) target_rcut=1.5;
  if((is_decr || target_rcut>rcut) && target_rcut>wcut)
  {
    rcut = target_rcut;
  }
}

// transform weight to Phi (=number density*Unit sphere volume)
double Graph_Construct::weight2Phi(double weight)
{
  return input.N / DNVCELL * pow(weight, 0.5*DIM) * Vunitsphere;
}

// transform Phi to weight
double Graph_Construct::Phi2weight(double Phi)
{
  return pow(DNVCELL * Phi / input.N / Vunitsphere, 2./DIM);
}

// file names
const char* Graph_Construct::pointfile = "points.dat";
const char* Graph_Construct::tdataprefix = "tdata";
const char* Graph_Construct::tdatasuffix = ".dat";


int Graph_Construct::savepoints() const
{
  std::ofstream ofs(pointfile, std::ios::out | std::ios::binary);
  if(!ofs) return -1;
  for(auto const &p: points_original)
  {
    ofs.write((char*)p.data(), sizeof(point));
  }
  return 0;
}

int Graph_Construct::loadpoints()
{
  std::ifstream ifs(pointfile, std::ios::in | std::ios::binary);
  if(!ifs) return -1;
  for(auto const &p: points_original)
  {
    ifs.read((char*)p.data(), sizeof(point));
    if(!ifs) return -2; // read error
  }
  return 0;
}

int Graph_Construct::savestate(const Triangulation_Data &TData, const char *midfix) const
{
  std::stringstream sst;
  sst << tdataprefix << "_" << midfix << ".dat";
  std::string fdname = sst.str();
  std::cout << "Save to " << fdname << std::endl;
  std::ofstream ofs(fdname, std::ios::out | std::ios::binary);
  uint64_t sizebuf;
  double dbbuf;
  // save size info
  sizebuf = TData.celllist.size();
  ofs.write((char*)&sizebuf, sizeof(uint64_t));
  sizebuf = TData.facetlist.size();
  ofs.write((char*)&sizebuf, sizeof(uint64_t));
  sizebuf = TData.cellsize;
  ofs.write((char*)&sizebuf, sizeof(uint64_t));
  sizebuf = TData.facetsize;
  ofs.write((char*)&sizebuf, sizeof(uint64_t));
  // save cells
  std::vector<Base_Cell> bcells(TData.celllist.size());
  for(const auto &val : TData.cellmap)
  {
    bcells[val.second] = val.first;
  }
  for(Cell_handle rp=0; rp<TData.celllist.size(); ++rp)
  {
    ofs.write((char*)(bcells[rp].get_idxs()), sizeof(int32_t)*(DIM+1));
    ofs.write((char*)(TData.celllist[rp].data()), sizeof(point));
  }
  // save facets
  for(const auto &facet: TData.facetlist)
  {
    ofs.write((char*)&facet.cell1, sizeof(Cell_handle));
    ofs.write((char*)&facet.cell2, sizeof(Cell_handle));
    ofs.write((char*)&facet.weight, sizeof(double));
  }
  // save counters
  sizebuf = count_neighbor;
  ofs.write((char*)&sizebuf, sizeof(uint64_t));
  sizebuf = count_neighbor2;
  ofs.write((char*)&sizebuf, sizeof(uint64_t));
  sizebuf = ev_point;
  ofs.write((char*)&sizebuf, sizeof(uint64_t));
  dbbuf = maxcellr;
  ofs.write((char*)&dbbuf, sizeof(double));
  dbbuf = maxnbrdis;
  ofs.write((char*)&dbbuf, sizeof(double));
  return 0;
}

int Graph_Construct::loadstate(Triangulation_Data &TData, const std::vector<std::string> &files)
{
  uint64_t cellsize=0, facetsize=0, dupcells=0;
  for(const auto &fdname: files)
  {
    std::cout << "Read from " << fdname << std::endl;
    std::ifstream ifs(fdname, std::ios::in | std::ios::binary);
    if(!ifs) return -1;
    uint64_t sizebuf, ncell, nfacet;
    double dbbuf;
    // load size variables
    ifs.read((char*)&ncell, sizeof(uint64_t));
    ifs.read((char*)&nfacet, sizeof(uint64_t));
    ifs.read((char*)&sizebuf, sizeof(uint64_t));
    cellsize += sizebuf;
    ifs.read((char*)&sizebuf, sizeof(uint64_t));
    facetsize += sizebuf;
    // load cell
    int32_t cellbuf[DIM+1];
    point pointbuf;
    Cell_handle cellidx;
    std::map<Cell_handle, Cell_handle> chmap;
    for(uint64_t rp=0; rp<ncell; ++rp)
    {
      ifs.read((char*)&cellbuf, sizeof(int32_t)*(DIM+1));
      ifs.read((char*)pointbuf.data(), sizeof(point));
      Base_Cell fc(cellbuf, 2); // option 2 do not sort index: already sorted
      TriMapType::iterator it1 = TData.cellmap.lower_bound(fc);
      if(it1== TData.cellmap.end() || it1->first != fc)
      {
        // not exist yet
        cellidx = TData.add_cell(it1, fc, pointbuf);
      }
      else
      {
        // already exist
        cellidx = it1->second;
        ++dupcells;
      }
      chmap.insert({rp, cellidx});
    }
    // load facet
    Cell_handle cell1, cell2;
    for(uint64_t rp=0; rp<nfacet; ++rp)
    {
      ifs.read((char*)&cell1, sizeof(Cell_handle));
      ifs.read((char*)&cell2, sizeof(Cell_handle));
      ifs.read((char*)&dbbuf, sizeof(double));
      cell1 = chmap[cell1];
      cell2 = chmap[cell2];
      TData.add_facet(dbbuf, cell1, cell2);
    }
    // load counter
    ifs.read((char*)&sizebuf, sizeof(uint64_t));
    count_neighbor += sizebuf;
    ifs.read((char*)&sizebuf, sizeof(uint64_t));
    count_neighbor2 += sizebuf;
    ifs.read((char*)&sizebuf, sizeof(uint64_t));
    ev_point += sizebuf;
    ifs.read((char*)&dbbuf, sizeof(double));
    if(maxcellr<dbbuf) maxcellr = dbbuf;
    ifs.read((char*)&dbbuf, sizeof(double));
    if(maxnbrdis<dbbuf) maxnbrdis = dbbuf;
    if(!ifs) return -2;
  }
  TData.cellsize = cellsize-dupcells;
  TData.facetsize = facetsize;
  std::sort(TData.facetlist.rbegin(), TData.facetlist.rend());
  std::cout << "Load states successfully" << std::endl;
  return 0;
}

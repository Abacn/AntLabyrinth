//
//  graph_construct.hpp
//  QhullPercolation
//
//  Created by Yi Hu on 3/21/19.

#ifndef graph_construct_hpp
#define graph_construct_hpp

#include <fstream>
#include <cstddef>
#include <vector>
#include <map>

// multithread header
#ifdef _MULTITHREAD
#include <atomic>
#endif

// user header
#include "dim.h"
#include "read_input.hpp"
#include "commondef.hpp"

// Optimization: if within this number of convex hull of points evaluated and no new facet is added,
// stop building convex hull. May save 10% total time. Now only used in DIM>=7
// Set to 0 to hide this optimization
#if DIM>=9
#define GRAPH_CONSTRUCT_QHULL_SKIPNUM 20
#elif DIM==8
#define GRAPH_CONSTRUCT_QHULL_SKIPNUM 50
#elif DIM==7
#define GRAPH_CONSTRUCT_QHULL_SKIPNUM 80
#elif DIM==6
#define GRAPH_CONSTRUCT_QHULL_SKIPNUM 100
#endif

struct Triangulation_Data
{
  TriArrType celllist;
  TriMapType cellmap;
  std::vector<Facet_decorate> facetlist;
  std::size_t cellsize, facetsize;
  Triangulation_Data(): cellsize(0), facetsize(0){}
  void add_facet(void); // formally add (small and never expected to visit) a facet
  void add_facet(double weight, Cell_handle cell1, Cell_handle cell2); // add true facet
  Cell_handle add_cell(void);
  Cell_handle add_cell(const Full_Cell &fc);
  Cell_handle add_cell(const Base_Cell &bc, const point &circcumcenter);
  Cell_handle add_cell(TriMapType::const_iterator hint, const Full_Cell &fc);
  Cell_handle add_cell(TriMapType::const_iterator hint, const Base_Cell &bc, const point &circcumcenter);
  void clear(){celllist.clear(); cellmap.clear(); facetlist.clear();}
};

class Graph_Construct
{
public:
  Graph_Construct(ReadInputRLG input);
  ~Graph_Construct();
  int graphrun(const bool rerunflag=false);
  int partialrun(int32_t startx);
  int graphdebug();
private:
  // methods
  // generate randomly distributed points
  void counterReset();
  void genrandomPoints();
  void buildTriangulation(Triangulation_Data &TData, int32_t startidx=0, int32_t stepidx=1);
  int eval_percolation(Triangulation_Data &TData);
  // calculate inverse vectors of point idx, const tag = thread safe
  int32_t calcinv(int32_t idx, points &points_rel_, points &points_inv_, int32_t &inv_cur_size_, std::vector<Point_handle> &inv_idx_) const;
  // using qhull to build the network
  void exec_qhull(int32_t currentidx, Triangulation_Data &TData, const points &points_rel_, const points &points_inv_, const int32_t cur_size, const std::vector<Point_handle> &inv_idx_, void *qhullhandle=nullptr);
  // adapt rcut if overflow
  void adapt_rcut(bool is_decr);
  // calculating facet weight
  double facet_weight(const int32_t point_arr[DIM], const points& points_rel_);
  double cell_weight(const int32_t point_arr[DIM+1], const points& points_rel_, point& circumcenter);
  int graphsum(const Triangulation_Data &TData, std::size_t vstfacet, std::size_t vstcell);
  double weight2Phi(double weight);
  double Phi2weight(double Phi);
  // extend output
  void graphsumextend(const Triangulation_Data &TData, const char *outputline);
  int savepoints() const;
  int loadpoints();
  int savestate(const Triangulation_Data &TData, const char *midfix) const;
  int loadstate(Triangulation_Data &TData, const std::vector<std::string> &files);
  // members
  ReadInputRLG input;
  std::ofstream fsummary;
  double Vunitsphere;
  double threshold;
  long count;
  points points_original;
#ifndef _MULTITHREAD
  points points_rel;        // shifted point vectors
  points points_inv;        // inverse point vectors
  std::vector<Point_handle> inv_idx;    // maintain original indices in the re-ordered points_inv
  int32_t inv_cur_size;
  // count the number of neighbors of a point, and its variance
  uint64_t count_neighbor, count_neighbor2, ev_point;
  double rcut, maxnbrdis, maxcellr;
#else
  int32_t *prs;
  std::vector<points> point_rel_omp;
  std::vector<points> point_inv_omp;
  std::vector<std::vector<Point_handle> > inv_idxs;
  std::vector<Triangulation_Data> tdatas;
  int32_t *cur_sizes;
  int async_fun(int32_t idx, int rp);
  // thread safe counter
  std::atomic<uint64_t> count_neighbor, count_neighbor2, ev_point;
  // maximum neighbor distance, maximum cell square circumradius
  std::atomic<double> rcut, maxnbrdis, maxcellr;
#endif
  double wcut;
  bool optflag, defaultoptflag;
  bool singletdataflag; // if the system is running on single TData
  // debug info
  std::ofstream fdbgsummary;
  // last element is used to prevent overflow only
  // const members
  const static char* pointfile;
  const static char* tdataprefix;
  const static char* tdatasuffix;
};

#endif /* graph_construct_hpp */

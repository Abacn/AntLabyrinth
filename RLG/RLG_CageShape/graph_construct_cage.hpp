//
//  graph_construct.hpp
//  RLG_Delaunay
//


#ifndef graph_construct_hpp
#define graph_construct_hpp

#include <fstream>
#include <iomanip>
#include <vector>
#include <iterator>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <stdint.h>

#include <Eigen/Core>
#include <Eigen/LU>

#include <CGAL/Dimension.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Triangulation.h>
#include <CGAL/Delaunay_triangulation.h>

#include "read_cageshape.hpp"
#include "cellsample.hpp"
#include "cagestat.hpp"
#include "dim.h"

typedef CGAL::Epick_d<CGAL::Dimension_tag<DIM> > DIM_tag;
typedef CGAL::Triangulation< DIM_tag > Triangulation;

typedef Triangulation::Point Point;
typedef Triangulation::Face Face;
typedef Triangulation::Facet Facet;
typedef Triangulation::Facet_iterator Facet_iterator;
typedef Triangulation::Vertex_handle Vertex_handle;
typedef Triangulation::Full_cell_iterator Cell_iterator;
typedef Triangulation::Full_cell_handle Full_cell_handle;
typedef std::vector<Full_cell_handle> Full_cells;

// store the status of percolation
// perco_mark[0]: touch the positive boundary;
// perco_mark[1]: touch the negative boundary;
// Each bit record if the void pathway approaches the boundary
// Percolation judgement: perco_mark[0] & perco_mark[1]
union Perco_flag{
  uint16_t perco_mark[2];   // is percolated in one direction
  uint32_t perco_all;       // convenient for merge two set
};

class Graph_Construct
{
public:
  Graph_Construct(ReadInputRLGCShape &input);
  ~Graph_Construct();
  int graphrun();
  void dumpall();
  double weight2Phi(double weight);
  double Phi2weight(double Phi);
private:
  // types
  using cellmap_type = std::unordered_map<Full_cell_handle, Cell_samples>;
  using cellmap_iter = cellmap_type::const_iterator;
  using celldata_type = std::vector<std::pair<Full_cell_handle, Cell_samples> >;
  // functions
  void construct_network();       // build facet network
  void insert_points();
  int graphdump(bool isperco);
  int graphsum(bool isperco);
  int process_facet(cellmap_iter cellit, const Full_cell_handle cellhdl);
  Cell_samples create_cell_sample(Full_cell_handle cell);
  // sample over one cluster
  void samplecluster(std::vector<Cell_single_sample> &allsample, double &clvolume, double &sqweight, double &simplexttlv,  double &Delta2, double &largest_cell_r, Vector &largest_cell_center);
  int samplecluster_v1(celldata_type &sampledata, std::vector<Cell_single_sample> &allsample, double &simplexttlv, double &clvolume, double &sqweight, Vector &possum, double &squaresum);
  int samplecluster_v2(celldata_type &sampledata, std::vector<Cell_single_sample> &allsample, double &simplexttlv, double &clvolume, double &sqweight, Vector &possum, double &squaresum);
  // members
  ReadInputRLGCShape &input;
  std::vector<Triangulation::Point> points;
  CGAL::Delaunay_triangulation<DIM_tag> tris;
  std::unordered_set<Full_cell_handle> boundary_cells;
  cellmap_type cellmap;
  std::uint32_t count;
  std::ofstream fout, fsummary, fdataA;
  double Vunitsphere;
  AllStat *allstats;
};

// decoration class containing additional information of a facet
class Facet_decorate
{
public:
  Facet_decorate();
  Facet_decorate(double weight, Full_cell_handle cell);
  ~Facet_decorate();
  Full_cell_handle cell;
  double weight;
  bool operator<(const Facet_decorate &facetB) const;
};


#endif /* graph_construct_hpp */

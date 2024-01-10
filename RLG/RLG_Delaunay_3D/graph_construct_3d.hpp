//
//  graph_construct_3d.hpp
//  RLG_Delaunay_3D
//


#ifndef graph_construct_3d_hpp
#define graph_construct_3d_hpp

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <iterator>
#include <map>
#include <unordered_map>
#include <stdint.h>
#include <math.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Random.h>
#include <CGAL/function_objects.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Timer.h>

#include "parameters.hpp"
#include "read_input.hpp"
#include "dim.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Triangulation_3<K>      Triangulation;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K>   Gt;
typedef CGAL::Periodic_3_Delaunay_triangulation_3<Gt>         P3DT3;
typedef CGAL::Vector_3<CGAL::Epick>    Vector;
typedef P3DT3::Point             Point;
typedef P3DT3::Offset            Offset;
typedef P3DT3::Iso_cuboid        Iso_cuboid;
typedef P3DT3::Vertex_handle     Vertex_handle;
typedef P3DT3::Cell_handle       Cell_handle;
typedef P3DT3::Cell_iterator     Cell_iterator;
typedef P3DT3::Locate_type       Locate_type;

typedef CGAL::Creator_uniform_3<double, Point> Creator;

// utility functions
Vector dif_points_PBC(const Point &point1, const Point &point2);
Vector dif_cells_PBC(const Cell_handle cell1, const Cell_handle cell2);
Vector dif_vec_PBC(Vector &vec);

// flag of percolation
class Perco_flag{
public:
  Perco_flag(): perco(0){}
  ~Perco_flag(){}
  // raw percolation status
  uint32_t getPerco()
  {
    return perco;
  }
  // percolation status at one specific dimension
  bool isPerco(int dim)
  {
    return (perco >> dim) & 1;
  }
  // percolation status for all dimensions
  bool isPerco_all()
  {
    return perco == perco_all;
  }
  // set percolation flag at a specific dimension
  void setPerco(uint32_t percobits)
  {
    perco |= percobits;
  }
  static const std::uint32_t perco_all=(1<<DIM)-1;
private:
  uint32_t perco;       // convenient for merge two set
};

// For recording the results
class Result_record{
private:
  double perco_phi;
  long n_allfacet;
  long n_checkedfacet;
  long n_allcell;
  long n_checkedcell;
public:
  Result_record();
  void set_value(double s_perco_phi, long s_n_allfacet, long s_n_checkedfacet, long s_n_allcell, long s_n_checkedcell);
  void print_value(std::ofstream &fout);
};

// decoration class containing additional information of a facet
class Facet_decorate
{
public:
  Facet_decorate(double weight, const Cell_handle cell1, const Cell_handle cell2);
  ~Facet_decorate();
  bool operator<(const Facet_decorate&fb) const;
  Cell_handle cell1, cell2;
  double weight;
private:
};

class Graph_Construct_3d
{
public:
  Graph_Construct_3d(ReadInputRLG input);
  ~Graph_Construct_3d();
  int graphrun();
  int graphtest();
  const Result_record get_results();
protected:
  double facet_weight(const std::vector<Vector> &offpoints, int covertex) const;
  double norm_squared_distance(const Vector &pa, const Vector &pb) const;
  double weight2rho(double weight) const;
  double weight2Phi(double weight) const;
private:
  // functions
  using Facetlist = typename std::vector<Facet_decorate>;
  void traverse_facet(Facetlist &facetlist);
  int graphdump();
  int graphdebug(const char *line);
  // members
  double Vunitsphere;
  ReadInputRLG input;
  std::vector<Point> points;
  P3DT3 tris;
  long count;
  std::ofstream fout, fdbgsummary;
  Result_record result_record;
};

// disjoint set containing full cell handles
class Cell_disjoint_set
{
public:
  uint32_t merge(const Facet_decorate &fc);
  std::size_t get_cell_size(){return cellset.size();}
  Cell_disjoint_set(): current_clu_size(0), max_clu_size(0){}
  ~Cell_disjoint_set(){}
  uint32_t getsysperco(){return syspercoflag.getPerco();}
  bool isPerco(int dim){return syspercoflag.isPerco(dim);}
  bool isPerco_all(){return syspercoflag.isPerco_all();}
  unsigned long getcurrent_clu_size(){return current_clu_size;}
  unsigned long getmax_clu_size(){return max_clu_size;}
  unsigned long getclustercount(){return clusterset.size();}
  bool is_current_max(){return current_clu_size==max_clu_size;}
private:
  // structure
  struct Cell_info{
    Cell_info();
    Cell_info(Cell_handle self, Cell_handle parent);
    Cell_handle parent;
    Vector vec_to_parent;
    long rank;
  };
  struct Cluster_info{
    Cluster_info();
    unsigned long size;
    Perco_flag percoflag;
  };
  using CellMapType = typename std::map<Cell_handle, Cell_info>;
  using ClusterMapType = typename std::unordered_map<Cell_handle, Cluster_info>;
  // member
  CellMapType cellset;
  ClusterMapType clusterset;
  Perco_flag syspercoflag;
  unsigned long current_clu_size, max_clu_size;
  // method
  CellMapType::iterator findset(CellMapType::iterator it);
  void mergecluster(Cell_handle cell1, Cell_handle cell2);
  void incrementcluster(Cell_handle cell1);
  void newcluster(Cell_handle cell1);
};

#endif /* graph_construct_3d_hpp */

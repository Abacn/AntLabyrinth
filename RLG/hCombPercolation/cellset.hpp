//
//  cellset.hpp
//  QhullPercolation
//
//  Created by Yi Hu on 3/25/19.

#ifndef cellset_hpp
#define cellset_hpp

#include <vector>
#include <map>
#include <unordered_map>
#include <cstdint>
#include "commondef.hpp"

using std::uint32_t;

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

// disjoint set containing full cell handles
class Cell_disjoint_set
{
public:
  std::uint32_t merge(const Facet_decorate &fc);
  std::size_t get_cell_size(){return cellset.size();}
  Cell_disjoint_set(const TriArrType &celllist): celllist(celllist), current_clu_size(0), max_clu_size(0){}
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
    Cell_info(Cell_handle self, Cell_handle parent, const Cell_disjoint_set& outer_inst);
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
  const TriArrType &celllist;
  Perco_flag syspercoflag;
  unsigned long current_clu_size, max_clu_size;
  // method
  CellMapType::iterator findset(CellMapType::iterator it);
  void mergecluster(Cell_handle cell1, Cell_handle cell2);
  void incrementcluster(Cell_handle cell1);
  void newcluster(Cell_handle cell1);
};

#endif /* cellset_hpp */

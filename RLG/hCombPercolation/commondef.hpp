//
//  commondef.hpp
//  QhullPercolation
//
//  Created by Yi Hu on 3/22/19.

#ifndef commondef_h
#define commondef_h

#include <cstdint>
#include <array>
#include <valarray>
#include <vector>
#include <map>
#include <cstring>
#include "dim.h"

#if DIM>=8
// Use Dn+ quantizer for d>=8 and even dimension, Dn for else.
#define DNPLUSFLAG
#if DIM%2!=0
// For odd dimension, the dual Dn is shifted by (1/2, 1/2, ..., 0)
#define DNLAMFLAG
#endif
#endif

using point = std::array< double, DIM >;
using points = std::vector< point >;
using Vector = std::valarray< double >;

// Base cell, only has indices stored
class Base_Cell
{
public:
  Base_Cell();
  // constructor. sort_arg_in: sort parse-in argument or not
  Base_Cell(int32_t idxs_unord[DIM+1], int sort_arg_in=0);
  bool get_common_facet(const Base_Cell &cellB, int32_t arr_out[DIM]) const;
  int32_t head() const; // first point index of the cell
  const int32_t* get_idxs() const{return (const int32_t*)idxs;}
  int32_t get_idx(int dim) const{return idxs[dim];}
  bool operator<(const Base_Cell &a) const;
  bool operator==(const Base_Cell &a) const;
  bool operator!=(const Base_Cell &a) const;
  // simple sorting for small array
  static void sort_indices(int32_t idxs[DIM+1])
  {
    int32_t ltmp;
    int i, j;
    for(i=0; i<DIM; ++i)
    {
      for(j=i+1; j<DIM+1; ++j)
      {
        if(idxs[i]>idxs[j])
        {
          ltmp = idxs[i];
          idxs[i] = idxs[j];
          idxs[j] = ltmp;
        }
      }
    }
  }
  // only sort facet
  static void sort_facet_indices(int32_t idxs[DIM])
  {
    int32_t ltmp;
    int i, j;
    for(i=0; i<DIM-1; ++i)
    {
      for(j=i+1; j<DIM; ++j)
      {
        if(idxs[i]>idxs[j])
        {
          ltmp = idxs[i];
          idxs[i] = idxs[j];
          idxs[j] = ltmp;
        }
      }
    }
  }
protected:
  int32_t idxs[DIM+1]; // indices of vertices, sorted
};

// full cell class
class Full_Cell: public Base_Cell
{
public:
  Full_Cell();
  // constructor. sort_arg_in: sort parse-in argument or not
  Full_Cell(int32_t idxs_unord[DIM+1], int sort_arg_in=0);
  Full_Cell(int32_t idxs_unord[DIM+1], const point &circumcenter, int sort_arg_in=0);
  void set_circumcenter(const point &circumcenter_){circumcenter = circumcenter_;}
  point get_circumcenter() const{return circumcenter;}
private:
  point circumcenter;
};

using Cell_handle = uint32_t; // Cell_handle is just an index
using Point_handle = int32_t; // Point_handle is just an index
extern const Cell_handle Null_cell_handle;

// Cell storage types
using TriArrType = typename std::vector<point>;
using TriMapType = typename std::map<Base_Cell, Cell_handle>;

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

// time consuming information
#ifdef DEBUG
#include <time.h>
extern double qhulltime, execqhulltime, calcinvtime, clustersampletime, clustercounttime;
#endif

// status information
#if DIM>=10
#define commondef_h_STATUS_INTERVAL 10
#elif DIM>=9
#define commondef_h_STATUS_INTERVAL 20
#else
#define commondef_h_STATUS_INTERVAL 100
#endif

#endif /* commondef_h */

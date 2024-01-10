//
//  commondef.cpp
//  QhullPercolation
//
//  Created by Yi Hu on 3/26/19.
//

#include "commondef.hpp"

const Cell_handle Null_cell_handle = ~(Cell_handle)0;

Base_Cell::Base_Cell(){}

// constructor. sort_arg_in: sort parse-in argument or not
Base_Cell::Base_Cell(int32_t idxs_unord[DIM+1], int sort_arg_in/*=0*/)
{
  if(0==sort_arg_in)
  {
    memcpy(idxs, idxs_unord, sizeof(int32_t)*(DIM+1));
    sort_indices(idxs);
  }
  else if(1==sort_arg_in)
  {
    sort_indices(idxs_unord);
    memcpy(idxs, idxs_unord, sizeof(int32_t)*(DIM+1));
  }
  else
  {
    memcpy(idxs, idxs_unord, sizeof(int32_t)*(DIM+1));
  }
}

Full_Cell::Full_Cell(): Base_Cell() {}

Full_Cell::Full_Cell(int32_t idxs_unord[DIM+1], int sort_arg_in/*=0*/):
Base_Cell(idxs_unord, sort_arg_in)
{}

Full_Cell::Full_Cell(int32_t idxs_unord[DIM+1], const point &circumcenter, int sort_arg_in/*=0*/):
Base_Cell(idxs_unord, sort_arg_in),
circumcenter(circumcenter)
{}

// return the first point index of the cell
int32_t Base_Cell::head() const
{
  return idxs[0];
}

// get the facet of two neighboring cell (if so) and store in arr_out, otherwise return false
// only first DIM values of arr_out is valid. Use DIM+1 in case of parse in two same cell
bool Base_Cell::get_common_facet(const Base_Cell &cellB, int32_t arr_out[DIM+1]) const
{
  int rp=0, rq=0, rr=0;
  while(rp<=DIM && rq<=DIM)
  {
    if(idxs[rp] == cellB.idxs[rq])
    {
      arr_out[rr] = idxs[rp];
      ++rp; ++rq; ++rr;
    }
    else if(idxs[rp] < cellB.idxs[rq])
    {
      ++rp;
    }
    else // idxs[rp] > idxs[rq]
    {
      ++rq;
    }
  }
  if(rr != DIM)
    return false;
  else
    return true;
}

bool Base_Cell::operator<(const Base_Cell &a) const
{
  for(int k=0; k<=DIM; k++)
  {
    if (idxs[k]==a.idxs[k])
      continue;
    else if (idxs[k]<a.idxs[k])
      return true;
    else
      return false;
  }
  return false;
  }

// ==
bool Base_Cell::operator==(const Base_Cell &a) const
{
  for(int k=0; k<=DIM; k++)
  {
    if (idxs[k]!=a.idxs[k])
      return false;
  }
  return true;
  // return std::equal(x, x+D, a.x); // slower
}

// != should be consistent with ==
bool Base_Cell::operator!=(const Base_Cell &a) const
{
  for(int k=0; k<=DIM; k++)
  {
    if (idxs[k]!=a.idxs[k])
      return true;
  }
  return false;
}

Facet_decorate::Facet_decorate(double weight, const Cell_handle cell1, const Cell_handle cell2):
cell1(cell1), cell2(cell2), weight(weight)
{
}

Facet_decorate::~Facet_decorate(){}

bool Facet_decorate::operator<(const Facet_decorate&fb) const
{
  return weight < fb.weight;
}

#ifdef DEBUG
#include <time.h>
double qhulltime, execqhulltime, calcinvtime, clustersampletime, clustercounttime;
#endif

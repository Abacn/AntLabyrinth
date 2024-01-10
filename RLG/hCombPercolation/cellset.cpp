//
//  cellset.cpp
//  QhullPercolation
//
//  Created by Yi Hu on 3/25/19.
#include <iostream>

#include "parameters.hpp"
#include "cellset.hpp"
#include "myutility.hpp"

// Created as parent
Cell_disjoint_set::Cell_info::Cell_info():
parent(Null_cell_handle), vec_to_parent(Vector(0.0, DIM)), rank(1)
{
}

// Created as child
Cell_disjoint_set::Cell_info::Cell_info(Cell_handle self, Cell_handle parent, const Cell_disjoint_set& outer_inst):
parent(parent), rank(0)
{
  const point &pa = outer_inst.celllist[self];
  const point &pb = outer_inst.celllist[parent];
  vec_to_parent = MyUtility::dif_points_PBC(pa, pb);
}

Cell_disjoint_set::Cluster_info::Cluster_info(): size(2) // link two cell, cluster size must be 2 at initial
{}

void Cell_disjoint_set::mergecluster(Cell_handle cell1, Cell_handle cell2)
{
  clusterset[cell1].size += clusterset[cell2].size;
  clusterset[cell1].percoflag.setPerco(clusterset[cell2].percoflag.getPerco());
  clusterset.erase(cell2);
}

void Cell_disjoint_set::incrementcluster(Cell_handle cell1)
{
  ++clusterset[cell1].size;
}

void Cell_disjoint_set::newcluster(Cell_handle cell1)
{
  clusterset.insert(ClusterMapType::value_type(cell1, Cluster_info()));
}

// link two cell
uint32_t Cell_disjoint_set::merge(const Facet_decorate &fc)
{
  uint32_t mgflag = 0;
  int rp;
  double dtmp;
  Cell_handle cell1 = fc.cell1;
  Cell_handle cell2 = fc.cell2;
  // link cell2 to cell1
  CellMapType::iterator it1 = cellset.lower_bound(cell1);
  CellMapType::iterator it2 = cellset.lower_bound(cell2);
  CellMapType::iterator root0 = cellset.end();
  if(it1== cellset.end() || it1->first != cell1)
  {
    if(it2 == cellset.end() || it2->first != cell2)
    {
      // both cells have not been added yet
      it1 = cellset.insert(it1, CellMapType::value_type(cell1, Cell_info())); // add cell1
      it2 = cellset.insert(it2, CellMapType::value_type(cell2, Cell_info(cell2, cell1, *this))); // add cell2 to cell1
      root0 = it1;
      newcluster(cell1);
    }
    else
    {
      // only cell 1 not been added yet. merge cell 1 to cell 2
      it1 = cellset.insert(it1, CellMapType::value_type(cell1, Cell_info(cell1, cell2, *this)));
      root0 = findset(it1); // route compression
      incrementcluster(root0->first);
    }
  }
  else
  {
    if(it2 == cellset.end() || it2->first != cell2)
    {
      // only cell 2 not been added yet. merge cell 2 to cell 1
      it2 = cellset.insert(it2, CellMapType::value_type(cell2, Cell_info(cell2, cell1, *this)));
      root0 = findset(it2); // route compression
      incrementcluster(root0->first);
    }
    else
    {
      // merge two existing cells
      CellMapType::iterator root1 = findset(it1);
      CellMapType::iterator root2 = findset(it2);
      if(root1 != root2)
      {
        Vector dis1to2 = MyUtility::dif_points_PBC(celllist[it1->first], celllist[it2->first]);
        if(root1->second.rank > root2->second.rank)
        {
          // add root 2 to root 1
          root2->second.parent = root1->first;
          root2->second.vec_to_parent = it1->second.vec_to_parent - dis1to2 - it2->second.vec_to_parent;
          root0 = root1;
          mergecluster(root1->first, root2->first);
        }
        else
        {
          // add root 1 to root 2
          root1->second.parent = root2->first;
          root1->second.vec_to_parent = it2->second.vec_to_parent + dis1to2 - it1->second.vec_to_parent;
          root0 = root2;
          mergecluster(root2->first, root1->first);
        }
        if(root1->second.rank == root2->second.rank)
        {
          // equal. Already added root 1 to root 2. Now icrement rank
          ++(root2->second.rank);
        }
      }
      else
      {
        // check percolation. 1 to 2
        Vector disA = MyUtility::dif_points_PBC(celllist[it1->first], celllist[it2->first]);
        // the findset always compressed the route, so
        Vector disB = it1->second.vec_to_parent - it2->second.vec_to_parent;
        Vector diffVec = disA - disB;
        //std::cout << it1->second.vec_to_parent[0] << "\t" <<  it1->second.vec_to_parent[1] << "|" << it2->second.vec_to_parent[0] << "\t" <<  it2->second.vec_to_parent[1] << std::endl;
        for(rp=0; rp<DIM; ++rp)
        {
          dtmp = diffVec[rp]*diffVec[rp];
          if(dtmp > 0.01)
          {
            // std::cout << diffVec[0] << "\t" <<  diffVec[1] << "\t" << diffVec[2] << std::endl;
            mgflag |= (1 << rp);
          }
        }
        root0 = root1;
        clusterset[root0->first].percoflag.setPerco(mgflag);
        syspercoflag.setPerco(mgflag);
      }
    }
  }
  current_clu_size = clusterset[root0->first].size;
  if(max_clu_size<current_clu_size) max_clu_size=current_clu_size;
  return clusterset[root0->first].percoflag.getPerco();
}

// find the root of a disjoint set
Cell_disjoint_set::CellMapType::iterator Cell_disjoint_set::findset(CellMapType::iterator it)
{
  CellMapType::iterator itp = it; // iterator to parent
  std::vector<CellMapType::iterator> route; // the route of finding root
  Cell_handle par = itp->second.parent;
  if(par == Null_cell_handle) return itp;
  do
  {
    route.push_back(itp);
    itp = cellset.find(par);
    par = itp->second.parent;
  }while(par != Null_cell_handle);
  // path compression
  if(route.size() > 1)
  {
    Vector dis = route.back()->second.vec_to_parent;
    for(auto itit = std::next(route.rbegin()); itit != route.rend(); ++itit)
    {
      (*itit)->second.parent = itp->first;
      dis = dis + (*itit)->second.vec_to_parent;
      (*itit)->second.vec_to_parent = dis;
    }
  }
  return itp;
}


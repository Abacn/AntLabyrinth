//
//  traverse_facet.cpp
//  RLG_Delaunay_3D
//
//  Created by Yi Hu on 4/9/19.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include "graph_construct_3d.hpp"

using namespace std;

void Graph_Construct_3d::traverse_facet(Facetlist &facetlist)
{
  
if(input.outmode == OUTPUT_MODE_ALL)
{
  stringstream sst;
  sst << input.datafile << "_facet_" << count << ".dat";
  string fname;
  sst >> fname;
  ofstream of(fname);
  of << "DIM: " << DIM << "\n";
  
  for(Cell_iterator fit=tris.cells_begin(); fit!=tris.cells_end(); ++fit)
  {
    Cell_handle cell1 = fit;
    vector<Vector> cellpoints;
    cellpoints.push_back(Vector(0.0, 0.0, 0.0));
    Point p0 = cell1->vertex(0)->point();
    for(int rp=1; rp<=DIM; ++rp)
    {
      Point p1 = cell1->vertex(rp)->point();
      Vector vec = dif_points_PBC(p1, p0);
      cellpoints.push_back(vec);
    }
    vector<double> vectmp(DIM+1);
    for(int rp=0; rp<=DIM; ++rp)
    {
      Cell_handle cell2 = cell1->neighbor(rp);
      // Output all facet weights
      double weight = facet_weight(cellpoints, rp);
      vectmp[rp] = weight2rho(weight);
      if(cell1 < cell2)
      {
        Facet_decorate fd(weight, cell1, cell2);
        facetlist.push_back(fd);
      }
    }
    sort(vectmp.begin(), vectmp.end());
    for(int rp=0; rp<=DIM; ++rp)
    {
      of << vectmp[rp] << "\n";
    }
  }
  of.close();
  }
  else
  {
    for(Cell_iterator fit=tris.cells_begin(); fit!=tris.cells_end(); ++fit)
    {
      Cell_handle cell1 = fit;
      vector<Vector> cellpoints;
      cellpoints.push_back(Vector(0.0, 0.0, 0.0));
      Point p0 = cell1->vertex(0)->point();
      for(int rp=1; rp<=DIM; ++rp)
      {
        Point p1 = cell1->vertex(rp)->point();
        Vector vec = dif_points_PBC(p1, p0);
        cellpoints.push_back(vec);
      }
      for(int rp=0; rp<=DIM; ++rp)
      {
        Cell_handle cell2 = cell1->neighbor(rp);
        if(cell1 < cell2)
        {
          double weight = facet_weight(cellpoints, rp);
          Facet_decorate fd(weight, cell1, cell2);
          facetlist.push_back(fd);
        }
      }
    }
    }
}

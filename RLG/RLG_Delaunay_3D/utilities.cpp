//
//  utilities.cpp
//  RLG_Delaunay_3D
//
//  Created by Yi Hu on 10/6/18.
//

#include "graph_construct_3d.hpp"

Vector dif_points_PBC(const Point &point1, const Point &point2)
{
  // point1 - point2
  Vector vec = point1 - point2;
  return dif_vec_PBC(vec);
}

Vector dif_cells_PBC(const Cell_handle cell1, const Cell_handle cell2)
{
  // self - parent
  const Point point1 = cell1->vertex(0)->point();
  const Point point2 = cell2->vertex(0)->point();
  return dif_points_PBC(point1, point2);
}

Vector dif_vec_PBC(Vector &vec)
{
  double dtmp[DIM];
  for(int rp=0; rp<DIM; ++rp)
  {
    dtmp[rp] = vec[rp];
    if(dtmp[rp] > Parameter::HalfBoxL) dtmp[rp] -= Parameter::BoxLength;
    if(dtmp[rp] < -Parameter::HalfBoxL) dtmp[rp] += Parameter::BoxLength;
  }
  return Vector(dtmp[0], dtmp[1], dtmp[2]);
}

//
//  graph_construct_3d.cpp
//  RLG_Delaunay_3D
//  Calculate the void percolation threshold of RLG
//  in 3D periodic box

#include <algorithm>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <CGAL/Random.h>
#include <Eigen/Core>
#include <Eigen/LU>

#include "graph_construct_3d.hpp"

// Convienient Macro
#define ttds tris.tds()  // get triangulation result

using namespace std;
// constructor
Graph_Construct_3d::Graph_Construct_3d(ReadInputRLG input): input(input)
{
  // overwrite the default random number generator by specified random seed
  CGAL::get_default_random() = CGAL::Random(input.seed);
  Vunitsphere = pow(M_PI, DIM*0.5)/tgamma(1.0+DIM*0.5);
  count = 0;
}

// deconstructor
Graph_Construct_3d::~Graph_Construct_3d()
{
  
}

// master working entry
int Graph_Construct_3d::graphrun()
{
  CGAL::Random_points_in_cube_3<Point, Creator> rand_it(0.5*Parameter::BoxLength); // point generator (-0.5, 0.5)
  if(!points.empty())
  {
    points.clear();
    tris.clear();
  }
#ifdef DEBUG
  cout << "  Insert points..." << endl;
#endif
  Vector halfvec(0.5, 0.5, 0.5);
  for (int rp=0 ; rp < input.N ; ++rp)
  {
    points.push_back(*rand_it + halfvec);
    ++rand_it;
  }
  tris.insert(points.begin(), points.end(), true); // insert into tris
#ifdef DEBUG
  cout << "  Traverse facets..." << endl;
#endif
  Facetlist facetlist;
  traverse_facet(facetlist);
  sort(facetlist.rbegin(), facetlist.rend()); // descend sorting
  // Perco_flag pcstatus = {.perco_all = 0};
  Cell_disjoint_set cdset;
  // merge graph
#ifdef DEBUG
  // cout << tris.number_of_vertices() << " " << tris.number_of_edges() << " " << tris.number_of_facets() << " "
  // << tris.number_of_cells() << ";" << " " << facetlist.size() << endl;
  cout << "  Linking cells..." << endl;
#endif
  vector<Facet_decorate>::iterator fit;
  long countfacet = 0;
  std::uint32_t percoflag;
  const int n_cri = 5;
  double weights[n_cri] = {0.0};
  bool pcflag[n_cri] = {false};
  for(fit = facetlist.begin(); !pcflag[n_cri-1] && fit != facetlist.end(); ++fit)
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
  //       double perco_rho = input.N / pow(1/perco_radius, DIM);
  double perco_phi;
  if(fit == facetlist.end())
  {
    perco_phi = NAN;
  }
  else
  { // tracer radius is 1
    perco_phi = weight2Phi(weights[n_cri-1]);
  }
  result_record.set_value(perco_phi, facetlist.size(), countfacet,
                          tris.number_of_cells(), cdset.get_cell_size() );
  ++count;
  graphdump();
  if(input.outmode >= OUTPUT_MODE_PRODUCT)
  {
    std::stringstream dbgoutstr;
    for(int rp=0; rp<n_cri; ++rp)
    {
      // percolation threshold in different critera
      dbgoutstr << weight2Phi(weights[rp]) << "\t";
    }
    // percolation cluster size, number of clusters
    dbgoutstr << cdset.getmax_clu_size() << "\t" << cdset.getclustercount();
    // maximum cell square circumradius
    graphdebug(dbgoutstr.str().c_str());
  }
  return 0;
}

// dump results
int Graph_Construct_3d::graphdump()
{
  if(1 == count)
  {
    fout.open(input.datafile + ".dat", std::ofstream::out);
    fout << "Dimension: " << DIM << "\n";
    fout << "N:         " << input.N << "\n";
    fout << "Phi\tNfacet\tNvfacet\tNcell\tNvcell" << endl;
  }
  result_record.print_value(fout);
  return 0;
}

int Graph_Construct_3d::graphdebug(const char *line)
{
  if(1 == count)
  {
    fdbgsummary.open("dbgsummary.dat");
    fdbgsummary << "Phi1\tPhi2\tPhi3\tPhi4\tPhi5";
    fdbgsummary << "\tpcclustersize\tNCluster";
    fdbgsummary << std::endl;
  }
  fdbgsummary << line << std::endl;
  return 0;
}

// return result class itself
const Result_record Graph_Construct_3d::get_results()
{
  return result_record;
}

// protected member methods

// Facet weight is the square of circumradius
double Graph_Construct_3d::facet_weight(const vector<Vector> &offpoints, int covertex) const
{
  Eigen::Matrix<double, DIM+1, DIM+1> matA;  // augmented matrix
  int rp, rq, t_rp, t_rq;
  double weight, dtmp;
  matA(0, 0) = 0.0;
  for(t_rp=1, rp=0; rp<=DIM; ++rp)
  {
    if(rp == covertex) continue;
    matA(t_rp, t_rp) = 0.0;
    matA(t_rp, 0) = matA(0, t_rp) = 1.0;
    for(rq=rp+1, t_rq=t_rp+1; rq<=DIM; ++rq)
    {
      if(rq == covertex) continue;
      dtmp = norm_squared_distance(offpoints[rp], offpoints[rq]);
      matA(t_rp, t_rq) = matA(t_rq, t_rp) = dtmp;
      ++t_rq;
    }
    ++t_rp;
  }
  auto matinv = matA.inverse();
  weight = -0.5*matinv(0, 0);
  return weight;
}

double Graph_Construct_3d::norm_squared_distance(const Vector &pa, const Vector &pb) const
{
  double result = 0.0, dtmp;
  for(int rp=0; rp<DIM; ++rp)
  {
    dtmp = pa[rp] - pb[rp];
    result += dtmp*dtmp;
  }
  return result;
}

double Graph_Construct_3d::weight2rho(double weight) const
{
  return input.N * pow(weight, 0.5*DIM);
}

// transform weight to Phi (=number density*Unit sphere volume)
double Graph_Construct_3d::weight2Phi(double weight) const
{
  return input.N * pow(weight, 0.5*DIM) * Vunitsphere;
}

// Test functions
int Graph_Construct_3d::graphtest()
{
  
  return 0;
}

#undef ttds

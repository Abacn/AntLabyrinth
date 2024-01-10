//
//  cellsample.cpp
//  RLG_Cageshape
//
//  Created by Yi Hu on 4/29/19.
//

#include <algorithm>
#include <random>
#include <iostream>
#include <cmath>
#include <vector>

#include <Eigen/Core>
#include <Eigen/LU>

#include "cellsample.hpp"
#include "dim.h"

namespace SampleUtility{
  std::mt19937 rg;
  // ----- namespace private -----
  std::mt19937 rgpoisson;  // random generator.
  // Use different generator for myrand() and poissonrand() for reproducibility
  std::uniform_real_distribution <> rrand; // uniform distribution in (0, 1]
  std::poisson_distribution<> pois;
  // -----
  // random number [0, 1)
  double myrand()
  {
    return rrand(rg);
  }
  
  // set random seed
  void setseed(unsigned int seed)
  {
    rg.seed(seed);
    rgpoisson.seed(seed);
  }
  
  // set poisson distribution parameter
  void setpoisson(double param)
  {
    pois = std::poisson_distribution<>(param);
  }
  
  int poissonrand()
  {
    return pois(rgpoisson);
  }
  
  // square distance of pa to the origin
  double norm_squared_distance(const Vector &pa)
  {
    double result = 0.0;
    for(int rp=0; rp<DIM; ++rp)
    {
      result += pa[rp]*pa[rp];
    }
    return result;
  }
  
  // square distance between pa
  double norm_squared_distance(const Vector &pa, const Vector &pb)
  {
    double result = 0.0, dtmp;
    for(int rp=0; rp<DIM; ++rp)
    {
      dtmp = pa[rp] - pb[rp];
      result += dtmp*dtmp;
    }
    return result;
  }
  
  // dot product
  double dot(const Vector &pa, const Vector &pb)
  {
    double result = 0.;
    for(int rp=0; rp<DIM; ++rp) result += pa[rp]*pb[rp];
    return result;
  }
  
  // translate pb so that pa is at the origin
  void translate(const Vector &pa, Vector &pb)
  {
    for(int rp=0; rp<DIM; ++rp) pb[rp] = pb[rp]-pa[rp];
  }
  
  // projection of point to hyperplane
  Vector projection(std::vector<Vector> sp, Vector p)
  {
    int rp, rq;
    int dim = sp.size()-1;
    Vector p2({0});
    if(0 == dim)
    {
      return sp.front();
    }
    else
    {
      // translate so that S0 is the origin
      for(rp=1; rp<=dim; ++rp) translate(sp.front(), sp[rp]);
      translate(sp.front(), p);
      // for(rp=0; rp<dim; ++rp) sp[0][rp] = 0.; // translate S0 not needed
      // find projection p2 of p to the linear subspace with the set of basis vectors sp
      Eigen::MatrixXd matA(dim, dim);
      Eigen::MatrixXd vecb(dim, 1);
      for(rp=1; rp<=dim; ++rp)
      {
        matA(rp-1, rp-1) = dot(sp[rp], sp[rp]);
        for(rq=rp+1; rq<=dim; ++rq)
        {
          matA(rq-1, rp-1) = matA(rp-1, rq-1) = dot(sp[rp], sp[rq]);
        }
        vecb(rp-1) = dot(p, sp[rp]);
      }
      Eigen::MatrixXd veca = matA.lu().solve(vecb);
      /*
       for(rp=0; rp<=dim; ++rp)
       std::cout << sp[rp][0] << " " << sp[rp][1] << "\n";
       std::cout << p[0] << " " << p[1] << "\n\n";
       */
      // get p2
      for(rp=0; rp<dim; ++rp)
      {
        for(rq=0; rq<DIM; ++rq)
        {
          p2[rq] += veca(rp)*sp[rp+1][rq];
        }
      }
      for(rp=0; rp<DIM; ++rp)
      {
        p2[rp] += sp[0][rp]; // translate back
      }
    }
    return p2;
  }
  
  SampleChoice::SampleChoice(std::vector<double> weights): cumuw(weights.size())
  {
    if(0 == weights.size()) return;
    double totalw = 0.;
    for(int rp=0; rp<weights.size(); ++rp)
    {
      totalw += weights[rp];
      cumuw[rp] = totalw;
    }
  }
  
  std::size_t SampleChoice::getSample()
  {
    // binary search
    double dtmp = myrand()*cumuw.back();
    std::size_t l = 0, r = cumuw.size()-1, mid;
    while(l<r)
    {
      mid = (l + r)/2;
      if(dtmp<cumuw[mid])
        r = mid;
      else
        l = mid+1;
    }
    return l;
  }
  
  double SampleChoice::getw()
  {
    return cumuw.back();
  }
}

using namespace SampleUtility;

// local help functions
using Matdim = typename Eigen::Matrix<double, DIM, DIM>;        // matrix of dimension dim
using Matdimp1 = typename Eigen::Matrix<double, DIM+1, DIM+1>;  // matrix of dimension dim+1
using Matdimp2 = typename Eigen::Matrix<double, DIM+2, DIM+2>;  // matrix of dimension dim+2
using Vecdim = typename Eigen::Matrix<double, DIM, 1>;          // vector of dimension dim
using Vecdimp2 = typename Eigen::Matrix<double, DIM+2, 1>;      // vector of dimension dim+2

// calculate sign of cell determinant
char calc_celld(const Matdim matA)
{
  double det = matA.determinant();
  return sgn(det);
}

// calculate both cell volume and store the sign of the determinant
char calc_celld(const Matdim matA, double &volume)
{
  double det = matA.determinant();
  volume = fabs(det)/tgamma(DIM+1);
  return sgn(det);
}

// Determine if the point A and the cofactor(i) are in the same side of hyperplane vecs\vecs[i]
char calc_halfplane(const Vector &vecA, int cofc, const Simplex_array<Vector> &vecs, double &subvolume)
{
  Matdim matA, matB;
  char sign1, sign2;
  int t_rp=0;
  for(int rp=0; rp<=DIM; ++rp)
  {
    if(rp!=cofc)
    {
      for(int rq=0; rq<DIM; ++rq)
      {
        matA(t_rp, rq) = vecs[rp][rq] - vecA[rq];
        matB(t_rp, rq) = vecs[rp][rq] - vecs[cofc][rq];
      }
      ++t_rp;
    }
  }
  sign1 = calc_celld(matA, subvolume);
  sign2 = calc_celld(matB);
  return sign1*sign2;
}

// simplex to vector square distance. would change the input
// Oleg Golubitsky, Vadim Mazalov and Stephen M. Watt, An Algorithm to Compute the Distance from a Point to a Simplex
double dis_simplex_to_point(std::vector<Vector> &sp, Vector &p)
{
  int rp, rq;
  double dtmp, result = 0.;
  int dim = sp.size()-1;
  if(0 == dim)
  {
    result = norm_squared_distance(sp.front(), p);
  }
  else
  {
    // translate so that S0 is the origin
    for(rp=1; rp<=dim; ++rp) translate(sp.front(), sp[rp]);
    translate(sp.front(), p);
    // for(rp=0; rp<dim; ++rp) sp[0][rp] = 0.; // translate S0 not needed
    // find projection p2 of p to the linear subspace with the set of basis vectors sp
    Eigen::MatrixXd matA(dim, dim);
    Eigen::MatrixXd vecb(dim, 1);
    for(rp=1; rp<=dim; ++rp)
    {
      matA(rp-1, rp-1) = dot(sp[rp], sp[rp]);
      for(rq=rp+1; rq<=dim; ++rq)
      {
        matA(rq-1, rp-1) = matA(rp-1, rq-1) = dot(sp[rp], sp[rq]);
      }
      vecb(rp-1) = dot(p, sp[rp]);
    }
    Eigen::MatrixXd veca = matA.lu().solve(vecb);
    /*
    for(rp=0; rp<=dim; ++rp)
    std::cout << sp[rp][0] << " " << sp[rp][1] << "\n";
    std::cout << p[0] << " " << p[1] << "\n\n";
    */
    // get p2
    Vector p2({0});
    for(rp=0; rp<dim; ++rp)
    {
      for(rq=0; rq<DIM; ++rq)
      {
        p2[rq] += veca(rp)*sp[rp+1][rq];
      }
    }
    // check if projection on the cell
    dtmp = 0.;
    for(rp=0; rp<dim; ++rp)
    {
      if(veca(rp)<0.) break;
      dtmp += veca(rp);
    }
    if(rp==dim)
    {
      if(dtmp<=1)
      {
        // case 1: projection is inside
        result = norm_squared_distance(p, p2);
      }
      else
      {
        // case 3: Distance(S\S0, P)
        std::vector<Vector> sq(std::next(sp.begin()), sp.end());
        result = dis_simplex_to_point(sq, p);
      }
    }
    else
    {
      // case 2: projection is outside
      std::vector<Vector> sq(1, Vector({0.}));
      for(rq=0; rq<dim; ++rq)
      {
        if(veca(rq)>0) sq.push_back(sp[rq+1]);
      }
      result = dis_simplex_to_point(sq, p);
    }
  }
  return result;
}

// calculate cell volume
double calc_cellvolume(const Simplex_array<Vector> &vecs)
{
  Matdim matA;
  for(int rp=0; rp<DIM; ++rp)
  {
    for(int rq=0; rq<DIM; ++rq)
    {
      matA(rp, rq) = vecs[rp+1][rq] - vecs[0][rq];
    }
  }
  return fabs(matA.determinant())/tgamma(DIM+1);
}

// calculate the circumcenter coordinate
Vector calc_cellcenter(const Simplex_array<Vector> &vecs, double &r)
{
  int rp, rq;
  Matdimp2 matA;
  Vecdimp2 vecb = Vecdimp2::Zero();
  double dtmp;
  matA(0, 0) = 0.;
  vecb(0) = 1.;
  for(int rp=0; rp<=DIM; ++rp)
  {
    matA(0, rp+1) = matA(rp+1, 0) = 1.;
    matA(rp+1, rp+1) = 0.;
    for(int rq=rp+1; rq<=DIM; ++rq)
    {
      dtmp = norm_squared_distance(vecs[rp], vecs[rq]);
      matA(rp+1, rq+1) = matA(rq+1, rp+1) = dtmp;
    }
  }
  Vecdimp2 vecx = matA.lu().solve(vecb);
  Vector cellcen({0.});
  for(rp=0; rp<=DIM; ++rp)
  {
    for(rq=0; rq<DIM; ++rq)
    {
      cellcen[rq] += vecx(rp+1)*vecs[rp][rq];
    }
  }
  r = -0.5*vecx(0); // r is square circumradius
  return cellcen;
}

Cell_single_sample::Cell_single_sample(double weight_, double r_, const Vector &vec): weight(weight_), r(r_)
{
  for(int rp=0; rp<DIM; ++rp)
  {
    coord[rp] = vec[rp];
  }
}

Vector Cell_single_sample::getvec()
{
  Vector vec({0.});
  for(int rp=0; rp<DIM; ++rp)
  {
    vec[rp] = coord[rp];
  }
  return vec;
}

// weighted square distance between two samples, return weight in ret_weighr
double Cell_single_sample::norm_square_distance(const Cell_single_sample &s2, double &ret_weight) const
{
  double result = 0., dtmp;
  for(int rp=0; rp<DIM; ++rp)
  {
    dtmp = coord[rp] - s2.coord[rp];
    result += dtmp*dtmp;
  }
  ret_weight = weight*s2.weight;
  return result;
}

// point to simplex distance
double Cell_samples::dis_to_point(const Vector p) const
{
  std::vector<Vector> simplex(vecs.begin(), vecs.end());
  Vector point(p);
  return dis_simplex_to_point(simplex, point);
}

// default cell sample
Cell_samples::Cell_samples(): volume(0) {}

// initialize a cell. Do not run sampling yet
Cell_samples::Cell_samples(const Simplex_array<Vector> &vertices): vecs(vertices)
{
  int rp;
  volume = calc_cellvolume(vecs);
  center = calc_cellcenter(vecs, r);
  // building subcells
  for(rp=0; rp<=DIM; ++rp)
  {
    subcellsign[rp] = calc_halfplane(center, rp, vecs, subvolume[rp]);
  }
}

// generate random numbers in Barycentric system of a simplex
void Cell_samples::genRandom(double arr[DIM+1]) const
{
  int rp;
  arr[0] = 0.;
  for(rp=1; rp<=DIM; ++rp)
  {
    arr[rp] = myrand();
  }
  std::sort(arr+1, arr+DIM+1);
  for(rp=0; rp<DIM; ++rp)
  {
    arr[rp] = arr[rp+1] - arr[rp];
  }
  arr[DIM] = 1.-arr[rp];
}

// sampling in a simplex. For space optimization, only the samples with distance to nearest vertex
// more than storethd is stored

void Cell_samples::resample(const long *samplestrat)
{
  std::vector<Vector> nulloutpoints;
  resample(samplestrat, nulloutpoints, true);
}

void Cell_samples::resample(const long *samplestrat, const std::vector<Vector> &outpoints, bool vertex_flag/*=true*/)
{
  double arr[DIM+1];
  long rp, rq, rr;
  bool has_outpoints = !outpoints.empty();
  clearsample();
  // divide subcell samples
  long dupnum = 0, maxdupnum = samplestrat[1]/samplestrat[0], voidsample = 0;
  Vector pos({0.});
  while(voidsample<samplestrat[2] && dupnum<maxdupnum)
  {
    for(rq = 0; rq<samplestrat[0]; ++rq) // each subcell samplepercell samples
    {
      genRandom(arr);
      for(rr=0; rr<DIM; ++rr) pos[rr] = 0.;
      for(rp=0; rp<=DIM; ++rp)
      {
        for(rr=0; rr<DIM; ++rr)
        {
          pos[rr] += vecs[rp][rr]*arr[rp];
        }
      }
      double minr = 1., nowr;
      // vertices
      if(vertex_flag)
      {
        for(rr=0; rr<=DIM; ++rr)
        {
          nowr = norm_squared_distance(pos, vecs[rr]);
          if(nowr<minr)
          {
            minr=nowr;
            if(minr<ONE) goto loop_end;
          }
        }
      }
      // oupoints
      if(has_outpoints)
      {
        for(const auto &outpoint: outpoints)
        {
          nowr = norm_squared_distance(pos, outpoint);
          if(nowr<minr)
          {
            minr=nowr;
            if(minr<ONE) goto loop_end;
          }
        }
      }
      samples.push_back(Cell_single_sample(volume/samplestrat[0], minr, pos)); // weight set after
      ++voidsample;
      loop_end:;
    }
    ++dupnum;
  }
  double reweight = 1.0/dupnum;
  for(Cell_single_sample &rf: samples)
  {
    rf.weight *= reweight;
  }
}

bool Cell_samples::onesample(const std::vector<Vector> &outpoints, bool vertex_flag/*=true*/)
{
  double arr[DIM+1];
  Vector pos({0.});
  int rp, rr;
  bool has_outpoints = !outpoints.empty();
  genRandom(arr);
  for(rr=0; rr<DIM; ++rr) pos[rr] = 0.;
  for(rp=0; rp<=DIM; ++rp)
  {
    for(rr=0; rr<DIM; ++rr)
    {
      pos[rr] += vecs[rp][rr]*arr[rp];
    }
  }
  double minr = 1., nowr;
  // vertices
  if(vertex_flag)
  {
    for(rr=0; rr<=DIM; ++rr)
    {
      nowr = norm_squared_distance(pos, vecs[rr]);
      if(nowr<minr)
      {
        minr=nowr;
        if(minr<ONE) return false;
      }
    }
  }
  // oupoints
  if(has_outpoints)
  {
    for(const auto &outpoint: outpoints)
    {
      nowr = norm_squared_distance(pos, outpoint);
      if(nowr<minr)
      {
        minr=nowr;
        if(minr<ONE) return false;
      }
    }
  }
  samples.push_back(Cell_single_sample(1., minr, pos)); // weight set after
  return true;
}

void Cell_samples::clearsample()
{
  samples.clear();
}

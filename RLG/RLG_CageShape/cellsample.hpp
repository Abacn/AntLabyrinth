//
//  cellsample.hpp
//  RLG_Cageshape
//
//  Created by Yi Hu on 4/29/19.
//

#ifndef cellsample_hpp
#define cellsample_hpp

#include <array>
#include <random>
#include <valarray>
#include <vector>
#include "dim.h"
using Vector = std::array<double, DIM>;
template<typename T> using Simplex_array = typename std::array<T, DIM+1>;


namespace SampleUtility{
  // get sign
  template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
  }
  extern std::mt19937 rg;
  void setseed(unsigned int seed);
  void setpoisson(double param);
  double myrand();
  int poissonrand();
  double norm_squared_distance(const Vector &pa, const Vector &pb);
  double norm_squared_distance(const Vector &pa);
  double dot(const Vector &pa, const Vector &pb);
  void translate(const Vector &pa, Vector &pb);
  Vector projection(std::vector<Vector> sp, Vector p);
  class SampleChoice{ // choose one sample
  public:
    SampleChoice(std::vector<double> weights);
    std::size_t getSample();
    double getw();
  private:
    std::vector<double> cumuw; // cumulative weights
  };
}
// a simple sample
struct Cell_single_sample{
  // Now use double
  double coord[DIM];   // coordinate to the cell
  double weight;   // the weight (dV) of this sample
  double r;        // square distance to the closest vertex
  Cell_single_sample(double weight_, double r_, const Vector &vec);
  Cell_single_sample(){};
  Vector getvec();
  double norm_square_distance(const Cell_single_sample &s2, double &ret_weight) const;
};

// sample the volume in a cell
struct Cell_samples{
  double volume, r;
  Vector center;
  Simplex_array<Vector> vecs;
  Simplex_array<double> subvolume;
  Simplex_array<char> subcellsign;
  std::vector<Cell_single_sample> samples;
  Cell_samples();
  Cell_samples(const Simplex_array<Vector> &vertices);
  void resample(const long *samplestrat);    // redo the sampling of N samples
  void resample(const long *samplestrat, const std::vector<Vector> &outpoints, bool vertex_flag=true);
  bool onesample(const std::vector<Vector> &outpoints, bool vertex_flag=true);
  void clearsample();
  // generate random numbers in Barycentric system
  void genRandom(double arr[DIM+1]) const;
  double dis_to_point(const Vector p) const;
  const double ONE = 1.;
};

#endif /* cellsample_hpp */

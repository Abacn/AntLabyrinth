//
//  cagestat.hpp
//  RLG_Cageshape
//
//  Created by Yi Hu on 4/29/19.
//

#ifndef cagestat_hpp
#define cagestat_hpp

#include "cellsample.hpp"

using SamplesType = typename std::vector<Cell_single_sample>;
// interface of histogram statistics
struct Histbins{
  double binstart, binend, logbinstart, logbinend;
  long ndata, logndata;
  double dr, invdr, logdr, loginvdr;
  std::vector<double> data, logdata;
  double dataleft, logdataleft;
  double denom;
  // virtual process(...)
  virtual void dump(const std::string &filename);
  Histbins(double binstart, double binend, int binperl, double logbinstart, double logbinend, int logbinperl);
  virtual ~Histbins() = default;
};

// self van hove function statistics
struct SelfVanHove: virtual Histbins{
  SelfVanHove(double binstart, double binend, int binperl, double logbinstart, double logbinend, int logbinperl);
  void process(const SamplesType &samples, double cluv, double sqweight);
  void dump(const std::string &filename) override;
};

// correlation function statistics
struct CorrFunc: virtual Histbins{
  CorrFunc(double binstart, double binend, int binperl, double logbinstart, double logbinend, int logbinperl);
  void process(const SamplesType &samples, double cluv, double sqweight);
  void dump(const std::string &filename) override;
};

// mass distribution around the furthest point
struct MassDist: virtual Histbins{
  MassDist(double binstart, double binend, int binperl, double logbinstart, double logbinend, int logbinperl);
  void process(const SamplesType &samples, const Vector &center);
  void dump(const std::string &filename) override;
};

// Mean Displacement distribution of tracers. Note that P(D^2) = P(D)/(2*D)
struct MSDDist: virtual Histbins{
  MSDDist(double binstart, double binend, int binperl, double logbinstart, double logbinend, int logbinperl);
  void process(const SamplesType &samples, double cluv, double sqweight);
  void dump(const std::string &filename) override;
};

class AllStat{
  Histbins* para;
  SelfVanHove* selfvhf;
  CorrFunc* corrfunc;
  MassDist* massdist;
  MSDDist* msddist;
  double chidyn; // variance of sd
public:
  AllStat(bool dostat[4], double binstart, double binend, int binperl, double logbinstart, double logbinend, int logbinperl);
  ~AllStat();
  int processAll(const SamplesType &samples, double cluv, double sqweight, const Vector &center);
  void dumpAll(const std::string &filename);
  inline double getchidyn(){return chidyn;}
};
#endif /* cagestat_hpp */

//
//  cagestat.cpp
//  RLG_Cageshape
//
//  Created by Yi Hu on 4/29/19.
//
#include <cmath>
#include <iostream>
#include <fstream>

#include "cagestat.hpp"
#include "cellsample.hpp"

// Base class
Histbins::Histbins(double binstart_, double binend_, int binperl, double logbinstart_, double logbinend_, int logbinperl):
binstart(binstart_), binend(binend_), logbinstart(log(logbinstart_)), logbinend(log(logbinend_))
{
  // linear hist
  dr = 1./binperl;
  invdr = binperl;
  ndata = ceil((binend-binstart)*invdr);
  data = std::vector<double>(ndata+1, 0.); // the last entry store all >= cutoff
  dataleft = 0.;
  // loghist
  logdr = log(10.)/logbinperl;
  loginvdr = 1./logdr;
  logndata = ceil((logbinend-logbinstart)*loginvdr);
  logdata = std::vector<double>(logndata+1, 0.); // the last entry store all >= cutoff
  logdataleft = 0.;
  // normalization factor
  denom = 0.;
}

void Histbins::dump(const std::string &filename)
{
  // linear
  std::ofstream ofs(filename+"_lin.dat");
  ofs << denom << "\n"; // the first line is the weight of this histogram
  double lftr=binstart+dr/2;
  for(unsigned long rp=0; rp<=ndata; ++rp)
  {
    double rval = lftr+dr*rp;
    double pval = data[rp]/(denom*dr);
    ofs << rval << "\t" << pval << "\n";
  }
  ofs.close();
  // log
  ofs.open(filename+"_log.dat");
  ofs << denom << "\n"; // the first line is the weight of this histogram
  double logoffset = logbinstart + logdr/2;
  for(unsigned long rp=0; rp<=logndata; ++rp)
  {
    double rval = logoffset+logdr*rp;
    double pval = logdata[rp]/(denom*(exp(rval+logdr/2) - exp(rval-logdr/2)));
    ofs << exp(rval) << "\t" << pval << "\n";
  }
}

// Self Van Hove function
SelfVanHove::SelfVanHove(double binstart, double binend, int binperl, double logbinstart, double logbinend, int logbinperl):
Histbins(binstart, binend+(binend-binstart)*0.5, binperl, logbinstart, logbinend+(logbinend-logbinstart)*0.5, logbinperl)
{
  
}

void SelfVanHove::process(const SamplesType &samples, double cluv, double sqweight)
{
  std::cerr << "Warning: function not implemented yet" << std::endl;
  // TODO
}

void SelfVanHove::dump(const std::string &filename)
{
  Histbins::dump(filename);
}

// Correlation function
CorrFunc::CorrFunc(double binstart, double binend, int binperl, double logbinstart, double logbinend, int logbinperl):
Histbins(binstart, binend, binperl, logbinstart, logbinend, logbinperl)
{
  
}

void CorrFunc::process(const SamplesType &samples, double cluv, double sqweight)
{
  std::cerr << "Warning: function not implemented yet" << std::endl;
  // TODO
}

void CorrFunc::dump(const std::string &filename)
{
  double Vunitsphere = pow(M_PI, DIM*0.5)/tgamma(1.0+DIM*0.5);
  // linear
  std::ofstream ofs(filename+"_lin.dat");
  ofs << denom << "\n"; // the first line is the weight of this histogram
  double lftr=binstart+dr/2;
  for(unsigned long rp=0; rp<=ndata; ++rp)
  {
    double rval = lftr+dr*rp;
    double pval = data[rp]/(denom*Vunitsphere*(pow(rval+dr/2, DIM) - pow(rval-dr/2, DIM)));
    ofs << rval << "\t" << pval << "\n";
  }
  ofs.close();
  // log
  ofs.open(filename+"_log.dat");
  ofs << denom << "\n"; // the first line is the weight of this histogram
  double logoffset = logbinstart + logdr/2;
  for(unsigned long rp=0; rp<=logndata; ++rp)
  {
    double logrval = logoffset+logdr*rp;
    double rval = exp(logrval);
    double pval = logdata[rp]/(denom*Vunitsphere*(exp((logrval+logdr/2)*DIM) - exp((logrval-logdr/2)*DIM)));
    ofs << rval << "\t" << pval << "\n";
  }
}

// Mass distribution with respect to the innerst site
MassDist::MassDist(double binstart, double binend, int binperl, double logbinstart, double logbinend, int logbinperl):
Histbins(binstart, binend, binperl, logbinstart, logbinend, logbinperl)
{
  
}

void MassDist::process(const SamplesType &samples, const Vector &center)
{
  std::cerr << "Warning: function not implemented yet" << std::endl;
  // TODO
}

void MassDist::dump(const std::string &filename)
{
  double Vunitsphere = pow(M_PI, DIM*0.5)/tgamma(1.0+DIM*0.5);
  // linear
  std::ofstream ofs(filename+"_lin.dat");
  ofs << denom << "\n"; // the first line is the weight of this histogram
  double lftr=binstart+dr/2;
  for(unsigned long rp=0; rp<=ndata; ++rp)
  {
    double rval = lftr+dr*rp;
    double pval = data[rp]/(denom*Vunitsphere*(pow(rval+dr/2, DIM) - pow(rval-dr/2, DIM)));
    ofs << rval << "\t" << pval << "\n";
  }
  ofs.close();
  // log
  ofs.open(filename+"_log.dat");
  ofs << denom << "\n"; // the first line is the weight of this histogram
  double logoffset = logbinstart + logdr/2;
  for(unsigned long rp=0; rp<=logndata; ++rp)
  {
    double logrval = logoffset+logdr*rp;
    double rval = exp(logrval);
    double pval = logdata[rp]/(denom*Vunitsphere*(exp((logrval+logdr/2)*DIM) - exp((logrval-logdr/2)*DIM)));
    ofs << rval << "\t" << pval << "\n";
  }
}

// Displacement distribution of tracers
MSDDist::MSDDist(double binstart, double binend, int binperl, double logbinstart, double logbinend, int logbinperl):
Histbins(binstart, binend, binperl, logbinstart, logbinend, logbinperl)
{
  
}

void MSDDist::process(const SamplesType &samples, double cluv, double sqweight)
{
  std::cerr << "Warning: function not implemented yet" << std::endl;
  // TODO
}

void MSDDist::dump(const std::string &filename)
{
  Histbins::dump(filename);
}

AllStat::AllStat(bool dostat[4], double binstart, double binend, int binperl, double logbinstart, double logbinend, int logbinperl): para(nullptr), selfvhf(nullptr), corrfunc(nullptr), massdist(nullptr)
{
  for(int rp=0; rp<4; ++rp)
  {
    if(dostat[rp])
    {
      switch(rp)
      {
        case 0: para = selfvhf = new SelfVanHove(binstart, binend, binperl, logbinstart, logbinend, logbinperl); break;
        case 1: para = corrfunc = new CorrFunc(binstart, binend, binperl, logbinstart, logbinend, logbinperl); break;
        case 2: para = massdist = new MassDist(binstart, binend, binperl, logbinstart, logbinend, logbinperl); break;
        case 3: para = msddist = new MSDDist(binstart, binend, binperl, logbinstart, logbinend, logbinperl);
            break;
      }
    }
  }
}

AllStat::~AllStat()
{
  if(selfvhf!=nullptr) delete selfvhf;
  if(corrfunc!=nullptr) delete corrfunc;
  if(massdist!=nullptr) delete massdist;
  if(msddist!=nullptr) delete msddist;
}

// process for all statistics. Use para as the parameter for histogram. All histogram should be in same division
int AllStat::processAll(const SamplesType &samples, double cluv, double sqweight, const Vector &center)
{
  int nth, lognth, rp, rq, rr;
  unsigned long spsize = samples.size();
  // check hist consistent
  double weight, dtmp, rij, rijsq;
  // normalization factor of self van hove function
  double scf = 1./(cluv*cluv-sqweight);
  double msd_ = 0., qsd_ = 0.;
  chidyn = 0.;
  if(selfvhf != nullptr)
  {
    selfvhf->denom += 1.;
  }
  // normalization factor of correlation length
  double scfi = 1./(cluv-sqweight/cluv);
  if(corrfunc != nullptr)
  {
    corrfunc->denom += 1.;
  }
  // normalization factor of mass distribution
  double scfm = 1.;
  if(massdist != nullptr)
  {
    massdist->denom += 1.;
  }
  // normalization for MSD distribution
  double scfms = 1./cluv;
  if(msddist != nullptr)
  {
    msddist->denom += 1.;
  }
  // check sample size, special case
  if(spsize<=1)
  {
    if(spsize==1)
    {
      // single sample in void, consider as a zero volume sample
      nth = -para->binstart*para->invdr;
      weight = samples[0].weight;
      if(selfvhf != nullptr)
      {
        dtmp = 1.;
        if(nth<0) selfvhf->dataleft += dtmp;
        else selfvhf->data[0] += dtmp;
        selfvhf->logdataleft += dtmp;
      }
      if(corrfunc != nullptr)
      {
        dtmp = weight;
        if(nth<0) corrfunc->dataleft += dtmp;
        else corrfunc->data[0] += dtmp;
        corrfunc->logdataleft += dtmp;
      }
      if(massdist != nullptr)
      {
        dtmp = weight;
        if(nth<0) massdist->dataleft += dtmp;
        else massdist->data[0] += dtmp;
        massdist->logdataleft += dtmp;
      }
      if(msddist != nullptr)
      {
        dtmp = 1.;
        if(nth<0) msddist->dataleft += dtmp;
        else msddist->data[0] += dtmp;
        msddist->logdataleft += dtmp;
      }
    }
    return 1;
  }
  for(rp=0; rp<spsize; ++rp)
  {
    const auto &s1 = samples[rp];
    double msd = 0., qsd = 0.;
    // O(n^2) part, self van hove and correlation
    for(rq=0; rq<spsize; ++rq)
    {
      if(rp==rq) continue;
      const auto &s2 = samples[rq];
      // lin hist add data
      rijsq = s1.norm_square_distance(s2, weight);
      msd += rijsq*s2.weight;
      qsd += rijsq*rijsq*s2.weight;
      rij = sqrt(rijsq);
      nth = floor((rij-para->binstart)*para->invdr);
      lognth = floor((log(rij)-para->logbinstart)*para->loginvdr);
      // self van hove
      if(selfvhf != nullptr)
      {
        dtmp = weight*scf;
        if(nth<0) selfvhf->dataleft += dtmp;
        else if(nth >= selfvhf->ndata) selfvhf->data[selfvhf->ndata] += dtmp;
        else selfvhf->data[nth] += dtmp;
        if(lognth<0) selfvhf->logdataleft += dtmp;
        else if(lognth >= selfvhf->logndata) selfvhf->logdata[selfvhf->logndata] += dtmp;
        else selfvhf->logdata[lognth] += dtmp;
      }
      if(corrfunc != nullptr)
      {
        dtmp = weight*scfi;
        if(nth<0) corrfunc->data[nth] += dtmp;
        else if(nth >= corrfunc->ndata) corrfunc->data[corrfunc->ndata] += dtmp;
        else corrfunc->data[nth] += dtmp;
        if(lognth<0) corrfunc->logdataleft += dtmp;
        else if(lognth >= corrfunc->logndata) corrfunc->logdata[corrfunc->logndata] += dtmp;
        else corrfunc->logdata[lognth] += dtmp;
      }
    }
    // O(n) part
    msd /= cluv-s1.weight;
    qsd /= cluv-s1.weight;
    msd_ += s1.weight*msd;
    qsd_ += s1.weight*qsd;
    // massdist
    if(massdist != nullptr)
    {
      rij = 0.;
      for(rr=0; rr<DIM; ++rr)
      {
        dtmp = s1.coord[rr] - center[rr];
        rij += dtmp*dtmp;
      }
      rij = sqrt(rij);
      nth = (rij-para->binstart)*para->invdr;
      lognth = (log(rij)-para->logbinstart)*para->loginvdr;
      
      dtmp = s1.weight*scfm;
      if(nth<0) massdist->dataleft += dtmp;
      else if(nth >= massdist->ndata) massdist->data[massdist->ndata] += dtmp;
      else massdist->data[nth] += dtmp;
      if(lognth<0) massdist->logdataleft += dtmp;
      else if(lognth >= massdist->logndata) massdist->logdata[massdist->logndata] += dtmp;
      else massdist->logdata[lognth] += dtmp;
    }
    if(msddist != nullptr)
    {
      dtmp = s1.weight*scfms;
      rij = sqrt(msd);
      nth = (rij-para->binstart)*para->invdr;
      lognth = (log(rij)-para->logbinstart)*para->loginvdr;
      if(nth<0) msddist->dataleft += dtmp;
      else if(nth >= para->ndata) msddist->data[para->ndata] += dtmp;
      else msddist->data[nth] += dtmp;
      if(lognth<0) msddist->logdataleft += dtmp;
      else if(lognth >= para->logndata) msddist->logdata[para->logndata] += dtmp;
      else msddist->logdata[lognth] += dtmp;
    }
  }
  msd_ /= cluv;
  qsd_ /= cluv;
  chidyn = qsd_ - msd_*msd_;
  return 0;
}

void AllStat::dumpAll(const std::string &filename)
{
  if(selfvhf!=nullptr) selfvhf->dump(filename+"_svh");
  if(corrfunc!=nullptr) corrfunc->dump(filename+"_corrf");
  if(massdist!=nullptr) massdist->dump(filename+"_massd");
  if(msddist!=nullptr) msddist->dump(filename+"_md");
}

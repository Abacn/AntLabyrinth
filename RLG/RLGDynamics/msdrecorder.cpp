//
//  msdrecorder.cpp
//  RLGDynamics
//
//  Created by Yi Hu on 8/6/19.
//

#include <sstream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <climits>
#include <numeric>

#include "msdrecorder.hpp"
#include "myutility.hpp"

// Record MSD on-the-fly int sleft, int sinterval_, int sright
// lefttime = pow(2.0, sleft)/sqrt((double)DIM)
// span = sright - sleft
MSDRecorder::MSDRecorder(double lefttime_, int sinterval_, int span, int nmult_/*=10*/)
: nmult(nmult_), nmultspl(1L << nmult_), sinterval(sinterval_), lefttime(lefttime_), lastsd(0.0)
{
  sqsumcount = 0;
  sqsums = qqsums = 0.0;
  sn = span; // those t>tmax/2 will not be recorded in this scheme
  tpref = new double[sinterval];
  r2s = new double*[sinterval];
  r4s = new double*[sinterval];
  rcov = new double*[sinterval];
  // logarithmic mean cage
  logr2s = new double*[sinterval];
  logr4s = new double*[sinterval];
  counts = new long*[sinterval];
  for(long i=0; i<sinterval; ++i)
  {
    tpref[i] = lefttime*pow(2.0, (double)i/sinterval);
    r2s[i] = new double[sn+1];
    r4s[i] = new double[sn+1];
    rcov[i] = new double[sn+1];
    logr2s[i] = new double[sn+1];
    logr4s[i] = new double[sn+1];
    counts[i] = new long[sn+1];
    for(long j=0; j<=sn; ++j)
    {
      r2s[i][j] = 0.0;
      r4s[i][j] = 0.0;
      rcov[i][j] = 0.0;
      logr2s[i][j] = 0.0;
      logr4s[i][j] = 0.0;
      counts[i][j] = 0;
    }
  }
}

double MSDRecorder::getfinalsd() const
{
  return lastsd;
}

void MSDRecorder::getfinalmsd(double *sqs, double *qqs) const
{
  if(sqsumcount>1)
  {
    // get cage size in single sample by averaging the final position.
    *sqs = sqsums/(double)sqsumcount;
    *qqs = qqsums/(double)sqsumcount;
  }
  else
  {
    *sqs = r2s[0][sn]/counts[0][sn];
    *qqs = r4s[0][sn]/counts[0][sn];
  }
}

int MSDRecorder::getinitfinalmsd(double *sqs, double *qqs) const
{
  // get cage size in single sample by averaging the MSD-t curve.
  // aka averaged over both initial and final positions of the tracer.
  // Only used in sphere geometry and in high density.
  int i=0, j, k=0, retval=3;
  double mintime = 20.0/sqrt((double)DIM);
  bool lefttimepassed = false;
  std::vector<double> r2sarray, r4sarray;
  while(i<sinterval*sn)
  {
    for(j=0; j<sinterval; ++j)
    {
      if(counts[j][k] > 8) // >O(1) to have enough samples calculating variance
      {
        ++i;
        if(lefttimepassed || lefttime*pow(2, k)*tpref[j] > mintime)
        {
          r2sarray.push_back(r2s[j][k]/counts[j][k]);
          r4sarray.push_back(r4s[j][k]/counts[j][k]);
        }
      }
      else
      {
        goto outloop;
      }
    }
    ++k;
  }
outloop:;
  // analysis
  if(r2sarray.size() >= 10)
  {
    int divide1 = (int)(r2sarray.size())/3, divide2 = (int)r2sarray.size()*2/3;
    double r2part1 = std::accumulate(r2sarray.begin(), std::next(r2sarray.begin(), divide1), 0.0) / divide1,
           r2part2 = std::accumulate(std::next(r2sarray.begin(), divide1), std::next(r2sarray.begin(), divide2), 0.0) / (divide2 - divide1),
           r2part3 = std::accumulate(std::next(r2sarray.begin(), divide2), r2sarray.end(), 0.0) / (r2sarray.size() - divide2);
    double r4part1 = std::accumulate(r4sarray.begin(), std::next(r4sarray.begin(), divide1), 0.0) / divide1,
           r4part2 = std::accumulate(std::next(r4sarray.begin(), divide1), std::next(r4sarray.begin(), divide2), 0.0) / (divide2 - divide1),
           r4part3 = std::accumulate(std::next(r4sarray.begin(), divide2), r4sarray.end(), 0.0) / (r4sarray.size() - divide2);
    if(r2part1/r2part2 < 0.9)
    {
      // going up, cut first 1/3
      *sqs = 0.5*(r2part1+r2part2);
      *qqs = 0.5*(r4part1+r4part2);
      retval = 1;
    }
    else if(r2part3/r2part2 < 0.9)
    {
      // going down, cut last 1/3
      *sqs = 0.5*(r2part3+r2part2);
      *qqs = 0.5*(r4part3+r4part2);
      retval = 1;
    }
    else
    {
      *sqs = (r2part1+r2part2+r2part3)/3;
      *qqs = (r4part1+r4part2+r4part3)/3;
      retval = 0;
    }
  }
  else
  {
    *sqs = std::accumulate(r2sarray.begin(), r2sarray.end(), 0.0)/r2sarray.size();
    *qqs = std::accumulate(r4sarray.begin(), r4sarray.end(), 0.0)/r2sarray.size();
  }
  return retval;
}

void MSDRecorder::dump(double *r2s_, double *r4s_, double *chidis_, double *rcov_, int *samplecount_, int sz)
{
  // add to overwritearr
  int i=0, j, k=0;
  while(i<sz)
  {
    for(j=0; j<sinterval; ++j)
    {
      if(counts[j][k] > 8) // >O(1) to have enough samples calculating variance
      {
        double dtmp = r2s[j][k]/counts[j][k];
        r2s_[i] += dtmp;
        r4s_[i] += r4s[j][k]/counts[j][k];
        rcov_[i] += rcov[j][k]/counts[j][k];
        chidis_[i] += dtmp*dtmp;
        ++samplecount_[i];
        ++i;
      }
      else
      {
        goto outloop;
      }
    }
    ++k;
  }
outloop:;
}

void MSDRecorder::dumplog(double *logr2s_, double *logr4s_, int sz)
{
  // add to overwritearr
  int i=0, j, k=0;
  while(i<sz)
  {
    for(j=0; j<sinterval; ++j)
    {
      if(counts[j][k] > 8) // >O(1) to have enough samples calculating variance
      {
        double dtmp = logr2s[j][k]/counts[j][k];
        logr2s_[i] += dtmp;
        logr4s_[i] += logr4s[j][k]/counts[j][k];
        ++i;
      }
      else
      {
        goto outloop;
      }
    }
    ++k;
  }
outloop:;
}

void MSDRecorder::dumpone(const char* fname)
{
  std::ofstream ofs(fname);
  int j, k=0;
  for(k=0;k<=sn;++k)
  {
    for(j=0; j<sinterval; ++j)
    {
      if(counts[j][k] > 0)
      {
        double tnow = tpref[j]*pow(2.0, k);
        double r2snow = r2s[j][k]/counts[j][k];
        double r4snow =  r4s[j][k]/counts[j][k];
        ofs << std::setprecision(10);
        ofs << tnow << '\t';
        ofs << std::setprecision(6);
        ofs << r2snow << '\t' << r4snow << '\t' << 1 << '\n';
      }
      else
      {
        goto outloop;
      }
    }
  }
outloop:;
  ofs.close();
}

MSDRecorder::~MSDRecorder()
{
  for(long i=0; i<sinterval; ++i)
  {
    delete [] r2s[i];
    delete [] r4s[i];
    delete [] rcov[i];
    delete [] logr2s[i];
    delete [] logr4s[i];
    delete [] counts[i];
  }
  delete [] tpref;
  delete [] r2s;
  delete [] r4s;
  delete [] rcov;
  delete [] logr2s;
  delete [] logr4s;
  delete [] counts;
}


AveMSDRecorder::AveMSDRecorder(double lefttime_, int sinterval_, int span, int nmult_/*=10*/)
: MSDRecorder(lefttime_, sinterval_, span, nmult_)
{
  // initialize snapshot queue and related paras
  Vector zerovec;
  zerovec.setzero();
  smark = new long[sinterval];
  nextrecs = new double[sinterval];
  inextrecs = new long[sinterval];
  snapshotsarr = new Vector*[sinterval];
  snapshotsque = new limited_queue<Vector>*[sinterval];
  nextarr = new long*[sinterval];
  if(sn <= nmult) // no short time scheme
  {
    spart1 = 0;
    spart2 = sn;
  }
  else
  {
    spart1 = sn - nmult;  // number of first part stat
    spart2 = nmult;       // number of second part stat = 2^spart2-1;
  }
  for(long i=0; i<sinterval; ++i)
  {
    smark[i] = 0;
    nextrecs[i] = tpref[i];
    inextrecs[i] = 0;
    if(spart1 > 0)
    {
      snapshotsarr[i] = new Vector[spart1];
      nextarr[i] = new long[spart1];
      for(long j=0; j<spart1; ++j)
      {
        snapshotsarr[i][j] = zerovec;
        nextarr[i][j] = 1L << j;
      }
    }
    if(spart2>0)
    {
      snapshotsque[i] = new limited_queue<Vector>(pow(2,spart2-1));
      snapshotsque[i]->push(zerovec);
    }
  }
  nextrec = nextrecs[0];
}

/**
 * record square displacement if needed.
 * return value: 0: not record; 1: record; 2: record last one;
 * ((not activated) -1: record sq>DIM)
 *
 * For each time interval, take average over different begin and end times.
 */
int AveMSDRecorder::record(const Tracer &tracer, double tnow, double tnext)
{
  long i, imark, timeoffset2;
  int ret_val = 0;
  double t_elapse, sd, sq, logsd;
  //only record when current time exceeds next record time
  if(nextrec < tnext)
  {
    nextrec = INFINITY;
    ret_val = 1;
    for(i=0; i<sinterval; ++i)
    {
      while(nextrecs[i] < tnext)
      {
        t_elapse = nextrecs[i] - tnow;
        assert(t_elapse >= 0.);
        Vector tx = tracer.getX() + tracer.v*t_elapse;  // tracer position at check point
        imark = smark[i];
        while(imark<spart1) // part 1 stat
        {
          auto tmpvec  = snapshotsarr[i][imark] - tx;
          tmpvec.quadsum(sd, sq);
          logsd = log(sd);
          snapshotsarr[i][imark] = tx;
          r2s[i][imark] += sd;
          r4s[i][imark] += sd*sd;
          rcov[i][imark] += sd*sd - sq;
          logr2s[i][imark] += logsd;
          logr4s[i][imark] += logsd*logsd;
          counts[i][imark]++;
          nextarr[i][imark] += 1L<<imark;
          if(counts[i][imark] == nmultspl)
          {
            ++smark[i];
          }
          if(!(nextarr[i][imark] & (1L<<imark))) break;
          ++imark;
        }
        if(imark == spart1) // part 2 stat
        {
          if(spart2>0)
          {
            for(timeoffset2 = 1L; timeoffset2 <= snapshotsque[i]->size(); timeoffset2 <<= 1)
            {
              auto tmpvec = (*snapshotsque[i])[snapshotsque[i]->size() - timeoffset2] - tx;
              tmpvec.quadsum(sd, sq);
              logsd = log(sd);
              /*
              if(sd>DIM)
              {
                ret_val = -1;
                break;
              }*/
              r2s[i][imark] += sd;
              r4s[i][imark] += sd*sd;
              rcov[i][imark] += sd*sd - sq;
              logr2s[i][imark] += logsd;
              logr4s[i][imark] += logsd*logsd;
              counts[i][imark]++;
              ++imark;
            }
            snapshotsque[i]->push(tx);
          }
          // sq displacement of tracer itself
          lastsd = tx.norm_squared();
          sqsums += lastsd;
          qqsums += lastsd*lastsd;
          sqsumcount++;
          // the tail
          if(inextrecs[i] == 1L<<sn)
          {
            logsd = log(lastsd);
            r2s[i][sn] += lastsd;
            r4s[i][sn] += lastsd*lastsd;
            logr2s[i][sn] += logsd;
            logr4s[i][sn] += logsd*logsd;
            counts[i][sn]++;
            ret_val = 2;
          }
        }
        if(smark[i] < spart1)// resolve case 1 nextrecs[i]
        {
          inextrecs[i] = nextarr[i][smark[i]];
          nextrecs[i] =tpref[i]*inextrecs[i];
        }
        else // resolve case 2 nextrecs[i]
        {
          inextrecs[i] = (1L << spart1)*(counts[i][spart1]+1);
          nextrecs[i] =tpref[i]*inextrecs[i];
        }
      }
      if(nextrecs[i] < nextrec)
      {
        nextrec = nextrecs[i];
      }
    }
    if(nextrecs[i] < nextrec)
    {
      nextrec = nextrecs[i];
    }
  }
  return ret_val;
}

AveMSDRecorder::~AveMSDRecorder()
{
  for(int i=0; i<sinterval; ++i)
  {
    if(spart1 > 0)
    {
      delete [] snapshotsarr[i];
      delete [] nextarr[i];
    }
    if(spart2>0)
    {
      delete snapshotsque[i];
    }
  }
  delete [] smark;
  delete [] nextrecs;
  delete [] inextrecs;
  delete [] snapshotsarr;
  delete [] snapshotsque;
  delete [] nextarr;
}


GaussianMSDRecorder::GaussianMSDRecorder(double lefttime_, int sinterval_, int span, int nmult_/*=10*/)
: MSDRecorder(lefttime_, sinterval_, span, nmult_)
{
  smark = new long[nmultspl];
  velocities = new double[nmultspl];
  accutimes = new double[nmultspl];
  nextrecs = new double[nmultspl];
  lastupdtimes = new double[nmultspl];
  for(int i=0; i<nmultspl; ++i)
  {
    smark[i] = 0;
    velocities[i] = 1.0;
    accutimes[i] = 0.0;
    lastupdtimes[i] = 0.0;
  }
  reassignvelocity(0.0);
}

/**
 * record square displacement if needed.
 * return value always 0.
 *
 * Averaged over velocity assignments (see reassignvelocity)
 */
int GaussianMSDRecorder::record(const Tracer &tracer, double tnow, double tnext)
{
  double t_elapse, sq, logsd;
  long i, j, jmark;
  for(i=0; i<nmultspl; ++i)
  {
    if(smark[i] > sinterval*sn) continue;
    while(nextrecs[i] < tnext)
    {
      t_elapse = nextrecs[i] - tnow;
      assert(t_elapse >= 0.);
      Vector tx = tracer.getX() + tracer.v*t_elapse;  // tracer position at check point
      tx.quadsum(lastsd, sq);
      logsd = log(lastsd);
      j = smark[i] % sinterval;
      jmark = smark[i] / sinterval;
      r2s[j][jmark] += lastsd;
      r4s[j][jmark] += lastsd*lastsd;
      rcov[j][jmark] += lastsd*lastsd - sq;
      logr2s[j][jmark] += logsd;
      logr4s[j][jmark] += logsd*logsd;
      counts[j][jmark]++;
      ++smark[i];
      double nextt = tpref[smark[i] % sinterval] * (1 << (smark[i] / sinterval));
      double delt = nextt - accutimes[i];
      nextrecs[i] = lastupdtimes[i] + delt * velocities[i];
    }
  }
  return 0;
}

/**
 * Rescale the velocity rate with a series of chi-squared distributed random
 * variable (in d=3 chi-squared is Maxwell-Boltzmann distribution).
 */
int GaussianMSDRecorder::reassignvelocity(double timestamp)
{
  for (int i=0; i<nmultspl; ++i)
  {
    // reset velocity with probability 1/d
    if(smark[i] > sinterval*sn || (timestamp > 0 && MyUtility::rand() > 1.0/DIM)) continue;
    double t_elapsed = timestamp - lastupdtimes[i];
    // increment accutimes
    accutimes[i] += t_elapsed / velocities[i];
    // lin (t->0) <\hat Delta / \hat t> = 1 in all dimension
    velocities[i] = sqrt(MyUtility::randc() / DIM);
    // recalculate nextrecs
    // next rec in gaussian-scaled scale
    double nextt = tpref[smark[i] % sinterval] * (1 << (smark[i] / sinterval));
    double delt = nextt - accutimes[i];
    assert(delt >= 0.0);
    // next rec in uniform scale
    nextrecs[i] = timestamp + delt * velocities[i];
    lastupdtimes[i] = timestamp;
  }
  return 0;
}

bool GaussianMSDRecorder::samplefinished()
{
  return counts[0][sn] == nmultspl;
}

GaussianMSDRecorder::~GaussianMSDRecorder()
{
  delete [] smark;
  delete [] velocities;
  delete [] accutimes;
  delete [] nextrecs;
  delete [] lastupdtimes;
}


EscapeRecorder::EscapeRecorder(double Dcutleft_, double Dcutright_, double dDelta_, int escapecountlimit_):
Dcutleft(Dcutleft_), Dcutright(Dcutright_), dDelta(dDelta_), DeltaNext(Dcutleft_), crtesc(0),
nextfntmsdidx(0), escapecountlimit(escapecountlimit_)
{
  if(Dcutright >= Dcutleft)
  {
    nesctime = (int)round((Dcutright-Dcutleft)/dDelta_)+1;
    esctimes = new double[nesctime];
    for(int rp=0; rp<nesctime; ++rp)
    {
      esctimes[rp] = 0.;
    }
  }
  else
  {
    // disabled
    nesctime = 0;
    esctimes = nullptr;
  }
}

/** Record first escape time at a shell radius, if needed. */
int EscapeRecorder::record(const Tracer &tracer, double tnow, double t_elapse)
{
  // escapecountlimit-0: do not record escape; 1-only rcord escape time; 2+, dump finite size scaling
  if(escapecountlimit == 0 || crtesc >= nesctime) return -1; // already full
  Vector endvec = tracer.coord + tracer.v*t_elapse;
  double nownorm = 0.0, xdotv = 0.0;
  double endnorm = endvec.norm_squared();
  if(DeltaNext <= endnorm)
  {
    // find t_elapse at DeltaNext
    nownorm = tracer.coord.norm_squared();
    xdotv = tracer.coord.dot(tracer.v);
  }
  while(DeltaNext <= endnorm)
  {
    double t_esc = MyUtility::rootsolverB(1.0, 2.0*xdotv, nownorm-DeltaNext);
    esctimes[crtesc++] = tnow + t_esc;
    DeltaNext += dDelta;
  }
  if(crtesc >= nextdumps[nextfntmsdidx] && crtesc <= escapecountlimit)
  {
    while(crtesc >= nextdumps[nextfntmsdidx+1]) nextfntmsdidx++;
    // finite size scaling msd
    return nextdumps[nextfntmsdidx++];
  }
  else
  {
    return 0;
  }
}

void EscapeRecorder::dump(int count)
{
  if(nesctime == 0) return;
  std::stringstream sst;
  sst << "escape_" << count << ".dat";
  std::ofstream ofs(sst.str());
  ofs << std::setprecision(10);
  for(int rp=0; rp<crtesc; ++rp)
  {
    ofs << Dcutleft+dDelta*rp << "\t" << esctimes[rp] << "\n";
  }
  ofs.close();
}

EscapeRecorder::~EscapeRecorder()
{
  if(nesctime > 0) delete [] esctimes;
}

const int EscapeRecorder::nextdumps[] = {4,8,10,12,14,16,20,24,30,40,50,INT_MAX};

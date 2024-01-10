//
//  hcomb.cpp
//  RLGDynamics
//
//  Created by Yi Hu on 8/5/19.
//

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdio>

#include "hcomb.hpp"
#include "myutility.hpp"

using namespace HCombPBC;

#ifndef HCOMBPBC_IS_LEECH
// inscribed radius of the pbc box
constexpr double SQRINSCRIBEDR=0.5;
constexpr double CIRCURMSCRIBER=1.0;
#else
constexpr double SQRINSCRIBEDR=1.0;
constexpr double CIRCURMSCRIBER=M_SQRT2;
#endif

HComb::HComb(ReadInputDy &input):
input(input), obstacles(input.N), count(0),
scale(pow(input.N/input.phi*MyUtility::getVunitsphere()/Vcell, -1.0/DIM)),
sqscale(scale*scale), quadscale(sqscale*sqscale), invscale(1./scale), invsqscale(invscale*invscale),
maxt(pow(2.0, input.texp)), maxsinglet(sqrt(SQRINSCRIBEDR)/scale - 1.),
ncollision(0L), nescape(0L), nlistsize(0L)
{
  std::cout << "     scale : " << scale << std::endl;
  if(input.rcut > sqrt(SQRINSCRIBEDR)*invscale)
  {
    std::cerr << "warning: initial rcut too large. Automatically adjusted to " << sqrt(SQRINSCRIBEDR)*invscale << std::endl;
    updatercut(sqrt(SQRINSCRIBEDR));
  }
  else
  {
    updatercut(input.rcut*scale);
  }
  // validity check
  assert(maxsinglet > 0.);
  // adjust maxt according to t scaling
  if(input.tscaleflag == 2 || input.tscaleflag == 3)
  {
    maxt /= sqrt((double)DIM);
  }
  else if(input.tscaleflag == 4)
  {
    maxt /= DIM;
  }
  // neighbor list initialization
  if(input.rcuttype != 2)
  {
    // nlists = new NeighborList<norm_squared_PBC>[input.N];
    std::cerr << "warning: rcuttype unsupported for hcomb box. Use 2 instead" << std::endl;
    input.rcuttype = 2;
  }
  if(sqscale > SQRINSCRIBEDR)
  {
    // obstacle radius greater than half box length, PBC fail
    std::cout << "PBC fail! exit..." << std::endl;
    exit(-1);
  }
  nsampletime = (input.texp - input.sleft)*input.sinterval + 1;
  timearray = new double[nsampletime];
  r2array = new double[nsampletime];
  r4array = new double[nsampletime];
  logr2array = new double[nsampletime];
  logr4array = new double[nsampletime];
  chidisarray = new double[nsampletime];
  rcovarray = new double[nsampletime];
  samplecount = new int[nsampletime];
  for(int i=0; i<nsampletime; ++i)
  {
    timearray[i] = pow(2., input.sleft + (double)i/input.sinterval);
    if(2 == input.tscaleflag || 3 == input.tscaleflag) timearray[i]/=sqrt((double)DIM);
    else if(4 == input.tscaleflag) timearray[i]/=DIM;
    r2array[i] = r4array[i] = chidisarray[i] = rcovarray[i] = 0.;
    logr2array[i] = logr4array[i] = 0.;
    samplecount[i] = 0;
  }
  nextrecord = 0;
  // ~ 5 brownian collision per record or 5 per \hat{t}
  // std::min(0.1, 1.0/(2 * input.phi)) / sqrt(DIM*1.0)
  br.setup(timearray[0], 0.5/sqrt(DIM*1.0));
}

HComb::~HComb()
{
  if(input.outmode >= OUTPUT_MODE_PRODUCT)
  {
    ofd.close();
  }
  delete []timearray;
  delete []r2array;
  delete []r4array;
  delete []logr2array;
  delete []logr4array;
  delete []chidisarray;
  delete []rcovarray;
  delete []samplecount;
}

// generate random distribution
void HComb::genrandomPoints()
{
  int32_t i=0, j=0, k=0, mirror[DIM];
// For Dn, covering radius is 1 (DIM<=4) or sqrt(DIM)/2 (DIM>4)
// For Dn+, covering radius is 1(d=8),
  double norm, dsum, dtmp;
  while(i<input.N)
  {
    dsum = 0.;
    for(j=0; j<DIM; ++j)
    {
      // [-sqrt(2), sqrt(2)]
      dtmp = (MyUtility::rand()*2.0 - 1.0)*CIRCURMSCRIBER;
      obstacles[i][j] = dtmp;
    }
    latticepoint(obstacles[i].x, mirror);
    ++k;
    norm = obstacles[i].norm_squared();
    if(sqscale < norm) ++i; // make a cavity at the origin
  }
//#ifdef DEBUG
  std::cout << " Random point accept rate: " << (double)input.N/k << std::endl;
//#endif
  if(input.outmode >= OUTPUT_MODE_DEBUG)
  {
    std::stringstream sst;
    sst << "obstacles_" << count << ".dat";
    std::ofstream ofs(sst.str());
    ofs << scale << "\n"; // scale
    for(i=0; i<input.N; ++i)
    {
      for(j=0; j<DIM; ++j)
      {
        ofs << obstacles[i][j] << "\t";
      }
      ofs << "\n";
    }
  }
  tracer.setzero();
  if(input.outmode >= OUTPUT_MODE_DEBUG)
  {
    std::stringstream sst;
    sst << "traj_" << count << ".dat";
    oftraj.open(sst.str());
    for(i=0; i<DIM; ++i)
    {
      oftraj << tracer.coord[i] << "\t";
    }
    oftraj << "\n";
  }
}

// random initial velocity
void HComb::genrandomVelocity()
{
  int i;
  double dtmp;
  if(input.tscaleflag == 4)
  {
    // Brownian, scaled by \hat{\Delta(\hat{t} = 1)} = 1 in small density limit
    dtmp = scale * 1.0 / sqrt(br.get_delt()*DIM);
    for(i=0; i<DIM; ++i)
    {
      tracer.v[i] = MyUtility::randn() * dtmp;
    }
  }
  else
  {
    double dsum = 0.;
    for(i=0; i<DIM; ++i)
    {
      dtmp = MyUtility::randn();
      tracer.v[i] = dtmp;
      dsum += dtmp*dtmp;
    }
    dsum = scale/sqrt(dsum);
    for(i=0; i<DIM; ++i)
    {
      tracer.v[i] *= dsum;
    }
  }
}

Point_handle HComb::getnextcollide(Point_handle lastcollide, double &t_elapse)
{
  double tmin = INFINITY, t;
  Point_handle idxmin = NULL_COLLISION;
  double A, B, C;
  int i, rp;
  A = tracer.v.norm_squared();
  Vector vec;
  
  for(const auto &pv: nlist.neighbors)
  {
    i = pv.first;
    const Vector &obi = obstacles[i], &pvs = pv.second;
    for(rp=0; rp<DIM; ++rp)
    {
      vec[rp] = obi[rp] - tracer.coord[rp] + pvs[rp];
    }
    B = vec.dot(tracer.v);
    if(B <= 0.) continue; // away
    C = vec.norm_squared();
    // if(0.5 < C) continue; // longer than half box
    B *= -2.0;
    C -= sqscale;
    if(C < 0.0)
    {
      std::cerr << "Assertion Fail: C: " << C << "<0. nowt: " << nowt << ", lastcollide: " << lastcollide << std::endl;
      abort();
    }
    t = MyUtility::rootsolverA(A, B, C);
    if(t < tmin)
    {
      tmin = t;
      idxmin = i;
    }
  }
  vec = tracer.coord - nlisttracer;
  if(std::isfinite(tmin))
  {
    // check if the collision point is outside the (neighbor shell - 1)
    Vector vecb = vec + tracer.v*tmin;
    if(vecb.norm_squared() > rtol)  // runs outside of the shell
    {
      tmin = INFINITY;
    }
  }

  if(tmin == INFINITY) // noncollide
  {
    B = 2.0*vec.dot(tracer.v);
    C = vec.norm_squared() - rtol;
    nlist.deconstruct(); // neighbor list expired
    if(C > 0.)
    {
      std::cerr << "Assertion Fail: C: " << C << ">0. nowt: " << nowt << ", rtol: " << rtol << std::endl;
      abort();
    }
    tmin = MyUtility::rootsolverB(A, B, C);
    idxmin = NULL_COLLISION;
  }
  
  t_elapse = tmin;
  return idxmin;
}

void HComb::record(std::vector<double> const &vals)
{
  if(input.outmode >= OUTPUT_MODE_PRODUCT)
  {
    if(0 == count)
    {
      ofd.open("dbgsummary.dat");
      ofd << "NCollision\tNescape\tFinalv\tFinalDelta\tFinalMSD\tFinalMQD\n";
    }
    ofd << ncollision << '\t' << nescape << '\t' << tracer.v.norm_squared()/sqscale << std::setprecision(8);
    for(int rp=0; rp<vals.size(); ++rp) ofd << '\t' << vals[rp];
    ofd << std::endl;
  }
  
}

void HComb::updatetracer(Point_handle thiscollide, double t_elapse)
{
  // update tracer position and velocity
  int i;
  tracer.coord += tracer.v*t_elapse;
  if(thiscollide == BROWNIAN_COLLISION)
  {
    // reset velocity at brownian time interval. Refer to De Michele's
    // event-driven Brownian dynamics algorithm.
    genrandomVelocity();
  }
  else if(thiscollide != NULL_COLLISION)
  {
    // collided with obstacle
    ++ncollision;
    // update v
    Vector n = obstacles[thiscollide] - tracer.coord;
    pbc(n);
    double f = -2.0*tracer.v.dot(n)/n.norm_squared();
    for(i=0; i<DIM; ++i)
    {
      tracer.v[i] += f*n[i];
    }
    if (input.tscaleflag==3)
    {
      // For gaussian distributed velocities in Newtonian dynamics, reassign velocity rate at collision
      dynamic_cast<GaussianMSDRecorder*>(msdrecorder)->reassignvelocity(nowt + t_elapse);
    }
  }
  if(input.outmode >= OUTPUT_MODE_DEBUG)
  {
    for(i=0; i<DIM; ++i)
    {
      oftraj << tracer.coord[i] + tracer.pbcvec[i] << "\t";
    }
    if(thiscollide == NULL_COLLISION)
    {
      oftraj << "\t" << -1 << std::endl;
    }
    else
    {
      oftraj << "\t" << thiscollide << std::endl;
    }
  }
}

void HComb::updatercut(double newrcut)
{
  rcut = newrcut;
  sqrcut = rcut*rcut;
  rtol = (rcut-scale)*(rcut-scale);
  rconstC = rcut*(2*scale-rcut);
  nlist.deconstruct();
}

void HComb::updatenlist()
{
  int32_t mirror[DIM];
  int i;
  if(nlist.status) nlist.deconstruct();

  // first move the tracer inside of the box
#ifdef DNPLUSFLAG
  int res = latticepoint(tracer.coord.x, mirror);
  if(res == 2)
  {
    for(i=0; i<DIM; ++i)
    {
      tracer.pbcvec[i] += 0.5; // Dn+/Dn0+ symmetry
    }
  }
#else
  latticepoint(tracer.coord.x, mirror);
#endif // DNPLUSFLAG
  for(i=0; i<DIM; ++i)
  {
#ifndef HCOMBPBC_IS_LEECH
    tracer.pbcvec[i] += mirror[i];
#else
    tracer.pbcvec[i] += mirror[i] / (2*M_SQRT2);
#endif // HCOMBPBC_IS_LEECH
  }
  nlist.construct(nlisttracer = tracer.coord, sqrcut, obstacles);
  // increment nlist build count
  ++nescape;
  nlistsize += nlist.neighbors.size();
}

bool HComb::samplefinished()
{
  if (input.tscaleflag==3)
  {
    // all gaussian velocities reached max time
    return dynamic_cast<GaussianMSDRecorder*>(msdrecorder)->samplefinished();
  }
  else
  {
    // uniform velocity
    return nowt > maxt;
  }
}

int HComb::process()
{
  using std::chrono::steady_clock;
  // the obstacle index of last collision. -1 means not hit
  double t_elapse;
  uint64_t nlastcollision = 0;
  Point_handle lastcollide = NULL_COLLISION, thiscollide;
  updatercut(input.rcut*scale);
  genrandomPoints();
  genrandomVelocity();
  double sleft = input.tscaleflag==1 ? pow(2.0, input.sleft) :
      input.tscaleflag==4 ? pow(2.0, input.sleft)/DIM :
      pow(2.0, input.sleft)/sqrt((double)DIM);
  if (input.tscaleflag==3)
  {
    msdrecorder = new GaussianMSDRecorder(sleft, input.sinterval, input.texp-input.sleft);
  }
  else
  {
    msdrecorder = new AveMSDRecorder(sleft, input.sinterval, input.texp-input.sleft);
  }
  nowt = 0.;
  nextrecord = 0;
  ncollision = nescape = nlistsize = 0L;
  if(input.tscaleflag == 4)
  {
    br.reset();
  }
 
  std::ofstream fstatus;
  bool recordstatus = false;
  clock_type t1b, t2b, t1c, t2c;
  duration_type tcollide, tnlist;
  if(input.repeatrun==1)
  {
    recordstatus = true;
    t1b = steady_clock::now();
  }
  t1c = steady_clock::now();
  // main loop for simulation. Process one event in each loop.
  while(!samplefinished())
  {
#ifdef DEBUG
    clock_type t1, t2, t3, t4, t5;
    t1 = steady_clock::now();
#endif
    if(!nlist.status)
    {
      t2c = steady_clock::now();
      updatenlist();
      tnlist += steady_clock::now() - t2c;
    }
#ifdef DEBUG
    t2 = steady_clock::now();
#endif
    thiscollide = getnextcollide(lastcollide, t_elapse);

    if(input.tscaleflag == 4 && br.reached(nowt + t_elapse))
    {
      t_elapse = br.inc() - nowt;
      thiscollide = BROWNIAN_COLLISION;
    }

#ifdef DEBUG
    t3 = steady_clock::now();
#endif
    // record MSD and other metrics if needed
    msdrecorder->record(tracer, nowt, nowt + t_elapse);
#ifdef DEBUG
    t4 = steady_clock::now();
#endif
    updatetracer(thiscollide, t_elapse);
#ifdef DEBUG
    t5 = steady_clock::now();
    nlisttime += t2-t1;
    collidetime += t3-t2;
    stattime += t4-t3;
    intgtime += t5-t4;
    if(!(ncollision&0xffff) && ncollision != nlastcollision) // 65536 per collision
#else
    if(!(ncollision&0xfffff) && ncollision != nlastcollision) // 1048576~10^6 per collision
#endif
    {
      // determine if update rcut or not
      t2c = steady_clock::now();
      tcollide = t2c - t1c - tnlist;
      double ecratio = tnlist/tcollide;
      nlastcollision = ncollision;
      if(ecratio > 1.0)
      {
        double newrcut = rcut*pow(1.05, 1.0/DIM);
        if(newrcut < sqrt(SQRINSCRIBEDR))
        {
          updatercut(newrcut); // increase 5%
#ifdef DEBUG
          std::cout << "ecratio " << ecratio << ", change rcut... " << rcut/scale << std::endl;
#endif
        }
      }
      else if(ecratio < 0.25)
      {
        double newrcut = rcut*pow(0.95, 1.0/DIM);
        if(newrcut > scale)
        {
          updatercut(newrcut); // shrink 5%
#ifdef DEBUG
          std::cout << "ecratio " << ecratio << ", change rcut... " << rcut/scale << std::endl;
#endif
        }
      }
#ifdef DEBUG
      else std::cout << "ecratio " << ecratio << std::endl;
#endif
      t1c = t2c;
      tnlist = steady_clock::now() - t2c;
      if(recordstatus)
      {
        t2b = steady_clock::now();
        if(static_cast<duration_type>(t2b-t1b).count()>3600.0)
        {
          // dump temporary status
          fstatus.open("status.txt");
          fstatus << "(" << nowt << " / " << maxt << ") " << ncollision << "\t" << nescape << "\t" << rcut*invscale << "\n";
          fstatus.close();
          t1b = t2b;
          // dump temporary msd
          msdrecorder->dump(r2array, r4array, chidisarray, rcovarray, samplecount, nsampletime);
          msdrecorder->dumplog(logr2array, logr4array, nsampletime);
          dump();
          // reset
          for(int i=0; i<nsampletime && samplecount[i]>0; ++i)
          {
            r2array[i] = r4array[i] = chidisarray[i] = rcovarray[i] = 0.;
            samplecount[i] = 0;
          }
        }
      }
    }
    lastcollide = thiscollide;
    nowt += t_elapse;
  }
  std::vector<double> vals(3);
  msdrecorder->getfinalmsd( &vals[1], &vals[2] );
  vals[0] = msdrecorder->getfinalsd();
  record(vals);
#ifndef DEBUG
  if(recordstatus) std::remove("status.txt");
#endif
  // clear neighbor list
  nlist.deconstruct();
  // close trajectory file
  if(input.outmode >= OUTPUT_MODE_DEBUG)
  {
    oftraj.close();
  }
  msdrecorder->dump(r2array, r4array, chidisarray, rcovarray, samplecount, nsampletime);
  msdrecorder->dumplog(logr2array, logr4array, nsampletime);
  delete msdrecorder;
  ++count;
  return 0;
}

void HComb::dump()
{
  std::ofstream ofs(input.datafile+".dat");
  ofs << std::setprecision(10);
  for(int i=0; i<nsampletime && samplecount[i]>0; ++i)
  {
    ofs << timearray[i] << "\t" << r2array[i]/samplecount[i]/sqscale
        << "\t" << r4array[i]/samplecount[i]/quadscale << "\t" << chidisarray[i]/samplecount[i]/quadscale
        << "\t" << rcovarray[i]/(samplecount[i]*quadscale*DIM*(DIM-1));
    ofs << "\t" << samplecount[i] << "\n";
  }
  ofs.close();
  if(input.logdatafile.size() > 0)
  {
    ofs.open(input.logdatafile+".dat");
    ofs << std::setprecision(10);
    for(int i=0; i<nsampletime && samplecount[i]>0; ++i)
    {
      ofs << timearray[i] << "\t" << logr2array[i]/samplecount[i] - log(sqscale)
          << "\t" << logr4array[i]/samplecount[i] - log(quadscale);
      ofs << "\t" << samplecount[i] << "\n";
    }
    ofs.close();
  }
}

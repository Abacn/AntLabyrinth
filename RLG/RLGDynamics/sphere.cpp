//
//  sphere.cpp
//  RLGDynamics
//
//  Created by Yi Hu on 7/16/19.
//

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <climits>
#include "sphere.hpp"
#include "myutility.hpp"

Sphere::Sphere(ReadInputDy &input):
input(input), count(0), expectN(input.phi*(pow(input.rmax, DIM) - 1.)), pois(expectN),
rconstA(pow(input.rmax, DIM) - 1.),
shelltol(input.rmax*input.rmax), // or (input.rmax-1.0)*(input.rmax-1.0)
maxt(pow(2.0, input.texp)),
nlistflag(false), nlists(nullptr), ncollision(0L), nescape(0L)
{
  updatercut(input.rcut);
  // validity check
  assert(input.rmax > 1.);
  assert(rtol > 0.);
  // adjust maxt if scale mode is 2
  if(input.tscaleflag != 1)
  {
    maxt /= sqrt((double)DIM);
  }
  // set poisson distribution
  std::cout << "   expectN : " << expectN << std::endl;
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
    if(1 != input.tscaleflag) timearray[i]/=sqrt((double)DIM);
    r2array[i] = 0.;
    r4array[i] = 0.;
    chidisarray[i] = 0.;
    rcovarray[i] = 0.;
    logr2array[i] = 0.;
    logr4array[i] = 0.;
    samplecount[i] = 0;
  }
  nextrecord = 0;
  // do not use neighbor list if the shell size is too small
  if(input.rmax*2 > rcut)
  {
    nlistflag = true;
    nlists = new NeighborList<MyUtility::norm_squared>;
  }
  // additional information output
}

Sphere::~Sphere()
{
  if(input.outmode >= OUTPUT_MODE_PRODUCT)
  {
    ofd.close();
  }
  if(nlistflag)
  {
    delete nlists;
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
void Sphere::genrandomPoints()
{
  int32_t i=0, j=0;
  double norm, dtmp;
  if(input.phi <= 0) // try to read obstacles_*.dat
  {
    std::stringstream sst;
    // sst << "obstacles_0.dat"; // debug usage only, otherwise
    sst << "obstacles_" << count << ".dat";
    std::ifstream ifs(sst.str());
    if(!ifs.is_open())
    {
      std::cerr << "Try to open " << sst.str() << " failed... exit" << std::endl;
      exit(-1);
    }
    N = 0;
    obstacles.clear();
    point obtmp;
    while(1)
    {
      for(j=0; j<DIM; ++j)
      {
        ifs >> obtmp[j];
      }
      if(!ifs.eof())
      {
        obstacles.push_back(obtmp);
        ++N;
      }
      else
      {
        break;
      }
    }
  }
  else
  {
    N = pois(MyUtility::rg);
    obstacles.resize(N);
    for(i=0; i<N; ++i)
    {
      norm = 0.;
      for(j=0; j<DIM; ++j)
      {
        dtmp = MyUtility::randn();
        obstacles[i][j] = dtmp;
        norm += dtmp*dtmp;
      }
      norm = pow(MyUtility::rand()*rconstA + 1.0, 1./DIM)/sqrt(norm);
      for(j=0; j<DIM; ++j)
      {
        obstacles[i][j] *= norm;
      }
    }
    if(input.outmode >= OUTPUT_MODE_DEBUG)
    {
      std::stringstream sst;
      sst << "obstacles_" << count << ".dat";
      std::ofstream ofs(sst.str());
      ofs << 1.0 << "\n"; // scale=1
      for(i=0; i<N; ++i)
      {
        for(j=0; j<DIM; ++j)
        {
          ofs << obstacles[i][j] << "\t";
        }
        ofs << "\n";
      }
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
void Sphere::genrandomVelocity()
{
  int i;
  double dtmp, dsum=0.;
  for(i=0; i<DIM; ++i)
  {
    dtmp = MyUtility::randn();
    tracer.v[i] = dtmp;
    dsum += dtmp*dtmp;
  }
  dsum = 1.0/sqrt(dsum);
  for(i=0; i<DIM; ++i)
  {
    tracer.v[i] *= dsum;
  }
}

Point_handle Sphere::getnextcollide(Point_handle lastcollide, double &t_elapse)
{
  double tmin = INFINITY, t;
  Point_handle idxmin = NULL_COLLISION, i;
  double A, B, C;
  A = tracer.v.norm_squared();
  if(!nlistflag)
  {
    for(i=0; i<N; ++i)
    {
      Vector vec = obstacles[i] - tracer.coord;
      B = vec.dot(tracer.v);
      if(B <= 0.) continue; // away
      B *= -2.0;
      C = vec.norm_squared() - 1.0;
      assert(C >= 0.);
      t = MyUtility::rootsolverA(A, B, C);
      if(t < tmin)
      {
        tmin = t;
        idxmin = i;
      }
    }
  }
  else
  {
    for(auto i: nlists->neighbors) // duplicate
    {
      Vector vec = obstacles[i] - tracer.coord;
      B = vec.dot(tracer.v);
      if(B <= 0.) continue; // away
      B *= -2.0;
      C = vec.norm_squared() - 1.0;
      assert(C >= 0.);
      t = MyUtility::rootsolverA(A, B, C);
      if(t < tmin)
      {
        tmin = t;
        idxmin = i;
      }
    }
    // check if the collision point is outside the (neighbor shell - 1)
    Vector vec = tracer.coord - nlisttracer;
    Vector vecb = vec + tracer.v*tmin;
    if(!std::isfinite(tmin) || vecb.norm_squared() > rtol)  // runs outside of the shell
    {
      B = 2.0*vec.dot(tracer.v);
      C = vec.norm_squared() - rtol;
      tmin = MyUtility::rootsolverB(A, B, C); // cut it then
      idxmin = NULL_COLLISION;
      nlists->deconstruct();
    }
  }
  // check percolation
  if(std::isfinite(tmin))
  {
    Vector vecb = tracer.coord + tracer.v*tmin;
    if( vecb.norm_squared() > shelltol)
    {
      tmin = INFINITY;
      idxmin = NULL_COLLISION;
    }
  }
  t_elapse = tmin;
  return idxmin;
}


void Sphere::record(std::vector<double> const &vals)
{
  // vals: finalt, finaldt, finalmsd, finalmqd, initfinalmsd, initfinalmqd
  if(0 == count)
  {
    ofd.open("dbgsummary.dat");
    ofd << "NCollision\tNescape\tFinalt\tFinalDelta\tFinalMSD\tFinalMQD\tInitFinalMSD\tInitFinalMQD\tNobstacle\tstatusCode\n";
  }
  ofd << ncollision << "\t" << nescape << "\t" << vals[0] << "\t" << std::setprecision(8);
  for(uint32_t rp=1; rp<vals.size(); ++rp) ofd << vals[rp] << "\t";
  ofd << N << "\t" << statuscode << std::endl;
}

void Sphere::updatetracer(Point_handle thiscollide, double t_elapse)
{
  // update tracer position
  int i;
  for(i=0; i<DIM; ++i)
  {
    tracer.coord[i] += t_elapse*tracer.v[i];
  }
  if(thiscollide != NULL_COLLISION)
  {
    // collide with obstacle
    ++ncollision;
    // update v
    Vector n = obstacles[thiscollide] - tracer.coord;
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
  else
  {
    // escape from one neighbor list
    ++nescape;
  }
  if(input.outmode >= OUTPUT_MODE_DEBUG)
  {
    for(i=0; i<DIM; ++i)
    {
      oftraj << tracer.coord[i] << "\t";
    }
    oftraj << thiscollide << std::endl; // thiscollide
  }
}

void Sphere::updatercut(double newrcut)
{
  rcut = newrcut;
  sqrcut = rcut*rcut;
  rtol = (rcut-1.0)*(rcut-1.0);
  rconstC = rcut*(2.0-rcut);
  if(nlistflag) nlists->deconstruct();
}

bool Sphere::samplefinished()
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

int Sphere::process()
{
  using std::chrono::steady_clock;
  
  int retval = 0;
  // the obstacle index of last collision. -1 means not hit
  double t_elapse;
  Point_handle lastcollide = NULL_COLLISION, thiscollide;
  double sleft = input.tscaleflag==1?pow(2.0, input.sleft):pow(2.0, input.sleft)/sqrt((double)DIM);
  if (input.tscaleflag==3)
  {
    msdrecorder = new GaussianMSDRecorder(sleft, input.sinterval, input.texp-input.sleft);
  }
  else
  {
    msdrecorder = new AveMSDRecorder(sleft, input.sinterval, input.texp-input.sleft);
  }
  // in sphere scheme input.rcuttype serves as a cutoff finite size scaling
  EscapeRecorder escrec(1.0/DIM, shelltol, 1.0/DIM, input.rcuttype);
  updatercut(input.rcut);
  genrandomPoints();
  genrandomVelocity();
  nowt = 0.;
  nextrecord = 0;
  ncollision = nescape = 0L;
  statuscode = 0;
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
    // update nlist
    if(nlistflag && !nlists->status)
    {
      t2c = steady_clock::now();
      nlists->construct(nlisttracer = tracer.coord, sqrcut, obstacles);
      tnlist += steady_clock::now() - t2c;
    }
#ifdef DEBUG
    t2 = steady_clock::now();
#endif
    thiscollide = getnextcollide(lastcollide, t_elapse);
    if(t_elapse == INFINITY) // percolated
    {
      // set t_elapse to when the tracer reaxh shelltol
      double b = tracer.coord.dot(tracer.v)*2.0;
      double c = tracer.coord.norm_squared() - shelltol;
      t_elapse = MyUtility::rootsolverB(1.0, b, c);
      msdrecorder->record(tracer, nowt, nowt + t_elapse);
      int rst = escrec.record(tracer, nowt, t_elapse);
      if(rst > 0)
      {
        std::stringstream sst;
        sst << "msd_" << count << "_" << rst << ".dat";
        msdrecorder->dumpone(sst.str().c_str());
      }
      updatetracer(thiscollide, t_elapse);
      statuscode = -1;
      nowt += t_elapse;
      break;
    }
#ifdef DEBUG
    t3 = steady_clock::now();
#endif
    // record MSD and other metrics if needed
    msdrecorder->record(tracer, nowt, nowt + t_elapse);
    int rst = escrec.record(tracer, nowt, t_elapse);
    if(rst > 0)
    {
      std::stringstream sst;
      sst << "msd_" << count << "_" << rst << ".dat";
      msdrecorder->dumpone(sst.str().c_str());
    }
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
    if(!(ncollision&0xffff)) // 65536 per collision
#else
    if(!(ncollision&0xfffff)) // 1048576~10^6 per collision
#endif
    {
      if(nlistflag)
      {
        // determine if update rcut or not
        t2c = steady_clock::now();
        tcollide = t2c - t1c - tnlist;
        double ecratio = tnlist/tcollide;
        // reset counter
        if(ecratio > 1.0)
        {
          double newrcut = rcut*pow(1.05, 1.0/DIM);
          if(newrcut < 2*input.rmax)
          {
            updatercut(newrcut); // increase 5%
#ifdef DEBUG
            std::cout << "ecratio " << ecratio << ", change rcut... " << rcut << std::endl;
#endif
          }
        }
        else if(ecratio < 0.25)
        {
          double newrcut = rcut*pow(0.95, 1.0/DIM);
          if(newrcut > 1.0)
          {
            updatercut(newrcut); // shrink 5%
#ifdef DEBUG
            std::cout << "ecratio " << ecratio << ", change rcut... " << rcut << std::endl;
#endif
          }
        }
#ifdef DEBUG
        else std::cout << "ecratio " << ecratio << std::endl;
#endif
        t1c = t2c;
        tnlist = steady_clock::now() - t2c;
      }
      if(recordstatus)
      {
        t2b = steady_clock::now();
        if(static_cast<duration_type>(t2b-t1b).count()>3600.0)
        {
          // dump temporary status
          fstatus.open("status.txt");
          fstatus << "(" << nowt << " / " << maxt << ") " << ncollision << "\t" << nescape << "\t" << rcut << "\n";
          fstatus.close();
          t1b = t2b;
          // dump temporary msd
          msdrecorder->dumpone("msd.dat");
        }
      }
    }
    lastcollide = thiscollide;
    nowt += t_elapse;
  }
  std::vector<double> vals(6);
  msdrecorder->getfinalmsd( &vals[2], &vals[3] );
  int rtcode = msdrecorder->getinitfinalmsd( &vals[4], &vals[5] );
  if(-1==statuscode)
  {
    vals[0] = nowt;
    vals[1] = shelltol;
    record(vals);
  }
  else
  {
    vals[0] = maxt;
    vals[1] = msdrecorder->getfinalsd();
    statuscode = rtcode;
    record(vals);
  }
#ifndef DEBUG
  if(recordstatus) std::remove("status.txt");
#endif
  // clear neighbor list
  nlists->deconstruct();
  // close trajectory file
  if(input.outmode >= OUTPUT_MODE_DEBUG)
  {
    oftraj.close();
  }
  msdrecorder->dump(r2array, r4array, chidisarray, rcovarray, samplecount, nsampletime);
  msdrecorder->dumplog(logr2array, logr4array, nsampletime);
  delete msdrecorder;
  if(input.rcuttype > 0)
  {
    escrec.dump(count);
  }
  ++count;
  return retval;
}

void Sphere::dump()
{
  std::ofstream ofs(input.datafile+".dat");
  ofs << std::setprecision(10);
  for(int i=0; i<nsampletime && samplecount[i] > 0; ++i)
  {
    ofs << timearray[i] << "\t" << r2array[i]/samplecount[i] << "\t"
        << r4array[i]/samplecount[i] << "\t" << chidisarray[i]/samplecount[i]
        << "\t" << rcovarray[i]/(samplecount[i]*DIM*(DIM-1));
    ofs << "\t" << samplecount[i] << "\n";
  }
  ofs.close();
  if(input.logdatafile.size() > 0)
  {
    ofs.open(input.logdatafile+".dat");
    ofs << std::setprecision(10);
    for(int i=0; i<nsampletime && samplecount[i]>0; ++i)
    {
      ofs << timearray[i] << "\t" << logr2array[i]/samplecount[i]
          << "\t" << logr4array[i]/samplecount[i];
      ofs << "\t" << samplecount[i] << "\n";
    }
    ofs.close();
  }
}

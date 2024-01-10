//
//  box.cpp
//  RLGDynamics
//  Periodic box
//
//  Created by Yi Hu on 6/13/19.
//

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include "box.hpp"
#include "myutility.hpp"

Box::Box(ReadInputDy &input):
input(input), obstacles(input.N), count(0),
scale(pow(input.N/input.phi*MyUtility::getVunitsphere(), -1.0/DIM)),
sqscale(scale*scale), quadscale(sqscale*sqscale), invscale(1./scale), invsqscale(invscale*invscale),
rcut(input.rcut*scale), sqrcut(rcut*rcut), rtol((rcut-scale)*(rcut-scale)), rconstC(rcut*(2*scale-rcut)),
maxt(pow(2.0, input.texp)), maxsinglet(0.5/scale - 1.), nlists(nullptr),
ncollision(0L), nescape(0L)
{
  // validity check
  assert(rtol > 0.);
  assert(maxsinglet > 0);
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
  if(input.rcuttype == 1)
  {
    nlists = new NeighborList<norm_squared_PBC>[input.N];
  }
  else
  {
    // only construct neighbor list w.r.t
    nlists = new NeighborList<norm_squared_PBC>;
  }
  std::cout << "     scale : " << scale << std::endl;
#ifdef GRIDLIST_ON
  grids.update_rcut(rcut);
  if(grids.is_valid())
  {
    // display gridsize=L*L*L
    std::cout << "  gridsize : " << grids.get_size() << std::endl;
    // set maxsinglet into proper value
    maxsinglet = rcut/scale - 1.;
  }
#endif
  if(scale > 0.5)
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
  // ~ 10 brownian collision per record or 10 per \hat{t}
  // std::min(0.1, 1.0/(2 * input.phi)) / sqrt(DIM*1.0)
  br.setup(timearray[0], 0.5/sqrt(DIM*1.0));
}

Box::~Box()
{
  if(input.outmode >= OUTPUT_MODE_PRODUCT)
  {
    ofd.close();
  }
  if(input.rcuttype == 1)
  {
    delete [] nlists;
  }
  else
  {
    delete nlists;
  }
  delete []timearray;
  delete []r2array;
  delete []r4array;
  delete []rcovarray;
  delete []logr2array;
  delete []logr4array;
  delete []chidisarray;
  delete []samplecount;
}

// generate random distribution
void Box::genrandomPoints()
{
  uint32_t i=0, j=0;
  double norm;
  while(i<input.N)
  {
    for(j=0; j<DIM; ++j)
    {
      obstacles[i][j] = MyUtility::rand() - 0.5; // [-0.5, 0.5)
    }
    norm = obstacles[i].norm_squared();
    if(sqscale < norm) ++i; // make a cavity at the origin
  }
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
void Box::genrandomVelocity()
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

Point_handle Box::getnextcollide(Point_handle lastcollide, double &t_elapse)
{
  double tmin = INFINITY, t;
  Point_handle idxmin = NULL_COLLISION, i;
  double A, B, C;
  std::vector<Point_handle> possible_neighbors;
  A = tracer.v.norm_squared();
  Vector vec;
  if(1 == input.rcuttype && lastcollide == NULL_COLLISION)
  {
    if(grids.is_valid())
    {
      grids.get_neighbor(tracer.coord, possible_neighbors);
      for(auto i: possible_neighbors)
      {
        vec = obstacles[i] - tracer.coord;
        pbc(vec);
        C = vec.norm_squared();
        if(sqrcut < C) continue; // longer than the shell radius
        B = vec.dot(tracer.v);
        if(B <= 0.) continue; // away
        B *= -2.0;
        C -= sqscale;
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
      for(i=0; i<input.N; ++i)
      {
        vec = obstacles[i] - tracer.coord;
        pbc(vec);
        C = vec.norm_squared();
        if(0.25 < C) continue; // longer than half box
        B = vec.dot(tracer.v);
        if(B <= 0.) continue; // away
        B *= -2.0;
        C -= sqscale;
        assert(C >= 0.);
        t = MyUtility::rootsolverA(A, B, C);
        if(t < tmin)
        {
          tmin = t;
          idxmin = i;
        }
      }
    }
    if(tmin > maxsinglet)
    {
      tmin = maxsinglet;
      idxmin = NULL_COLLISION;
    }
  }
  else
  {
    // set lastcollide = 0 if input.rcuttype == 2
    if(1 == input.rcuttype)
    {
      if(!nlists[lastcollide].status) // lazy construct
      {
        if(grids.is_valid())
        {
          grids.get_neighbor(obstacles[lastcollide], possible_neighbors);
          nlists[lastcollide].construct(lastcollide, sqrcut, obstacles, possible_neighbors);
        }
        else
        {
          nlists[lastcollide].construct(lastcollide, sqrcut, obstacles);
        }
      }
    }
    else
    {
      lastcollide = 0;
      if(!nlists->status)
      {
        if(grids.is_valid())
        {
          grids.get_neighbor(tracer.coord, possible_neighbors);
          nlists->construct(nlisttracer = tracer.coord, sqrcut, obstacles, possible_neighbors);
        }
        else
        {
          nlists->construct(nlisttracer = tracer.coord, sqrcut, obstacles);
        }
      }
    }
    for(auto i: nlists[lastcollide].neighbors) // duplicate
    {
      vec = obstacles[i] - tracer.coord;
      pbc(vec);
      B = vec.dot(tracer.v);
      if(B <= 0.) continue; // away
      C = vec.norm_squared();
      if(0.25 < C) continue; // longer than half box
      B *= -2.0;
      C -= sqscale;
      assert(C >= 0.);
      t = MyUtility::rootsolverA(A, B, C);
      if(t < tmin)
      {
        tmin = t;
        idxmin = i;
      }
    }
    if(std::isfinite(tmin))
    {
      // check if the collision point is outside the (neighbor shell - 1)
      if(1 == input.rcuttype)
      {
        vec = tracer.coord - obstacles[lastcollide];
      }
      else
      {
        vec = tracer.coord - nlisttracer;
      }
      pbc(vec);
      Vector vecb = vec + tracer.v*tmin;
      if(vecb.norm_squared() > rtol)  // runs outside of the shell
      {
        tmin = INFINITY;
      }
    }
    else
    {
      if(1 == input.rcuttype)
      {
        vec = tracer.coord - obstacles[lastcollide];
      }
      else
      {
        vec = tracer.coord - nlisttracer;
      }
      pbc(vec);
    }
    if(tmin == INFINITY) // noncollide
    {
      B = 2.0*vec.dot(tracer.v);
      if(1 == input.rcuttype)
      {
        C = rconstC;
      }
      else
      {
        C = vec.norm_squared() - rtol;
        nlists->deconstruct(); // neighbor list expired
      }
      assert(C <= 0.);
      tmin = MyUtility::rootsolverB(A, B, C);
      idxmin = NULL_COLLISION;
    }
  }
  t_elapse = tmin;
  return idxmin;
}

void Box::record(std::vector<double> const &vals)
{
  if(input.outmode >= OUTPUT_MODE_PRODUCT)
  {
    if(0 == count)
    {
      ofd.open("dbgsummary.dat");
      ofd << "NCollision\tNescape\tFinalv\tFinalDelta\tFinalMSD\tFinalMQD\n";
    }
    ofd << ncollision << '\t' << nescape << '\t' << tracer.v.norm_squared()/sqscale << std::setprecision(8);
    for(uint32_t rp=0; rp<vals.size(); ++rp) ofd << '\t' << vals[rp];
    ofd << std::endl;
  }
}

void Box::updatetracer(Point_handle thiscollide, double t_elapse)
{
  // update tracer position
  int i;
  for(i=0; i<DIM; ++i)
  {
    tracer.coord[i] += t_elapse*tracer.v[i];
    if(0.5 <= tracer.coord[i])
    {
      tracer.coord[i] -= 1.;
      tracer.pbcvec[i] += 1.;
    }
    else if(tracer.coord[i] < -0.5)
    {
      tracer.coord[i] += 1.;
      tracer.pbcvec[i] -= 1.;
    }
  }

  if(thiscollide == BROWNIAN_COLLISION)
  {
    // reset velocity at brownian time interval. Refer to De Michele's
    // event-driven Brownian dynamics algorithm.
    genrandomVelocity();
  }
  else if(thiscollide != NULL_COLLISION)
  {
    // collide with obstacle
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
  else
  {
    // escape from one neighbor list
    ++nescape;
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

bool Box::samplefinished()
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

int Box::process()
{
  using std::chrono::steady_clock;
  // the obstacle index of last collision. -1 means not hit
  double t_elapse;
  uint64_t nlastcollision = 0;
  Point_handle lastcollide = NULL_COLLISION, thiscollide;
  genrandomPoints();
#ifdef GRIDLIST_ON
  if(grids.is_valid())
  {
    grids.set_points(obstacles);
  }
#endif
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
  ncollision = nescape = 0L;
  if(input.tscaleflag == 4)
  {
    br.reset();
  }

  std::ofstream fstatus;
  bool recordstatus = false;
  clock_type t1b, t2b;
  if(input.repeatrun==1)
  {
    recordstatus = true;
    t1b = steady_clock::now();
  }
  // main loop for simulation. Process one event in each loop.
  while(!samplefinished())
  {
#ifdef DEBUG
    clock_type t1, t2, t3, t4;
    t1 = steady_clock::now();
#endif

    thiscollide = getnextcollide(lastcollide, t_elapse);
    
    if(input.tscaleflag == 4 && br.reached(nowt + t_elapse))
    {
      t_elapse = br.inc() - nowt;
      thiscollide = BROWNIAN_COLLISION;
    }

#ifdef DEBUG
    t2 = steady_clock::now();
#endif
    msdrecorder->record(tracer, nowt, nowt + t_elapse);
#ifdef DEBUG
    t3 = steady_clock::now();
#endif
    updatetracer(thiscollide, t_elapse);
#ifdef DEBUG
    t4 = steady_clock::now();
    stattime += t3-t2;
    collidetime += t2-t1;
    intgtime += t4-t3;
    if(!(ncollision&0xffff) && ncollision != nlastcollision) // 65536 per collision
#else
    if(!(ncollision&0xfffff) && ncollision != nlastcollision) // 1048576~10^6 per collision
#endif
    {
      nlastcollision = ncollision;
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
  // record MSD and other metrics if needed
  std::vector<double> vals(3);
  msdrecorder->getfinalmsd( &vals[1], &vals[2] );
  vals[0] = msdrecorder->getfinalsd();
  record(vals);
#ifndef DEBUG
  if(recordstatus) std::remove("status.txt");
#endif
  // clear neighbor list
  if(1 == input.rcuttype)
  {
    for(uint32_t i=0; i<input.N; ++i)
      nlists[i].deconstruct();
  }
  else
  {
    nlists->deconstruct();
  }
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

void Box::dump()
{
  std::ofstream ofs(input.datafile+".dat");
  ofs << std::setprecision(10);
  for(int i=0; i<nsampletime && samplecount[i]>0; ++i)
  {
    ofs << timearray[i] << "\t" << r2array[i]/samplecount[i]/sqscale
        << "\t" << r4array[i]/samplecount[i]/quadscale << "\t" << chidisarray[i]/samplecount[i]/quadscale
        << "\t" << rcovarray[i]/(samplecount[i] * quadscale * DIM * (DIM-1) ) ;
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

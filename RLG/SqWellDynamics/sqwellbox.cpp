//
//  sqwellbox.cpp
//  SqWellDynamics
//
//  Created by Yi Hu on 5/1/20.
//

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <cstdio>
#include <boost/algorithm/string.hpp>

#include "sqwellbox.hpp"
#include "myutility.hpp"

int ReadInputDySq::read(const char* inputf)
{
  int error = 0;
  std::ifstream infile(inputf);
  std::string outstr;
  if(!infile)
  {
    std::cout << "Can't open " << inputf << " for input." << std::endl;
    error = 2;
    return error;
  }
  else
  {
    std::cout << "Reading input from file " << inputf << std::endl;
  }
  char buf[300],c;
  infile.get(buf,300,'='); infile.get(c); infile >> outstr; // system type
  systype = getSystemType(outstr);
  if(SYSTEM_TYPE_ERROR == systype)
  {
    error = 4;
    return error;
  }
  infile.get(buf,300,'='); infile.get(c); infile >> repeatrun;
  if(systype == SYSTEM_TYPE_BOX || systype == SYSTEM_TYPE_HCOMB)
  {
    infile.get(buf,300,'='); infile.get(c); infile >> N;
  }
  else if(boost::iequals(outstr, "sphere"))
  {
    infile.get(buf,300,'='); infile.get(c); infile >> rmax;
  }
  infile.get(buf,300,'='); infile.get(c); infile >> phi;
  infile.get(buf,300,'='); infile.get(c); infile >> texp;
  infile.get(buf,300,'='); infile.get(c);
  if(!(infile >> sleft >> sinterval >> tscaleflag))
  {
    tscaleflag = 2; // backward compactibility
    infile.clear();
  }
  infile.get(buf,300,'='); infile.get(c); infile >> seed;
  infile.get(buf,300,'='); infile.get(c); infile >> datafile;
  infile.get(buf,300,'='); infile.get(c); infile >> outstr;
  infile.get(buf,300,'='); infile.get(c); infile >> rcuttype;
  infile.get(buf,300,'='); infile.get(c); infile >> rcut;
  infile.get(buf,300,'='); infile.get(c); infile >> xi;
  infile.get(buf,300,'='); infile.get(c); infile >> lambda;
  outmode = getOutputMode(outstr);
  if(infile.eof())
  {
    std::cout << "Error reading input file " << inputf << std::endl;
    error = 3;
  }
  std::cout << "Repeat run : " << repeatrun << "\n";
  if(systype == SYSTEM_TYPE_BOX || systype == SYSTEM_TYPE_HCOMB)
  {
    std::cout << "         N : " << N << "\n";
  }
  else if(boost::iequals(outstr, "sphere"))
  {
    std::cout << "      rmax : " << rmax << "\n";
  }
  std::cout << "       phi : " << phi << "\n";
  std::cout << "       2^t : " << texp << "\n";
  std::cout << "Sample interval (2^t";
  if(tscaleflag == 1)
  {
    std::cout << ")";
  }
  else
  {
    std::cout << " scaled)"; // time now scaled as t/sqrt(dim)
  }
  std::cout << " : " << sleft << " 1/" << sinterval << "\n";
  
  std::cout << "      Seed : " << seed << "\n";
  std::cout << "  datafile : " << datafile << "\n";
  std::cout << "  rcuttype : " << rcuttype << "\n";
  std::cout << "      rcut : " << rcut << "\n";
  std::cout << " xi/lambda : " << xi << " " << lambda << std::endl;
  return error;
}


SqwellBox::SqwellBox(ReadInputDySq &input):
input(input), obstacles(input.N), count(0),
sconstA(pow(1+input.lambda, DIM)-1.0),
scale(pow((input.N/input.phi - (exp(input.xi)-1.0)*sconstA)*MyUtility::getVunitsphere(), -1.0/DIM)),
sqscale(scale*scale), quadscale(sqscale*sqscale), invscale(1./scale), invsqscale(invscale*invscale),
scalepdim(pow(scale, DIM)), rconstA(sconstA*scalepdim),
maxt(pow(2.0, input.texp)), maxsinglet(0.5/scale - 1.), nlists(nullptr),
ncollision(0L), nescape(0L), nlistsize(0L)
{
  updatercut(input.rcut*scale);
  // validity check
  assert(rtol > 0.);
  assert(maxsinglet > 0);
  // adjust maxt if scale mode is 2
  if(input.tscaleflag == 2)
  {
    maxt /= sqrt((double)DIM);
  }
  // neighbor list initialization
  if(input.rcuttype == 1)
  {
    std::cerr << "rcuttype 1 no longer unsupported. Use 2 instead" << std::endl;
    exit(-1);
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
  chidisarray = new double[nsampletime];
  rcovarray = new double[nsampletime];
  samplecount = new int[nsampletime];
  for(int i=0; i<nsampletime; ++i)
  {
    timearray[i] = pow(2., input.sleft + (double)i/input.sinterval);
    if(2==input.tscaleflag) timearray[i]/=sqrt((double)DIM);
    r2array[i] = 0.;
    r4array[i] = 0.;
    chidisarray[i] = 0.;
    rcovarray[i] = 0.;
    samplecount[i] = 0;
  }
  nextrecord = 0;
  // additional information output
}

SqwellBox::~SqwellBox()
{
  if(input.outmode >= OUTPUT_MODE_PRODUCT)
  {
    ofd.close();
  }
  delete nlists;
  delete []timearray;
  delete []r2array;
  delete []r4array;
  delete []chidisarray;
  delete []rcovarray;
  delete []samplecount;
}

// generate random distribution
void SqwellBox::genrandomPoints()
{
  uint32_t i=0, j=0;
  double norm, dtmp;
  double factorA = (exp(input.xi)-1.0)*sconstA;
  double factorB = 1.0/(MyUtility::getVunitsphere()*scalepdim);
  double propA = factorA/(factorA+factorB); // probability of extra obstacle in square well
  for(i=0; i<input.N; ++i)
  {
    if(MyUtility::rand() < propA)
    {
      norm = 0.;
      for(j=0; j<DIM; ++j)
      {
        dtmp = MyUtility::randn();
        obstacles[i][j] = dtmp;
        norm += dtmp*dtmp;
      }
      norm = pow(MyUtility::rand()*rconstA + scalepdim, 1./DIM)/sqrt(norm);
      for(j=0; j<DIM; ++j)
      {
        obstacles[i][j] *= norm;
      }
    }
    else
    {
      do{
        for(j=0; j<DIM; ++j)
        {
          obstacles[i][j] = MyUtility::rand() - 0.5; // [-0.5, 0.5)
        }
        norm = obstacles[i].norm_squared();
        // make a cavity at the origin
      }while(norm < sqscale);
    }
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
void SqwellBox::genrandomVelocity()
{
  int i;
  double dtmp, dsum=0.;
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

Point_handle SqwellBox::getnextcollide(Point_handle lastcollide, double &t_elapse)
{
  double tmin = INFINITY, t;
  Point_handle idxmin = NULL_COLLISION;
  double A, B, C;
  A = tracer.v.norm_squared();
  Vector vec;

  for(auto i: nlists->neighbors) // duplicate
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
    vec = tracer.coord - nlisttracer;
    pbc(vec);
    Vector vecb = vec + tracer.v*tmin;
    if(vecb.norm_squared() > rtol)  // runs outside of the shell
    {
      tmin = INFINITY;
    }
  }
  else
  {
    vec = tracer.coord - nlisttracer;
    pbc(vec);
  }
  if(tmin == INFINITY) // noncollide
  {
    B = 2.0*vec.dot(tracer.v);
    C = vec.norm_squared() - rtol;
    nlists->deconstruct(); // neighbor list expired
    assert(C <= 0.);
    tmin = MyUtility::rootsolverB(A, B, C);
    idxmin = NULL_COLLISION;
  }
  
  t_elapse = tmin;
  return idxmin;
}

void SqwellBox::record(double t_elapse)
{
  int ret_val = 0;
  ret_val = msdrecorder->record(tracer, nowt, nowt + t_elapse);
  if(2 == ret_val && input.outmode >= OUTPUT_MODE_PRODUCT)
  {
    if(0 == count)
    {
      ofd.open("dbgsummary.dat");
      ofd << "NCollision\tNescape\tFinalv\tFinalDelta\tNlistSize\n";
    }
    ofd << ncollision << "\t" << nescape << "\t" << tracer.v.norm_squared()/sqscale << "\t"
    << std::setprecision(8) <<  msdrecorder->getfinalsd()/sqscale << "\t" << nlistsize/nescape << std::endl;
  }
}

void SqwellBox::updatetracer(Point_handle thiscollide, double t_elapse)
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
  if(thiscollide != NULL_COLLISION)
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

void SqwellBox::updatercut(double newrcut)
{
  rcut = newrcut;
  sqrcut = rcut*rcut;
  rtol = (rcut-scale)*(rcut-scale);
  rconstC = rcut*(2*scale-rcut);
  if(nlists) nlists->deconstruct();
#ifdef GRIDLIST_ON
  int retval = grids.update_rcut(rcut);
  if(1==retval)
  {
    grids.set_points(obstacles);
  }
#endif
}

void SqwellBox::updatenlist()
{
  std::vector<Point_handle> possible_neighbors;
  if(nlists->status) nlists->deconstruct();
  if(grids.is_valid())
  {
    grids.get_neighbor(tracer.coord, possible_neighbors);
    nlists->construct(nlisttracer = tracer.coord, sqrcut, obstacles, possible_neighbors);
  }
  else
  {
    nlists->construct(nlisttracer = tracer.coord, sqrcut, obstacles);
  }
  ++nescape;
  nlistsize += nlists->neighbors.size();
}

bool SqwellBox::samplefinished()
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

int SqwellBox::process()
{
  using std::chrono::steady_clock;
  // the obstacle index of last collision. -1 means not hit
  double t_elapse;
  Point_handle lastcollide = NULL_COLLISION, thiscollide;
  genrandomPoints();
#ifdef GRIDLIST_ON
  if(grids.is_valid())
  {
    grids.set_points(obstacles);
  }
#endif
  genrandomVelocity();
  double sleft = input.tscaleflag==1?pow(2.0, input.sleft):pow(2.0, input.sleft)/sqrt((double)DIM);
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
  std::ofstream fstatus;
  bool recordstatus = false;
  clock_type t1b, t2b, t1c, t2c;
  std::chrono::duration<double> tcollide, tnlist;
  if(input.repeatrun==1)
  {
    recordstatus = true;
    t1b = steady_clock::now();
  }
  t1c = steady_clock::now();
  while(!samplefinished())
  {
#ifdef DEBUG
    clock_type t1, t2, t3, t4, t5;
    t1 = steady_clock::now();
#endif
    if(!nlists->status)
    {
      t2c = steady_clock::now();
      updatenlist();
      tnlist += steady_clock::now() - t2c;
    }
#ifdef DEBUG
    t2 = steady_clock::now();
#endif
    thiscollide = getnextcollide(lastcollide, t_elapse);
#ifdef DEBUG
    t3 = steady_clock::now();
#endif
    record(t_elapse);
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
        // determine if update rcut or not
        t2c = steady_clock::now();
        tcollide = t2c - t1c - tnlist;
        double ecratio = tnlist.count()/tcollide.count();
        // reset counter
        if(ecratio > 1.0)
        {
          double newrcut = rcut*pow(1.05, 1.0/DIM);
          if(newrcut < invscale*0.5)
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
  delete msdrecorder;
  ++count;
  return 0;
}

void SqwellBox::dump()
{
  std::ofstream ofs("msd.dat");
  ofs << std::setprecision(10);
  for(int i=0; i<nsampletime && samplecount[i]>0; ++i)
  {
    ofs << timearray[i] << "\t" << r2array[i]/samplecount[i]/sqscale << "\t"
        << r4array[i]/samplecount[i]/quadscale << "\t" << chidisarray[i]/samplecount[i]/quadscale
        << "\t" << rcovarray[i]/(samplecount[i]*quadscale*DIM*(DIM-1));
    ofs << "\t" << samplecount[i] << "\n";
  }
  ofs.close();
}

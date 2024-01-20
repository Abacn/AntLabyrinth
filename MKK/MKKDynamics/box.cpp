/*
 Modified 2022 by Yi Hu, Newtonian dynamics of MKK and MK model

 Updated July 24, 2009 to ~include hardwall boundary condition option~,
 as well as polydisperse spheres.

 Packing of hard spheres via molecular dynamics
 Developed by Monica Skoge, 2006, Princeton University
 Contact: Monica Skoge (mskoge.princeton.edu) with questions
 This code may be used, modified and distributed freely.
 Please cite:

 "Packing Hyperspheres in High-Dimensional Euclidean Spaces"
 M. Skoge, A. Donev, F. H. Stillinger and S. Torquato, 2006

 if you use these codes.
 */

#include "box.h"
#include "utility.h"
#include <iostream>
#include <stdlib.h>
#include <iomanip>

using std::chrono::steady_clock;
// unfortunate bad programming
// need some raw functions
Box* tmpbox;
double normfunc_(const int i, const int j, vector<> &vec)
{
  vector<> minimg;
  tmpbox->getshift(i, j, minimg, vec);
  return minimg.norm_squared();
}

void vectoridx_(vector<>& vec, const int i)
{
  return tmpbox->vectoridx(vec, i);
}

const uint64_t Box::COUNT_LIMIT;

//==============================================================
//==============================================================
//  Class Box: Fills box with hardspheres to given packing fraction
//  and evolves spheres using molecular dynamics!
//==============================================================
//==============================================================


//==============================================================
// Constructor
//==============================================================
Box::Box(read_input const& input):
N(input.N),
boxvolume(1==BOX_TYPE ? 1.0 : 2.0),
maxpf(input.maxpf),
h(N + 1),
neighs(nullptr),
tstart(steady_clock::now()),
tnlist(duration_type::zero()),
calcpmsd(input.calcpmsd),
create_counter(0),
pmsd_counter(0)
{
  // intialize parameter values
  rx = diameterRange(input.polydispersity);
  rscale = gerRfrompfrac(input.initialpf * boxvolume, rx, N);
  rinitial = rscale;
  rfinal = gerRfrompfrac(input.maxpf * boxvolume, rx, N);
  rmeanfin = rfinal * (1.0 - rx * rx);
  growthrate = input.growthrate;
  grsq = rinitial * rinitial * growthrate * growthrate;
  shiftscale = input.shiftscale;

  // claim memory
  s = new Sphere[N];
  x0 = new vector<>[N];
  h.s = s;

  // set zero to status values
  gtime = 0.;
  rtime = 0.;
  ncollisions = 0;
  ntransfers = 0;
  nchecks = 0;
  toteventscheck = toteventscoll = toteventstransfer = 0;
  ncycles = 0;
  xmomentum = 0.;
  collisionrate = 0.;
  nlastbuildnlist = 0;
  nextsampletime = sampletimedelt = 1e20; // set to very large value to dismiss sample before measure switch on

  if (calcpmsd)
  {
    sums_pmsd = new double[N];
    sum_pmqd = 0.;
    for(int i = 0; i < N; i++)
      sums_pmsd[i] = 0.0;
    if (calcpmsd == 2)
    {
      dist_sd = new uint64_t[sd_nbin];
      for(int i = 0; i < sd_nbin; i++)
      dist_sd[i] = 0;
    }
  }
  else
  {
    sums_pmsd = nullptr;
  }
  std::cout << "   rx : " << rx << '\n';
  std::cout << "   rinitial : " << rinitial << '\n';
  std::cout << "   rfinal : " << rfinal << '\n';

  neighs = new nlist(N);
  double setrv = 2 * rfinal * (1 + rx);
  if (setrv > sqrt(0.5))
  {
    std::cout << "Error: final sphere radius " << setrv << " too small. Cannot construct neighbor list." << std::endl;
    exit(1);
  }
  if (input.rfactor <= 1.0)
  {
    std::cout << "Error: invalid neighor list shell size: " << input.rfactor << std::endl;
    exit(1);
  }
  setrv *= input.rfactor;
  if (setrv > sqrt(0.5))
  {
    std::cout << "Attempted rfactor too large, changed to " << input.rfactor * sqrt(0.5) / setrv << std::endl;
    setrv = sqrt(0.5);
  }

  if(input.rfactor > 1.0)
  {
    // use neighbor list
    neighs->setrv(setrv, rfinal * (1 + rx));
  }
  else
  {
    std::cout << "Error: invalid neighor list shell size: " << input.rfactor << std::endl;
    exit(-1);
  }

  allshifts = new vector<> * [N - 1];
  for (int i = 0; i < N; i++)
  {
    if (i + 1 < N)
      allshifts[i] = new vector<>[N - i - 1];
  }
  tmpbox = this;
}

//==============================================================
// Destructor
//==============================================================
Box::~Box()
{
  delete[] s;
  delete[] x0;
  if (calcpmsd)
  {
    delete [] sums_pmsd;
    if (calcpmsd == 2)
    {
      delete [] dist_sd;
    }
  }
  for (int i=0; i<N-1; i++)
  {
    delete [] allshifts[i];
  }
  delete [] allshifts;
  if (neighs != nullptr) delete neighs;
}


//==============================================================
// ReadFile
//==============================================================
void Box::ReadPositions(const std::string filename, const std::string shiftname/*=""*/)
{
  if(shiftscale != 0.0)
  {
    std::ifstream inshift(shiftname, std::ios::binary);
    if(!inshift.is_open())
    {
      std::cout << "Error: cannot read from " << shiftname << "." << std::endl;
      exit(-1);
    }
    for(int rp=0; rp<N-1; ++rp)
    {
      for(int rq=0; rq<N-rp-1; ++rq)
      {
        allshifts[rp][rq].read(inshift);
      }
    }
  }
  else
  {
    // set to zero
    for(int rp=0; rp<N-1; ++rp)
    {
      for(int rq=0; rq<N-rp-1; ++rq)
      {
        for(int rr=0; rr<DIM; ++rr)
        {
          allshifts[rp][rq].x[rr] = 0.0;
        }
      }
    }
  }
  // open file to read in arrays
  std::ifstream infile(filename);
  if (!infile)
  {
    std::cout << "error, can't open " << filename << std::endl;
    exit(-1);
  }
  else
  {
    int dim;
    // Two format is supported:
    // (i) Format described in https://cims.nyu.edu/~donev/Packing/C++/
    // (ii) raw coordinates
    // first check what kind of file it is
    std::string line;
    std::getline(infile, line);
    if (line.size() > 3)
    {
      // read raw r
      infile.seekg(0, std::ios::beg);

      for (int i = 0; i < N; i++)
      {
        s[i].m = 1;
        for (int k = 0; k < DIM; k++)
        {
          infile >> s[i].x[k]; // read in position
        }
        if (rx > 0.0)
        {
          infile >> s[i].r; // read r
        }
        else
        {
          s[i].r = 1.0; // monodisperse, just assign as rinitial
          infile.ignore(256, '\n');
        }
        s[i].i = i;
        s[i].lutime = gtime;
        s[i].species = 1;
      }
      return;
    }
    else
    {
      dim = atoi(line.c_str());
      infile.ignore(256, '\n');
    }

    if (dim != DIM)  // quit if dimensions don't match
    {
      std::cout << "error, dimensions don't match" << std::endl;
      exit(-1);
    }
    int tempN;
    infile >> tempN; infile.ignore(256, '\n'); // ignore the N1 N2 line
    if (tempN != N)
    {
      std::cout << "error, system sizes don't match" << std::endl;
      exit(-1);
    }
    infile.ignore(256, '\n'); // N
    double tmprx;
    infile >> rscale >> tmprx; infile.ignore(256, '\n'); // r rx
    if (abs(tmprx - rx) > 1e-10)
    {
      std::cout << "error, polydispersity don't match" << std::endl;
      exit(-1);
    }
    infile.ignore(256, '\n'); // 1 0 0 ...
    infile.ignore(256, '\n'); // T T T
  }
  double r_mean=0;


  for (int i=0; i<N; i++)
  {
    s[i].m=1;
    for (int k=0; k<DIM; k++)
    {
      infile >> s[i].x[k]; // read in position
    }
    infile >> s[i].r;      // read in radius
    //s[i].gr=growthrate*s[i].r;
    r_mean = r_mean + s[i].r;
    //s[Ncurrent] = Sphere(Ncurrent, xrand, gtime, radius, growth_rate, mass, species);
    s[i].i = i;
    s[i].lutime = gtime;
    s[i].species = 1; //
  }
  r_mean = r_mean/N;
  if(fabs(r_mean - (1-rx*rx)) > 0.01)
  {
    std::cout << "error: the global radius does not equal to the mean radius." << std::endl;
    exit(-1);
  }
  //for (int i=0; i<N; i++)
  //  s[i].gr=growthrate*s[i].r;
  infile.close();
}

//==============================================================
// ReadVelocities
//==============================================================
int Box::ReadVelocities(const std::string filename)
{
  // open file to read in arrays
  std::ifstream infile(filename);
  if (!infile.is_open()) return -1;
  infile.ignore(256, '\n');  // ignore the dim time line
  infile.ignore(256, '\n');  // ignore the #sphere 1 line
  infile.ignore(256, '\n');  // ignore the #sphere line
  infile.ignore(256, '\n');  // ignore the diameter line
  infile.ignore(1000, '\n'); // ignore the 100 010 001 line
  infile.ignore(256, '\n');  // ignore the T T T line

  for (int i=0; i<N; i++)
  for (int k=0; k<DIM; k++)
  infile >> s[i].v[k];

  infile.close();
  return 0;
}

void Box::CopyPositions(const Box &b, double temp)
{
  for (int i=0; i<N; i++)
  {
    for (int k=0; k<DIM; k++)
    {
      s[i].r = b.s[i].r;
      //s[i].gr = b.s[i].gr;
      s[i].species = b.s[i].species;
      s[i].m = b.s[i].m;
    }
  }
  VelocityGiver(temp);
  SetInitialEvents();
}

//==============================================================
// Creates pair shift
//==============================================================
void Box::CreateShift()
{
  if (shiftscale == 0.0)
  {
    // set to zero
    for (int rp = 0; rp < N - 1; ++rp)
    {
      for (int rq = 0; rq < N - rp - 1; ++rq)
      {
        for (int rr = 0; rr < DIM; ++rr)
        {
          allshifts[rp][rq].x[rr] = 0.0;
        }
      }
    }
  }
  else if(shiftscale > 0)
  {
    // shifts r are in final r mean scale
    // Final radius set the length scale
    // variance = shiftscale^2 / DIM
    std::normal_distribution<double> gsdis(0.0, shiftscale / sqrt((double)DIM));
    for (int rp = 0; rp < N - 1; ++rp)
    {
      for (int rq = 0; rq < N - rp - 1; ++rq)
      {
        for (int rr = 0; rr < DIM; ++rr)
        {
          allshifts[rp][rq].x[rr] = gsdis(rg)*(2*rmeanfin); // length scale: 2 r
        }
      }
    }
  }
  else
  {
    // MK model -- flat shift. particle coordinate created first
    vector<> result;
    for (int rp = 0; rp < N - 1; ++rp) // particle rp
    {
      for (int rq = 0; rq < N - rp - 1; ++rq) // particle rp+rq+1
      {
        int idxj = rp + rq + 1;
        double contactdis = rscale * rscale * (s[rp].r + s[idxj].r) * (s[rp].r + s[idxj].r);
        do
        {
          for (int rr = 0; rr < DIM; ++rr) allshifts[rp][rq].x[rr] = 2 * rrand(rg) - 1.0;
          setminimg(allshifts[rp][rq].x);
          getshift(rp, rp + rq + 1, result);
        } while (result.norm_squared() <= contactdis);
      }
    }
  }
}


//==============================================================
// Recreates all N spheres at random positions
//==============================================================
void Box::RecreateSpheres(const std::string fileprefix, double temp)
{
  if(shiftscale == 0.0)
    ReadPositions(fileprefix + "_config.dat");  // reads in positions of spheres
  else
    ReadPositions(fileprefix + "_config.dat", fileprefix + "_shifts.dat");
  if(ReadVelocities(fileprefix + "_v.dat") != 0)
  {
    std::cout << "Cannot read from " << fileprefix + "_v.dat" << ". generate new" << std::endl;
    VelocityGiver(temp);      // gives spheres initial velocities
  }
  SetInitialEvents();
}


//==============================================================
// Creates all N spheres at random positions
//==============================================================
void Box::CreateSpheres(double temp)
{
  double rtom2min = 1.0 / ((1.0 + rx) * (1.0 + rx));
  double rtom2step = (1.0 / ((1.0 - rx) * (1.0 - rx)) - 1.0 / ((1.0 + rx) * (1.0 + rx))) / (N - 1);
  if (shiftscale < 0)
  {
    // create spheres first, assign shift later
    for (int i = 0; i < N; i++)
    {
      CreateSphere(i, 1.0 / sqrt(rtom2min + i * rtom2step));
    }
    CreateShift();
  }
  else
  {
    int Ncurrent = 0;
    CreateShift();
    create_counter = 0;
    for (int i = 0; i < N; i++)
    {
      CreateSphere(Ncurrent, 1.0 / sqrt(rtom2min + i * rtom2step));
      Ncurrent++;
    }
    if (Ncurrent != N)
    {
      std::cout << "problem! only made " << Ncurrent << " out of " << N << " desired spheres" << std::endl;
      exit(1);
    }
  }
  VelocityGiver(temp);
  SetInitialEvents();
}


//==============================================================
// Creates a sphere of random radius r at a random unoccupied position
//==============================================================
void Box::CreateSphere(int Ncurrent, double givenR/* = 0.*/)
{
  int keeper;    // boolean variable: 1 means ok, 0 means sphere already there
  vector<> xrand;  // random new position vector
  double d = 0.;
  double radius;
  double rsqrscale = rscale * rscale;
  double growth_rate;
  double mass;
  int species;

  if (givenR == 0.)
  {
    // generate random radius between r*(1-x) and r*(1+x)
    double rtom2 = 1.0/((1.0+rx)*(1.0 + rx))
    + rrand(rg)*(1.0/((1.0-rx)*(1.0-rx)) - 1.0 /((1.0 + rx)*(1.0 + rx)));  // r^(-2)
    radius = 1.0 / sqrt(rtom2);
  }
  else
  {
    radius = givenR;
  }
  growth_rate = growthrate;
  mass = 1.;
  species = 1;

  while (create_counter < COUNT_LIMIT * N)
  {
    keeper = 1;

    int32_t itmp[DIM];
    for(int k=0; k<DIM; k++) xrand[k] = 2*rrand(rg) - 1;
    DnLatticePoint(xrand.x, itmp);

    if (shiftscale >= 0)
    {
      // check overlap in finite-range MKK (shift created first)
      for (int i = 0; i < Ncurrent; i++)  // check if overlapping other spheres
      {
        vector<> dx = s[i].x - xrand + allshifts[i][Ncurrent - i - 1]; // i<j
        setminimg(dx.x);
        d = dx.norm_squared();
        if (d <= rsqrscale * (radius + s[i].r) * (radius + s[i].r)) // overlapping!
        {
          keeper = 0;
          create_counter++;
          break;
        }
      }
    }

    if (keeper == 1)
      break;
  }

  if (create_counter >= COUNT_LIMIT * N)
  {
    std::cout << "counter >= " << COUNT_LIMIT * N << std::endl;
    exit(-1);
  }
  s[Ncurrent] = Sphere(Ncurrent, xrand, gtime, radius, growth_rate,
                       mass, species);
  if (shiftscale >= 0)
  {
    if ((Ncurrent + 1) % 10 == 0) std::cout << '.' << std::flush;
    if ((Ncurrent + 1) % 500 == 0) std::cout << std::endl;
  }
}

//==============================================================
// Velocity Giver, assigns initial velocities from Max/Boltz dist.
// Modification by Patrick Charbonneau: make sure that the
// center of mass velocity is 0,
// and that the total kinetic energy is equal to the corresponding T.
//==============================================================
void Box::VelocityGiver(double T)
{
  double sumv[DIM],sumv2,fs;
  int i,k;
  for(k=0;k<DIM;k++)
  {
    sumv[k]=0;
  }

  for (i=0; i<N; i++)
  {
    for (k=0; k<DIM; k++)
    {
      if (T==0.)
      {
        s[i].v[k] = 0.;
      }
      else
      {
        s[i].v[k] = Velocity(T);
      }
      sumv[k]+=s[i].v[k];
    }
  }

  //Eliminate the center of mass displacement
  sumv2=0;
  for(k=0;k<DIM;k++)
  {
    sumv[k]/=N;
    for(i=0;i<N;i++)
    {
      s[i].v[k]-=sumv[k];
      sumv2+=s[i].v[k]*s[i].v[k];
    }
  }

  //Rescale the velocities to the set temperature
  fs=sqrt(DIM*T*(N-1)/sumv2);
  for(k=0;k<DIM;k++)
  {
    for(i=0;i<N;i++)
    {
      s[i].v[k]*=fs;
    }
  }
}


//==============================================================
// Velocity, gives a single velocity from Max/Boltz dist.
//==============================================================
double Box::Velocity(double T)
{
  double rand;                       // random number between -0.5 and 0.5
  double sigmasquared = T;    // Assumes M = mass of sphere = 1
  double sigma = sqrt(sigmasquared); // variance of Gaussian
  double stepsize = 1000.;           // stepsize for discretization of integral
  double vel = 0.0;                  // velocity
  double dv=sigma/stepsize;
  double p=0.0;

  rand = rrand(rg) - 0.5;
  if(rand < 0)
  {
    rand = -rand;
    dv = -dv;
  }

  while(fabs(p) < rand) // integrate until the integral equals rand
  {
    p += dv * 0.39894228 * exp(-vel*vel/(2.*sigmasquared))/sigma;
    vel += dv;
  }
  return vel;
}


//==============================================================
// Finds next events for all spheres..do this once at beginning
//==============================================================
void Box::SetInitialEvents()
{
  h.clear();
  for (int i = 0; i < N; i++)  // set all events to checks
  {
    Event e(gtime, i, INF);
    s[i].nextevent = e;
    h.insert(i);
  }
  neighs->construct(vectoridx_,  normfunc_);
}


//==============================================================
// Finds next event for sphere i
//==============================================================
Event Box::FindNextEvent(int i)
{
  Event t = FindNextTransfer(i);
  Event c = FindNextCollision(i);
  if ((c.time < t.time)&&(c.j == INF)) // next event is check at DBL infinity
    return c;
  else if (c.time < t.time) // next event is collision!
  {
    CollisionChecker(c);
    return c;
  }
  else // next event is transfer!
    return t;
}


//==============================================================
// Checks events of predicted collision partner to keep collisions
// symmetric
//==============================================================
void Box::CollisionChecker(Event c)
{
  int i = c.i;
  int j = c.j;
  Event cj(c.time,j,i,-c.shift);

  // j should have NO event before collision with i!
  if (!(c.time  < s[j].nextevent.time))
    std::cout << i << " " <<  j << " error collchecker, s[j].nextevent.time= " << s[j].nextevent.time << " " << s[j].nextevent.j << ", c.time= " << c.time << std::endl;

  int k = s[j].nextevent.j;
  if ((k < N) && (k!=i)) // j's next event was collision so give k a check
    s[k].nextevent.j = INF;

  // give collision cj to j
  s[j].nextevent = cj;
  h.upheap(h.index[j]);
}


//==============================================================
// Find next transfer for sphere i
//==============================================================
Event Box::FindNextTransfer(int i)
{
  double ttime = DBL_LARGE;

  vector<> xi = s[i].x + s[i].v*(gtime - s[i].lutime);
  vector<> vi = s[i].v;

  vector<> dx = s[i].x - neighs->getstartpos(i);
  setminimg(dx.x);
  double A = s[i].v.norm_squared(),
  B = vector<>::dot(dx, s[i].v),
  C = dx.norm_squared() - neighs->getsqrhfsh();
  double delt = B*B - A * C;
  if (delt < 0.0 || C > 1e-10)
  {
    // already escaped!
    //#ifdef DEBUG
    //      std::cout << i << "escaped (delt, dr)" << delt << " " << -C << std::endl;
    //#endif
    ttime = 0;
  }
  else
  {
    ttime = (sqrt(delt) - B) / A;
  }
  if (ttime < 0)
    ttime = 0;
  // make the event and return it
  Event e = Event(ttime + gtime, i, INF-1);
  return e;
}

// neighbor list version
void Box::ForAllNeighbors(int i, Neighbor& operation)
{
  auto p = neighs->getnlist(i);
  for (int rp = 0; rp < p.first; ++rp)
  {
    operation.Operation(p.second[rp], neighs->getshift(i, rp));
  }
}


//==============================================================
// PredictCollision
//==============================================================
void Box::PredictCollision(int i, int j, const vector<> &shift, collision* ccollision)
{
  double ctimej;

  if (i != j)
  {
    // calculate updated position and velocity of i and j
    vector<> dv = s[i].v - s[j].v;
    vector<> dx = s[i].x - s[j].x - shift + s[i].v * (gtime - s[i].lutime) - s[j].v * (gtime - s[j].lutime);
    double rxsum = s[i].r + s[j].r;
    double r_now = rxsum * (rscale + rinitial * gtime * growthrate);

    double A, B, C;
    A = dv.norm_squared() - grsq * rxsum * rxsum;
    B = vector<>::dot(dx, dv) - r_now * rinitial * growthrate * rxsum;
    C =dx.norm_squared() - r_now * r_now;

    ctimej = rootsolverA(A, B, C) + gtime;

    if (isnan(ctimej))
    {
      std::cout << "error, " << i << " and " << j <<
        " are overlapping at time " << gtime << std::endl;
      std::cout << "A, B, C = " << A << " " << " " << B <<
        " " << " " << C << std::endl;
      std::cout << "unknown error" << std::endl;
      // WriteLastConfiguration("config_onerror.dat");
      PrintStatistics();
    }

    if (ctimej < gtime)
      std::cout << "error in find collision ctimej < 0" << std::endl;

    if ((ctimej < ccollision->ctime) && (ctimej < s[j].nextevent.time))
    {
      ccollision->ctime = ctimej;
      ccollision->cpartner = j;
      ccollision->cpartnerpboffset = shift;
    }
  }
}

// solve for the positive root in A*x^2 + 2*B*x + C = 0 assumes C >= 0
double Box::rootsolverA(double a, double b, double c)
{
  double x = DBL_LARGE;

  if (c < -1E-12)
  {
    return NAN;
  }
  double det = b * b - a * c;
  if (det > -10. * DBL_EPSILON)
  {
    if (det < 0.)  // determinant can be very small for double roots
      det = 0.;
    if (b < 0.)
      x = c / (-b + sqrt(det));
    else if ((a < 0.) && (b > 0.))
      x = -(b + sqrt(det)) / a;
    else
      x = DBL_LARGE;
  }
  return x;
}

//==============================================================
// Find next collision
//==============================================================
Event Box::FindNextCollision(int i)
{
  collision cc(i, this);

  ForAllNeighbors(i, cc);
  Event e;
  if (cc.cpartner == i)  // found no collisions in neighboring cells
  {
    if (cc.ctime != DBL_LARGE)
      std::cout << "ctime != DBL_LARGE" << std::endl;
    e = Event(DBL_LARGE,i,INF);  // give check at double INF
  }
  else
    e = Event(cc.ctime,i,cc.cpartner,cc.cpartnerpboffset);

  return e;
}


//==============================================================
// Returns first event
//==============================================================
int Box::ProcessEvent()
{
  // Extract first event from heap
  int i = h.extractmax();
  Event e = s[i].nextevent; // current event
  Event f;                  // replacement event
  int retval = -1;

  Statistics(e.time);

  if ((e.j>=0)&&(e.j<N))  // collision!
  {
    retval = 0;
    ncollisions++;
    //std::cout << "collision between " << e.i << " and " << e.j << " at time " << e.time << std::endl;
    Collision(e);
    f = FindNextEvent(i);
    s[i].nextevent = f;
    h.downheap(1);
    if (f.time < e.time)
    {
      std::cout << "error, replacing event with < time" << std::endl;
      exit(-1);
    }

    // make sure collision was symmetric and give j a check
    if ((s[e.j].nextevent.j != i)||(s[e.j].nextevent.time != gtime))
    {
      std::cout << "error collisions not symmetric" << std::endl;
      std::cout << "collision between " << e.i << " and " << e.j << " at time " << e.time << std::endl;
      std::cout << "but " << e.j << " thinks it has " << s[e.j].nextevent.j<< " "  << s[e.j].nextevent.time << std::endl;
      exit(-1);
    }
    else  // give j a check
      s[e.j].nextevent.j = INF;
  }
  else if (e.j==INF)      // check!
  {
    retval = 1;
    nchecks++;
    //std::cout << "check for " << e.i << " at time " << e.time << std::endl;
    f = FindNextEvent(i);
    s[i].nextevent = f;
    h.downheap(1);
  }
  else if (e.j==INF-1)
  {
    clock_type event_start = steady_clock::now();
    retval = 2;
    // for neighbor list: sphere escaped from neighbor list possible, reconstruct neighborlist!
    gtime = e.time;
    Synchronize(false);

    // check if need bump/shrink nlist rv
    if (nlastbuildnlist > 99)
    {
      nlastbuildnlist = 0;
      duration_type t_all = steady_clock::now() - tlastresize;
      // fradction of time used to construct neighbor list
      double nlist_fraction = tnlist.count() / t_all.count();
      // std::cout << "nlist fraction: " << tnlist.count() << '/' << t_all.count() << std::endl;
      if (nlist_fraction > 0.65)
      {
        neighs->tryadjustrv(1);
      }
      else if (nlist_fraction < 0.3)
      {
        neighs->tryadjustrv(2);
      }
    }
    if (nlastbuildnlist == 0)
    {
      tlastresize = event_start;
      tnlist = duration_type::zero();
    }

    // then set all event to check
    SetInitialEvents();
    ntransfers++;
    nlastbuildnlist++;

    tnlist += steady_clock::now() - event_start;
  }
  else
  {
    std::cout << "Unknown event!" << e.j << std::endl;
    exit(-1);
  }

  return retval;
}


//==============================================================
// Processes a collision
//=============================================================
void Box::Collision(Event e)
{
  double ctime = e.time;
  int i = e.i;
  int j = e.j;
  gtime = ctime;

  // Update positions and cells of i and j to ctime
  s[i].x += s[i].v*(gtime-s[i].lutime);
  s[j].x += s[j].v*(gtime-s[j].lutime);

  // Check to see if a diameter apart
  double r_sum = (s[i].r + s[j].r)*(rscale + rinitial*growthrate*gtime);

  // make unit vector out of displacement vector
  vector<> dhat = s[i].x - s[j].x - e.shift;
  double dhatsqr = dhat.norm_squared();
  double distance = dhatsqr - r_sum*r_sum;
  if (distance * distance > 10. * DBL_EPSILON)
  {
    std::cout << "box::Collision\n";
    std::cout << s[i].x << "\n";
    std::cout << s[j].x << "\n";
    std::cout << i << " " << j << " overlap " << distance << std::endl;
  }
  s[i].lutime = gtime;
  s[j].lutime = gtime;

  double dhatmagnitude = sqrt(dhatsqr);
  dhat /= dhatmagnitude; // normalize

  vector<> vipar = dhat*vector<>::dot(s[i].v, dhat);  // parallel comp. vi
  vector<> vjpar = dhat*vector<>::dot(s[j].v, dhat);  // parallel comp. vj
  vector<> viperp = s[i].v - vipar;                   // perpendicular comp. vi
  vector<> vjperp = s[j].v - vjpar;                   // perpendicular comp. vj

  // half of reletive velocity due to growthrate
  double irgrate = rinitial * s[i].r * growthrate;
  double jrgrate = rinitial * s[j].r * growthrate;

  // particles has same mass ?
  s[i].v = vjpar + dhat*jrgrate + viperp;
  s[j].v = vipar - dhat*irgrate + vjperp;

  // momentum exchange
  double xvelocity;   // exchanged velocity
  xvelocity = vector<DIM>::dot(s[i].v - s[j].v, dhat) - (irgrate + jrgrate);
  xmomentum += xvelocity*dhatmagnitude*s[i].m*s[j].m*2/(s[i].m+s[j].m);
}


//==============================================================
// Output event heap...purely used for debugging
//==============================================================
void Box::OutputEvents()
{
  h.print();
}

//==============================================================
// Computes the total energy
//==============================================================
double Box::Energy()
{
  double E=0;
  for (int i=0; i<N; i++)
  E += s[i].m*s[i].v.norm_squared();
  // YI: changed N-1 to N
  return 0.5*E/N;
}

//==============================================================
// Calliborate velociity such that bulk momentum is zero
//==============================================================
void Box::CallibVelocity()
{
  vector<> vtmp;
  for(int rp=0; rp<N; ++rp)
  {
    vtmp += s[rp].v;
  }
  vtmp /= N;
  for(int rp=0; rp<N; ++rp)
  {
    s[rp].v -= vtmp;
  }
  energy = Energy();
  double vavg = sqrt(2.*M*energy);

  for (int i=0; i<N; i++)
  s[i].v /= vavg/sqrt(DIM);
}

//==============================================================
// Computes the pressure
// return pressure and energy
//==============================================================
std::pair<double, double> Box::Pressure()
{
  double E = Energy();
  double tnow = rtime + gtime;
  double pressure = 1 + xmomentum / (2. * E * N * (tnow-plasttime));
  plasttime = tnow;
  xmomentum = 0;
  return { pressure, E };
}

//==============================================================
// Calculates the packing fraction
//==============================================================
double Box::PackingFraction()
{
  // lazily compute packing fraction
  double rfactor = 0.;
  for (int i=0; i<N; i++)
  {
    rfactor += pow(s[i].r*(rscale + rinitial*gtime*growthrate), DIM);
  }
  double v = rfactor*VOLUMESPHERE;
  return pf = v / boxvolume;
}

//==============================================================
// Calculates the optimal ngrids, assuming ngrids is not updated
// automatically and is very conservative
//==============================================================
//int box::Optimalngrids()
//{
//  double maxr;

//  maxr = pow(exp(lgamma(1.+((double)(DIM))/2.))*maxpf/
//       (N*(bidispersityfraction + (1.-bidispersityfraction)*
//     pow(bidispersityratio, DIM))),
//       1./DIM)/sqrt(PI);

//  return (int)(1./(2.*maxr));
//  return 1;
//}


//==============================================================
// Processes n events
// option: 0: increment for every event; 1: increment only for
// collision
// nextt: if not 0, process event until time
//==============================================================
void Box::Process(int n, int option, double nextt)
{
  int retval;
  double deltat = gtime;
  if (nextt == 0.) {
    int i=0;
    while (i<n)
    {
      retval = ProcessEvent();
      if (!option || 0 == retval) ++i; // is collision
    }
  } else {
    while((gtime + rtime)/(rmeanfin* 2.0) < nextt) {
      retval = ProcessEvent();
    }
  }
  deltat = gtime - deltat;
  double oldenergy = energy;
  energy = Energy();        // kinetic energy

  energychange = ((oldenergy - energy)/oldenergy)*100; // percent change in energy

  if (deltat != 0.)
  {
    collisionrate = ((double)(ncollisions))/deltat;
  }

  // reset to 0
  toteventscheck += nchecks;
  toteventscoll += ncollisions;
  toteventstransfer += ntransfers;
  ncollisions = 0;
  ntransfers = 0;
  nchecks = 0;
  ncycles++;
}

//==============================================================
// Maybe Run statistics
//==============================================================
void Box::Statistics(double ctime)
{
  if (calcpmsd)
  {
    double currentt = ctime + rtime;
    if(currentt >= nextsampletime)
    {
      double deltt = nextsampletime - rtime;
      double msd, mqd = 0.;
      double scalefactor = rmeanfin*rmeanfin*4;
      for (int i=0; i<N; ++i)
      {
        msd = vector<>::norm_squared(s[i].x + s[i].v*(deltt-s[i].lutime), x0[i]);
        if (calcpmsd == 2)
        {
          // collect distribution of square displacements
          int use_n = (int)((log2(msd * DIM / scalefactor) - sd_left) * sd_delta) + 1;
          if (use_n < 0) use_n = 0;
          else if (use_n >= sd_nbin) use_n = sd_nbin - 1;
          ++dist_sd[use_n];
        }
        sums_pmsd[i] += msd;
        sum_pmqd += msd*msd;
      }
      pmsd_counter++;
      do
      {
        nextsampletime += sampletimedelt;
      }
      while (nextsampletime < currentt);
    }
  }
}

//==============================================================
// Prints statistics for n events
// mode: 0 - print final and instant statistics. Useful for debug
//       1 - print final statistics
//==============================================================
void Box::PrintStatistics(int mode/*=0*/)
{
  if (calcpmsd && pmsd_counter>10)
  {
    // summary: MSD, MQD, particle MSD
    std::ofstream ofs("dbgsummary.dat");
    ofs.precision(10);
    double scalefactor = rmeanfin*rmeanfin*4;
    // MSD, MQD, <PMSD^2>
    double pmsd = 0., pmqd = sum_pmqd / (N*pmsd_counter*scalefactor*scalefactor), psmsd=0.;
    for(int i=0; i<N; ++i)
    {
      pmsd += sums_pmsd[i];
      psmsd += sums_pmsd[i]*sums_pmsd[i];
    }
    pmsd /= scalefactor*N*pmsd_counter;
    psmsd /= N*scalefactor*scalefactor*pmsd_counter*pmsd_counter;
    ofs << pmsd << '\t' << pmqd << '\t' << psmsd << '\n';
    ofs.close();
    if (calcpmsd == 2)
    {
      // distribution of SD
      std::ofstream ofs2("sddist.dat");
      double cutoff_left, cutoff_right, p_deltad;
      for (int i=0; i<sd_nbin; ++i)
      {
        // right cutoff of distribution point, except for the last one, which is actually inf
        if (i == 0) cutoff_left = 0.;
        else cutoff_left = cutoff_right;
        cutoff_right = exp2(sd_left + ((double) i) / sd_delta);
        p_deltad = (double) dist_sd[i] / (N*pmsd_counter) / (cutoff_right - cutoff_left);
        ofs2 << cutoff_right << '\t' << p_deltad << '\n';
      }
      ofs2.close();
      // distribution of cage MSD
      std::ofstream ofs3("psddist.dat");
      for(int i=0; i<N; ++i)
      {
        ofs3 << DIM * sums_pmsd[i] / (scalefactor * pmsd_counter) << '\n';
      }
      ofs3.close();
    }
  }
  auto pressurepair = Pressure();
  std::cout << "packing fraction = " << PackingFraction() << std::endl;
  if(0==mode)
  {
    std::cout << "gtime = " << gtime << std::endl;
    std::cout << "processing # events = " << ncollisions+ntransfers+nchecks << ", # collisions = " << ncollisions << ", # transfers = " << ntransfers << ", # checks = " << nchecks << std::endl;
    std::cout << "growthrate = " << growthrate << std::endl;
    std::cout << "reduced pressure = " << pressurepair.first << std::endl;
  }
  std::cout << "total time = " << rtime+gtime << std::endl;
  std::cout << "kinetic energy = " << pressurepair.second << std::endl;
  std::cout << "total # events = " << toteventscheck+toteventscoll+toteventstransfer << ", # collisions = " << toteventscoll << ", # transfers = " << toteventstransfer << ", # checks = " << toteventscheck << std::endl;
  std::cout << "collisionrate = " << collisionrate << std::endl;
  std::cout << "-----------------" << std::endl;
}


//==============================================================
// Updates spheres to gtime, reset gtime to zero, synchronizes, and can change growth rate
//==============================================================
void Box::Synchronize(bool rescale)
{
  for (int i=0; i<N; i++)
  {
    s[i].x = s[i].x + s[i].v*(gtime-s[i].lutime);
    s[i].nextevent.time -= gtime;

    if (s[i].nextevent.time < 0.)
      std::cout << "error, event times negative after synchronization" << std::endl;
    if (rescale == true)   // give everyone checks
    {
      s[i].nextevent = Event(0., i, INF);
    }
    s[i].lutime = 0.;
  }
  rscale += rinitial*gtime*growthrate;       // r defined at gtime = 0
  rtime += gtime;
  gtime = 0.;

  // calliborate velocity
  CallibVelocity();
}


//reset growthrate and time to zero
void Box::Reset()
{
  CallibVelocity();

  for (int i=0; i<N; i++)
  {
    s[i].x = s[i].x + s[i].v*(gtime-s[i].lutime);
    // move back to box
    setminimg(s[i].x.x);
    x0[i] = s[i].x;
    s[i].nextevent.time -= gtime;

    if (s[i].nextevent.time < 0.)
      std::cout << "error, event times negative after synchronization" << std::endl;
    s[i].nextevent = Event(0., i, INF);
    s[i].lutime = 0.;
  }
  rscale += rinitial*gtime*growthrate;       // r defined at gtime = 0
  gtime = 0.;
  growthrate = 0.0;
  grsq = 0.0;
  rtime = 0.0;
  xmomentum = 0.0;
  plasttime = 0.0;
  SetInitialEvents();
  //check the mean radius
  //  double r_mean=0, rmin=1e6, rmax=0;
  //  for (int i=0; i<N; i++){
  //    r_mean = r_mean + s[i].r;
  //    if (s[i].r>rmax)
  //     rmax = s[i].r;
  //   if (s[i].r<rmin)
  //     rmin = s[i].r;
  //  }
  //  r_mean = r_mean/N;

  //  if(fabs(r_mean-r)>100.*DBL_EPSILON)
  //  {
  //    std::cout << "error: the global radius does not equal to the mean radius." << std::endl;
  //    exit(-1);
  //  }
  //  else
  //    r = r_mean;
  Process(N);
}

void Box::StartMeasure(double nextsampletime, double sampletimedelt)
{
  Reset();
  this->nextsampletime = nextsampletime * rmeanfin * 2.0;
  this->sampletimedelt = sampletimedelt * rmeanfin * 2.0;
}

//==============================================================
// Run time
//==============================================================
void Box::RunTime()
{
  duration_type elapsed = steady_clock::now() - tstart;
  std::cout << "run time = " << elapsed.count() << '\n';
  std::cout << "nlist rfactor = " << neighs->getRFactor() << std::endl;
}


//==============================================================
// Write configuration
//==============================================================
void Box::WriteConfiguration(const std::string wconfigfile)
{
  if (gtime != 0.)   // synchronize spheres if not currently synchronized
    Synchronize(false);

  std::ofstream output(wconfigfile,std::ios::out | std::ios::app);

  // make header
  output << DIM << " " << rtime << "\n";
  output << N << " " << 1 << "\n";
  output << N << "\n";
  output << std::setprecision(16) << rscale << " " <<  rx << "\n";

  // Lattice vectors:
  for (int i=0; i<DIM; i++)
  for (int j=0; j<DIM; j++)
  {
    if (i==j)
      output << 1 << " ";
    else
      output << 0 << " ";
  }
  output << "\n";

  // Boundary conditions:
  for (int i=0; i<DIM; i++)
  output << "T" << " ";  // T = periodic BC
  output << "\n";

  for (int i=0; i<N; i++)  // output diameter and positions
  {
    for (int k=0; k<DIM; k++)
    output << std::setprecision(16) << s[i].x[k] << " ";
    output << std::setprecision(16) << s[i].r << "";
    //output << std::setprecision(16) << s[i].gr << " ";
    //output << std::setprecision(16) << s[i].m << " ";
    output << "\n";
  }

  output.close();
}

//==============================================================
// Write shifts
//==============================================================
void Box::WriteShifts(const std::string shiftfile) const
{
  std::ofstream outshift(shiftfile, std::ios::binary);
  for(int rp=0; rp<N-1; ++rp)
  {
    for(int rq=0; rq<N-rp-1; ++rq)
    {
      allshifts[rp][rq].write(outshift);
    }
  }
}


//==============================================================
// Write LAST configuration
//==============================================================
void Box::WriteLastConfiguration(const std::string wconfigfile)
{
  if (gtime != 0.)   // synchronize spheres if not currently synchronized
    Synchronize(false);

  std::ofstream output(wconfigfile,std::ios::out);

  // make header
  output << DIM << " " << rtime << "\n";
  output << N << " " << 1 << "\n";
  output << N << "\n";
  // output << std::setprecision(16) << 2*r  << "\n";
  output << std::setprecision(16) << rscale << " " << rx << "\n";
  // Lattice vectors:
  for (int i=0; i<DIM; i++)
  for (int j=0; j<DIM; j++)
  {
    if (i==j)
      output << 1 << " ";
    else
      output << 0 << " ";
  }
  output << "\n";

  // Boundary conditions:
  for (int i=0; i<DIM; i++)
  output << "T" << " ";  // T = periodic BC
  output << "\n";


  for (int i=0; i<N; i++)  // output diameter and positions
  {
    for (int k=0; k<DIM; k++)
    output << std::setprecision(16) << s[i].x[k] << " ";
    output << std::setprecision(16) << s[i].r << "";
    //output << std::setprecision(16) << s[i].gr << " ";
    //output << std::setprecision(16) << s[i].m << " ";
    output << "\n";
  }

  output.close();
}


//==============================================================
// Write velocitities
//==============================================================
void Box::WriteVelocities(const std::string wconfigfile)
{
  //Already synchronized from above
  std::ofstream output(wconfigfile,std::ios::out);

  // make header
  output << DIM << " " << rtime << "\n";
  output << N << " " << 1 << "\n";
  output << N << "\n";
  //output << std::setprecision(16) << 2*r  << "\n";
  output << std::setprecision(16) << rscale << " " << rx << "\n";
  // Lattice vectors:
  for (int i=0; i<DIM; i++)
  for (int j=0; j<DIM; j++)
  {
    if (i==j)
      output << 1 << " ";
    else
      output << 0 << " ";
  }
  output << "\n";

  // Boundary conditions:
  for (int i=0; i<DIM; i++)
  output << "T" << " ";  // T = periodic BC
  output << "\n";


  // remember to get 16 digits of accuracy
  for (int i=0; i<N; i++)
  {
    for (int k=0; k<DIM; k++)
    output << std::setprecision(16) << s[i].v[k] << " ";
    output << "\n";
  }

  output.close();
}

//==============================================================
// get shift coordinates
//==============================================================
void Box::getshift(int i, int j, vector<>& result)
{
  if (i > j)
  {
    result = s[i].x - s[j].x - allshifts[j][i - j - 1];
  }
  else if (i == j)
  {
    result.setzero();
    return;
  }
  else
  {
    result = s[i].x - s[j].x + allshifts[i][j - i - 1];
  }
  setminimg(result.x);
}

void Box::getshift(int i, int j, vector<> &result, vector<> &shift)
{
  getshift(i, j, result);
  shift = s[i].x - s[j].x - result;
}

// neighbor displacement under shift
vector<> Box::getshift(const int i, const int j)
{
  vector<> x;
  getshift(i, j, x);
  return x;
}

// input index, output particle coordinates
void Box::vectoridx(vector<>& x, const int i)
{
  x = s[i].x;
}

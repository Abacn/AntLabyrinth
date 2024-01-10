//
//  graph_construct.cpp
//  RLG_Delaunay
// grow the cage between Poisson distributed particles from the origin of a dD-sphere

#include <algorithm>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <unordered_set>
#include <map>
#include <queue>
#include <vector>
#include <math.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/exceptions.h>

#include "graph_construct_cage.hpp"
#include "cellsample.hpp"
#include "cagestat.hpp"
#include "myutility.hpp"

// Convienient Macro
#define ttds tris.tds()  // get triangulation result

// constructor
Graph_Construct::Graph_Construct(ReadInputRLGCShape &input): input(input), tris(DIM), count(0), allstats(nullptr)
{
  // overwrite the default random number generator by specified random seed
  Vunitsphere = input.Vunitsphere; // volume of unit sphere
  input.N = input.phi*input.rconstA;
  std::cout << "expected N : " << input.N << std::endl;
  SampleUtility::setpoisson(input.N);
  CGAL::get_default_random() = CGAL::Random(input.seed);
  if(input.dostat)
  {
    allstats = new AllStat(input.dostats, input.binstart, input.rmax-1., input.nbinperl, input.logbinstart, input.rmax-1, input.nbinperlogl);
  }
}

// deconstructor
Graph_Construct::~Graph_Construct()
{
  fout.close();
  fsummary.close();
  if(allstats != nullptr) delete allstats;
}

// Given cell1 belongs to the cavity, add cell2 to the cavity or not
// facet weight = circumradius of the facet
// otherwise, weight = min(circumradius of the cell
// return value: 0 - facet not connected;
// -1 - cell already exist;
// 1 - cell added, circumcenter at different side of common facet;
// 2 - cell added, circumcenter at same side of common facet;
int Graph_Construct::process_facet(cellmap_iter cellit, const Full_cell_handle cellhdlB)
{
  // if cellhdl is already in visited cells
  int ret_val = 0;
  if(cellmap.find(cellhdlB) != cellmap.end())
  {
    ret_val = -1;
  }
  else
  {
    const Full_cell_handle cellhdlA = cellit->first;
    const Cell_samples &cellA = cellit->second;
    Cell_samples cellB = create_cell_sample(cellhdlB);
    if(cellA.subcellsign[cellhdlA->index(cellhdlB)] == -1)
    {
      if(cellA.r > 1.0)
      {
        ret_val = 2;
      }
    }
    else if(cellB.subcellsign[cellhdlB->index(cellhdlA)] == -1)
    {
      if(cellB.r > 1.0)
      {
        ret_val = 2;
      }
    }
    else // different side, then calculate facet weight
    {
      double dist = SampleUtility::norm_squared_distance(cellA.center, cellB.center);
      double weight = (0.5*(cellA.r*cellB.r + cellB.r*dist + dist*cellA.r) - 0.25*(cellA.r*cellA.r + cellB.r*cellB.r + dist*dist))/dist;
      if(weight > 1.0)
      {
        ret_val = 1;
      }
    }
    if(0 != ret_val)
    {
      cellmap.insert({cellhdlB, cellB});
    }
  }
  return ret_val;
}

// create Cell_samples object
Cell_samples Graph_Construct::create_cell_sample(Full_cell_handle cell)
{
  Simplex_array<Vector> vertices;
  for(int rp=0; rp<=DIM; ++rp)
  {
    const Point &pa = cell->vertex(rp)->point();
    for(int rq=0; rq<DIM; ++rq)
    {
      vertices[rp][rq] = pa[rq];
    }
  }
  return Cell_samples(vertices);
}

// master working entry
int Graph_Construct::graphrun()
{
  int rp;
  if(input.outmode >= OUTPUT_MODE_DEBUG)
  {
    std::cout << '[' << MyUtility::timer.get() << ']' << " Insert points..." << std::endl;
  }
  try
  {
    insert_points();
  }
  catch(const CGAL::Assertion_exception &e)
  {
    std::cerr << e.what() << std::endl;
    return -1;
  }
  // build network
  if(input.outmode >= OUTPUT_MODE_DEBUG)
  {
    std::cout << '[' << MyUtility::timer.get() << ']' << " Locate the origin..." << std::endl;
  }
  // Full cell at origin
  std::vector<double> zerovs(DIM, 0.0);
  Full_cell_handle fc_trv = tris.locate(Point(zerovs.begin(), zerovs.end()));
  cellmap.insert({fc_trv, create_cell_sample(fc_trv)});
  if(input.outmode >= OUTPUT_MODE_DEBUG)
  {
    std::cout << '[' << MyUtility::timer.get() << ']' << " Grow voids..." << std::endl;
  }
  std::queue<Full_cell_handle> cavitycells;
  cavitycells.push(fc_trv);
  bool isperco = false;
  while(!isperco && !cavitycells.empty())
  {
    // get a cavity cell
    fc_trv = cavitycells.front();
    auto cit = cellmap.find(fc_trv);
    // Traverse neighbor cells
    for(rp=0; rp<=DIM; ++rp)
    {
      Full_cell_handle nb = fc_trv->neighbor(rp);
      if(tris.is_infinite(nb))
      {
        isperco = true;
        break;
      }
      int ret_val = process_facet(cit, nb);
      if(ret_val > 0)
      {
        // new cell
        cavitycells.push(nb);
      }
    }
    cavitycells.pop();
  }
  if(input.outmode >= OUTPUT_MODE_DEBUG)
  {
    std::cout << '[' << MyUtility::timer.get() << ']' << " Random sampling in void space..." << std::endl;
  }
  if(input.samplestrat[0] > 0)
  {
    graphdump(isperco);
  }
  else
  {
    graphsum(isperco);
  }
  ++count;
  return 0;
}

// insert random point in shell (r=1 to rmax)
void Graph_Construct::Graph_Construct::insert_points()
{
  double disti, ra, rb;
  std::vector<double> modifiedpoint(DIM);
  int i, j, N;
  CGAL::Random_points_in_ball_d<Point> rand_it(DIM, input.rmax); // point generator (-1, 1)
  if(!points.empty())
  {
    points.clear();
    tris.clear();
    boundary_cells.clear();
    cellmap.clear();
  }
  N = SampleUtility::poissonrand();
  // std::cout << N << std::endl;
  CGAL::cpp11::copy_n(rand_it, N, std::back_inserter(points));
  // modify radius to the center
  for(i=0; i<N; ++i)
  {
    disti = MyUtility::norm_squared_distance(points[i]);
    if(disti <= 1.0)
    {
      // move it into outer shell
      ra = sqrt(disti);
      rb = pow(pow(ra, DIM)*input.rconstA + 1., 1.0/DIM);
      rb /= ra;
      for(j=0; j<DIM; ++j)
      {
        modifiedpoint[j] = points[i][j]*rb;
      }
      // replace the original point
      points[i] = Point(modifiedpoint.begin(), modifiedpoint.end());
    }
  }
  tris.insert(points.begin(), points.end()); // insert into tris
  if(input.outmode == OUTPUT_MODE_ALL)
  {
    std::stringstream sst;
    sst << "obstacles_" << count << ".dat";
    std::ofstream of(sst.str());
    of << std::setprecision(16);
    for(int i=0; i<N; ++i)
    {
      for(int j=0; j<DIM; ++j)
      {
        of << points[i][j] << " ";
      }
      of << "\n";
    }
  }
}

// sample over cluster
void Graph_Construct::samplecluster(std::vector<Cell_single_sample> &allsample, double &clvolume, double &sqweight, double &simplexttlv,  double &Delta2, double &largest_cell_r, Vector &largest_cell_center)
{
  double squaresum = 0.;
  double dtmp, furthest_r;
  celldata_type sampledata;
  Vector possum({0});
  clvolume = simplexttlv = squaresum = Delta2 = largest_cell_r = furthest_r = 0.;
  // first build sample data structures
  for(auto mit: cellmap)
  {
    sampledata.push_back(mit);
    auto &retpair = sampledata.back();
    const Cell_samples &sp = retpair.second;
    if(sp.r > largest_cell_r)
    {
      largest_cell_r = sp.r;
      largest_cell_center = sp.center;
    }
    // get the furthest circumcenter from the origin
    // use this to determine which sampling method should use
    dtmp = SampleUtility::norm_squared_distance(sp.center);
    if(dtmp>furthest_r) furthest_r = dtmp;
  }
  // release memory for cellmap
  cellmap.clear();
  // deside which scheme to use
  // v2: if expected not many void vertices and small cavity
  // if too few samples, try alternate method
  if(sampledata.size() <= (DIM-1)*10 || (furthest_r < 0.04 && DIM <= 6))
  {
    samplecluster_v2(sampledata, allsample, simplexttlv, clvolume, sqweight, possum, squaresum);
  }
  else
  {
    int ret_val = samplecluster_v1(sampledata, allsample, simplexttlv, clvolume, sqweight, possum, squaresum);
    if(ret_val)
    {
      if(input.outmode >= OUTPUT_MODE_DEBUG)
      {
          std::cout << '[' << MyUtility::timer.get() << ']' << " method 1 fail for " << allsample.size() << " samples. Try method 2." << std::endl;
          //std::cout << '[' << MyUtility::timer.get() << ']' << " method 1 fail for " << allsample.size() << std::endl;
      }
      samplecluster_v2(sampledata, allsample, simplexttlv, clvolume, sqweight, possum, squaresum);
    }
  }
  // calcualte MSD over cluster
  if(allsample.size() <= 1)
  {
    Delta2 = 0.;
  }
  else
  {
    Delta2 = 2.*(squaresum/clvolume - SampleUtility::norm_squared_distance(possum)/(clvolume*clvolume));
  }
  double reweight;
  // if too many samples, cut it and reweight the samples
  // This should not happen after the scheme is changed and if the input the parameter is good.
  if(allsample.size() > input.samplestrat[3])
  {
    // truncate allsample
    std::shuffle(std::begin(allsample), std::end(allsample), SampleUtility::rg);
    allsample.resize(input.samplestrat[3]);
    reweight = 0.;
    sqweight = 0.;
    for(const auto &sit: allsample)
    {
      reweight += sit.weight;
      sqweight += sit.weight*sit.weight;
    }
    reweight = clvolume/reweight;
    sqweight *= reweight*reweight;
    for(auto &sit: allsample)
    {
      sit.weight *= reweight;
    }
  }
}

int Graph_Construct::samplecluster_v1(celldata_type &sampledata, std::vector<Cell_single_sample> &allsample, double &simplexttlv, double &clvolume, double &sqweight, Vector &possum, double &squaresum)
{
  int rp, rq;
  double dtmp, dtmp2;
  long dupnum = 0, maxdupnum = input.samplestrat[1]/input.samplestrat[0], voidsample = 0, cellidx;
  long chkdupnum = maxdupnum*10/input.samplestrat[2];
  Vector vectmp;
  std::set<Vertex_handle> outvecset;
  std::vector<double> simplexvs;
  allsample.clear();
  clvolume = sqweight = squaresum = 0.;
  for(rp=0; rp<DIM; ++rp) possum[rp] = 0.;
  for(auto &mit: sampledata)
  {
    auto fc_trv = mit.first;
    auto &samples = mit.second;
    // filter out those simplices where inner space is fully occupied
    double furthest_dist = samples.r - samples.dis_to_point(samples.center);
    if(furthest_dist > 1.0)
    {
      simplexvs.push_back(samples.volume);
    }
    else
    {
      simplexvs.push_back(0.0);
    }
    // check if center in cell
    for(rp=0; rp<=DIM; ++rp)
    {
      if(samples.subcellsign[rp] == -1)
      {
        auto vit = outvecset.insert(fc_trv->vertex(rp));
      }
    }
  }
  // distribute outpoints
  std::vector<std::vector<Vector> > alloutpoints(sampledata.size());
  for(cellidx=0; cellidx<sampledata.size(); ++cellidx)
  {
    auto fc_trv = sampledata[cellidx].first;
    auto &samples = sampledata[cellidx].second;
    // consider reduce the number of alloutpoints
    for(auto outpoint: outvecset)
    {
      if(!fc_trv->has_vertex(outpoint)) // Full cell does not contain this vertex
      {
        // point to vector
        for(rp=0; rp<DIM; ++rp) vectmp[rp] = outpoint->point()[rp];
        dtmp = SampleUtility::norm_squared_distance(samples.center, vectmp);
        dtmp2 = sqrt(samples.r)+1.;
        if(dtmp < dtmp2*dtmp2) // overlap possible
        {
          // check cell to point distance
          dtmp = samples.dis_to_point(vectmp);
          if(dtmp<1.)
          {
            // obstacle (radius=1) not belonging to the cell overlaps with the cell
            alloutpoints[cellidx].push_back(vectmp);
          }
        }
      }
    }
  }
  // Then do sampling
  SampleUtility::SampleChoice schoice(simplexvs);
  simplexttlv = schoice.getw();
  while(voidsample<input.samplestrat[2] && dupnum<maxdupnum)
  {
    for(rq = 0; rq<input.samplestrat[0]; ++rq) // each subcell samplepercell samples
    {
      cellidx = schoice.getSample();
      bool is_void = sampledata[cellidx].second.onesample(alloutpoints[cellidx], true);
      if(is_void) ++voidsample;
    }
    if(dupnum == chkdupnum && voidsample == 0)
    {
      // v1 not efficiency, abort;
      return -1;
    }
    ++dupnum;
  }
  
  double sampleweight = simplexttlv/(dupnum*input.samplestrat[0]);
  for(auto &mit: sampledata)
  {
    auto &spb = mit.second;
    for(auto &sit: spb.samples) // iteration over one cell
    {
      sit.weight = sampleweight;
      clvolume += sit.weight;
      sqweight += sit.weight*sit.weight;
      for(rp=0; rp<DIM; ++rp)
      {
        dtmp = sit.coord[rp];
        squaresum += sit.weight*dtmp*dtmp;
        possum[rp] += sit.weight*dtmp;
      }
    }
    allsample.insert(allsample.end(), spb.samples.begin(), spb.samples.end());
    spb.samples.clear();
  }
  if(allsample.size() < input.samplestrat[2])
    return 1;
  else
    return 0;
}

int Graph_Construct::samplecluster_v2(celldata_type &sampledata, std::vector<Cell_single_sample> &allsample, double &simplexttlv, double &clvolume, double &sqweight, Vector &possum, double &squaresum)
{
  int rp, rq;
  long dupnum = 0, maxdupnum = input.samplestrat[1]/input.samplestrat[0], voidsample = 0, cellidx;
  double sqdis, dtmp;
  allsample.clear();
  clvolume = sqweight = squaresum = 0.;
  for(rp=0; rp<DIM; ++rp) possum[rp] = 0.;
  std::vector<Triangulation::Point> auxvertices; // vertices that covers the void space
  Simplex_array<Vector> vertices;
  std::set<Vertex_handle> invecset;
  std::vector<Vector> invecvecs;
  // data used to check if void is valid
  std::set<Full_cell_handle> cellshasvoid;
  std::vector<double> simplexvs;
  std::vector<Cell_samples> sampleslist;
  for(auto mit: sampledata)
  {
    auto fc_trv = mit.first;
    auto &samples = mit.second;
    cellshasvoid.insert(fc_trv);
    for(rp=0; rp<=DIM; ++rp)
    {
      invecset.insert(fc_trv->vertex(rp));
      std::vector<Vector> sp;
      for(rq=0; rq<=DIM; ++rq)
      {
        if(rp==rq) continue;
        sp.push_back(samples.vecs[rq]);
      }
      Vector vec = SampleUtility::projection(sp, samples.center);
      sqdis = SampleUtility::norm_squared_distance(vec, sp.front());
      if(sqdis<=1.)
      {
        dtmp = sqrt((1. - sqdis)/(samples.r - sqdis));
        for(rq=0; rq<DIM; ++rq)
        {
          vec[rq] = vec[rq]*(1-dtmp) + samples.center[rq]*dtmp;
        }
        // check if auxvertex is overlapped by the opposite vertex
        double dist = SampleUtility::norm_squared_distance(samples.vecs[rp], vec);
        if(dist > 1.0)
        {
          // std::cout << vec[0] << "\t" << vec[1] << "\n";
          auxvertices.push_back(Point(vec.begin(), vec.end()));
        }
      }
    }
  }
  assert(auxvertices.size() > DIM);
  for(auto vhd: invecset)
  {
    invecvecs.push_back(MyUtility::point2vec(vhd->point()));
  }
  if(input.outmode >= OUTPUT_MODE_DEBUG)
  {
    std::cout << '[' << MyUtility::timer.get() << ']' << " " << auxvertices.size() << " void shape vertices" << std::endl;
  }
  CGAL::Delaunay_triangulation<DIM_tag> auxtris(DIM);
  auxtris.insert(auxvertices.begin(), auxvertices.end());
  if(input.outmode >= OUTPUT_MODE_DEBUG)
  {
    std::vector<Face> edges;
    std::back_insert_iterator<std::vector<Face> > out(edges);
    auxtris.tds().incident_faces(auxtris.infinite_vertex(), 1, out);
    std::cout << '[' << MyUtility::timer.get() << ']' << " " << edges.size() << " on the convex hull" << std::endl;
  }
  int n_finite_cells = 0, n_cell_in_void = 0, n_cell_in_other = 0, n_cell_ocpy = 0;
  for(auto fc_trv=auxtris.finite_full_cells_begin(); fc_trv!=auxtris.finite_full_cells_end(); ++fc_trv)
  {
    ++n_finite_cells;
    // get vector
    for(rp=0; rp<=DIM; ++rp)
    {
      const Point &pa = fc_trv->vertex(rp)->point();
      for(rq=0; rq<DIM; ++rq)
      {
        vertices[rp][rq] = pa[rq];
      }
    }
    Cell_samples spb(vertices);
    // first try to get one sample
    for(rp=0; rp<input.samplestrat[0]; ++rp)
    {
      if(spb.onesample(invecvecs, false)) break;
    }
    if(rp==input.samplestrat[0])
    {
      // cell in occupied space, just skip
      ++n_cell_ocpy;
      continue;
    }
    else
    {
      // std::cout << rp << " " << spb.volume << "\n";
      const Vector &vecfront = spb.samples.back().getvec();
      Full_cell_handle fhd = tris.locate(Point(vecfront.begin(), vecfront.end()));
      if(cellshasvoid.find(fhd) == cellshasvoid.end())
      {
        // This case should be rare
        ++n_cell_in_other;
        continue; // this cell is not in some other cavity
      }
      else
      {
        ++n_cell_in_void;
        spb.clearsample();
      }
    }
    simplexvs.push_back(spb.volume);
    sampleslist.push_back(spb);
  }
  if(input.outmode >= OUTPUT_MODE_DEBUG)
  {
    std::cout << '[' << MyUtility::timer.get() << ']' << " n_finite_cells/n_cell_in_void/n_cell_in_other/n_cell_ocpy:"
    << n_finite_cells << " " << n_cell_in_void << " " << n_cell_in_other << " " << n_cell_ocpy << std::endl;
  }
  SampleUtility::SampleChoice schoice(simplexvs);
  while(voidsample<input.samplestrat[2] && dupnum<maxdupnum)
  {
    for(rq = 0; rq<input.samplestrat[0]; ++rq) // each subcell samplepercell samples
    {
      cellidx = schoice.getSample();
      bool is_void = sampleslist[cellidx].onesample(invecvecs, false);
      if(is_void) ++voidsample;
    }
    ++dupnum;
  }
  simplexttlv = schoice.getw();
  double sampleweight = simplexttlv/(dupnum*input.samplestrat[0]);
  for(auto &spb: sampleslist)
  {
    for(auto &sit: spb.samples) // iteration over one cell
    {
      sit.weight = sampleweight;
      clvolume += sit.weight;
      sqweight += sit.weight*sit.weight;
      for(rp=0; rp<DIM; ++rp)
      {
        dtmp = sit.coord[rp];
        squaresum += sit.weight*dtmp*dtmp;
        possum[rp] += sit.weight*dtmp;
      }
    }
    allsample.insert(allsample.end(), spb.samples.begin(), spb.samples.end());
    spb.samples.clear();
  }
  if(allsample.size() < input.samplestrat[2])
    return 1;
  else
    return 0;
}

// dump results
int Graph_Construct::graphdump(bool isperco)
{
  if(input.samplestrat[0]<=0) return -1; // do not do sample
  size_t cellsize = cellmap.size();
  if(0==count)
  {
    if(input.outmode>=OUTPUT_MODE_SAVESPL)
    {
      fout.open(input.datafile+"_cage.dat", std::ios::binary);
    }
    if(input.outmode>=OUTPUT_MODE_PRODUCT)
    {
      fdataA.open(input.datafile+"_fth.dat"); // point furthest to the boundary
      fdataA << "rmax" << "\t" << "centerCoords" << "\n";
    }
    if(input.outmode>=OUTPUT_MODE_SUMMARY)
    {
      fsummary.open(input.datafile+"_out.dat");
      fsummary << "Dimension: " << DIM << "\n";
      fsummary << "rmax:      " << input.rmax << "\n";
      // 1. number of cell the cluster included
      // 2. is the cluster percolated
      // 3. cluster volume
      // 4. Infinite time MSD (Delta^2)
      // 5. Coordinate of the point furthest from the cluster boundary
      // 6. Variance of SD: <D^2>-<D>^2
      fsummary << "NSimplex\tIsPerco\tNsample\tVolume\tDelta2\tVoidRatio\tChidyn" << std::endl;
    }
  }
  int rp, rq;
  double dtmp;
  double clvolume = 0., sqweight = 0., simplexttlv = 0., Delta2 = 0., largest_cell_r = 0, chidyn = NAN;
  Vector largest_cell_center({0});
  std::vector<Cell_single_sample> allsample;
  if(!isperco)
  {
    samplecluster(allsample, clvolume, sqweight, simplexttlv, Delta2, largest_cell_r, largest_cell_center);
    // process samples
    if(input.dostat)
    {
      allstats->processAll(allsample, clvolume, sqweight, largest_cell_center);
      chidyn = allstats->getchidyn();
    }
    // save samples
    if(input.outmode>=OUTPUT_MODE_SAVESPL)
    {
      for(const auto &sit: allsample)
      {
        fout.write((char*)&count, sizeof(count));
        dtmp = sit.weight;
        fout.write((char*)&dtmp, sizeof(dtmp));
        fout.write((char*)sit.coord, sizeof(sit.coord));
      }
    }
    //if(clvolume<0.) clvolume = 0.;
    //if(Delta2<0.) Delta2 = 0.;
  }
  else
  {
    // percolated
    clvolume = Delta2 = largest_cell_r = NAN;
  }

  if(input.outmode>=OUTPUT_MODE_SUMMARY)
  {
    fsummary << cellsize << "\t" << isperco << "\t" << allsample.size() << "\t" << std::setprecision(8) << clvolume << "\t" << Delta2
             << std::setprecision(6) << "\t" << clvolume/simplexttlv << "\t" << chidyn << std::endl;
  }
  if(input.outmode>=OUTPUT_MODE_PRODUCT)
  {
    fdataA << std::setprecision(8) << largest_cell_r;
    for(rp=0; rp<DIM; ++rp)
    {
      fdataA << "\t" << std::setprecision(8) << largest_cell_center[rp];
    }
    fdataA << std::endl;
  }
  return 0;
}

// dump all statistics
void Graph_Construct::dumpall()
{
  if(input.dostat && count>0)
  {
    allstats->dumpAll(input.datafile);
  }
}

// dump summary
int Graph_Construct::graphsum(bool isperco)
{
  size_t cellsize = cellmap.size();
  if(0==count)
  {
    fsummary.open(input.datafile+"_out.dat");
    fsummary << "Dimension: " << DIM << "\n";
    fsummary << "rmax:      " << input.rmax << "\n";
    // 1. number of cell the cluster included
    // 2. is the cluster percolated
    fsummary << "NSimplex\tIsPerco" << std::endl;
  }
  fsummary << cellsize << isperco << std::endl;
  return 0; // not used yet
}

// transform weight to Phi (=number density*Unit sphere volume)
double Graph_Construct::weight2Phi(double weight)
{
  return input.N / (pow(weight, -0.5*DIM)*input.rconstA);
}

// transform Phi to weight
double Graph_Construct::Phi2weight(double Phi)
{
  return pow(input.phi*input.rconstA/input.N, 2./DIM);
}

// private member methods
void Graph_Construct::construct_network()
{
  ;
}


#undef ttds

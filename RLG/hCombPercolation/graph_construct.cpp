//
//  graph_construct.cpp
//  QhullPercolation
//
//  Created by Yi Hu on 3/21/19.

#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "read_input.hpp"
#include "graph_construct.hpp"
#include "myutility.hpp"
#include "cellset.hpp"
#include "parameters.hpp"

// Constructor
Graph_Construct::Graph_Construct(ReadInputRLG input):
input(input), threshold(0.0), count(0),
points_original(input.N),
#ifndef _MULTITHREAD
points_rel(input.N),
points_inv(input.N),
inv_idx(input.N),
#else
point_rel_omp(input.n_thread, points(input.N)),
point_inv_omp(input.n_thread, points(input.N)),
inv_idxs(input.n_thread, std::vector<Point_handle>(input.N)),
#endif
rcut(input.rcut),
defaultoptflag(input.defaultoptflag)
{
  double phicut, rcut_fix;
  MyUtility::setseed(input.seed);  // set random seed
  Vunitsphere = pow(M_PI, DIM*0.5)/tgamma(1.0+DIM*0.5); // volume of unit sphere
#ifndef _MULTITHREAD
  inv_cur_size = 0;
#else
  prs = new int32_t[input.n_thread];
  cur_sizes = new int32_t[input.n_thread];
#endif
  phicut = Parameter::getPhiCut(input.N);
  wcut = Phi2weight(phicut);
  std::cout << "      wcut : " << wcut << std::endl;
  if(rcut < wcut)
  {
    rcut_fix = wcut*pow(1.2, 2.0/DIM);
    std::cout << "Error: rcut(" << rcut << ") too small! change to " << rcut_fix << "..." << std::endl;
  }
  if(rcut>1.5) // cutoff sphere volume greater than box size, bad
  {
    rcut_fix=1.5;
    std::cout << "Error: rcut(" << rcut << ") too large! change to " << rcut_fix << "..." << std::endl;
    rcut=rcut_fix;
    if(rcut < wcut)
    {
      std::cout << "Error: adjusted rcut too small! Imply system size too small... Quit" << std::endl;
      exit(-1);
    }
  }
}

// Destructor
Graph_Construct::~Graph_Construct()
{
#ifdef _MULTITHREAD
  delete [] prs;
  delete [] cur_sizes;
#endif
}

// Run master
int Graph_Construct::graphrun(const bool rerunflag)
{
  Triangulation_Data TData;
  counterReset();
  if(!rerunflag)
  {
    optflag = defaultoptflag;
#ifdef DEBUG
    std::cout << "  Generate random points..." << std::endl;
#endif
    genrandomPoints();
  }
  else
  {
    // re-run
    // two possible fix option;
    // (1) decrease wcut. May not work if wcut still too large
    wcut = wcut*pow(0.9, 2.0/DIM);
    optflag = defaultoptflag;
  }
#ifdef DEBUG
  std::cout << "  Build triangulation..." << std::endl;
#endif
  buildTriangulation(TData);
#ifdef DEBUG
  std::cout << "  Linking cells..." << std::endl;
#endif
  bool perco_ret_val = eval_percolation(TData);
  if(perco_ret_val) return -1;
  else return 0;
}

// Partial Run master
int Graph_Construct::partialrun(int32_t startx)
{
  counterReset();
  optflag = defaultoptflag;
  std::cout << "    startx : " << startx << std::endl;
  if(startx == -1)
  {
    // generate and save random configurations
    genrandomPoints();
    savepoints();
  }
  else if(startx>=0)
  {
    // change to startidx
    int32_t stepidx = (int)ceil((double)(input.N-DIM+1)/input.repeatrun);
    if(startx < stepidx)
    {
      // Build triangulation and save
      if(loadpoints())
      {
        std::cerr << "[" << startx << "] " << "Error! Cannot read points.dat" << std::endl;
        return -2;
      }
      Triangulation_Data TData;
      buildTriangulation(TData, startx, stepidx);
      savestate(TData, std::to_string(startx).c_str());
    }
    else if(startx==stepidx)
    {
      // last one, evaluate triangulation
      Triangulation_Data TData;
      std::vector<std::string> filelist;
      // check validity
      bool badflag = false;
      for(int32_t nowidx = 0; nowidx<stepidx; ++nowidx)
      {
        std::stringstream sst;
        sst << tdataprefix << "_" << nowidx << tdatasuffix;
        std::string st = sst.str();
        if (FILE *f = fopen(st.c_str(), "r"))
        {
          fclose(f);
          filelist.push_back(st);
        }
        else
        {
          if(!badflag)
          {
            std::cerr << "[" << startx << "] " << "Error! Following files not found:\n";
            badflag = true;
          }
          std::cerr << st << "\n";
        }
      }
      if(badflag)
      {
        std::cerr << std::flush;
        return -3;
      }
      int load_ret_Val = loadstate(TData, filelist);
      if(load_ret_Val)
      {
        std::cerr << "[" << startx << "] Loadstate error\n";
        return -1;
      }
      bool perco_ret_val = eval_percolation(TData);
      if(perco_ret_val)
      {
        std::cerr << "[" << startx << "] Percolation detection error\n";
        return -1;
      }
      // percolation is successfully evaluated
#ifndef DEBUG
      remove(pointfile);
      for(const auto &fdname: filelist)
      {
        remove(fdname.c_str());
      }
#endif
    }
    else
    {
      return -4; // invalid startx
    }
  }
  return 0;
}

// unit tests
int Graph_Construct::graphdebug()
{
  // last one, evaluate triangulation
  Triangulation_Data TData;
  std::vector<std::string> filelist = {"tdata_0.dat"};
  std::cout << "Load from tdata_0.dat..." << std::endl;
  int load_ret_Val = loadstate(TData, filelist);
  if(load_ret_Val)
  {
    std::cout << "Load state fail... Generate new..." << std::endl;
    genrandomPoints();
    TData.clear();
    buildTriangulation(TData);
    savestate(TData, "0");
  }
  std::vector<Facet_decorate> wholefctlist;
  wholefctlist.swap(TData.facetlist);
  int32_t atmp[DIM+1]; // store common facet
  // arrange cells
  std::vector<Base_Cell> bcells(TData.celllist.size());
  for(const auto &val : TData.cellmap)
  {
    bcells[val.second] = val.first;
  }
  for(int rp=0; rp<100; ++rp)
  {
    // test the percolation threshold if include only first ev_point points
    ev_point = round(0.01*(rp+1)*input.N);
    TData.facetlist.clear();
    for(auto const &fct: wholefctlist)
    {
      bcells[fct.cell1].get_common_facet(bcells[fct.cell2], atmp);
      if(atmp[0] < ev_point) TData.facetlist.push_back(fct);
    }
    eval_percolation(TData);
  }
  return 0;
}

// dump summary
int Graph_Construct::graphsum(const Triangulation_Data &TData, std::size_t vstfacet, std::size_t vstcell)
{
  std::size_t totalfacet = TData.facetsize, totalcell=TData.cellsize;
  if(1==count)
  {
    fsummary.open(input.datafile+"_out.dat", std::ofstream::out);
    fsummary << "Dimension: " << DIM << "\n";
    fsummary << "N:         " << input.N << "\n";
    fsummary << "Phi\tNFacet\tNvFacet\tNCell\tNvCell" << std::endl;
  }
  else
  {
    fsummary.open(input.datafile+"_out.dat", std::ios_base::app);
  }
  fsummary << std::setprecision(9) << weight2Phi(threshold) << "\t";
  fsummary << totalfacet << "\t" << vstfacet << "\t" << totalcell << "\t" << vstcell << std::endl;
  fsummary.close();
  if(fsummary.bad())
  {
    return -1;
  }
  else
  {
    return 0;
  }
}

void Graph_Construct::graphsumextend(const Triangulation_Data &TData, const char *outputline)
{
  if(input.outmode >= OUTPUT_MODE_PRODUCT)
  {
    // count overflow facet
    if(count == 1)
    {
      fdbgsummary.open("dbgsummary.dat");
      fdbgsummary << "NFarr\tNCarr\tNevp";
      fdbgsummary << "\tPhi1\tPhi2\tPhi3\tPhi4\tPhi5";
      fdbgsummary << "\tpcclsize\tNClu";
      fdbgsummary << "\tNneighbor\tVar\tmaxcellr\tmaxnbrdis" << std::endl;
    }
    else
    {
      fdbgsummary.open("dbgsummary.dat", std::ios_base::app);
    }
    fdbgsummary << outputline << std::endl;
    fdbgsummary.close();
  }
}

void Triangulation_Data::add_facet(void)
{
  ++facetsize;
}

void Triangulation_Data::add_facet(double weight, Cell_handle cell1, Cell_handle cell2)
{
  facetlist.push_back(Facet_decorate(weight, cell1, cell2));
  ++facetsize;
}

Cell_handle Triangulation_Data::add_cell(void)
{
  ++cellsize;
  return Null_cell_handle;
}

Cell_handle Triangulation_Data::add_cell(const Full_Cell &fc)
{
  Cell_handle cellidx = (Cell_handle)celllist.size();
  cellmap.insert(TriMapType::value_type(fc, cellidx));
  celllist.push_back(fc.get_circumcenter());
  ++cellsize;
  return cellidx;
}

Cell_handle Triangulation_Data::add_cell(const Base_Cell &bc, const point &circcumcenter)
{
  Cell_handle cellidx = (Cell_handle)celllist.size();
  cellmap.insert(TriMapType::value_type(bc, cellidx));
  celllist.push_back(circcumcenter);
  ++cellsize;
  return cellidx;
}

Cell_handle Triangulation_Data::add_cell(TriMapType::const_iterator hint, const Full_Cell &fc)
{
  Cell_handle cellidx = (Cell_handle)celllist.size();
  cellmap.insert(hint, TriMapType::value_type(fc, cellidx));
  celllist.push_back(fc.get_circumcenter());
  ++cellsize;
  return cellidx;
}

Cell_handle Triangulation_Data::add_cell(TriMapType::const_iterator hint, const Base_Cell &bc, const point &circcumcenter)
{
  Cell_handle cellidx = (Cell_handle)celllist.size();
  cellmap.insert(hint, TriMapType::value_type(bc, cellidx));
  celllist.push_back(circcumcenter);
  ++cellsize;
  return cellidx;
}

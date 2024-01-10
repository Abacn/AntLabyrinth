//
//  main.cpp
//  RLG_Delaunay
//

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>

#include "cagestat.hpp"
#include "cellsample.hpp"
#include "graph_construct_cage.hpp"
#include "myutility.hpp"
#include "read_cageshape.hpp"

int triangulation_main(ReadInputRLGCShape &input)
{
  long rp, stepview, nextview;
  MyUtility::timer.set();
  stepview = input.repeatrun / 50;
  nextview = stepview;
  Graph_Construct graph(input);
  std::cout << '[' << MyUtility::timer.get() << ']' << " Start calculation..." << std::endl;
  rp = 0;
  for(rp = 0; rp<input.repeatrun; ++rp)
  {
    if(rp >= nextview)
    {
      std::cout << "(" << rp <<  ")" << std::endl;
      nextview = rp+stepview;
      graph.dumpall();
    }
    while(graph.graphrun());
  }
  graph.dumpall();
  std::cout << "(" << rp <<  ")" << std::endl;
  //#ifdef DEBUG
  std::cout << "Time per run (s):" << MyUtility::timer.get()/input.repeatrun << std::endl;
  //#endif
  return 0;
}

int processold_main(ReadInputRLGCShape &input, std::vector<std::string> filelists, std::vector<std::string> centerlists)
{
  // process all
  int rp, rq;
  long stepview, nextview;
  if(filelists.size() != centerlists.size())
  {
    return -1;
  }
  stepview = filelists.size() / 50;
  nextview = stepview;
  AllStat allstat(input.dostats, input.binstart, input.rmax-1., input.nbinperl, input.logbinstart, input.rmax-1, input.nbinperlogl);
  Vector center, coord;
  double r_center;
  std::uint32_t cluidx, oldclidx;
  float weight;
  std::vector<Cell_single_sample> samples;
  double ttlweight=0.;
  float fcoords[DIM];
  double cluv = 0., sqweight = 0; // cluster volume
  for(rp=0; rp<filelists.size(); ++rp)
  {
    if(rp >= nextview)
    {
      std::cout << filelists[rp] << std::endl;
      nextview = rp+stepview;
      allstat.dumpAll(input.datafile);
    }
    std::ifstream ifs(filelists[rp], std::ios::binary);
    std::ifstream ifc(centerlists[rp]);
    if(!ifs)
    {
      std::cout << "Warning: open \"" << filelists[rp] << "\" failed" << std::endl;
      continue;
    }
    if(!ifc)
    {
      std::cout << "Warning: open \"" << centerlists[rp] << "\" failed" << std::endl;
      continue;
    }
    ifs.read((char*)&cluidx, sizeof(cluidx));    // read the first sample cluster index
    for(oldclidx = cluidx; !ifs.eof(); ifs.read((char*)&cluidx, sizeof(cluidx)))
    {
      if(cluidx != oldclidx) // previous cluster ends
      {
#ifdef DEBUG
        if(samples.size()>=1000) std::cout << "Cluster " << cluidx << " has " << samples.size() << " samples." << std::endl;
#endif
        // read center
        ifc >> r_center;
        for(rq=0; rq<DIM; ++rq)
        {
          ifc >> center[rq];
        }
        // process a cluster
        allstat.processAll(samples, cluv, sqweight, center);
        samples.clear();
        oldclidx = cluidx;
        cluv = sqweight = 0.;
      }
      ifs.read((char*)&weight, sizeof(weight));
      ifs.read((char*)fcoords, sizeof(fcoords));
      for(rq=0; rq<DIM; ++rq) coord[rq] = (double)fcoords[rq];
      samples.push_back(Cell_single_sample(weight, 0.f, coord)); // singlesample.r not used yet
      cluv += weight;
      sqweight += weight*weight;
    }
    // read the last center
    ifc >> r_center;
    for(rq=0; rq<DIM; ++rq)
    {
      ifc >> center[rq];
    }
    // process a the last cluster
    allstat.processAll(samples, cluv, sqweight, center);
    ttlweight += cluv;
  }
  // dump
  allstat.dumpAll(input.datafile);
  return 0;
}

// Main function
int main(int argc, const char * argv[])
{
  int retval = 0;
  const char *inputf, *flistf, *dftinputf = "input_CageShape.dat", *dftfilelistf = "samplefilelist.dat";
  if(argc>=2)
  {
    if(!strcmp(argv[1], "--help"))
    {
      std::cout << "Cavity construction in " << DIM << "D. \nUsage: cavconst [parameterfile] [samplelistfile]" << std::endl;
    }
  }
  if(argc < 3)
  {
    if (argc < 2)
    {
      inputf = dftinputf;
    }
    else
    {
      inputf = argv[1];
    }
    flistf = dftfilelistf;
  }
  else
  {
    inputf = argv[1];
    flistf = argv[2];
  }
  ReadInputRLGCShape input;
  int error = input.read(inputf);
  if (error) return error;
  SampleUtility::setseed(input.seed);
  if(input.dostat != 2)
  {
    retval = triangulation_main(input);
  }
  else
  {
    // read filelist from  and do stat
    std::ifstream infile(flistf);
    std::vector<std::string> filelists, centerlists;
    std::string toReplace("_cage.dat");
    if(!infile)
    {
      std::cout << "Can't open " << flistf << " for sample list." << std::endl;
      return 2;
    }
    while(!infile.eof())
    {
      std::string line;
      std::getline(infile, line);
      if(!line.empty())
      {
        std::size_t pos = line.find(toReplace);
        if(pos == std::string::npos) continue; // not valid input name
        filelists.push_back(line);
        line.replace(pos, toReplace.length(), "_fth.dat");
        centerlists.push_back(line);
      }
    }
    if(filelists.empty())
    {
      std::cout << "Empty file list." << std::endl;
      return 4;
    }
    retval = processold_main(input, filelists, centerlists);
  }
  return retval;
}

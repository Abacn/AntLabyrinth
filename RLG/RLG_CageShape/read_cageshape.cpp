//
//  read_input_shell.cpp
//  RLG_Delaunay
//
//  Created by Yi Hu on 3/19/19.
//

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/algorithm/string.hpp>
#include "read_cageshape.hpp"
#include "dim.h"

int ReadInputRLGCShape::read(const char* inputf)
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
  infile.get(buf,300,'='); infile.get(c); infile >> repeatrun;
  infile.get(buf,300,'='); infile.get(c); infile >> phi;
  infile.get(buf,300,'='); infile.get(c); infile >> rmax;
  infile.get(buf,300,'='); infile.get(c); infile >> seed;
  infile.get(buf,300,'='); infile.get(c); infile >> datafile;
  infile.get(buf,300,'='); infile.get(c); infile >> outstr;
  infile.get(buf,300,'='); infile.get(c);
  // default sample numbers; maximum sample numbers; smallest number of void samples in a cell;
  // largest number of void samples dumped in a cluster
  // if samplestrat[0]>0 then do sampling
  infile >> samplestrat[0] >> samplestrat[1] >> samplestrat[2] >> samplestrat[3];
  // statistics parameter
  infile.get(buf,300,'='); infile.get(c); infile >> dostat;
  if(dostat)
  {
    // number of bins per length
    infile.get(buf,300,'='); infile.get(c); infile >> dostats[0];
    infile.get(buf,300,'='); infile.get(c); infile >> dostats[1];
    infile.get(buf,300,'='); infile.get(c); infile >> dostats[2];
    infile.get(buf,300,'='); infile.get(c); infile >> dostats[3];
    infile.get(buf,300,'='); infile.get(c); infile >> binstart >> nbinperl;
    infile.get(buf,300,'='); infile.get(c); infile >> logbinstart >> nbinperlogl;
  }
  else
  {
    dostats[0] = dostats[1] = dostats[2] = false;
  }
  // some useful constants
  Vunitsphere = pow(M_PI, DIM*0.5)/tgamma(1.0+DIM*0.5);
  rconstA = pow(rmax, DIM) - 1.;
  if(boost::iequals(outstr, "debug"))
  {
    outmode = OUTPUT_MODE_DEBUG;
  }
  else if(boost::iequals(outstr, "product"))
  {
    outmode = OUTPUT_MODE_PRODUCT;
  }
  else if(boost::iequals(outstr, "savespl"))
  {
    outmode = OUTPUT_MODE_SAVESPL;
  }
  else if(boost::iequals(outstr, "all"))
  {
    outmode = OUTPUT_MODE_ALL;
  }
  else if(boost::iequals(outstr, "summary"))
  {
    outmode = OUTPUT_MODE_SUMMARY;
  }
  else
  {
    outmode = OUTPUT_MODE_PRODUCT;
  }
  if(infile.eof())
  {
    std::cout << "Error reading input file " << inputf << std::endl;
    error = 3;
  }
  std::cout << "       DIM : " << DIM << "\n";
  std::cout << "Repeat run : " << repeatrun << "\n";
  std::cout << "       Phi : " << phi << "\n";
  std::cout << "      rmax : " << rmax << "\n";
  std::cout << "      Seed : " << seed<< "\n";
  std::cout << "  datafile : " << datafile << "\n";
  std::cout << "    dostat : " << dostats[0] << " " <<  dostats[1] << " " << dostats[2] << " " << dostats[3] << "\n";
  std::cout << "cellsample : " << samplestrat[0] << " " << samplestrat[1] << " " << samplestrat[2] << " " << samplestrat[3] << "\n";
  return error;
}

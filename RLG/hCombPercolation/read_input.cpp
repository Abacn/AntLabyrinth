//
//  read_input.cpp
//  RLG_Delaunay
//
//  Created by Yi Hu on 8/6/18.

#include "read_input.hpp"
#include <string>

int ReadInputRLG::read(const char* inputf)
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
  infile.get(buf,300,'='); infile.get(c); infile >> N;
  infile.get(buf,300,'='); infile.get(c); infile >> seed;
  infile.get(buf,300,'='); infile.get(c); infile >> datafile;
  infile.get(buf,300,'='); infile.get(c); infile >> outstr;
  infile.get(buf,300,'='); infile.get(c); infile >> defaultoptflag;
  infile.get(buf,300,'='); infile.get(c); infile >> rcut;
  infile.get(buf,300,'='); infile.get(c); infile >> n_thread;
  if(!outstr.compare("debug"))
  {
    outmode = OUTPUT_MODE_DEBUG;
  }
  else if(!outstr.compare("savegraph"))
  {
    outmode = OUTPUT_MODE_SAVEGRAPH;
  }
  else if(!outstr.compare("product"))
  {
    outmode = OUTPUT_MODE_PRODUCT;
  }
  else if(!outstr.compare("all"))
  {
    outmode = OUTPUT_MODE_ALL;
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
  std::cout << "Repeat run : " << repeatrun << "\n";
  std::cout << "         N : " << N << "\n";
  std::cout << "      Seed : " << seed<< "\n";
  std::cout << "  datafile : " << datafile << "\n";
  std::cout << "dftoptflag : " << defaultoptflag << "\n";
  std::cout << "      rcut : " << rcut << std::endl;
  std::cout << "   nthread : " << n_thread << std::endl;
  return error;
}

#include <iomanip>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cstring>

#include "read_input.h"
#include "coordinate.h"

//================================================================
//
// Source File for input
//
//================================================================
int read_input_lattice::read(const char *inputf)
{
  int error = 0;
  

  std::ifstream infile(inputf);
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
  char buf[100],c;
  infile.get(buf,100,'='); infile.get(c); infile >> L;
  infile.get(buf,100,'='); infile.get(c); infile >> pf;
  infile.get(buf,100,'='); infile.get(c); infile >> maxtime;
  infile.get(buf,100,'='); infile.get(c);
  infile.width(NAME_LEN-1); infile >> datafile;
  infile.get(buf,100,'='); infile.get(c); infile >> seed;
  infile.get(buf,100,'='); infile.get(c); infile >> repeat;
  if(infile.eof())
  {
    std::cout << "Error reading input file " << inputf << std::endl;
    error = 3;
  }
  std::cout << "   DIM: " << DIM << std::endl;
  std::cout << "   L : " << L << std::endl;
  std::cout << "   N : " << pow(L,DIM) << std::endl;
  std::cout << "   pf : " << pf << std::endl;
  std::cout << "   maxtime : " << maxtime << std::endl;
  std::cout << "   datafile : " << datafile << std::endl;
  std::cout << "   seed : " << seed << std::endl;
  std::cout << "   repeat : " << repeat << std::endl;
  if(seed==-1)
  {
    std::cout << "   starting time : " << tprint << std::endl;
  }


  return error;
}

int read_input_leath::read(const char *inputf)
{
  int error = 0;
  
  
  std::ifstream infile(inputf);
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
  char buf[100],c;
  infile.get(buf,100,'='); infile.get(c); infile >> L;
  infile.get(buf,100,'='); infile.get(c); infile >> pf;
  infile.get(buf,100,'='); infile.get(c); infile >> maxcluster;
  infile.get(buf,100,'='); infile.get(c);
  infile.width(NAME_LEN-1); infile >> datafile;
  infile.get(buf,100,'='); infile.get(c); infile >> seed;
  infile.get(buf,100,'='); infile.get(c); infile >> maxsite;
  if(infile.eof())
  {
    std::cout << "Error reading input file " << inputf << std::endl;
    error = 3;
  }
  
  // std::cout << "   L : " << L << std::endl;
  std::cout << "   DIM: " << DIM << std::endl;
  std::cout << "   pf : " << pf << std::endl;
  std::cout << "   maxcluster : " << maxcluster << std::endl;
  std::cout << "   datafile : " << datafile << std::endl;
  std::cout << "   seed : " << seed << std::endl;
  std::cout << "   maxsite : " << maxsite << std::endl;
  if(seed==-1)
  {
    std::cout << "   starting time : " << tprint << std::endl;
  }
  
  
  return error;
}

int read_input_ant::read(const char *inputf)
{
  int error = 0;
  
  
  std::ifstream infile(inputf);
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
  char buf[100],c;
  infile.get(buf,100,'='); infile.get(c); infile >> L;
  infile.get(buf,100,'='); infile.get(c); infile >> pf;
  infile.get(buf,100,'='); infile.get(c); infile >> maxcluster;
  infile.get(buf,100,'='); infile.get(c);
  infile.width(NAME_LEN-1); infile >> datafile;
  infile.get(buf,100,'='); infile.get(c); infile >> seed;
  infile.get(buf,100,'='); infile.get(c); infile >> arraylen;
  if(infile.eof())
  {
    std::cout << "Error reading input file " << inputf << std::endl;
    error = 3;
  }
  
  std::cout << "   DIM: " << DIM << std::endl;
  std::cout << "   pf : " << pf << std::endl;
  std::cout << "   maxcluster : " << maxcluster << std::endl;
  std::cout << "   datafile : " << datafile << std::endl;
  std::cout << "   seed : " << seed << std::endl;
  std::cout << "   Arraylen : " << arraylen << std::endl;
  if(seed==-1)
  {
    std::cout << "   starting time : " << tprint << std::endl;
  }
  
  
  return error;
}

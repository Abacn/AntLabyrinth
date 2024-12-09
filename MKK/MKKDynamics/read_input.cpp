#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cstring>

#include "read_input.h"

#define BUFSIZE 1024

//================================================================
//
// Source File for input
//
//================================================================
int read_input::read(const char *fname)
{
  int error = 0;

  std::ifstream infile;
  infile.open(fname);
  if(!infile)
  {
      std::cout << "Can't open " << fname << " for input." << std::endl;
      error = 2;
      return error;
  }
  else
  {
      std::cout << "Reading input from file " << fname << std::endl;
  }
  char buf[BUFSIZE],c;
  infile.get(buf, BUFSIZE,'='); infile.get(c); infile >> eventspercycle;
  infile.get(buf, BUFSIZE,'='); infile.get(c); infile >> N;
  infile.get(buf, BUFSIZE, '='); infile.get(c); infile >> polydispersity;
  infile.get(buf, BUFSIZE,'='); infile.get(c); infile >> initialpf;
  infile.get(buf, BUFSIZE,'='); infile.get(c); infile >> maxpf;
  infile.get(buf, BUFSIZE,'='); infile.get(c); infile >> temp;
  infile.get(buf, BUFSIZE,'='); infile.get(c); infile >> growthrate;
  infile.get(buf, BUFSIZE, '='); infile.get(c); infile >> rfactor;
  infile.get(buf, BUFSIZE,'='); infile.get(c); infile >> maxpressure;
  infile.get(buf, BUFSIZE,'='); infile.get(c); infile >> mintime;
  infile.get(buf, BUFSIZE,'='); infile.get(c); infile >> maxtime;
  infile.get(buf, BUFSIZE,'='); infile.get(c); infile >> maxcollisionrate;
  infile.get(buf, BUFSIZE,'='); infile.get(c); infile >> shiftscale;
  infile.get(buf, BUFSIZE,'='); infile.get(c); infile >> calcpmsd;
  infile.get(buf, BUFSIZE, '=');
  if (strstr(buf, "t_dist_start") != NULL)
  {
    // read t_dist_start
    infile.get(c); infile >> t_dist_start;
    infile.get(buf, BUFSIZE, '=');
  }
  infile.get(c);
  infile.width(NAME_LEN-1); infile >> readfile;
  infile.get(buf, BUFSIZE,'='); infile.get(c);
  infile.width(NAME_LEN-1); infile >> writefile;
  infile.get(buf, BUFSIZE,'='); infile.get(c);
  infile.width(NAME_LEN-1); infile >> configprefixes;
  infile.get(buf, BUFSIZE,'='); infile.get(c);
  infile.width(NAME_LEN-1); infile >> datafile;
  //If the seed is -1, then read the velocities from the velocityin.dat file, otherwise generate new velocities
  infile.get(buf, BUFSIZE,'='); infile.get(c); infile >> seed;
  if(seed==-1)
  {
      infile.get(buf, BUFSIZE,'='); infile.get(c); infile >> tprint;
  }
  if(infile.eof())
  {
      std::cout << "Error reading input file " << fname << std::endl;
      error = 3;
  }
  if (0 == strcmp(configprefixes, "null")) outconfig_switch = 0;
  else outconfig_switch = 1;

  std::cout << "   eventspercycle : " << eventspercycle << '\n';
  std::cout << "   N : " << N << '\n';
  std::cout << "   polydispersity : " << polydispersity << '\n';
  std::cout << "   initialpf : " << initialpf << '\n';
  std::cout << "   maxpf : " << maxpf << '\n';
  std::cout << "   temp : " << temp << '\n';
  std::cout << "   growthrate : " << growthrate << '\n';
  std::cout << "   maxpressure : " << maxpressure << '\n';
  std::cout << "   mintime : " << mintime << '\n';
  std::cout << "   maxtime : " << maxtime << '\n';
  std::cout << "   maxcollisionrate : " << maxcollisionrate << '\n';
  std::cout << "   shiftscale : " << shiftscale << '\n';
  std::cout << "   readfile : " << readfile << '\n';
  std::cout << "   writefile : " << writefile << '\n';
  std::cout << "   configfileprefix : " << configprefixes << '\n';
  std::cout << "   datafile : " << datafile << '\n';
  std::cout << "   seed : " << seed << std::endl;
  //   std::cout << "   waiting time: " << tw << std::endl;
  if(seed==-1)
  {
      std::cout << "   starting time : " << tprint << std::endl;
  }


  std::cout << std::endl;
  return error;
}

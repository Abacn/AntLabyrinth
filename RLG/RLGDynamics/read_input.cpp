//
//  read_input.cpp
//  RLG_Dynamics
//
//  Created by Yi Hu on 8/6/18.

#include <iostream>
#include <fstream>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "read_input.hpp"

int ReadInputDy::read(const char* inputf)
{
  std::string inputfs(inputf);
  int errcode = 0;
  if (inputfs.size()>4 && boost::iequals(inputfs.substr(inputfs.length() - 4), ".ini"))
  {
    errcode = read_ini(inputf);
  }
  else
  {
    errcode = read_dat(inputf);
  }
  if (errcode != 0)
  {
    return errcode;
  }
  std::cout << "Repeat run : " << repeatrun << "\n";
  if(systype == SYSTEM_TYPE_BOX || systype == SYSTEM_TYPE_HCOMB)
  {
    std::cout << "         N : " << N << "\n";
  }
  else if(systype == SYSTEM_TYPE_SPHERE)
  {
    std::cout << "      rmax : " << rmax << "\n";
  }
  std::cout << "       phi : " << phi << "\n";
  std::cout << "       2^t : " << texp << "\n";
  std::cout << "Sample interval (2^t";
  switch(tscaleflag)
  {
    case(1): std::cout << ")"; break;
    case(2): std::cout << " scaled)"; break; // time now scaled as t/sqrt(dim)
    case(3): std::cout << " gaussian)"; break;
    case(4): std::cout << " brownian)"; break;
    default: std::cerr << "\nInvalid tscaleflag " << tscaleflag << std::endl;
      return errcode = 5;
  }
  std::cout << " : " << sleft << " 1/" << sinterval << "\n";
  std::cout << "      Seed : " << seed << "\n";
  std::cout << "  datafile : " << datafile << "\n";
  if(logdatafile.size() > 0)
    std::cout << " logdata   : " << logdatafile << "\n";
  std::cout << "  rcuttype : " << rcuttype << "\n";
  std::cout << "      rcut : " << rcut << std::endl;
  return errcode;
}

int ReadInputDy::read_ini(const char* inputf)
{
  // read configurations in a .ini file
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(inputf, pt);
  std::cout << "Reading input from file " << inputf << std::endl;
  std::string outstr = pt.get<std::string>("runconfig.systemtype");
  systype = getSystemType(outstr);
  if(SYSTEM_TYPE_ERROR == systype)
  {
    return 4;
  }
  else if(SYSTEM_TYPE_SPHERE == systype)
  {
    rmax = pt.get<double>("obstacles.rmax");
    // in sphere geometry, rcuttype is used to set maximal escape shell in tracer escape statistics, instead.
    // TODO: dynamical planting the obstacles when rcuttype set to 0
    rcuttype = pt.get<int>("runconfig.rcuttype", 0);
  }
  else
  {
    // number of obstacles will be generated in sphere geometry
    N = pt.get<unsigned int>("obstacles.N");
    // rcuttype is basically a deprecated parameter and set to 2 by default
    rcuttype = pt.get<int>("runconfig.rcuttype", 2);
  }
  repeatrun = pt.get<int>("runconfig.repeatrun");
  phi = pt.get<double>("obstacles.Phi");
  sleft = pt.get<int>("runconfig.sampleLeft");
  texp = pt.get<int>("runconfig.texp");
  sinterval = pt.get<int>("runconfig.sampleInterval");
  tscaleflag = pt.get<int>("runconfig.tscaleFlag");
  seed = pt.get<unsigned int>("runconfig.seed");
  datafile = pt.get<std::string>("outputconfig.datafilePrefix", "msd");
  logdatafile = pt.get<std::string>("outputconfig.logdatafilePrefix", "-");
  outstr = pt.get<std::string>("outputconfig.outputLevel");
  outmode = getOutputMode(outstr);
  rcut = pt.get<double>("runconfig.rcut");
  return 0;
}

int ReadInputDy::read_dat(const char* inputf)
{
  // Legacy version of read configuration from file
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
  char buf[1000],c;
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
    // tscaleflag: 1 - do not rescale time
    // 2 - rescale time with sqrt(DIM)
    // 3 - use gaussian velocity
    // 4 - brownian velocity
    tscaleflag = 2; // backward compactibility
    infile.clear();
  }
  infile.get(buf,300,'='); infile.get(c); infile >> seed;
  infile.get(buf,300,'='); infile.get(c); infile >> datafile;
  infile.get(buf,300,'='); infile.get(c); infile >> logdatafile;
  if(logdatafile == "-")
  {
    // not output log msd
    logdatafile = "";
  }
  infile.get(buf,300,'='); infile.get(c); infile >> outstr;
  // rcuttype: 1 (only supported in box) neighbor lists constructed on obstacle
  // 2: neighbor list constructed on tracer
  // For sphere geometry, this parameter is used to set maximal escape shell in
  // tracer escape statistics, instead.
  infile.get(buf,300,'='); infile.get(c); infile >> rcuttype;
  infile.get(buf,300,'='); infile.get(c); infile >> rcut;
  outmode = getOutputMode(outstr);
  if(infile.eof())
  {
    std::cout << "Error reading input file " << inputf << std::endl;
    error = 3;
  }
  return error;
}

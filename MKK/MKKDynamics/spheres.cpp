//===========================================================
//===========================================================
//===========================================================
//
//  Molecular dynamics simulation of Mk/MKK model
//
//===========================================================
//===========================================================
//===========================================================

#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <iomanip>

#include "box.h"
#include "sphere.h"
#include "event.h"
#include "heap.h"
#include "read_input.h"
#include "utility.h"

int main(int argc, char **argv)
{
  read_input input;
  int mult;
  //double r, rmin, rmax;
  const char *dftinput = "inputfile.dat", *infname;

  std::cout << "DIM = " << DIM << std::endl;

  if (argc < 2) infname = dftinput;
  else infname = argv[1];

  int error = input.read(infname);

  if (error) return error;
  if (input.seed == -1)
    rg.seed(time(0));
  else
    rg.seed(input.seed);   // initialize the random number generator

  Box b(input);

  if(strcmp(input.readfile, "new")==0)
    input.readfile[0]=0;

  if (input.readfile[0]) // read in existing configuration
  {
    // read the header
    std::cout << "Reading in positions of spheres" << std::endl;
    b.RecreateSpheres(input.readfile, input.temp);
  }
  else
  {
    std::cout << "Creating new positions of spheres" << std::endl;
    b.CreateSpheres(input.temp);
    if (input.outconfig_switch == 1)
    {
      b.WriteLastConfiguration(std::string(input.configprefixes) + "_config_initial.dat");
      b.WriteVelocities(std::string(input.configprefixes) + "_v_initial.dat");
      if (0.0 != input.shiftscale)
        b.WriteShifts(std::string(input.configprefixes) + "_shifts.dat");
    }
  }
  std::cout << "nlist tol = " << sqrt(b.neighs->getsqrhfsh()) << std::endl;
  //  Clear output file
  //    std::ofstream output(input.writefile);
  //    output.close();

  std::ofstream output(input.datafile);
  output.precision(10);
  double tprint = 20 * input.mintime / DIM;
  double maxtime = input.maxtime / DIM;  // scaled maxtime
  mult=input.N;
  double sigma;

  // process all checks
  b.Process(mult);

  //procedure 1: compress to maxpf
  if (b.PackingFraction() < input.maxpf)
  {
    std::cout << "Procedure 1: compress to maxpf ..." << std::endl;
    while (b.PackingFraction() < input.maxpf) {
      b.Process(mult);
      sigma = b.rscale * 2.0;
      if ((b.gtime + b.rtime) / sigma > tprint)
      {
        b.Synchronize(true);
        auto thermos = b.Pressure();
        output << (b.gtime + b.rtime) / sigma << '\t' << b.PackingFraction() << '\t' << thermos[0] << std::endl;
        tprint = tprint * 1.1892071150027211;
      }
    }
    if (input.outconfig_switch == 1)
    {
      b.WriteLastConfiguration(std::string(input.configprefixes) + "_config_compressed.dat");
      b.WriteVelocities(std::string(input.configprefixes) + "_v_compressed.dat");
    }
  }

  //reset compression rate to zero;
  b.Reset();
  sigma = b.rmeanfin* 2.0;
  tprint = input.mintime / DIM;

  //procedure 2: equilibrate
  double MSD=0., rtyp2=0., MQD=0., dtmp;
  std::cout << "Procedure 2: equilibrate ..." << std::endl;
  double eqtime = std::min(1e3 / DIM, maxtime / 4.0);
  while ((b.gtime + b.rtime)/sigma < eqtime && MSD < 20.0 / DIM)
  {
    b.Process(mult*input.eventspercycle, 1);
    if((b.gtime + b.rtime)/sigma > tprint)
    {
      b.Synchronize(false);
      MSD = 0.;
      MQD = 0.;
      for (int i=0; i<b.N; i++){
        double dr2 = vector<>::norm_squared(b.s[i].x, b.x0[i]);
        MSD += dr2;
        MQD += dr2*dr2;
      }
      dtmp = sigma*sigma;
      MSD = MSD/input.N/dtmp;
      MQD = MQD/input.N/(dtmp*dtmp);
      auto thermos = b.Pressure();
      output << (b.gtime + b.rtime) / sigma << '\t' << thermos[0];
      output << '\t' << MSD << '\t' << MQD << std::endl;
      tprint = tprint * 1.1892071150027211;
    }
  }

  tprint = input.mintime / DIM;
  b.StartMeasure(maxtime/4, maxtime*3/4000, tprint, maxtime);
  std::ofstream outMSD("MSD.dat");
  outMSD.precision(10);

  //procedure 3: measure MSD and rtyp
  std::cout << "Procedure 3: measure MSD and rtyp ..." << std::endl;

  while ((b.gtime + b.rtime)/sigma < maxtime)
  {
    b.Process(0, 1, tprint);

    b.Synchronize(false);
    // time evolution statistics
    MSD = 0.;
    rtyp2 = 0.;
    MQD = 0.;
    for (int i=0; i<b.N; i++) {
      double dr2 = vector<>::norm_squared(b.s[i].x, b.x0[i]);
      MSD = MSD + dr2;
      rtyp2 = rtyp2 + log(dr2);
      MQD += dr2*dr2;
    }
    dtmp = sigma*sigma;
    MSD = MSD/input.N/dtmp;
    rtyp2 = exp(rtyp2/input.N)/dtmp;
    MQD = MQD/input.N/(dtmp*dtmp);
    auto thermos = b.Pressure();
    outMSD << (b.gtime + b.rtime) / sigma << '\t' << MSD << '\t' << rtyp2 << '\t' << MQD << '\t' << thermos[0];
    outMSD << std::endl;
    tprint = tprint * 1.1892071150027211; // 2^(1/4)
  }
  output.close();
  outMSD.close();
  b.PrintStatistics(1);
  b.RunTime();

  if (input.outconfig_switch == 1)
  {
    b.WriteLastConfiguration(std::string(input.configprefixes) + "_config_last.dat");
    b.WriteVelocities(std::string(input.configprefixes) + "_v_last.dat");
  }
  return 0;
}

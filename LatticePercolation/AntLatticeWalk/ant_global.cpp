//
//  ant_global.cpp
//  AntLatticeWalk
//
//  Dump the data. Call this single function to uniformize the output format
//

#include <iostream>
#include <fstream>
#include <iomanip>

#include "ant.hpp"

char randstate[256];

int AntWalker::antdump()
{
  std::ofstream output(input.datafile, std::ios::out);
  output << "count: " << count << "\n\n2^n\tMSD\n";
  for(int i=0; i<arraylen; ++i)
  {
    output << i << "\t" << std::setprecision(10) << msdarray[i]/count << "\n";
  }
  return 0;
}


//
//  ant_main.cpp
//  AntLatticeWalk
//
//  Created by Yi Hu on 9/30/18.
//

#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "read_input.h"
#include "ant_hbd.hpp"

using namespace std;
// Main function
int main(int argc, const char * argv[])
{
  const char *inputf, *dftinputf = "inputfile_ant.dat";
  long stepview, nextview, nowclu;
#ifdef DEBUG
  clock_t t1, t2;
  t1 = clock();
#endif
  if (argc < 2)
  {
    inputf = dftinputf;
  }
  else
  {
    inputf = argv[1];
  }
  read_input_ant input;
  int error = input.read(inputf);
  if (error) return error;
  HybridAntWalker antwalker(input);
  stepview = input.maxcluster / 50;
  nextview = stepview;
  cout << "Start calculation..." << endl;
  for(nowclu=0; nowclu<input.maxcluster; ++nowclu)
  {
    antwalker.antrun();
    if(nowclu >= nextview)
    {
      cout << antwalker.getCountStr() << endl;
      antwalker.antdump();
      nextview += stepview;
    }
  }
  cout << antwalker.getCountStr() << "\nEnd of calculation" << endl;
  antwalker.antdump();
#ifdef DEBUG
  t2 = clock();
  cout << (float)(t2-t1)/CLOCKS_PER_SEC << endl;
#endif
}

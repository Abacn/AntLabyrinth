//
//  main.cpp
//  MSDLimit
//
//  Created by Yi Hu on 5/3/18.
//

#include <iostream>
#include "cluster_counter.hpp"

int main(int argc, const char * argv[]) {
  const char *inputf, *dftinputf = "inputfile.dat";
  if (argc < 2)
  {
    inputf = dftinputf;
  }
  else
  {
    inputf = argv[1];
  }
  read_input_lattice input;
  int error = input.read(inputf);
  if (error) return error;
  
  ClusterCounter clcounter(input);
  for(int i=1; i<=input.repeat; ++i)
  {
    clcounter.clugen();
#if DIM<4
    if(i%100==0)
#else
      if(i%10==0)
#endif
      {
        std::cout << "("<< i << "/" << input.repeat << ") " << clcounter.getMSD() << std::endl;
      }
  }
  clcounter.cludump();
  return 0;
}

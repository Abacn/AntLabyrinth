#include <iostream>
#include <time.h>
#include "leath_counter.hpp"

// main function for the Leath algorithm to generate (finite) clusters

int main(int argc, const char * argv[]) {
  const char *inputf, *dftinputf = "inputfile_leath.dat";
  long stepview, nextview, nowsites, nowclu;
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
  read_input_leath input;
  int error = input.read(inputf);
  if (error) return error;
  
  LeathCounter clcounter(input);
  stepview = input.maxcluster / 50;
  nextview = stepview;
  std::cout << "Start calculation...\nNow (sites, clusters)" << std::endl;
  for(nowsites = nowclu = 0;
      nowsites < input.maxsite || nowclu < input.maxcluster;
      nowsites = clcounter.getCountSites(), nowclu = clcounter.getCount())
  {
    if(nowclu >= nextview)
    {
      std::cout << "(" << nowsites << " , " << nowclu << ") MSD: " <<  clcounter.getMSD() << std::endl;
      nextview = nowclu+stepview;
    }
    clcounter.clugen();
  }
  std::cout << "(" << nowsites << " , " << nowclu << ") MSD: " <<  clcounter.getMSD() << std::endl;
  std::cout << "End of calculation" << std::endl;
#ifdef DEBUG
  t2 = clock();
  std::cout << (float)(t2-t1)/CLOCKS_PER_SEC << std::endl;
#endif
  clcounter.cludump();
  return 0;
}

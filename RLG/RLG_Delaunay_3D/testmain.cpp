//
//  main.cpp
//  RLG_Delaunay
//

#include <iostream>
#include "read_input.hpp"
#include "graph_construct_3d.hpp"

// Main function
int main(int argc, const char * argv[])
{
  const char *inputf, *dftinputf = "input_RLG_Delaunay_3D.dat";
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
  ReadInputRLG input;
  int error = input.read(inputf);
  if (error) return error;
  std::cout << "Start calculation..." << std::endl;
  
  Graph_Construct graph(input);
  graph.graphtest(); 
  // graph.graphdump();
  
#ifdef DEBUG
  t2 = clock();
  std::cout << (float)(t2-t1)/CLOCKS_PER_SEC << std::endl;
#endif
}

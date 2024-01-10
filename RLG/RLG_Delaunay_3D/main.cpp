//
//  main.cpp
//  RLG_Delaunay_3D
//  Calculate the void percolation threshold of RLG
//  in 3D periodic box

#include <iostream>
#include "read_input.hpp"
#include "graph_construct_3d.hpp"

// Main function
int main(int argc, const char * argv[])
{
  const char *inputf, *dftinputf = "input_RLG_Delaunay_3D.dat";
  long rp, stepview, nextview;
//#ifdef DEBUG
  clock_t t1, t2;
  t1 = clock();
//#endif
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
  stepview = input.repeatrun / 50;
  nextview = stepview;
  Graph_Construct_3d graph(input);
  for(rp = 0; rp<input.repeatrun; ++rp)
  {
    graph.graphrun();
    if(rp >= nextview)
    {
      std::cout << "(" << rp <<  ")" << std::endl;
      nextview = rp+stepview;
    }
  }
  std::cout << "(" << rp <<  ")" << std::endl;
//#ifdef DEBUG
  t2 = clock();
  std::cout << "Time per run (s):" << (float)(t2-t1)/CLOCKS_PER_SEC/input.repeatrun << std::endl;
//#endif
}

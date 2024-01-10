//
//  main.cpp
//  D4Percolation
//  Compute the void percolation on D4/E8 honeycomb
//
//  Created by Yi Hu on 5/30/19.
//

#include <cstdio>
#include <iostream>
#include <fstream>
#include <chrono>
#include "commondef.hpp"
#include "read_input.hpp"
#include "graph_construct.hpp"

/*
 #include "myutility.hpp"
 // test main
 int main(int argc, const char * argv[])
 {
 #if DIM==4
 point pa({-0.1, -0.1, 0.1, 0.1});
 #elif DIM==8
 point pa({0.1, 0.1, 0.8, 1.3, 2.2, -0.6, -0.7, 0.9});
 #endif
 bool result = MyUtility::pbc(pa);
 std::cout << result << " ";
 for(int rp=0; rp<DIM; ++rp)
 std::cout << pa[rp] << " ";
 std::cout << std::endl;
 }
 */

// debug working function
int debug_main(ReadInputRLG&);

// Main function
int main(int argc, const char * argv[])
{
  const char *inputf, *dftinputf = "input_RLG_hcombhull.dat";
  long rp, stepview, nextview;
  int run_ret_flag;
  auto t1 = std::chrono::system_clock::now();
#ifdef DEBUG
  qhulltime = execqhulltime = calcinvtime = 0.0;
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
  if(input.outmode == OUTPUT_MODE_DEBUG)
  {
    return debug_main(input);
  }
  std::cout << "Start calculation..." << std::endl;
  stepview = input.repeatrun / 50;
  nextview = stepview;
  Graph_Construct graph(input);
  //#ifdef _OPENMP
  //  omp_set_num_threads((int)ceil(input.n_thread*0.5));
  //#endif
  if(input.outmode != OUTPUT_MODE_SAVEGRAPH)
  {
    for(rp = 0; rp<input.repeatrun; ++rp)
    {
      if(rp >= nextview)
      {
        std::cout << "(" << rp <<  ")" << std::endl;
        nextview = rp+stepview;
      }
      run_ret_flag = graph.graphrun();
      if(run_ret_flag < 0)
      {
        run_ret_flag = graph.graphrun(true); // run failure, re-run
        // still failure, optflag disable re-run
        if(run_ret_flag < 0)
        {
          run_ret_flag = graph.graphrun(true);
        }
      }
      // check stop signal
      if (FILE *f = fopen("perco.stop", "r"))
      {
        fclose(f);
        remove("perco.stop");
        break;
      }
    }
    std::cout << "(" << rp <<  ")" << std::endl;
    auto t2 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t2-t1;
    std::cout << "Time per run: " << elapsed_seconds.count()/input.repeatrun << " s\n";
#ifdef DEBUG
    std::cout << "  calcinv: " << calcinvtime/input.repeatrun << " s\n"
#ifdef _MULTITHREAD
    << "execqhull: " << (execqhulltime)/input.repeatrun << " s\n"
#else
    << "execqhull: " << (execqhulltime-qhulltime)/input.repeatrun << " s\n"
#endif
    << "    qhull: " << qhulltime/input.repeatrun << " s" << std::endl;
#endif
  }
  else
  {
    // save graph mode, take second parameter as starting N
    if(argc<3) exit(-1);
    int startx = std::atoi(argv[2]);
    graph.partialrun(startx);
    auto t2 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t2-t1;
    double tpr = elapsed_seconds.count();
    if(tpr>0.1)
    {
      std::cout << "Time per run: " << elapsed_seconds.count() << " s\n";
#ifdef DEBUG
      std::cout << "  calcinv: " << calcinvtime << " s\n"
#ifdef _MULTITHREAD
      << "execqhull: " << (execqhulltime) << " s\n"
#else
      << "execqhull: " << (execqhulltime-qhulltime) << " s\n"
#endif
      << "    qhull: " << qhulltime << " s" << std::endl;
#endif
    }
  }
  return 0;
}

int debug_main(ReadInputRLG &input)
{
  Graph_Construct graph(input);
  graph.graphdebug();
  return 0;
}

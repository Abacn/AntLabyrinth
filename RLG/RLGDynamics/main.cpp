//
//  main.cpp
//  QhullPercolation
//
//  Created by Yi Hu on 3/19/19.
//

#include <cstdio>
#include <iostream>
#include "commondef.hpp"
#include "myutility.hpp"
#include "read_input.hpp"
#include "box.hpp"
#include "hcomb.hpp"
#include "sphere.hpp"

// Main function
int main(int argc, const char * argv[])
{
  const char *inputf, *dftinputf = "input_RLG_dynamics.dat";
  long rp, stepview, nextview;
  int run_ret_flag;
  clock_type t1, t2;
  t1 = std::chrono::steady_clock::now();
#ifdef DEBUG
  stattime = collidetime = intgtime = nlisttime = duration_type::zero();
#endif
  if (argc < 2)
  {
    inputf = dftinputf;
  }
  else
  {
    inputf = argv[1];
  }
  ReadInputDy input;
  int error = input.read(inputf);
  if (error) return error;
  std::cout << "Start calculation..." << std::endl;
  stepview = input.repeatrun / 50;
  nextview = stepview;
  MyUtility::init(input.seed);
  BaseBox *boxhandle = nullptr;
  if(input.systype == SYSTEM_TYPE_BOX)
  {
    boxhandle = new Box(input);
  }
  else if(input.systype == SYSTEM_TYPE_HCOMB)
  {
    boxhandle = new HComb(input);
  }
  else if(input.systype == SYSTEM_TYPE_SPHERE)
  {
    boxhandle = new Sphere(input);
  }
  else
  {
    std::cerr << "Unknown system type" << std::endl;
    exit(-1);
  }
  for(rp = 0; rp<input.repeatrun; ++rp)
  {
    run_ret_flag = boxhandle->process();
    if(rp >= nextview)
    {
      std::cout << "(" << rp <<  ")" << std::endl;
      nextview = rp+stepview;
      boxhandle->dump();
    }
    // check stop signal
    if (FILE *f = fopen("rlgdynamics.stop", "r"))
    {
      fclose(f);
      remove("rlgdynamics.stop");
      break;
    }
  }
  boxhandle->dump();
  delete boxhandle;
  std::cout << "(" << rp <<  ")" << std::endl;
  t2 = std::chrono::steady_clock::now();
  std::cout << "Time per run:" << static_cast<duration_type>(t2-t1).count()/input.repeatrun << " s\n";
#ifdef DEBUG
  std::cout << "Update nlist: " << nlisttime.count()/input.repeatrun << " s\n"
  << " Searching collision: " << collidetime.count()/input.repeatrun << " s\n"
  << " Integration: " << intgtime.count()/input.repeatrun << " s\n"
  << " Statistics: " << stattime.count()/input.repeatrun << " s\n"
  << std::endl;
#endif
  return 0;
}

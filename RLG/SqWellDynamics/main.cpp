//
//  main.cpp
//  QhullPercolation
//
//  Created by Yi Hu on 3/19/19.
//

#include <iostream>

#include "commondef.hpp"
#include "myutility.hpp"
#include "sqwellbox.hpp"

// Main function
int main(int argc, const char * argv[])
{
  const char *inputf, *dftinputf = "input_RLG_dynamics_sqwell.dat";
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
  ReadInputDySq input;
  int error = input.read(inputf);
  if (error) return error;
  std::cout << "Start calculation..." << std::endl;
  stepview = input.repeatrun / 50;
  nextview = stepview;
  MyUtility::init(input.seed);
  SqwellBox *boxhandle = new SqwellBox(input);
  for(rp = 0; rp<input.repeatrun; ++rp)
  {
    run_ret_flag = boxhandle->process();
    if(rp >= nextview)
    {
      std::cout << "(" << rp <<  ")" << std::endl;
      nextview = rp+stepview;
      boxhandle->dump();
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

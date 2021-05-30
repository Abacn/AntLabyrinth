//
//  main.cpp
//  DnBondAnimal
//
//  Created by Yi Hu on 1/6/21.
//

#include <iostream>
#include "animal_bond.hpp"

//using clock_type = std::chrono::time_point<std::chrono::steady_clock>;
//using duration_type = std::chrono::duration<double>;

int main(int argc, const char * argv[]) {
  // insert code here...
  if(argc < 3)
  {
    std::cout << "Usage: DnBondAnimal D N" << std::endl;
    return 1;
  }
  int D, N;
  D = atoi(argv[1]);
  N = atoi(argv[2]);
  std::vector<unsigned int> preassignedidxs;
  for(int rp=3; rp<argc; ++rp)
  {
    preassignedidxs.push_back(atoi(argv[rp]));
  }
  std::cout << "DIM: " << D << ", " << "N: " << N;
  if(preassignedidxs.size() > 0)
  {
    std::cout << ", Firstidxs: ";
    for(auto rp: preassignedidxs)
    {
      std::cout << rp << ' ';
    }
  }
  std::cout << std::endl;
  //clock_type t1 = std::chrono::steady_clock::now();
#define RUND(DIM, N) case(DIM):{BondNeighborDn<DIM> dn; CountBondAnimal<DIM> an(dn); auto result=an.getCount(N, preassignedidxs); printResult(result);} break;
  switch(D)
  {
      RUND(2, N)
      RUND(3, N)
      RUND(4, N)
      RUND(5, N)
      RUND(6, N)
      RUND(7, N)
      RUND(8, N)
      RUND(9, N)
      RUND(10, N)
      RUND(11, N)
      RUND(12, N)
      RUND(13, N)
      RUND(14, N)
      RUND(15, N)
      RUND(16, N)
  }
#undef RUND
  //clock_type t2 = std::chrono::steady_clock::now();
  //float time_consumed = static_cast<duration_type>(t2-t1).count();
  //if(time_consumed>0.1f)
  //    std::cout << "Time consumed: " << time_consumed << " s" << std::endl;
  return 0;
}

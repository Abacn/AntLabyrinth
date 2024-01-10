//
//  main.cpp
//  DnAnimalPerimeter
//
//  Created by Yi Hu on 12/10/20.
//

#include <iostream>
#include "animal_site.hpp"

int main(int argc, const char * argv[]) {
  if(argc < 2)
  {
    std::cout << "Usage: DnLatticeAnimal D [filename]" << std::endl;
    return 1;
  }
  const char *fname, *dftname = "animal.txt";
  if(argc==3) fname = argv[2];
  else fname = dftname;
  int D;
  D = atoi(argv[1]);
  std::cout << "DIM: " << D << std::endl;
#define RUND(DIM) case(DIM):{PointNeighborDn<DIM> dn; LatticeAnimal<DIM> an; an.setAnimal(fname); int perimeter=getPerimeter(an, dn); std::cout << perimeter << std::endl;} break;
  switch(D)
  {
      RUND(2)
      RUND(3)
      RUND(4)
      RUND(5)
      RUND(6)
      RUND(7)
      RUND(8)
      RUND(9)
      RUND(10)
      RUND(11)
      RUND(12)
      RUND(13)
      RUND(14)
  }
#undef RUND
  return 0;
}

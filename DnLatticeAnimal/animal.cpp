//
//  animal.cpp
//  DnLatticeAnimal
//
//  Created by Yi Hu on 12/10/20.
//

#include <vector>
#include <utility>
#include <iostream>

void printResult(std::vector<std::pair<int, uint64_t> > const &ilist)
{
  for(auto const &hdl: ilist) std::cout << hdl.first << "\t" << hdl.second << std::endl;
}

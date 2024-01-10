//
//  read_input.hpp
//  RLG_Delaunay
//
//  Created by Yi Hu on 8/6/18.


#ifndef read_input_hpp
#define read_input_hpp
#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#define NAME_LEN 256

enum Output_mode{
  OUTPUT_MODE_SUMMARY,
  OUTPUT_MODE_PRODUCT,
  OUTPUT_MODE_SAVEGRAPH,
  OUTPUT_MODE_DEBUG,
  OUTPUT_MODE_ALL
};

class ReadInput {
public:
  virtual int read(const char *inputf) = 0;
};

class ReadInputRLG: public ReadInput {
public:
  // member method
  int read(const char *inputf);
  // member variable
  int repeatrun;
  int32_t N;
  unsigned int seed;
  std::string datafile;
  enum Output_mode outmode;
  bool defaultoptflag;
  double rcut;
  int32_t n_thread;
};
#endif /* read_input_hpp */

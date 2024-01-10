//
//  read_input.hpp
//  RLG_Delaunay
//
//  Created by Yi Hu on 8/6/18.


#ifndef read_input_hpp
#define read_input_hpp
#include <cstdint>
#include <string>

#include "basereadinput.hpp"
#define NAME_LEN 256

/** Read parameters input file. */
class ReadInputDy: public ReadInput {
public:
  // member method
  int read(const char *inputf);
  // member variable
  System_type systype;
  int repeatrun;
  uint32_t N;
  double phi, rmax;
  int sleft, texp, sinterval, tscaleflag;
  unsigned int seed;
  std::string datafile, logdatafile;
  enum Output_mode outmode;
  int rcuttype;
  double rcut;
private:
  // legacy reader. Read from .dat file
  int read_dat(const char *inputf);
  // ini reader. Read from .ini file
  int read_ini(const char *inputf);
};
#endif /* read_input_hpp */

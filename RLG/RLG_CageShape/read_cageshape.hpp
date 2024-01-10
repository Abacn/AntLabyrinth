//
//  read_input_shell.hpp
//  RLG_Delaunay
//
//  Created by Yi Hu on 3/19/19.
//

#ifndef read_input_shell_hpp
#define read_input_shell_hpp

enum Output_mode{
  OUTPUT_MODE_SUMMARY,
  OUTPUT_MODE_PRODUCT,
  OUTPUT_MODE_SAVESPL,
  OUTPUT_MODE_DEBUG,
  OUTPUT_MODE_ALL
};

class ReadInputRLGCShape{
public:
  // member method
  int read(const char *inputf);
  // member variable
  int repeatrun;
  double phi;
  double N;
  double rmax, rconstA, rconstB, Vunitsphere;
  unsigned int seed;
  std::string datafile;
  enum Output_mode outmode;
  int dostat; // 0: do not do stat on-the-fly; 1: do stat; 2: read samplefilelist.dat and do stats
  bool dostats[4];
  long samplestrat[4];
  double binstart, logbinstart;
  int nbinperl, nbinperlogl;
};

#endif /* read_input_shell_hpp */

#ifndef READ_INPUT_H
#define READ_INPUT_H

#define NAME_LEN 256

class read_input {
public:
  virtual int read(const char *inputf) = 0;
};

class read_input_lattice : public read_input {

public:

  int L;                         // # spheres
  double pf;                     // initial packing fraction
  double maxtime;                // deprecated in msdlim
  int seed;
  double tprint;
  char datafile[NAME_LEN];       // file to write statistics
  int read(const char *inputf);
  long repeat;
};

class read_input_leath : public read_input {
public:
  int L;                         // # deprecated in Leath algorithm
  double pf;                     // initial packing fraction
  long maxcluster;               // max number of cluster generated
  int seed;
  double tprint;
  char datafile[NAME_LEN];       // file to write statistics
  int read(const char *inputf);
  long maxsite;                  // max number of sites generated
};

class read_input_ant : public read_input {
public:
  // und and ord: useless; hybrid: below L use transfer matrix method
  int L;                         
  double pf;                     // initial packing fraction
  long maxcluster;               // max number of cluster generated
  int seed;
  double tprint;
  char datafile[NAME_LEN];       // file to write statistics
  int read(const char *inputf);
  int arraylen;                  // length of array that records the MSD
                                 // from time 2^1, ..., 2^n
};
#endif

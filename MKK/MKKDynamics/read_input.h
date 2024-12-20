#ifndef READ_INPUT_H
#define READ_INPUT_H

#define NAME_LEN 1024


class read_input {

public:

  double eventspercycle;          // # events per particle per cycle
  int N;                          // # spheres
  double polydispersity;          // polydispersity
  double initialpf;               // initial packing fraction
  double maxpf;            // maximum packing fraction
  double temp;                      // initial temperature (temp=0 means v=0)
  double growthrate;
  double rfactor;
  double maxpressure;
  double mintime;
  double maxtime;
  double maxcollisionrate;
  double shiftscale;                   // random shift scale
  long seed;
  double tprint;
  double t_dist_start;          // smallest time interval if run time dependent displacement distribution
  char readfile[NAME_LEN];     // file with configuration; if new, creates new
  char writefile[NAME_LEN];    // file to write configuration
  char configprefixes[NAME_LEN]; // file to write the configurations
  char datafile[NAME_LEN];       // file to write statistics
  int outconfig_switch;        // save configuration or not
  // calculate particle MSD/MQD or not
  // 0 - disabled
  // 1 - enable particle MSD/MQD
  // 2 - enabled 1 + distribution of MSD/MQD
  // 3 - preserved
  // a bit at 4 - enable viscosity
  int calcpmsd;
  int read(const char* fname);
};


#endif

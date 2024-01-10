//
//  msdrecorder.hpp
//  RLGDynamics
//
//  Created by Yi Hu on 8/6/19.
//

#ifndef msdrecorder_hpp
#define msdrecorder_hpp

#include "commondef.hpp"
#include "limited_queue.hpp"
#include "basebox.hpp"

/**
 * MSD over time base class. record method is abstract and left implementation
 * in derived class.
 */
class MSDRecorder{
protected:
  const int nmult; // 2^10 samples
  const long nmultspl;
  long sn, sinterval;          // sn: number of 2^n points
  long sqsumcount;
  double *tpref;               // prefactor term of t: t = tpref*2^nextt;
  double **r2s, **r4s, **rcov; // <Delta>, <Delta^2>, covariance <Delta_1 Delta_2>
  double **logr2s, **logr4s;
  double sqsums, qqsums;       // square displacement from origin
  double lefttime;
  double lastsd;               // last recorded tracer square displacement
  long **counts;
public:
  MSDRecorder(double lefttime_, int sinterval_, int span, int nmult_=10);
  virtual int record(const Tracer &tracer, double tnow, double tnext) = 0;
  void dump(double *r2s_, double *r4s_, double *chidis_, double *rcov_, int *samplecount_, int sz);
  void dumplog(double *r2s_, double *r4s_, int sz);
  void dumpone(const char* fname);
  // get final square displacement to the origin
  double getfinalsd() const;
  // get square and quadruple displacement with final positions averaged
  void getfinalmsd(double *sqs, double *qqs) const;
  // get square and quadruple displacement with both initial and final positions averaged
  // returns: 0-plateaued; 1-going up; 2-going down; 3-not enough data to average
  int getinitfinalmsd(double *sqs, double *qqs) const;
  virtual ~MSDRecorder();
};

/** Record MSD averaged by tracer initial and final location. */
class AveMSDRecorder: public MSDRecorder{
private:
  long spart1, spart2;
  long *smark;                  // part 1 array position mark
  double nextrec, *nextrecs;
  long *inextrecs;
  Vector **snapshotsarr;
  limited_queue<Vector> **snapshotsque;
  long **nextarr;
public:
  AveMSDRecorder(double lefttime_, int sinterval_, int span, int nmult_=10);
  virtual int record(const Tracer &tracer, double tnow, double tnext) override;
  virtual ~AveMSDRecorder();
};

/** Record MSD velocity scaled in gaussian for each segment. */
class GaussianMSDRecorder: public MSDRecorder{
private:
  // marking index of next record
  long *smark;
  // time at next record
  double *nextrecs;
  // current assigned velocities
  double *velocities;
  // accumulated reassigned time elapse when last velocity update
  double *accutimes;
  // traj time stamp at last velocity update
  double *lastupdtimes;
public:
  GaussianMSDRecorder(double lefttime_, int sinterval_, int span, int nmult_=10);
  virtual int record(const Tracer &tracer, double tnow, double tnext) override;
  int reassignvelocity(double timestamp);
  bool samplefinished();
  virtual ~GaussianMSDRecorder();
};

/** Record first escape times. Used in sphere geometry. */
class EscapeRecorder{
private:
  double Dcutleft, Dcutright, dDelta, DeltaNext;
  double *esctimes;
  int nesctime, crtesc;
  int nextfntmsdidx;
  int escapecountlimit;
  static const int nextdumps[]; // dump finite size shell
public:
  EscapeRecorder(double Dcutleft_, double Dcutright_, double dDelta_, int escapecountlimit_);
  int record(const Tracer &tracer, double tnow, double t_elapse);
  void dump(int count);
  ~EscapeRecorder();
};
#endif /* msdrecorder_hpp */

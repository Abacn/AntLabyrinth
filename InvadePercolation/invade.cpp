#include <string.h>
#include <math.h>
#include <set>
#include <queue>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "leathsitenode.hpp"
#include "invade.hpp"

/* Invade percolation
 * efficient way to calculate percolation threshold
 */

InvadePclt::InvadePclt(const read_input_ant &inp) :
input(inp), posmin{1, -1}
{
  count=0;
  pcarray = new double[input.arraylen];
  sqpcarray = new double[input.arraylen];
  tmparrcount = new long[input.arraylen];
  for(int i=0; i<DIM; ++i)
  {
    unitvc[i][i] = 1;
  }
  for(int i=0; i<input.arraylen; ++i)
  {
    tmparrcount[i] = (long)(pow(2, i/4.0)*100L);
  }
  maxtime = tmparrcount[input.arraylen-1];
#ifdef DEBUG
  std::cout << "Maxcluster size: " << maxtime << std::endl;
#endif
  arraylen = input.arraylen;
  memset(pcarray, 0, sizeof(double)*input.arraylen);
  memset(sqpcarray, 0, sizeof(double)*input.arraylen);
  // srandom(input.seed);
#if APPLY_RAND_TYPE==1
  rg.seed(input.seed);
#else
  initstate(input.seed, randstate, 256);
#endif
}

InvadePclt::~InvadePclt()
{
  delete [] pcarray;
  delete [] sqpcarray;
  delete [] tmparrcount;
}

InvadePclt_Zn::InvadePclt_Zn(const read_input_ant &inp): InvadePclt(inp) {}

InvadePclt_Dn::InvadePclt_Dn(const read_input_ant &inp): InvadePclt(inp) {}

InvadePclt_DS::InvadePclt_DS(const read_input_ant &inp): InvadePclt(inp) {}

InvadePclt_Znb::InvadePclt_Znb(const read_input_ant &inp): InvadePclt(inp) {}

InvadePclt_Dnb::InvadePclt_Dnb(const read_input_ant &inp): InvadePclt(inp)
{
  int nowidx = 0;
  for(int rp=0; rp<DIM-1; ++rp)
    for(int rq=rp+1; rq<DIM; ++rq)
    {
      dnbase[nowidx][rp] = dnbase[nowidx][rq] = 1;
      ++nowidx;
      dnbase[nowidx][rp] = 1; dnbase[nowidx][rq] = -1;
      ++nowidx;
    }
}

InvadePclt_DSb::InvadePclt_DSb(const read_input_ant &inp): InvadePclt(inp) {}

InvadePclt_Zn::~InvadePclt_Zn() {}

InvadePclt_Dn::~InvadePclt_Dn() {}

InvadePclt_DS::~InvadePclt_DS() {}

InvadePclt_Znb::~InvadePclt_Znb() {}

InvadePclt_Dnb::~InvadePclt_Dnb() {}

InvadePclt_DSb::~InvadePclt_DSb() {}

int InvadePclt_Zn::ivdrun()
{
  long nextrectime = 0;
  long csize=0;
  double weight, dtmp;
  int rp;
  LeathSiteNode posvc;
  std::set<LeathSiteNode> siterec = { {zerovc} };
  std::pair<std::set<LeathSiteNode>::iterator,bool> retval;
  std::priority_queue<PQueueNode> mypque;
  mypque.push(PQueueNode(siterec.begin(), 0.0));
  
  for(csize = 0; csize<=maxtime;)
  {
    const LeathSiteNode &popNode = mypque.top().getSiteNode();
    // pop the front element
    mypque.pop();
    // scan nearest neighbors
    for(rp=0; rp<DIM; ++rp)
    {
      // +1
      posvc = popNode + unitvc[rp];
      retval = siterec.insert(posvc);
      if(retval.second)
      {
        // Does not have the key
#if APPLY_RAND_TYPE==1
        weight = rrand(rg);
#else
        weight = (double)random()/(double)RAND_MAX;
#endif
        mypque.push(PQueueNode(retval.first, weight));
      }
      // -1
      posvc = popNode - unitvc[rp];
      retval = siterec.insert(posvc);
      if(retval.second)
      {
        // Does not have the key
#if APPLY_RAND_TYPE==1
        weight = rrand(rg);
#else
        weight = (double)random()/(double)RAND_MAX;
#endif
        mypque.push(PQueueNode(retval.first, weight));
      }
    }
    ++csize;
    // cout << csize << "\t" << siterec.size() << "\t" << double(csize)/siterec.size() << endl;
    if(csize == tmparrcount[nextrectime])
    {
      dtmp = (double)csize / siterec.size();
      pcarray[nextrectime] += dtmp;
      sqpcarray[nextrectime] += dtmp*dtmp;
      ++nextrectime;
    }
  }

  count++;
  return 0;
}

int InvadePclt_Dn::ivdrun()
{
  long nextrectime = 0;
  long csize=0;
  double weight, dtmp;
  int rp, rq;
  LeathSiteNode posvc;
  std::set<LeathSiteNode> siterec = { {zerovc} };
  std::pair<std::set<LeathSiteNode>::iterator,bool> retval;
  std::priority_queue<PQueueNode> mypque;
  mypque.push(PQueueNode(siterec.begin(), 0.0));
  
  for(csize = 0; csize<=maxtime;)
  {
    const LeathSiteNode &popNode = mypque.top().getSiteNode();
    // pop the front element
    mypque.pop();
    // scan nearest neighbors
    for(rp=0; rp<DIM; ++rp)
    {
      for(rq=rp+1; rq<DIM; ++rq)
      {
        for(int rr=0; rr<4;++rr)
        {
          posvc = popNode;
          if(rr&1) ++posvc[rp];
          else --posvc[rp];
          if(rr&2) ++posvc[rq];
          else --posvc[rq];
          retval = siterec.insert(posvc);
          if(retval.second)
          {
            // Does not have the key
#if APPLY_RAND_TYPE==1
            weight = rrand(rg);
#else
            weight = (double)random()/(double)RAND_MAX;
#endif
            mypque.push(PQueueNode(retval.first, weight));
          }
        }
      }
    }
    ++csize;
    // cout << csize << "\t" << siterec.size() << "\t" << double(csize)/siterec.size() << endl;
    if(csize == tmparrcount[nextrectime])
    {
      dtmp = (double)csize / siterec.size();
      pcarray[nextrectime] += dtmp;
      sqpcarray[nextrectime] += dtmp*dtmp;
      ++nextrectime;
    }
  }
  count++;
  return 0;
}

int InvadePclt_DS::ivdrun()
{
#if DIM>=6 && DIM<=9
  long nextrectime = 0;
  long csize=0;
  double weight, dtmp;
  int rp, rq;

  extern const int size_densebase;
  extern const SiteNode densebase[];
  SiteNode zerovec_D, posvc;
  std::set<SiteNode> siterec = { {zerovec_D} };
  std::pair<std::set<SiteNode>::iterator,bool> retval;
  std::priority_queue<PQueueNode > mypque;
  mypque.push(PQueueNode(siterec.begin(), 0.0));
  
  for(csize = 0; csize<=maxtime;)
  {
    const SiteNode &popNode = mypque.top().getSiteNode();
    // pop the front element
    mypque.pop();
    // scan nearest neighbors
    for(rp=0; rp<size_densebase; ++rp)
    {
      // +1
      posvc = popNode;
      for(rq=0; rq<ndim; ++rq)
      {
        posvc[rq] += densebase[rp][rq];
      }
      retval = siterec.insert(posvc);
      if(retval.second)
      {
        // Does not have the key
#if APPLY_RAND_TYPE==1
        weight = rrand(rg);
#else
        weight = (double)random()/(double)RAND_MAX;
#endif
        mypque.push(PQueueNode(retval.first, weight));
      }
      // -1
      posvc = popNode;
      for(rq=0; rq<ndim; ++rq)
      {
        posvc[rq] -= densebase[rp][rq];
      }
      retval = siterec.insert(posvc);
      if(retval.second)
      {
        // Does not have the key
#if APPLY_RAND_TYPE==1
        weight = rrand(rg);
#else
        weight = (double)random()/(double)RAND_MAX;
#endif
        mypque.push(PQueueNode(retval.first, weight));
      }
    }
    ++csize;
    // cout << csize << "\t" << siterec.size() << "\t" << double(csize)/siterec.size() << endl;
    if(csize == tmparrcount[nextrectime])
    {
      dtmp = (double)csize / siterec.size();
      pcarray[nextrectime] += dtmp;
      sqpcarray[nextrectime] += dtmp*dtmp;
      ++nextrectime;
    }
  }
  count++;
  return 0;
#else
  return -1;
#endif
}

int InvadePclt_Znb::ivdrun()
{
  long nextrectime = 0;
  long csize=0;
  double weight, dtmp;
  LeathSiteNode posvc;
  std::set<LeathSiteNode> siterec = { {zerovc} };
  std::set<LeathBondNode> bondrec;
  std::pair<std::set<LeathBondNode>::iterator,bool> retval;
  std::priority_queue<BQueueNode> mypque;
  auto zeroneighbors = getNeighbors(zerovc); // bond start from origin
  for(auto const &nd: zeroneighbors)
  {
    auto retval = bondrec.insert(nd);
#if APPLY_RAND_TYPE==1
        weight = rrand(rg);
#else
        weight = (double)random()/(double)RAND_MAX;
#endif
    mypque.push(BQueueNode(retval.first, weight));
  }
  // now add sites
  for(csize = 0; csize<=maxtime;)
  {
    const LeathBondNode &popNode = mypque.top().getSiteNode();
    // pop the front element
    mypque.pop();
    // scan nearest neighbors

    // edge start site
    posvc = popNode.first;
    bool addflag = true;
    if(siterec.find(posvc) != siterec.end())
    {
      ++posvc[popNode.second];
      if(siterec.find(posvc) != siterec.end())
      {
        // both side added
        addflag = false;
      }
    }
    if(addflag)
    {
      // add new sites
      auto zeroneighbors = getNeighbors(posvc); // bond start from origin
      for(auto const &nd: zeroneighbors)
      {
        auto retval = bondrec.insert(nd);
        if(retval.second)
        {
#if APPLY_RAND_TYPE==1
          weight = rrand(rg);
#else
          weight = (double)random()/(double)RAND_MAX;
#endif
          mypque.push(BQueueNode(retval.first, weight));
        }
      }
      siterec.insert(posvc);
    }
    ++csize;
    // cout << csize << "\t" << siterec.size() << "\t" << double(csize)/siterec.size() << endl;
    if(csize == tmparrcount[nextrectime])
    {
      dtmp = (double)csize / bondrec.size();
      pcarray[nextrectime] += dtmp;
      sqpcarray[nextrectime] += dtmp*dtmp;
      ++nextrectime;
    }
  }
  count++;
  return 0;
}

int InvadePclt_Dnb::ivdrun()
{
  long nextrectime = 0;
  long csize=0;
  double weight, dtmp;
  LeathSiteNode posvc;
  std::set<LeathSiteNode> siterec = { {zerovc} };
  std::set<LeathBondNode> bondrec;
  std::pair<std::set<LeathBondNode>::iterator,bool> retval;
  std::priority_queue<BQueueNode> mypque;
  auto zeroneighbors = getNeighbors(zerovc); // bond start from origin
  for(auto const &nd: zeroneighbors)
  {
    auto retval = bondrec.insert(nd);
#if APPLY_RAND_TYPE==1
        weight = rrand(rg);
#else
        weight = (double)random()/(double)RAND_MAX;
#endif
    mypque.push(BQueueNode(retval.first, weight));
  }
  // now add sites
  for(csize = 0; csize<=maxtime;)
  {
    const LeathBondNode &popNode = mypque.top().getSiteNode();
    // pop the front element
    mypque.pop();
    // scan nearest neighbors

    // edge start site
    posvc = popNode.first;
    bool addflag = true;
    if(siterec.find(posvc) != siterec.end())
    {
      posvc += dnbase[popNode.second];
      if(siterec.find(posvc) != siterec.end())
      {
        // both side added
        addflag = false;
      }
    }
    if(addflag)
    {
      // add new sites
      auto zeroneighbors = getNeighbors(posvc); // bond start from origin
      for(auto const &nd: zeroneighbors)
      {
        auto retval = bondrec.insert(nd);
        if(retval.second)
        {
#if APPLY_RAND_TYPE==1
          weight = rrand(rg);
#else
          weight = (double)random()/(double)RAND_MAX;
#endif
          mypque.push(BQueueNode(retval.first, weight));
        }
      }
      siterec.insert(posvc);
    }
    ++csize;
    // cout << csize << "\t" << siterec.size() << "\t" << double(csize)/siterec.size() << endl;
    if(csize == tmparrcount[nextrectime])
    {
      dtmp = (double)csize / bondrec.size();
      pcarray[nextrectime] += dtmp;
      sqpcarray[nextrectime] += dtmp*dtmp;
      ++nextrectime;
    }
  }
  count++;
  return 0;
}

int InvadePclt_DSb::ivdrun()
{
#if DIM>=6 && DIM<=9
  long nextrectime = 0;
  long csize=0;
  double weight, dtmp;
  
  extern const int size_densebase;
  extern const SiteNode densebase[];
  SiteNode zerovec_D, posvc;
  std::set<SiteNode> siterec = { {zerovec_D} };
  std::set<BondNode> bondrec;
  std::pair<std::set<BondNode>::iterator,bool> retval;
  std::priority_queue<BQueueNode> mypque;
  auto zeroneighbors = getNeighbors(zerovec_D); // bond start from origin
  for(auto const &nd: zeroneighbors)
  {
    auto retval = bondrec.insert(nd);
#if APPLY_RAND_TYPE==1
        weight = rrand(rg);
#else
        weight = (double)random()/(double)RAND_MAX;
#endif
    mypque.push(BQueueNode(retval.first, weight));
  }
  // now add sites
  for(csize = 0; csize<=maxtime;)
  {
    const BondNode &popNode = mypque.top().getSiteNode();
    // pop the front element
    mypque.pop();
    // scan nearest neighbors

    // edge start site
    posvc = popNode.first;
    bool addflag = true;
    if(siterec.find(posvc) != siterec.end())
    {
      posvc += densebase[popNode.second];
      if(siterec.find(posvc) != siterec.end())
      {
        // both side added
        addflag = false;
      }
    }
    if(addflag)
    {
      // add new sites
      auto zeroneighbors = getNeighbors(posvc); // bond start from origin
      for(auto const &nd: zeroneighbors)
      {
        auto retval = bondrec.insert(nd);
        if(retval.second)
        {
#if APPLY_RAND_TYPE==1
          weight = rrand(rg);
#else
          weight = (double)random()/(double)RAND_MAX;
#endif
          mypque.push(BQueueNode(retval.first, weight));
        }
      }
      siterec.insert(posvc);
    }
    ++csize;
    // cout << csize << "\t" << siterec.size() << "\t" << double(csize)/siterec.size() << endl;
    if(csize == tmparrcount[nextrectime])
    {
      dtmp = (double)csize / bondrec.size();
      pcarray[nextrectime] += dtmp;
      sqpcarray[nextrectime] += dtmp*dtmp;
      ++nextrectime;
    }
  }
  count++;
  return 0;
#else
  return -1;
#endif
}

std::vector<LeathBondNode> InvadePclt_Znb::getNeighbors(LeathSiteNode const &node)
{
  std::vector<LeathBondNode> result(2*DIM);
  for(int rp=0; rp<DIM; ++rp)
  {
    result[rp].first = node;
    result[rp].second = rp;
    result[DIM+rp].first = node;
    --result[DIM+rp].first.x[rp];
    result[DIM+rp].second = rp;
  }
  return result;
}

std::vector<LeathBondNode> InvadePclt_Dnb::getNeighbors(LeathSiteNode const &node)
{
  std::vector<LeathBondNode> result(2*size_dnbase);
  for(int rp=0; rp<size_dnbase; ++rp)
  {
    result[rp].first = node;
    result[rp].second = rp;
    result[size_dnbase+rp].first = node - dnbase[rp];
    result[size_dnbase+rp].second = rp;
  }
  return result;
}

std::vector<InvadePclt_DSb::BondNode> InvadePclt_DSb::getNeighbors(SiteNode const &node)
{
#if DIM>=6 && DIM<=9
  extern const int size_densebase;
  extern const SiteNode densebase[];
  std::vector<BondNode> result(2*size_densebase);
  for(int rp=0; rp<size_densebase; ++rp)
  {
    result[rp].first = node;
    result[rp].second = rp;
    result[size_densebase+rp].first = node - densebase[rp];
    result[size_densebase+rp].second = rp;
  }
#else
  std::vector<BondNode> result;
#endif
  return result;
}

int InvadePclt::ivddump()
{
  std::ofstream output(input.datafile, std::ios::out);
  long rp;
  output << "count: " << count << "\n\nsize\tpc\tVar\n";
  for(int i=0; i<arraylen; ++i)
  {
    rp = tmparrcount[i];
    output << rp << "\t" << std::setprecision(10) << pcarray[i]/count << "\t" << (sqpcarray[i]-pcarray[i]*pcarray[i]/count)/count << "\n";
  }
  return 0;
}

// Main function
int main(int argc, const char * argv[])
{
  const char *inputf, *dftinputf = "inputfile_ivd.dat";
  long stepview, nextview, nowclu;
#ifdef DEBUG
  clock_t t1, t2;
  t1 = clock();
#endif
  if (argc < 2)
  {
    inputf = dftinputf;
  }
  else
  {
    inputf = argv[1];
  }
  
  read_input_ant input;
  int error = input.read(inputf);
  if (error) return error;
  InvadePclt *ivdwalker;
  switch(input.L)
  {
    case(0):
      ivdwalker = new InvadePclt_Zn(input);
      std::cout << "Hypercubic lattice site percolation." << std::endl;
      break;
    case(1):
      ivdwalker = new InvadePclt_Dn(input);
      std::cout << "Dn lattice site percolation." << std::endl;
      break;
    case(2): // dense packing
      std::cout << "Dense lattice site percolation." << std::endl;
      switch(DIM)
      {
        // 2D analytical solvable
        case(2): std::cout << "Triangular lattice analytical result available: p_c = 1/2" << std::endl; return 0;
        // D3 D4 D5. Note that in 3D there are two dense packings where threshold are pretty similar but not identical
        case(3): case(4): case(5): ivdwalker = new InvadePclt_Dn(input); break;
        case(6): case(7): case(8): case(9): ivdwalker = new InvadePclt_DS(input); break;
      }
      break;
    case(10):
      ivdwalker = new InvadePclt_Znb(input);
      std::cout << "Hypercubic lattice bond percolation." << std::endl;
      break;
    case(11):
      ivdwalker = new InvadePclt_Dnb(input);
      std::cout << "Dn lattice bond percolation." << std::endl;
      break;
    case(12):
      std::cout << "Dense lattice bond percolation." << std::endl;
      switch(DIM)
      {
        // 2D analytical solvable
        case(2): std::cout << "Triangular lattice analytical result available: p_c = 2 sin (pi/18)" << std::endl; return 0;
        // D3 D4 D5. Note that in 3D there are two dense packings where threshold are pretty similar but not identical
        case(3): case(4): case(5): ivdwalker = new InvadePclt_Dnb(input); break;
        case(6): case(7): case(8): case(9): ivdwalker = new InvadePclt_DSb(input); break;
      }
      break;
    default:
      std::cout << "Unknown PBC type: " << input.L << std::endl; return(-1);
  }
  stepview = input.maxcluster / 50;
  nextview = stepview;
  std::cout << "Start calculation..." << std::endl;
  for(nowclu=0; nowclu<input.maxcluster; ++nowclu)
  {
    ivdwalker->ivdrun();
    if(nowclu >= nextview)
    {
      std::cout << "(" << nowclu << ") " << std::endl;
      ivdwalker->ivddump();
      nextview += stepview;
    }
  }
  std::cout << "(" << nowclu << ") \nEnd of calculation" << std::endl;
  ivdwalker->ivddump();
  delete ivdwalker;
#ifdef DEBUG
  t2 = clock();
  std::cout << (float)(t2-t1)/CLOCKS_PER_SEC << std::endl;
#endif
}

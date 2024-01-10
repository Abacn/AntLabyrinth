//
//  main.cpp
//  MDSDumpReader
//
//  Created by Yi Hu on 7/11/18.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <iterator>
#include <cmath>

using namespace std;

// Extrapolating the MSD by a biased Ensemble average. Refer to:
// PRE 67, 036101 (2003)

int msdExtrapolator(const string &inFileStr)
{
  double pf, p0, msd, weight, sumw, summsd;
  long perimeter, csize, size_limit, nrec;
  int nfile=0, npf;
  char buf[100], c;
  vector<string> filelist;
  vector<double> pflist;
  ifstream lstFile(inFileStr, ios::in);
  if(!lstFile.is_open())
  {
    cout << "Open file " << inFileStr << " failed\n";
    return 2;
  }
  string line;
  getline(lstFile, line);
  if(!line.compare(0, 4, "[pf]"))
  {
    getline(lstFile, line); // skip
    istringstream iss(line);
    pflist = vector<double>{istream_iterator<double>(iss), istream_iterator<double>()};
    // read pf list
  }
  while(!lstFile.eof())
  {
    // Until find [file] mark
    getline(lstFile, line);
    if(!line.compare(0, 6, "[file]")) break;
  }
  while(!lstFile.eof())
  {
    getline(lstFile, line);
    if(line.empty()) continue;
    filelist.push_back(line);
    nfile++;
  }
  lstFile.close();
  if(filelist.empty())
  {
    nfile = 1;
    filelist.push_back("msddump.dat"); // default input file
  }
  npf = pflist.size();
  vector<double> msdlist(npf, 0.0);
  vector<double> weightlist(npf, 0.0);
  cout << "File\tAverageMSD\n";
  /* Data file format
   4D
   p:  0.196
   Site limit:  606060606
   size  parameter  msd
   12  65  4.979166667
   */
  for(vector<string>::iterator it = filelist.begin(); it != filelist.end(); ++it)
  {
    summsd = 0.0;
    nrec = 0; // number of record
    ifstream inFile(*it, ios::in);
    string line;
    inFile.getline(buf, 256);
    inFile.get(buf,100,':'); inFile.get(c); inFile >> p0;
    inFile.get(buf,100,':'); inFile.get(c); inFile >> size_limit;
    inFile.getline(buf, 256);  // skip "\n"
    inFile.getline(buf, 256);  // skip head
    // Get p
    vector<double> pfdp0(npf, 0.0); // p/p0
    vector<double> qfdq0(npf, 0.0); // (1-p)/(1-p0)
    for(int i=0; i<npf; ++i)
    {
      pfdp0[i] = log(pflist[i]/p0);
      qfdq0[i] = log((1-pflist[i])/(1-p0));
    }
    while(!inFile.eof())
    {
      inFile >> csize >> perimeter >> msd;
      if(0==csize)
        continue; /// eof meet
      for(int i=0; i<npf; ++i)
      {
        weight = exp(csize*pfdp0[i] + perimeter*qfdq0[i]); // use pow directly will cause overflow
        weightlist[i] += weight;
        msdlist[i] += weight*msd;
      }
      ++nrec;
      summsd += msd;
    }
    cout << *it << "\t" << summsd/nrec << endl;
  }
  cout << "\npf\tmsd\n" << endl;
  for(int i=0; i<npf; ++i)
  {
    cout << pflist[i] << "\t" << msdlist[i]/weightlist[i] << "\n";
  }
  return 0;
}

// Calculating d ln<MSD>/d ln|p-p_c|
int msdExponentGuesser(const string &inFileStr)
{
  double pf, p0, msd, pc, dtmp;
  double termA, termB, summsd;       // all
  double termA_a, termB_a, summsd_a; // one file
  double variance_est, mean_est;
  long perimeter, csize, size_limit, nrec, nrec_a;
  int nfile=0;
  char buf[100], c;
  string dimstr, line;
  vector<string> filelist;
  ifstream lstFile(inFileStr, ios::in);
  bool pcflag = false;
  long nperco = 0; // number of clusters that reached the maximum
  if(lstFile.is_open())
  {
    /*
    while(!lstFile.eof())
    {
      // Until find [file] mark
      getline(lstFile, line);
      if(!line.compare(0, 6, "[file]")) break;
    } */
    while(!lstFile.eof())
    {
      getline(lstFile, line);
      if(line.empty()) continue;
      filelist.push_back(line);
    }
    lstFile.close();
  }
  if(filelist.empty())
  {
    filelist.push_back("msddump.dat"); // default input file
  }
  summsd = termA = termB = 0.0;
  variance_est = mean_est = 0.0;
  nrec = 0;
  // cout << "File\tAverageMSD\n";
  /* Data file format
   4D
   p:  0.196
   Site limit:  606060606
   size  parameter  msd
   12  65  4.979166667
   */
  for(vector<string>::iterator it = filelist.begin(); it != filelist.end(); ++it)
  {
    ifstream inFile(*it, ios::in);
    if(!inFile.is_open()) continue;
    if(!pcflag)
    {
      inFile.getline(buf, 256);
      if(!strcmp(buf, "2D"))
        pc = 0.59274605079210;
      else if(!strcmp(buf, "3D"))
        pc = 0.3116077;
      else if(!strcmp(buf, "4D"))
        pc = 0.19688561;
      else if(!strcmp(buf, "5D"))
        pc = 0.14079633;
      else if(!strcmp(buf, "6D"))
        pc = 0.109016661;
      else if(!strcmp(buf, "7D"))
        pc = 0.088951121;
      else if(!strcmp(buf, "8D"))
        pc = 0.075210128;
      else if(!strcmp(buf, "9D"))
        pc = 0.0652095348;
      else if(!strcmp(buf, "10D"))
        pc = 0.0575929488;
      else if(!strcmp(buf, "11D"))
        pc = 0.0515896843;
      else if(!strcmp(buf, "12D"))
        pc = 0.0467309755;
      else if(!strcmp(buf, "13D"))
        pc = 0.04271507960;
      else break; // Head damaged
      // get p0
      inFile.get(buf,100,':'); inFile.get(c); inFile >> p0;
      inFile.get(buf,100,':'); inFile.get(c); inFile >> size_limit;
      inFile.getline(buf, 256);  // skip "\n"
      inFile.getline(buf, 256);
      if(p0 <=0 || p0 >=1) continue; // damaged
      pcflag = true;
    }
    else
    {
      // skip head
      inFile.getline(buf, 256);
      inFile.get(buf,100,':'); inFile.get(c); inFile >> pf;
      if(fabs(pf-p0) > 1e-10)
      {
        cerr << "File " << *it << ": p0 not the same... Skip" << endl;
        continue;
      }
      inFile.getline(buf, 256); // skip "\n"
      inFile.get(buf,100,':'); inFile.get(c); inFile >> size_limit;
      inFile.getline(buf, 256); // skip "\n"
      inFile.getline(buf, 256);
    }
    // Get p
    summsd_a = termA_a = termB_a = 0.0;
    nrec_a = 0;
    while(!inFile.eof() && nrec_a<100000000) // Truncate the count at 1e8
    {
      if(!(inFile >> csize >> perimeter >> msd))
        continue; // eof meet
      else if(csize+perimeter >= size_limit)
        nperco++;
      dtmp = csize/p0 - perimeter/(1.0-p0);
      termA_a += dtmp*msd;
      termB_a += dtmp;
      summsd_a += msd;
      ++nrec_a;
    }
    dtmp = (termA_a/summsd_a - termB_a/nrec_a)*(pc-p0);
    variance_est += dtmp*dtmp*nrec_a;
    mean_est += dtmp*nrec_a;
    termA += termA_a;
    termB += termB_a;
    summsd += summsd_a;
    nrec += nrec_a;
    nfile++;
  }
  if(nrec == 0) return 1; // no record
  dtmp = (termA/summsd - termB/nrec)*(pc-p0);
  if (nfile>1)
    variance_est = 1.96*sqrt((variance_est - mean_est*mean_est/nrec)/(nrec*(nfile-1)));
  else
    variance_est = 0;
  // result and uncertainty: 95% confidence ofo interval
  cout << p0 << "\t" << dtmp << "\t" << setprecision(2) << variance_est;
  if(nperco)
    cout << "\t" << double(nperco)/double(nrec);
  cout << endl;
  return 0;
}

int main(int argc, const char * argv[])
{
  string inFileStr;
  /* Read filelist from the file indicated by argv[1] or "msdfilelist.txt" by default.
   If msdfilelist.txt also does not exist, read directly msddump.txt */
  if(argc>1)
  {
    inFileStr = argv[1];
  }
  else
  {
    inFileStr = "msdfilelist.txt";
  }
  // msdExtrapolator(inFileStr);
  msdExponentGuesser(inFileStr);
}

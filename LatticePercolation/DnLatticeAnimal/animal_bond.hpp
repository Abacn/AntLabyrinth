//
//  animal_bound.hpp
//  DnBondAnimal
//
//  Created by Yi Hu on 1/6/21.
//

#ifndef animal_bound_hpp
#define animal_bound_hpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <string>

#include "leathsitenode.h"

template <int D>
class BondNeighborBase{
public:
  using PointType = LeathSiteNode<D>;
  using BondType = std::pair<PointType, int>;
  BondNeighborBase(): hasnext(false) {}
  bool hasNext() const {return hasnext;}
  // set bond coordinate
  virtual void setPoint(const PointType &po) = 0;
  // get point incident to the bond
  // option: 0 - bond derive from this vertex
  // 1 - bond derive from opposite vertex
  virtual PointType getPoint(const BondType &bo, int option) const = 0;
  virtual const BondType &getNext() = 0;
  virtual bool isNeighbor(const BondType &ba, const BondType &bb) const = 0;
protected:
  bool hasnext;
  PointType vertex;
  BondType nextbond;
};

template <int D>
class BondNeighborDn: public BondNeighborBase<D>{
public:
  using PointType = typename BondNeighborBase<D>::PointType;
  using BondType = typename BondNeighborBase<D>::BondType;
  // get coordinate index from index idx
  inline static void getcoordidxs(int vertexidxs, int &idx1, int &idx2, int &sign)
  {
    constexpr int bytemask = 0xff;
    sign = vertexidxs & 1;
    vertexidxs >>= 1;
    idx2 = vertexidxs & bytemask;
    idx1 = vertexidxs >> 8;
  }
  // get bond idx from coords
  inline static int getbondidxs(int idx1, int idx2, int sign)
  {
    int result = idx1 << 9;
    result |= (idx2 << 1);
    result |= sign;
    return result;
  }
  // get opposite site coordinate
  virtual PointType getPoint(const BondType &bo, int option) const
  {
    if(!option) return bo.first;
    else
    {
      int idx1, idx2, idxsign;
      getcoordidxs(bo.second, idx1, idx2, idxsign);
      PointType result = bo.first;
      ++result[idx1];
      if(!idxsign) ++result[idx2];
      else --result[idx2];
      return result;
    }
  }
  // set bond
  virtual void setPoint(const PointType &po)
  {
    this->vertex = po;
    this->hasnext = true;
    now1 = 0;
    now2 = 1;
    nowsign = noworent = 0;
  }
  
  // get next neighboring bond
  virtual const BondType &getNext()
  {
    if(!this->hasnext) // do not have next, return current point
    {
      return this->nextbond;
    }
    this->nextbond.first = this->vertex;
    if(noworent)
    {
      --(this->nextbond.first)[now1];
      if(!nowsign) --(this->nextbond.first)[now2];
      else ++(this->nextbond.first)[now2];
    }
    this->nextbond.second = getbondidxs(now1, now2, nowsign);
    
    if(!(++nowsign & 1))
    {
      nowsign = 0;
      if(++now2 == D)
      {
        if(++now1 == D-1)
        {
          if(noworent)
          {
            this->hasnext = false;
          }
          else
          {
            noworent = 1;
            now1 = now2 = nowsign = 0;
          }
        }
        now2 = now1+1;
      }
    }
    return this->nextbond;
  }
  virtual bool isNeighbor(const BondType &ba, const BondType &bb) const
  {
    if(ba.first == bb.first) return true;
    else
    {
      const PointType paa = getPoint(ba, 1);
      if(paa == bb.first) return true;
      else
      {
        const PointType pbb = getPoint(bb, 1);
        if(paa == pbb || ba.first == pbb) return true;
        else return false;
      }
    }
  }
private:
  // now1, now2: coordinate assign 1s
  // nowsign, assign 1 or -1
  // noworent, modify vertex or not
  int now1, now2, nowsign, noworent;
};

// Count animal generically
template <int D>
class CountBondAnimal{
public:
  using PointType = typename BondNeighborBase<D>::PointType;
  using BondType = typename BondNeighborBase<D>::BondType;
  using BondSetType = std::set< BondType >;
  using BondVectorType = std::vector< BondType >;
  // neighboring bond list type
  using NBondsType = std::pair<PointType, BondSetType>;
  // constructor
  CountBondAnimal(BondNeighborBase<D> &pb): nbrhdl(pb)
  {
  }
  std::vector<std::pair<int, uint64_t> > getCount(int N, const std::vector<unsigned int> &asgdidxs={})
  {
    // get count
    std::unordered_map<int, uint64_t> gmap;
    std::vector<NBondsType> nbrsets; // point-by-point neighbor set
    std::vector<BondIteratorType> itvec; // iterators for BFS
    std::vector<std::pair<int, uint64_t> > result;
    if(N <= 0) return result;
    if(asgdidxs.size() >= N)
    {
      std::cerr << "Number of pre-assigned indices (size " << asgdidxs.size() << ") not less than N (" << N << ")" << std::endl;
      return result;
    }
    // set origin
    nbrsets.reserve(N);
    addtomap(nbrsets, itvec);
    // set starting idxs.size() bonds
    for(int rp=0; rp<asgdidxs.size(); ++rp)
    {
      for(int rq=0; rq<asgdidxs[rp]; ++rq)
      {
        itvec.back().next(nbrsets);
        if(itvec.back().it == nbrsets.back().second.end() )
        {
          // touched end
          std::cout << "0\t0" << std::endl;
          return result;
        }
      }
      addtomap(nbrsets, itvec);
    }
    while(itvec.size() > asgdidxs.size())
    {
      // dump
      if(itvec.back().it == nbrsets.back().second.end() )
      {
        // reach end
        nbrsets.pop_back();
        itvec.pop_back();
        if(!itvec.empty()) itvec.back().next(nbrsets);
      }
      else
      {
        if(nbrsets.size() == N)
        {
          countIncrement(nbrsets, itvec, gmap);
          itvec.back().next(nbrsets);
        }
        else
        {
          // add
          addtomap(nbrsets, itvec);
        }
      }
    }
    for(auto const &hdl: gmap) result.push_back(hdl);
    std::sort(result.begin(), result.end());
    return result;
  }
  static const BondType zerobond;
private:
  struct BondIteratorType{
    typename BondSetType::iterator it; // iterator handle
    int now; // index of set the iterator now locates
    typename BondSetType::iterator next(std::vector<NBondsType> const &addmap)
    {
      if(it!=addmap[now].second.end())
        ++it;
      while(it==addmap[now].second.end())
      {
        ++now;
        if(now >= addmap.size())
        {
          --now;
          break;
        }
        it = addmap[now].second.lower_bound(zerobond);
      }
      return it;
    }
  };
  
  // count perimeter
  void countIncrement(std::vector<NBondsType> &addmapvec, const std::vector<BondIteratorType> &itvec, std::unordered_map<int, uint64_t> &gmap)
  {
    int t = 0;
    // first increment 1-[N-1]
    for(auto const &hdl: addmapvec)
    {
      t += hdl.second.size();
    }
    // then last vertex
    bool addbondflag = true;
    PointType addpoint;
    const BondType &lastb = *itvec.back().it;
    // check if starting point has already been added
    int chkvertex = (int)(addmapvec.size());
    while(--chkvertex >= 0)
    {
      if(addmapvec[chkvertex].first == lastb.first) break;
    }
    if(chkvertex < 0)
    {
      addpoint = lastb.first;
    }
    else
    {
      // check if end point has already been added
      addpoint = nbrhdl.getPoint(lastb, 1);
      chkvertex = (int)(addmapvec.size());
      while(--chkvertex >= 0)
      {
        if(addmapvec[chkvertex].first == addpoint) break;
      }
      if(chkvertex >= 0)
      {
        // both sides of the bonded has been added, do not add more
        addbondflag = false;
      }
    }
    // add bonds incident to the point
    if(addbondflag)
    {
      nbrhdl.setPoint(addpoint);
      while(nbrhdl.hasNext())
      {
        
        bool existed = false;
        auto const &nbrb = nbrhdl.getNext();
        // check if nbrp, the other vertex than addpoint, has already been added
        PointType nbrp;
        if(nbrb.first == addpoint)
        {
          nbrp = nbrhdl.getPoint(nbrb, 1);
        }
        else
        {
          nbrp = nbrb.first;
        }
        for(int rp=0; rp<addmapvec.size(); ++rp)
        {
          if(nbrp == addmapvec[rp].first)
          {
            existed = true; break;
          }
        }
        if(!existed)
        {
          ++t;
        }
      }
    }
    t -= addmapvec.size();
    auto it = gmap.find(t);
    if(it == gmap.end()) gmap.insert({t, 1});
    else ++gmap[t];
  }
  // add point and neighbors to map
  void addtomap(std::vector<NBondsType> &addmapvec, std::vector<BondIteratorType> &itvec)
  {
    if(addmapvec.empty())
    {
      // initialization
      nbrhdl.setPoint(zerobond.first);
      addmapvec.push_back({zerobond.first, BondSetType() });
      while(nbrhdl.hasNext())
      {
        addmapvec.back().second.insert(nbrhdl.getNext());
      }
      itvec.push_back({addmapvec.back().second.lower_bound(zerobond), 0});
    }
    else
    {
      BondSetType newset;
      bool addbondflag = true;
      PointType addpoint;
      const BondType &lastb = *itvec.back().it;
      // check if starting point has already been added
      int chkvertex = (int)(addmapvec.size());
      while(--chkvertex >= 0)
      {
        if(addmapvec[chkvertex].first == lastb.first) break;
      }
      if(chkvertex < 0)
      {
        addpoint = lastb.first;
      }
      else
      {
        // check if end point has already been added
        addpoint = nbrhdl.getPoint(lastb, 1);
        chkvertex = (int)(addmapvec.size());
        while(--chkvertex >= 0)
        {
          if(addmapvec[chkvertex].first == addpoint) break;
        }
        if(chkvertex >= 0)
        {
          // both sides of the bonded has been added, do not add more
          addbondflag = false;
        }
      }
      // add bonds incident to the point
      if(addbondflag)
      {
        nbrhdl.setPoint(addpoint);
        
        while(nbrhdl.hasNext())
        {
          
          bool existed = false;
          auto const &nbrb = nbrhdl.getNext();
          // check if nbrp, the other vertex than addpoint, has already been added
          PointType nbrp;
          if(nbrb.first == addpoint)
          {
            nbrp = nbrhdl.getPoint(nbrb, 1);
          }
          else
          {
            nbrp = nbrb.first;
          }
          for(int rp=0; rp<addmapvec.size(); ++rp)
          {
            if(nbrp == addmapvec[rp].first)
            {
              existed = true; break;
            }
          }
          if(!existed)
          {
            newset.insert(nbrb);
          }
        }
      }
      addmapvec.push_back({addpoint, std::move(newset)});
      auto it2 = itvec.back();
      it2.next(addmapvec);
      itvec.push_back( it2 );
    }
  }
  
  BondNeighborBase<D> &nbrhdl;
};


template <int D>
const typename CountBondAnimal<D>::BondType CountBondAnimal<D>::zerobond{ CountBondAnimal<D>::PointType(), 0 };

void printResult(std::vector<std::pair<int, uint64_t> > const &ilist);


#endif /* animal_bound_hpp */

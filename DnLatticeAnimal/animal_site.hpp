//
//  animal.hpp
//  DnLatticeAnimal
//
//  Created by Yi Hu on 12/10/20.
//

#ifndef animal_hpp
#define animal_hpp

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
class PointNeighborBase{
public:
  using PointType = LeathSiteNode<D>;
  PointNeighborBase(): hasnext(false) {}
  bool hasNext() const {return hasnext;}
  virtual void setPoint(const PointType &po) = 0;
  virtual const PointType &getNext() = 0;
  virtual bool isNeighbor(const PointType &pa, const PointType &pb) const = 0;
protected:
  bool hasnext;
  PointType point, nextpoint;
};

template <int D>
class PointNeighborZn: public PointNeighborBase<D>{
public:
  using PointType = typename PointNeighborBase<D>::PointType;
  virtual void setPoint(const PointType &po)
  {
    this->point = po;
    this->hasnext = true;
    now = 0;
  }
  virtual const LeathSiteNode<D> &getNext()
  {
    if(!this->hasnext) // do not have next, return current point
    {
      return this->point;
    }
    this->nextpoint = this->point;
    if(now & 1) ++(this->nextpoint)[now >> 1];
    else --(this->nextpoint)[now >> 1];
    if(++now==2*D)
    {
      this->hasnext = false;
    }
    return this->nextpoint;
  }
  virtual bool isNeighbor(const PointType &pa, const PointType &pb) const
  {
    return pa.distance_squared_T(pb) <= 1;
  }
private:
  int now;
};

template <int D>
class PointNeighborDn: public PointNeighborBase<D>{
public:
  using PointType = typename PointNeighborBase<D>::PointType;
  virtual void setPoint(const PointType &po)
  {
    this->point = po;
    this->hasnext = true;
    now1 = 0;
    now2 = 1;
    nowsign = 0;
  }
  virtual const LeathSiteNode<D> &getNext()
  {
    if(!this->hasnext) // do not have next, return current point
    {
      return this->point;
    }
    this->nextpoint = this->point;
    if(nowsign & 1) ++(this->nextpoint)[now1];
    else --(this->nextpoint)[now1];
    if(nowsign & 2) ++(this->nextpoint)[now2];
    else --(this->nextpoint)[now2];
    if(++nowsign == 4)
    {
      nowsign = 0;
      if(++now2 == D)
      {
        if(++now1 == D-1)
        {
          this->hasnext = false;
        }
        now2 = now1+1;
      }
    }
    return this->nextpoint;
  }
  virtual bool isNeighbor(const PointType &pa, const PointType &pb) const
  {
    return pa.distance_squared_T(pb) <= 2;
  }
private:
  int now1, now2, nowsign;
};

template <int D>
class LatticeAnimal{
public:
  using PointType = LeathSiteNode<D>;
  using AnimalSetType = std::set<PointType >;
  using AnimalVectorType = std::vector<PointType >;
  // set animal content
  void setAnimal(const AnimalSetType &itp)
  {
    points = itp;
  }
  void setAnimal(const AnimalVectorType &itp)
  {
    points.clear();
    for(auto const &hdl: itp) points.insert(itp);
    setMin();
  }
  int setAnimal(const char *filename)
  {
    // read animal
    std::ifstream ifs(filename);
    std::string line;
    points.clear();
    while (std::getline(ifs, line))
    {
      if (std::all_of(line.begin(), line.end(), isspace)) continue;
      std::istringstream iss(line);
      LeathSiteNode<D> newpoint;
      for (int a, rp=0; !(iss >> a).fail() && rp<D; ++rp)
      {
        newpoint[rp] = a;
      }
      points.insert(newpoint);
    }
    return 0;
  }
  const AnimalSetType &getAnimal() const
  {
    return points;
  }
  // move the animal sorted
  void setMin()
  {
    if(points.size()<1) return;
    AnimalSetType newpoints;
    for(auto it=points.begin(); it != points.end(); ++it)
    {
      newpoints.insert(*it - *points.begin());
    }
    points = std::move(newpoints);
  }
  // operators
  bool operator==( const LatticeAnimal& rhs) { return points==rhs.points; }
  bool operator!=( const LatticeAnimal& rhs) { return points!=rhs.points; }
  bool operator<( const LatticeAnimal& rhs) { return points<rhs.points; }
  bool operator<=( const LatticeAnimal& rhs) { return points<=rhs.points; }
  bool operator>( const LatticeAnimal& rhs) { return points>rhs.points; }
  bool operator>=( const LatticeAnimal& rhs) { return points>=rhs.points; }
private:
  AnimalSetType points;
};

// Count animal generically
template <int D>
class CountLatticeAnimal{
public:
  using PointType = typename LatticeAnimal<D>::PointType;
  using AnimalSetType = typename LatticeAnimal<D>::AnimalSetType;
  using AnimalVectorType = typename LatticeAnimal<D>::AnimalVectorType;
  // constructor
  CountLatticeAnimal(PointNeighborBase<D> &pb): nbrhdl(pb)
  {
  }
  std::vector<std::pair<int, uint64_t> > getCount(int N, const std::vector<unsigned int> &asgdidxs={})
  {
    // get count
    std::unordered_map<int, uint64_t> gmap;
    std::vector<AnimalSetType> nbrsets; // point-by-point neighbor set
    std::vector<AnimalIteratorType> itvec; // iterators for BFS
    std::vector<std::pair<int, uint64_t> > result;
    if(N <= 0) return result;
    if(asgdidxs.size() >= N-1)
    {
      std::cerr << "Number of pre-assigned indices (size " << asgdidxs.size() << ") not less than N-1 (" << N-1 << ")" << std::endl;
      return result;
    }
    // set origin
    nbrsets.reserve(N);
    addtomap(nbrsets, itvec);
    // set starting idxs.size() bonds
    for(int rp=0; rp<asgdidxs.size(); ++rp)
    {
      addtomap(nbrsets, itvec);
      for(int rq=0; rq<asgdidxs[rp]; ++rq)
      {
        itvec.back().next(nbrsets);
        if(itvec.back().it == nbrsets.back().end() )
        {
          // touched end
          std::cout << "0\t0" << std::endl;
          return result;
        }
      }
    }
    while(itvec.size() > asgdidxs.size())
    {
      // dump
      if(itvec.back().it == nbrsets.back().end() )
      {
        // reach end
        nbrsets.pop_back();
        itvec.pop_back();
        if(itvec.size() > 1+asgdidxs.size())
          itvec.back().next(nbrsets);
        else
          break;
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
  static const PointType zerovec;
private:
  struct AnimalIteratorType{
    typename AnimalSetType::iterator it; // iterator handle
    int now; // index of set the iterator now locates
    typename AnimalSetType::iterator next(std::vector<AnimalSetType> const &addmap)
    {
      if(it!=addmap[now].end()) ++it;
      while(it==addmap[now].end())
      {
        ++now;
        if(now >= addmap.size())
        {
          --now;
          break;
        }
        it = addmap[now].upper_bound(zerovec);
      }
      return it;
    }
    
  };
  
  // count perimeter
  void countIncrement(std::vector<AnimalSetType> &addmapvec, const std::vector<AnimalIteratorType> &itvec, std::unordered_map<int, uint64_t> &gmap)
  {
    int t = 0;
    // first increment 1-[N-1]
    for(auto const &hdl: addmapvec)
    {
      t += hdl.size();
    }
    // then last vertex
    nbrhdl.setPoint(*itvec.back().it);
    while(nbrhdl.hasNext())
    {
      // check if nbrp is already a neighbor
      bool existed = false;
      auto const &nbrp = nbrhdl.getNext();
      for(int rp=0; rp<itvec.size()-1; ++rp)
      {
        if(nbrhdl.isNeighbor(nbrp, *itvec[rp].it))
        {
          existed = true; break;
        }
      }
      if(!existed)
      {
        ++t;
      }
    }
    t -= addmapvec.size();
    auto it = gmap.find(t);
    if(it == gmap.end()) gmap.insert({t, 1});
    else ++gmap[t];
  }
  // add point and neighbors to map
  void addtomap(std::vector<AnimalSetType> &addmapvec, std::vector<AnimalIteratorType> &itvec)
  {
    if(addmapvec.empty())
    {
      // initialization
      addmapvec.push_back({zerovec});
      itvec.push_back({addmapvec.back().begin(), 0});
    }
    else
    {
      nbrhdl.setPoint(*itvec.back().it);
      AnimalSetType newset;
      while(nbrhdl.hasNext())
      {
        // check if nbrp is already a neighbor
        bool existed = false;
        auto const &nbrp = nbrhdl.getNext();
        for(int rp=0; rp<itvec.size()-1; ++rp)
        {
          if(nbrhdl.isNeighbor(nbrp, *itvec[rp].it))
          {
            existed = true; break;
          }
        }
        if(!existed)
        {
          newset.insert(nbrp);
        }
      }
      addmapvec.push_back(std::move(newset));
      auto it2 = itvec.back();
      it2.next(addmapvec);
      itvec.push_back( it2 );
    }
  }
  
  PointNeighborBase<D> &nbrhdl;
};


template <int D>
const typename CountLatticeAnimal<D>::PointType CountLatticeAnimal<D>::zerovec;

void printResult(std::vector<std::pair<int, uint64_t> > const &ilist);

template<int D>
int getPerimeter(LatticeAnimal<D> const &la, PointNeighborBase<D> &nbrhdl)
{
  auto animals = la.getAnimal();
  typename LatticeAnimal<D>::AnimalSetType peris;
  for(auto const &p: animals)
  {
    nbrhdl.setPoint(p);
    peris.insert(p);
    while(nbrhdl.hasNext())
    {
      peris.insert(nbrhdl.getNext());
    }
  }
  return peris.size()-animals.size();
}

#endif /* animal_hpp */

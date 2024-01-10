/* 
   Packing of hard spheres via molecular dynamics
   Developed by Monica Skoge, 2006, Princeton University
   Contact: Aleksandar Donev (adonev@math.princeton.edu) with questions
   This code may be used, modified and distributed freely.
   Please cite:
   
   "Packing Hyperspheres in High-Dimensional Euclidean Spaces"
   	M. Skoge, A. Donev, F. H. Stillinger and S. Torquato, 2006
	
   if you use these codes.	
*/


#ifndef GRID_FIELD_H
#define GRID_FIELD_H

#include <string.h>
#include <random>
#include "coordinate.h"
#include "read_input.h"
// ======================================================================
// grid_field 
// ======================================================================

// A field of V-Coordinates on a D dimensional manifold




template<int D, class T>
class grid_field {

 protected:
  T* f;
  long elements, occupied;
  Coordinate<D, int> size;           // number of grid points for each dimension
  Coordinate<D, long> offset;
  std::mt19937_64 rg;
 public:

  grid_field();
  grid_field(const Coordinate<D, int>&);
  ~grid_field();

  T& get(const Coordinate<D, int>&);
  void set(const Coordinate<D,int>&,const int i);
  T& get(const long p);
  void set(const long p,const int i);
  Coordinate<D, int> get_size() const;
  void set_size(const Coordinate<D, int>&);
  void set_size(const int);
  void set_seed(const int);
  void initialize();
  void initialize(const int);
  int initialize(const read_input_lattice &input);
  long getElements(){return elements;}
  long getOccupied(){return occupied;}
  long coord2int(const Coordinate<D, int>&);
  Coordinate<D, int> int2coord(long) const;
  void applypbc(Coordinate<D, int>& pos);
  void applypbc(Coordinate<D, int>& pos, int rp);
};

// coordinate convert
template<int D, class T>
long grid_field<D, T>::coord2int(const Coordinate<D, int>& pos)
{
  long p=0;
  for(int i=0; i<D; i++)
    p += pos.x[i]*offset[i];
  return p;
}

template<int D, class T>
Coordinate<D, int> grid_field<D, T>::int2coord(long p) const
{
  Coordinate<D, int> coord;
  int i;
  for(i=0; i<D; ++i)
  {
    coord.x[i] = p%size.x[i];
    p /= size.x[i];
  }
  return coord;
}

// apply pbc
template<int D, class T>
void grid_field<D, T>::applypbc(Coordinate<D, int>& pos)
{
  for(int rp=0; rp<D; ++rp)
  {
    if(pos[rp]>=size[rp])
    {
      pos[rp]-=size[rp];
    }
    else if(pos[rp]<0)
    {
      pos[rp]+=size[rp];
    }
  }
}

template<int D, class T>
void grid_field<D, T>::applypbc(Coordinate<D, int>& pos, int rp)
{
  if(pos[rp]>=size[rp])
  {
    pos[rp]-=size[rp];
  }
  if(pos[rp]<0)
  {
    pos[rp]+=size[rp];
  }
}

// grid_field
// ~~~~~~~~~~~~
template<int D, class T>
grid_field<D, T>::grid_field()
: f(0), elements(0), occupied(0)
{
}

// grid_field
// ~~~~~~~~~~~~
template<int D, class T>
grid_field<D, T>::grid_field(const Coordinate<D, int>& s)
: f(0)
{
  set_size(s);
}


// ~grid_field
// ~~~~~~~~~~~~~
template <int D, class T>
grid_field<D, T>::~grid_field()
{
  if(f != 0)
    delete[] f;
}


// get_size
// ~~~~~~~~
template<int D, class T>
inline Coordinate<D, int> grid_field<D, T>::get_size() const
{
  return size;
}


// set_size
// ~~~~~~~~
template<int D, class T>
void grid_field<D, T>::set_size(const Coordinate<D, int>& s)
{
  if(f != 0)
    delete[] f;
  
  size = s;
  
  elements = 1;
  for(int i=0; i<D; i++) {
    offset[i] = elements;
    elements *= size[i];
  }
  
  f = new T[elements];
}


// set_size
// ~~~~~~~~
template<int D, class T>
void grid_field<D, T>::set_size(const int s)
{
  Coordinate<D, int> square;
  
  for(int k=0; k<D; k++)
    square[k] = s;
  
  set_size(square);
}


// set_seed
// ~~~~~~~~
template<int D, class T>
void grid_field<D, T>::set_seed(const int seed)
{
  ;
  
  rg.seed(seed);
}


// get
// ~~~
template<int D, class T>
inline T& grid_field<D, T>::get(const Coordinate<D, int>& pos)
{
  long p=0;
  for(int i=0; i<D; i++)
    p += pos.x[i]*offset[i];
  
  return f[p];
}

template<int D, class T>
inline T& grid_field<D, T>::get(const long p)
{
  return f[p];
}

//set
//~~
template<int D, class T>
void grid_field<D, T>::set(const Coordinate<D, int>& pos,const int value)
{
  long p=0;
  for(int i=0; i<D; i++)
    p += pos[i]*offset[i];
  
  f[p]=value;
}

template<int D, class T>
void grid_field<D, T>::set(const long p,const int value)
{
  f[p]=value;
}

// initialize
// ~~~
template<int D, class T>
void grid_field<D, T>::initialize()
{
  memset(f, 0, sizeof(T)*elements);
}

template<int D, class T>
void grid_field<D, T>::initialize(const int value)
{
  for(long i=0; i<elements; ++i)
  {
    f[i] = value;
  }
}

template<int D, class T>
int grid_field<D, T>::initialize(const read_input_lattice &input)
{
  double pf = input.pf;
  int vacant;
  long xrand=0;
  std::uniform_int_distribution<long> randrange(0L,elements);
  if(pf <= 0.5)
  {
    this->initialize();
    vacant = 0;
  }
  else
  {
    pf = 1 - pf;
    this->initialize(1);
    vacant = 1;
  }
  
  //Randomly place the other blocked sites on the lattice
  occupied = (long)(pf*elements);
  for(long i=0;i<occupied;)
  {
    xrand = randrange(rg);
    if(f[xrand]==vacant)
    {
      f[xrand] = !vacant;
      ++i;
    }
  }
  
  // origin is 0
  if(f[0]==0)
  {
    if(vacant!=0)
    {
     // In the case of pf>0.5, randoming find a filled space to switch
      do{
        xrand = randrange(rg);
      }while(f[xrand]==0);
    }
    f[xrand] = 0;
    f[0] = 1;
  }
  return 0;
}

#endif

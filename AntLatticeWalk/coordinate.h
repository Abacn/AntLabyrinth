#ifndef COORDINATE_H
#define COORDINATE_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include "dim.h"

template <int D, typename T=double>
class Coordinate {
  
  public:
  T x[D];
  
  public:
  
  Coordinate();
  Coordinate(const T[D]);
  Coordinate(const Coordinate&);
  ~Coordinate();

  Coordinate<D, T>& operator+=(const Coordinate<D, T>&);
  Coordinate<D, T>& operator-=(const Coordinate<D, T>&);
  Coordinate<D, T>& operator*=(const T);
  Coordinate<D, T>& operator/=(const T);
  Coordinate<D, T> operator+(const Coordinate<D, T>&) const;
  Coordinate<D, T> operator-(const Coordinate<D, T>&) const;
  Coordinate<D, T> operator*(const T) const;
  Coordinate<D, T> operator/(const T) const;
  Coordinate<D, T> operator%(const T) const;
  bool operator==(const Coordinate<D, T> &a) const;
  bool operator<(const Coordinate<D, T> &a) const;
  Coordinate<D, int> integer() const;
  Coordinate<D, long> Long() const;
  Coordinate<D, double> Double() const;
  static Coordinate<D, int> integer(const Coordinate<D, T>&);
  static Coordinate<D, long> Long(const Coordinate<D, T>&);
  static Coordinate<D, double> Double(const Coordinate<D, T>&);
  
  T& operator[](const unsigned int);
 
  double dot(const Coordinate<D, T>&) const;
  static double dot(const Coordinate<D, T>&, const Coordinate<D, T>&);
  Coordinate<D, T> multiply(const Coordinate<D, T>&) const;
  Coordinate<D, T> squared() const;
  T sum() const;
  double norm_squared() const;
  static double norm_squared(const Coordinate<D, T>&);

  void read(std::ifstream&);
  void write(std::ofstream&) const;
};


template <int D, typename T>
std::ostream& operator<<(std::ostream&, const Coordinate<D, T>&);


// constructor
// ~~~~~~~~~~~
template <int D, typename T>
Coordinate<D, T>::Coordinate()
{
  for(int k=0; k<D; k++)
    x[k] = 0;
}

template <int D, typename T>
Coordinate<D, T>::Coordinate(const T x_i[D])
{
  for(int k=0; k<D; k++)
    x[k] = x_i[k];
}

template <int D, typename T>
Coordinate<D, T>::Coordinate(const Coordinate<D, T> &v)
{
  for(int k=0; k<D; k++)
    x[k] = v.x[k];
}


// destructor
// ~~~~~~~~~~
template <int D, typename T>
Coordinate<D, T>::~Coordinate()
{
}


// +=
// ~~
template <int D, typename T>
inline Coordinate<D, T>& Coordinate<D, T>::operator+=(const Coordinate<D, T> &v)
{
  for(int k=0; k<D; k++)
    x[k] += v.x[k];

  return *this;
}


// -=
// ~~
template <int D, typename T>
inline Coordinate<D, T>& Coordinate<D, T>::operator-=(const Coordinate<D, T> &v)
{
  for(int k=0; k<D; k++)
    x[k] -= v.x[k];

  return *this;
}


// *=
// ~~
template <int D, typename T>
inline Coordinate<D, T>& Coordinate<D, T>::operator*=(const T s)
{
  for(int k=0; k<D; k++)
    x[k] *= s;

  return *this;
}

// /=
// ~~
template <int D, typename T>
inline Coordinate<D, T>& Coordinate<D, T>::operator/=(const T s)
{
  for(int k=0; k<D; k++)
    x[k] /= s;

  return *this;
}


// +
// ~ 
template <int D, typename T>
inline Coordinate<D, T> Coordinate<D, T>::operator+(const Coordinate<D, T> &a) const
{
  Coordinate<D, T> c;
  
  for(int k=0; k<D; k++)
    c.x[k] = x[k] + a.x[k];

  return c;
}


// -
// ~
template <int D, typename T>
inline Coordinate<D, T> Coordinate<D, T>::operator-(const Coordinate<D, T> &a) const
{
  Coordinate<D, T> c;

  for(int k=0; k<D; k++)
    c[k] = x[k] - a.x[k];

  return c;
}


// *
// ~
template <int D, typename T>
inline Coordinate<D, T> Coordinate<D, T>::operator*(const T s) const
{
  Coordinate<D, T> c;
  
  for(int k=0; k<D; k++)
    c[k] = x[k] * s;

  return c;
}


// /
// ~
template <int D, typename T>
inline Coordinate<D, T> Coordinate<D, T>::operator/(const T s) const
{
  Coordinate<D, T> c;

  for(int k=0; k<D; k++)
    c[k] = x[k] / s;

  return c;
}

// ==
// ~
template <int D, typename T>
inline bool Coordinate<D, T>::operator==(const Coordinate<D, T> &a) const
{
  for(int k=0; k<D; k++)
    {
      if (x[k]!=a.x[k])
	  return false;
    }
  return true;
  // return std::equal(x, x+D, a.x); // slower
}

// <
// ~
template <int D, typename T>
inline bool Coordinate<D, T>::operator<(const Coordinate<D, T> &a) const
{
  /* // an aggressive optimization
#if DIM==4
  long *pa = (long*)x, *pb = (long*)a.x;
  if ((*pa) > (*pb)) return false;
  else if ((*pa) < (*pb)) return true;
  else if ((*(pa+1)) < (*(pb+1))) return true;
  else return false;
#else */
  for(int k=0; k<D; k++)
  {
    if (x[k]==a.x[k])
      continue;
    else if (x[k]<a.x[k])
      return true;
    else
      return false;
  }
  return false;
}

// %
// ~
template <int D, typename T>
inline Coordinate<D, T> Coordinate<D, T>::operator%(const T s) const
{
  Coordinate<D, T> c;

  for(int k=0; k<D; k++)
    c[k] = x[k] % s;

  return c;
}


// integer
// ~~~~~~~
template <int D, typename T>
inline Coordinate<D, int> Coordinate<D, T>::integer() const
{
  Coordinate<D, int> c;

  for(int k=0; k<D; k++)
    c[k] = (int)x[k];

  return c;
}

template <int D, typename T>
inline Coordinate<D, int> Coordinate<D, T>::integer(const Coordinate<D, T>& v)
{
  return v.integer();
}


// long
// ~~~~~~~
template <int D, typename T>
inline Coordinate<D, long> Coordinate<D, T>::Long() const
{
  Coordinate<D, long> c;
  
  for(int k=0; k<D; k++)
    c[k] = (long)x[k];
  
  return c;
}

template <int D, typename T>
inline Coordinate<D, long> Coordinate<D, T>::Long(const Coordinate<D, T>& v)
{
  return v.Long();
}

// double
// ~~~~~~~
template <int D, typename T>
inline Coordinate<D, double> Coordinate<D, T>::Double() const
{
  Coordinate<D, double> c;

  for(int k=0; k<D; k++)
    c[k] = (double)x[k];

  return c;
}

template <int D, typename T>
inline Coordinate<D, double> Coordinate<D, T>::Double(const Coordinate<D, T>& v)
{
  return v.Double();
}



// []
// ~~
template <int D, typename T>
inline T& Coordinate<D, T>::operator[](const unsigned int i)
{
  return x[i];
}


// Dot
// ~~~
template <int D, typename T>
inline double Coordinate<D, T>::dot(const Coordinate<D, T> &a) const
{
  double d=0;

  for(int k=0; k<D; k++)
    d += x[k] * a.x[k];

  return d;
}

template <int D, typename T>
inline double Coordinate<D, T>::dot(const Coordinate<D, T> &a, const Coordinate<D, T> &b)
{
  return a.dot(b);
}


// NormSquared
// ~~~~~~~~~~~
template <int D, typename T>
inline double Coordinate<D, T>::norm_squared() const
{
  double d=0, dk;
  // Attention. long*long may overflow!
  for(int k=0; k<D; k++)
  {
    dk = (double)x[k];
    d += dk*dk;
  }
  return d;
}

template <int D, typename T>
inline double Coordinate<D, T>::norm_squared(const Coordinate<D, T>& v)
{
  return v.norm_squared();
}


// .*
// ~
template <int D, typename T>
Coordinate<D, T> Coordinate<D, T>::multiply(const Coordinate<D, T>& a) const
{
  Coordinate<D, T> c;
  
  for(int k=0; k<D; k++)
    c[k] = x[k] * a[k];
  
  return c;
}


// ^2
// ~
template <int D, typename T>
Coordinate<D, T> Coordinate<D, T>::squared() const
{
  Coordinate<D, T> c;
  
  for(int k=0; k<D; k++)
    c[k] = x[k] * x[k];
  
  return c;
}

// sum
// ~
template <int D, typename T>
T Coordinate<D, T>::sum() const
{
  T c = 0;
  
  for(int k=0; k<D; k++)
    c += x[k];
  
  return c;
}

// read
// ~~~~
template <int D, typename T>
void Coordinate<D, T>::read(std::ifstream& in)
{
  in.read((char*)x, sizeof(T)*D);
}

// write
// ~~~~~
template <int D, typename T>
void Coordinate<D, T>::write(std::ofstream& out) const
{
  out.write((const char*)x, sizeof(T)*D);
}



// Insertion
// ~~~~~~~~~
template <int D, typename T>
std::ostream& operator<<(std::ostream& os, const Coordinate<D, T>& v)
{
  os << "(";

  for(int k=0; k<D-1; k++)
    os << v.x[k] << ", ";

  os << v.x[D-1] << ")";

  return os;
}


// ======================================================================
// Coordinate_field
// ======================================================================

// A field of V-Coordinates on a D dimensional manifold

template<int V, int D, typename T=double>
class Coordinate_field {
  
 public:
  int elements;

 private:
  Coordinate<V, T>* f;
  Coordinate<D, int> size;           // number of grid points for each dimension
  Coordinate<D, int> offset;
 
 public:

  Coordinate_field();
  Coordinate_field(const Coordinate<D, int>&);
  ~Coordinate_field();

  Coordinate<D, int> get_size() const;
  void set_size(const Coordinate<D, int>&);

  Coordinate<V, T>& get(const Coordinate<D, int>&);

  void read(std::ifstream&);
  void write(std::ofstream&) const;

  static void swap(Coordinate_field<V, D, T>&, Coordinate_field<V, D, T>&);
};


// Coordinate_field
// ~~~~~~~~~~~~
template<int V, int D, typename T>
Coordinate_field<V, D, T>::Coordinate_field()
  : f(0), elements(0)
{
}


// Coordinate_field
// ~~~~~~~~~~~~
template<int V, int D, typename T>
Coordinate_field<V, D, T>::Coordinate_field(const Coordinate<D, int>& s)
  : f(0)
{
  set_size(s);
}

// ~Coordinate_field
// ~~~~~~~~~~~~~
template <int V, int D, typename T>
Coordinate_field<V, D, T>::~Coordinate_field()
{
  if(f != 0)
    delete[] f;
}

// get_size
// ~~~~~~~~
template<int V, int D, typename T>
inline Coordinate<D, int> Coordinate_field<V, D, T>::get_size() const
{
  return size;
}

// set_size
// ~~~~~~~~
template<int V, int D, typename T>
void Coordinate_field<V, D, T>::set_size(const Coordinate<D, int>& s)
{
  if(f != 0)
    delete[] f;

  size = s;

  elements = 1;
  for(int i=0; i<D; i++) {
    offset[i] = elements;
    elements *= size.x[i];
  }

  f = new Coordinate<V, T>[elements];
}

// get
// ~~~
template<int V, int D, typename T>
inline Coordinate<V, T>& Coordinate_field<V, D, T>::get(const Coordinate<D, int>& pos)
{
  int p=0;
  for(int i=0; i<D; i++)
    p += pos.x[i]*offset[i];

  return f[p];
}


// read
// ~~~~
template<int V, int D, typename T>
void Coordinate_field<V, D, T>::read(std::ifstream& in)
{
  in.read((char*)f, elements*sizeof(T)*V);
}


// write
// ~~~~~
template<int V, int D, typename T>
void Coordinate_field<V, D, T>::write(std::ofstream& out) const
{
  out.write((const char*)f, elements*sizeof(T)*V);
}

// swap
// ~~~~
template<int V, int D, typename T>
void Coordinate_field<V, D, T>::swap(Coordinate_field<V, D, T>& v1,
				 Coordinate_field<V, D, T>& v2)
{
  Coordinate<V, T>* f;

  f = v1.f;
  v1.f = v2.f;
  v2.f = f;
}



#endif 

#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <fstream>

template <int D, typename T=double>
struct vector {
  T x[D];
  
  vector();
  vector(double entry);
  vector(const T[D]);           
  vector(const vector&);
  ~vector();

  vector<D, T>& operator+=(const vector<D, T>&);
  vector<D, T>& operator-=(const vector<D, T>&);
  vector<D, T>& operator*=(const T);
  vector<D, T>& operator/=(const T);
  vector<D, T> operator+(const vector<D, T>&) const;
  vector<D, T> operator-(const vector<D, T>&) const;
  vector<D, T> operator*(const T) const;
  vector<D, T> operator/(const T) const;
  vector<D, T> operator%(const T) const;
  bool operator==(const vector<D, T> &a) const;

  void assum(const vector<D, T>&, const vector<D, T>&);
  void assub(const vector<D, T>&, const vector<D, T>&);
  
  vector<D, int> integer() const;
  vector<D, double> Double() const;
  static vector<D, int> integer(const vector<D, T>&); 
  static vector<D, double> Double(const vector<D, T>&); 
  
  T& operator[](const unsigned int);
  const T& operator[](const unsigned int ) const;
  double dot(const vector<D, T>&) const;
  static double dot(const vector<D, T>&, const vector<D, T>&);
  
  double norm_squared() const;
  void quadsum(double &, double &) const;
  static double norm_squared(const vector<D, T>&);

  void read(std::ifstream&);
  void write(std::ofstream&) const;
  void setzero();
};


template <int D, typename T>
std::ostream& operator<<(std::ostream&, const vector<D, T>&);


// constructor
// ~~~~~~~~~~~
template <int D, typename T>
vector<D, T>::vector() = default;

template <int D, typename T>
vector<D, T>::vector(double entry)
{
  for(int k=0; k<D; ++k)
    x[k] = entry;
}

template <int D, typename T>
vector<D, T>::vector(const T x_i[D])
{
  for(int k=0; k<D; ++k)
    x[k] = x_i[k];
}

template <int D, typename T>
vector<D, T>::vector(const vector<D, T> &v)
{
  for(int k=0; k<D; ++k)
    x[k] = v.x[k];
}


// destructor
// ~~~~~~~~~~
template <int D, typename T>
vector<D, T>::~vector()
{
}


// +=
// ~~
template <int D, typename T>
inline vector<D, T>& vector<D, T>::operator+=(const vector<D, T> &v)
{
  for(int k=0; k<D; k++)
    x[k] += v.x[k];

  return *this;
}


// -=
// ~~
template <int D, typename T>
inline vector<D, T>& vector<D, T>::operator-=(const vector<D, T> &v)
{
  for(int k=0; k<D; ++k)
    x[k] -= v.x[k];

  return *this;
}


// *=
// ~~
template <int D, typename T>
inline vector<D, T>& vector<D, T>::operator*=(const T s)
{
  for(int k=0; k<D; ++k)
    x[k] *= s;

  return *this;
}

// /=
// ~~
template <int D, typename T>
inline vector<D, T>& vector<D, T>::operator/=(const T s)
{
  for(int k=0; k<D; ++k)
    x[k] /= s;

  return *this;
}


// +
// ~ 
template <int D, typename T>
inline vector<D, T> vector<D, T>::operator+(const vector<D, T> &a) const
{
  vector<D, T> c;
  
  for(int k=0; k<D; ++k)
    c.x[k] = x[k] + a.x[k];

  return c;
}


// -
// ~
template <int D, typename T>
inline vector<D, T> vector<D, T>::operator-(const vector<D, T> &a) const
{
  vector<D, T> c;

  for(int k=0; k<D; ++k)
    c[k] = x[k] - a.x[k];

  return c;
}


// *
// ~
template <int D, typename T>
inline vector<D, T> vector<D, T>::operator*(const T s) const
{
  vector<D, T> c;
  
  for(int k=0; k<D; ++k)
    c[k] = x[k] * s;

  return c;
}


// /
// ~
template <int D, typename T>
inline vector<D, T> vector<D, T>::operator/(const T s) const
{
  vector<D, T> c;

  for(int k=0; k<D; ++k)
    c[k] = x[k] / s;

  return c;
}


// ==
// ~
template <int D, typename T>
inline bool vector<D, T>::operator==(const vector<D, T> &a) const
{
  for(int k=0; k<D; ++k)
    {
      if (!(x[k]==a.x[k]))
	return false;
    }
  return true;
}


// %
// ~
template <int D, typename T>
inline vector<D, T> vector<D, T>::operator%(const T s) const
{
  vector<D, T> c;

  for(int k=0; k<D; ++k)
    c[k] = x[k] % s;

  return c;
}

template <int D, typename T>
inline void vector<D, T>::assum(const vector<D, T>& a, const vector<D, T>& b)
{
  for(int k=0; k<D; ++k)
    x[k] = a[k] + b[k];
}

template <int D, typename T>
inline void vector<D, T>::assub(const vector<D, T>& a, const vector<D, T>& b)
{
  for(int k=0; k<D; ++k)
    x[k] = a[k] - b[k];
}

// integer
// ~~~~~~~
template <int D, typename T>
inline vector<D, int> vector<D, T>::integer() const
{
  vector<D, int> c;

  for(int k=0; k<D; ++k)
    c[k] = (int)x[k];

  return c;
}

template <int D, typename T>
inline vector<D, int> vector<D, T>::integer(const vector<D, T>& v)
{
  return v.integer();
}


// double
// ~~~~~~~
template <int D, typename T>
inline vector<D, double> vector<D, T>::Double() const
{
  vector<D, double> c;

  for(int k=0; k<D; ++k)
    c[k] = (double)x[k];

  return c;
}

template <int D, typename T>
inline vector<D, double> vector<D, T>::Double(const vector<D, T>& v)
{
  return v.Double();
}



// []
// ~~
template <int D, typename T>
inline T& vector<D, T>::operator[](const unsigned int i)
{
  return x[i];
}

// [] const
// ~~
template <int D, typename T>
const T& vector<D, T>::operator[](const unsigned int i) const
{
  return x[i];
}

// Dot
// ~~~
template <int D, typename T>
inline double vector<D, T>::dot(const vector<D, T> &a) const
{
  double d=0.0;

  for(int k=0; k<D; ++k)
    d += x[k] * a.x[k];

  return d;
}

template <int D, typename T>
inline double vector<D, T>::dot(const vector<D, T> &a, const vector<D, T> &b)
{
  return a.dot(b);
}


// NormSquared
// ~~~~~~~~~~~
template <int D, typename T>
inline double vector<D, T>::norm_squared() const
{
  double d=0.0;
  
  for(int k=0; k<D; ++k)
    d += x[k] * x[k];
  return d;
}

// return both sum x^2 and sum x^4
template <int D, typename T>
inline void vector<D, T>::quadsum(double &squared, double &quad) const
{
  squared = quad = 0.;
  double dtmp;
  for(int k=0; k<D; ++k)
  {
    dtmp = x[k] * x[k];
    squared += dtmp;
    quad += dtmp*dtmp;
  }
}

template <int D, typename T>
inline double vector<D, T>::norm_squared(const vector<D, T>& v)
{
  return v.norm_squared();
}


// read
// ~~~~
template <int D, typename T>
void vector<D, T>::read(std::ifstream& in)
{
  in.read((char*)x, sizeof(T)*D);
}

// write
// ~~~~~
template <int D, typename T>
void vector<D, T>::write(std::ofstream& out) const
{
  out.write((const char*)x, sizeof(T)*D);
}



// Insertion
// ~~~~~~~~~
template <int D, typename T>
std::ostream& operator<<(std::ostream& os, const vector<D, T>& v)
{
  os << "(";

  for(int k=0; k<D-1; ++k)
    os << v.x[k] << ", ";

  os << v.x[D-1] << ")";

  return os;
}

// Set zero
template <int D, typename T>
inline void vector<D, T>::setzero()
{
  for(int k=0; k<D; ++k)
    x[k] = 0;
}
#endif /*VECTOR_H*/

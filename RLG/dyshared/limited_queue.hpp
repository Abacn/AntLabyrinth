//
//  limited_queue.hpp
//  AntLatticeWalk
//
//  Queue with Limited size
#ifndef limited_queue_hpp
#define limited_queue_hpp

#include <stdio.h>
template <typename T>
class limited_queue {
private:
  unsigned long maxsize, now_size, now_start;
  T* tarray;
public:
  limited_queue(const unsigned long max_size);
  ~limited_queue();
  T& operator[](const unsigned long idx);
  int push(T& item);
  inline unsigned long size(){return now_size;}
  void pop();
};

template<typename T>
inline limited_queue<T>::limited_queue(const unsigned long max_size):
maxsize(max_size)
{
  tarray = new T[max_size];
  now_size = 0;
  now_start = 0;
}

template<typename T>
inline limited_queue<T>::~limited_queue()
{
  delete [] tarray;
}

// index operator
template<typename T>
inline T& limited_queue<T>::operator[](const unsigned long idx)
{
  unsigned long true_idx = now_start + idx;
  if(true_idx >= now_size) true_idx -= now_size;
  return tarray[true_idx];
}

// push one element at the back
template<typename T>
inline int limited_queue<T>::push(T& item)
{
  unsigned long end_idx = now_start + now_size;
  int push_flag;
  if(end_idx >= maxsize) end_idx -= maxsize;
  tarray[end_idx] = item;
  if(now_size == maxsize)
  {
    if(++now_start == maxsize) now_start = 0;
    push_flag = 1;
  }
  else
  {
    ++now_size;
    push_flag = 0;
  }
  return push_flag;
}

// pop the front element
template<typename T>
inline void limited_queue<T>::pop()
{
  if(++now_start == maxsize) now_start = 0;
  --now_size;
}

#endif /* limited_queue_hpp */

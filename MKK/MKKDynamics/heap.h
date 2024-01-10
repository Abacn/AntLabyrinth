//---------------------------------------------------------------------------
// Event heap maker
//---------------------------------------------------------------------------

#ifndef  HEAP_H
#define  HEAP_H

#include "event.h"
#include "sphere.h"

class EventHeap {

 public:

  // constructor and destructor
  EventHeap(int maxsize);
  EventHeap(const EventHeap &h);
  ~EventHeap();

  // variables
  int maxsize;   // max allowed number of events
  int N;         // current number of events
  int *a;
  Sphere *s;
  int *index;     // array of indices for each sphere
  //event minevent;


  // functions which operate on a binary heap

  void upheap(int k);
  void downheap(int k);
  void clear();
  void insert(int i);
  void replace(int i);
  int search(int j);
  void change(int i);
  int extractmax();
  void print();
  void checkindex();

};
#endif

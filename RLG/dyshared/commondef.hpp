//
//  commondef.hpp
//  RLGDynamics
//
//  Created by Yi Hu on 6/13/19.
//
// Common definitions

#ifndef commondef_hpp
#define commondef_hpp

#include <cmath>
#include <vector>
#include <chrono>

#include "dim.h"
#include "vector.h"

using point = vector< DIM >;         // point coordinates typw
using points = std::vector< point >; // collection of points
using Vector = vector< DIM >;        // vector type
using Point_handle = uint32_t;       // point handle aka index
using clock_type = std::chrono::time_point<std::chrono::steady_clock>;
using duration_type = std::chrono::duration<double>;
#ifdef DEBUG
// time of building neighbor list; searching collide; do the integration, do statistics
extern duration_type nlisttime, collidetime, intgtime, stattime;
#endif

#endif /* commondef_hpp */

// -*- C++ -*-
//----------------------------------------------------------------------
//
//                    William S. Klug & Feng Feng
//                University of California Los Angeles
//                 (C) 2004-2008 All Rights Reserved
//
//----------------------------------------------------------------------
//
#ifndef _VOOM_H_
#define _VOOM_H_

/* #ifdef __APPLE__ */
/* #undef _X */
/* #undef _F */
/* #undef _P */
/* #endif */

/* using namespace std; */
#include "mpi.h"
#include <iostream>
#include <iomanip>

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

#include <tvmet/Vector.h>
#include <tvmet/Matrix.h>

namespace voom {
  typedef double Real;
  typedef blitz::Array<double,1>          Array1D;
  typedef blitz::Array<double,2>          Array2D;
  typedef blitz::Array<double,3>          Array3D;
  typedef blitz::Array<double,4>          Array4D;
  
  typedef tvmet::Vector<double,3>     Vector3D;
  typedef tvmet::Matrix<double,3,3>   Tensor3D;

  typedef tvmet::Vector<double,2>     Vector2D;
  typedef tvmet::Matrix<double,2,2>   Tensor2D;
  
  using tvmet::dot;
  using tvmet::cross;
  using tvmet::norm2;
  using tvmet::trans;
  using tvmet::trace;
  
  enum ComputeRequest {
    nothing 	= 0,
    energy 	= 1,
    force 	= 2,
    stiffness 	= 4,
    geometry	= 8
  };

  enum GlobalConstraint {
    noConstraint = 0,
    multiplier = 1,
    penalty = 2,
    augmented = 3
  };

  enum OutputSelection {
    noOutput = 0,
    paraview = 1 
  };


}; // namespace voom

#endif // _VOOM_H_

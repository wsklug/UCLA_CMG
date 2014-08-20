// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2008 All Rights Reserved
//
//----------------------------------------------------------------------

#include "Hermite.h"

namespace voom {

  void Hermite::compute(double xi) {
    
    // standard form for Hermitian polynomials is with coordinate
    // ranging 0<=t<=1.  Need to transform derivatives by multiplying
    // by jacobian dt/dxi.

    double t = 0.5*(xi+1.0);
    double t2 = t*t;
    double t3 = t2*t;
    _functions[0] = 1.0 - 3.0*t2 + 2.0*t3;
    _functions[1] = t - 2.0*t2 + t3;
    _functions[2] = 3.0*t2 - 2.0*t3;
    _functions[3] = -t2 + t3;

    _derivatives[0] = 6.0*(t2-t);
    _derivatives[1] = 1.0 - 4.0*t + 3.0*t2;
    _derivatives[2] = 6.0*(t-t2);
    _derivatives[3] = -2.0*t + 3.0*t2;

    for(int i=0; i<4; i++) _derivatives[i] *= 0.5;

    _2derivatives[0] = 12.0*t - 6.0;
    _2derivatives[1] = -4.0 + 6.0*t;
    _2derivatives[2] = 6.0 - 12.0*t;
    _2derivatives[3] = -2.0 + 6.0*t;

    for(int i=0; i<4; i++) _2derivatives[i] *= 0.25;

    return;
  }

}

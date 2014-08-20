// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2005 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.2  2005/04/21 01:47:42  klug
// Added cvs loggin.
//
//----------------------------------------------------------------------

#include "ShapeTri3.h"

namespace voom {  
  void ShapeTri3::compute(const CoordinateArray & s) {
    _functions[0] = s(0);
    _functions[1] = s(1);
    _functions[2] = 1.0-s(0)-s(1);

    _derivatives[0] =  1.0,  0.0;
    _derivatives[1] =  0.0,  1.0;
    _derivatives[2] = -1.0, -1.0;
  }
};

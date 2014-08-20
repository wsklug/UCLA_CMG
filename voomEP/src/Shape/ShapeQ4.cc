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

#include "ShapeQ4.h"

namespace voom {  
  void ShapeQ4::compute(const CoordinateArray & s) {
    _functions[0] = 0.25*(1-s(0))*(1-s(1));
    _functions[1] = 0.25*(1+s(0))*(1-s(1));
    _functions[2] = 0.25*(1+s(0))*(1+s(1));
    _functions[3] = 0.25*(1-s(0))*(1+s(1));

    _derivatives[0] =  -0.25*(1-s(1)),  -0.25*(1-s(0));
    _derivatives[1] =   0.25*(1-s(1)),  -0.25*(1+s(0));
    _derivatives[2] =   0.25*(1+s(1)),   0.25*(1+s(0));
    _derivatives[3] =  -0.25*(1+s(1)),   0.25*(1-s(0));   

  }
};

// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2006 All Rights Reserved
//
//----------------------------------------------------------------------

#include "ShapeTri6.h"

namespace voom {  
  void ShapeTri6::compute(const CoordinateArray & s) {
    /* Nodal arrangement
    //! Changed the numbering of nodes to comply with triangle (http://www.cs.cmu.edu/~quake/triangle.highorder.html)


     * 
     * ^ s2 
     * |
     * 1
     * | \ 
     * 4  3
     * |    \
     * 2--5--0--> s1
     */
    _functions[3] = 4.0*s(0)*s(1);
    _functions[4] = 4.0*s(1)*(1.0-s(0)-s(1));
    _functions[5] = 4.0*s(0)*(1.0-s(0)-s(1));
    _functions[0] = s(0)          - 0.5*_functions[3] - 0.5*_functions[5];
    _functions[1] = s(1)          - 0.5*_functions[3] - 0.5*_functions[4];
    _functions[2] = 1.0-s(0)-s(1) - 0.5*_functions[5] - 0.5*_functions[4];
//     _functions(0) = s(0)*(2.0*s(0)-1.0);
//     _functions(1) = s(1)*(2.0*s(1)-1.0);
//     _functions(2) = (1.0-s(0)-s(1))*(2.0*(1.0-s(0)-s(1))-1.0);

    _derivatives[0] =	4.0*s(0)-1.0	   	,  0.0;
    _derivatives[1] =	0.0         	   	,  4.0*s(1)-1.0;
    _derivatives[2] =	4.0*(s(0)+s(1))-3.0	,  4.0*(s(0)+s(1))-3.0;
    _derivatives[3] =	4.0*s(1)		,  4.0*s(0);
    _derivatives[4] =  -4.0*s(1)		,  4.0*(1.0-s(0)-2.0*s(1));
    _derivatives[5] =  	4.0*(1.0-2.0*s(0)-s(1))	, -4.0*s(0);

  }
};

// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file ShapeTri3.h

  \brief Class for Linear 3-Node Triangular Finite Element shape functions.

*/

#if !defined(__ShapeTri3_h__)
#define __ShapeTri3_h__

#include "Shape.h"

namespace voom
{
  //! Linear 3-node triangular isoparametric shape functions.
  /*! These 3-node triangular shape functions are to be used for 2-D
    planar problems.  Two of the three barycentric coordinates are
    used as curveliear coordinates.  The nodal numbering arrangement
    and the geometry of the standard domain are determined by the
    method <tt>nodalCoordinates()</tt>.
  */
  class ShapeTri3 : public Shape<2> {
    
  public:

    //! constructor
    ShapeTri3(const CoordinateArray & s) {
      const unsigned int n = 3;
      _functions.resize(n);
      _derivatives.resize(n);
      _positions.resize(n);

      _positions[0] = 1.0, 0.0;
      _positions[1] = 0.0, 1.0;
      _positions[2] = 0.0, 0.0;

      compute(s);
    }

    //! compute the shape functions and derivatives
    void compute(const CoordinateArray & s);

    //! Return parametric coordinates of nodes.
    /*! Nodal arrangement   	\verbatim

    ^ s2        	
    |		
    2		
    | \		
    3--1--> s1			\endverbatim
    
    */
    PositionContainer nodalCoordinates() {
      return _positions;
    }
  };
}


#endif //#define __ShapeTri3_h__

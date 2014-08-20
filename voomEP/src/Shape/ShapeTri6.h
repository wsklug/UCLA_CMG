// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file ShapeTri6.h

  \brief Class for Quadratic 6-Node Triangular Finite Element shape functions.

*/

#if !defined(__ShapeTri6_h__)
#define __ShapeTri6_h__

#include "Shape.h"

namespace voom
{
  //! Quadratic 6-node triangular isoparametric shape functions.
  /*! These 6-node triangular shape functions are to be used for 2-D
    planar problems.  Two of the three barycentric coordinates are
    used as curveliear coordinates.  The nodal numbering arrangement
    and the geometry of the standard domain are determined by the
    method <tt>nodalCoordinates()</tt>.
  */
  class ShapeTri6 : public Shape<2> {
    
  public:

    //! constructor
    ShapeTri6(const CoordinateArray & s) {
      const unsigned int n = 6;
      _functions.resize(n);
      _derivatives.resize(n);
      _positions.resize(n);
      
      //! Return parametric coordinates of nodes.
      /*! Nodal arrangement    \verbatim
	
      ^ s2 		
      |		
      1		
      | \ 		
      4  3		
      |    \		
      2--5--0--> s1	    \endverbatim
      */
      _positions[0] = 1.0, 0.0;
      _positions[1] = 0.0, 1.0;
      _positions[2] = 0.0, 0.0;
      _positions[3] = 0.5, 0.5;
      _positions[4] = 0.0, 0.5;
      _positions[5] = 0.5, 0.0;

      compute(s);
    }

    //! compute the shape functions and derivatives
    void compute(const CoordinateArray & s);

    PositionContainer nodalCoordinates() {
      return _positions;
    }
  };
}


#endif //#define __ShapeTri6_h__

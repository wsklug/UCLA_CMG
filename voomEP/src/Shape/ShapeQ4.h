// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file ShapeQ4.h

  \brief Class for Linear 4-Node Quadrilateral Finite Element shape functions.

*/

#if !defined(__ShapeQ4_h__)
#define __ShapeQ4_h__

#include "Shape.h"

namespace voom
{
  //! Linear 4-node quadrilateral isoparametric shape functions.
  /*! These 4-node quadrilateral shape functions are to be used for 2-D
    planar problems.  The nodal numbering arrangement
    and the geometry of the standard domain are determined by the
    method <tt>nodalCoordinates()</tt>.
  */
  class ShapeQ4 : public Shape<2> {
    
  public:

    //! constructor
    ShapeQ4(const CoordinateArray & s) {
      const unsigned int n = 4;
      _functions.resize(n);
      _derivatives.resize(n);
      _positions.resize(n);

      _positions[0] = -1.0, -1.0;
      _positions[1] = 1.0, -1.0;
      _positions[2] = 1.0, 1.0;
      _positions[3] = -1.0, 1.0;

      compute(s);
    }

    //! compute the shape functions and derivatives
    void compute(const CoordinateArray & s);

    //! Return parametric coordinates of nodes.
    /*! Nodal arrangement   	\verbatim

    ^ s2        	
    |		
    4--3		
    |  |		
    1--2--> s1			\endverbatim
    
    */
    PositionContainer nodalCoordinates() {
      return _positions;
    }
  };
}


#endif //#define __ShapeTri3_h__

// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
//
//----------------------------------------------------------------------

/*! 
  \file ShapeTet4.h

  \brief Class for Linear 3-Node Triangular Finite Element shape functions.

*/

#if !defined(__ShapeTet4_h__)
#define __ShapeTet4_h__

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
  class ShapeTet4 : public Shape<3> {
    
  public:

    //! Constructor.  Calls compute function.
    ShapeTet4(const CoordinateArray & s) {
      const unsigned int n=4;
      _functions.resize(n);
      _derivatives.resize(n);
      _positions.resize(n);
    
      _positions[0] = 0.0, 0.0, 0.0;
      _positions[1] = 1.0, 0.0, 0.0;
      _positions[2] = 0.0, 1.0, 0.0;
      _positions[3] = 0.0, 0.0, 1.0;

      compute(s);
    }

    //! Compute the shape functions and derivatives.
    void compute(const CoordinateArray & s);

    //! Return parametric coordinates of nodes.
    /*! Nodal arrangement   	\verbatim
      Node 1 = 0.0, 0.0, 0.0;
      Node 2 = 1.0, 0.0, 0.0;
      Node 3 = 0.0, 1.0, 0.0;
      Node 4 = 0.0, 0.0, 1.0;    \endverbatim
    
    */
    PositionContainer nodalCoordinates() {
       return _positions;
       }
  };
}


#endif //#define __ShapeTet4_h__

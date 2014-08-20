// -*- C++ -*-
//----------------------------------------------------------------------
//
//                 William S. Klug, Melissa M. Gibbons
//                University of California Los Angeles
//                   (C) 2004-2007 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
//
//----------------------------------------------------------------------

/*! 
  \file ShapeTet10.h

  \brief Class for Quadratic 10-Node Tetrahedral Finite Element shape functions.

*/

#if !defined(__ShapeTet10_h__)
#define __ShapeTet10_h__

#include "Shape.h"

namespace voom
{
  //! Quadratic 10-node tetrahedral isoparametric shape functions.
  /*! These 10-node tetrahedral shape functions are to be used for 3-D
    problems.  The nodal numbering arrangement
    and the geometry of the standard domain are determined by the
    method <tt>nodalCoordinates()</tt>.
  */
  class ShapeTet10 : public Shape<3> {
    
  public:

    //! constructor
    ShapeTet10(const CoordinateArray & s) {
      const unsigned int n=10;
      _functions.resize(n);
      _derivatives.resize(n);
      _positions.resize(n);
    
      // positions are consistent with the node numbering of the input
      // TetGen generated quadratic tet mesh.

      _positions[0] = 0.0, 0.0, 0.0;
      _positions[1] = 1.0, 0.0, 0.0;
      _positions[2] = 0.0, 1.0, 0.0;
      _positions[3] = 0.0, 0.0, 1.0;

      _positions[4] = 0.5, 0.0, 0.0;
      _positions[5] = 0.5, 0.5, 0.0;
      _positions[6] = 0.0, 0.5, 0.0;

      _positions[7] = 0.0, 0.0, 0.5;
      _positions[8] = 0.5, 0.0, 0.5;
      _positions[9] = 0.0, 0.5, 0.5;


      compute(s);
    }

    //! compute the shape functions and derivatives
    void compute(const CoordinateArray & s);

    //! Return parametric coordinates of nodes.
    /*! Nodal arrangement   	\verbatim
      Node 1 = 0.0, 0.0, 0.0;
      Node 2 = 1.0, 0.0, 0.0;
      Node 3 = 0.0, 1.0, 0.0;
      Node 4 = 0.0, 0.0, 1.0;  
      Node 5 = 0.5, 0.0, 0.0;    
      Node 6 = 0.5, 0.5, 0.0;    
      Node 7 = 0.0, 0.5, 0.0;
      Node 8 = 0.0, 0.0, 0.5;    
      Node 9 = 0.5, 0.0, 0.5;    
      Node 10 = 0.0, 0.5, 0.5;  \endverbatim
    
    */
    PositionContainer nodalCoordinates() {
       return _positions;
       }
  };
}


#endif //#define __ShapeTet10_h__

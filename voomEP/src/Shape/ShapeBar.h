// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2008 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file ShapeBar.h

  \brief Class for Linear 2-node Bar Element

*/

#if !defined(__ShapeBar_h__)
#define __ShapeBar_h__

#include <vector>
#include "voom.h"
#include "Shape.h"


namespace voom
{
  //! 1-D Linear Lagrangian shape functions
  /*! These linear lagrangian are to be used for 1-D bar-like problems.  
    The single parametric coordinate has the range -1 <= xi <= 1.
  */
  class ShapeBar: public Shape<1> {
    
  public:

    //! constructor
    ShapeBar(const CoordinateArray & s) {
      _positions.resize(2);
      _functions.resize(2);
      _derivatives.resize(2);
      compute(s);
    }

    //! Destructor
    ~ShapeBar() {;}

    //! compute the shape functions and derivatives
    void compute(const CoordinateArray & s);
    
    PositionContainer nodalCoordinates() {
      return _positions;
    }

  };
}


#endif //#define __Hermite_h__

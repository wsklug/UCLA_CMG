
#if !defined(__ShapeHex8_h__)
#define __ShapeHex8_h__

#include "Shape.h"

namespace voom
{
  //! Linear 8-node hexahedral isoparametric shape functions.
  /*! These 8-node hexahedral shape functions are to be used for 3-D
    problems.  The nodal numbering arrangement and the geometry of the 
    standard domain are determined by the method <tt>nodalCoordinates()</tt>.
  */
  class ShapeHex8 : public Shape<3> {
    
  public:

    //! constructor
    ShapeHex8(const CoordinateArray & s) {
      const unsigned int n = 8;
      _functions.resize(n);
      _derivatives.resize(n);
      _positions.resize(n);

      _positions[0] = -1.0, -1.0,  1.0;
      _positions[1] =  1.0, -1.0,  1.0;
      _positions[2] =  1.0, -1.0, -1.0;
      _positions[3] = -1.0, -1.0, -1.0;
      _positions[4] = -1.0,  1.0,  1.0;
      _positions[5] =  1.0,  1.0,  1.0;
      _positions[6] =  1.0,  1.0, -1.0;
      _positions[7] = -1.0,  1.0, -1.0;       

      compute(s);
    }

    //! compute the shape functions and derivatives
    void compute(const CoordinateArray & s);
   
    PositionContainer nodalCoordinates() {
      return _positions;
    }

  };
}


#endif 

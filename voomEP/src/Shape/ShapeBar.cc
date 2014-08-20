#include "ShapeBar.h"

namespace voom {

  void ShapeBar::compute(const CoordinateArray & s) {
    _positions[0] = -1.0;
    _positions[1] =  1.0;

    _functions[0] = 0.5*(1. - s[0]);
    _functions[1] = 0.5*(1. + s[0]);
    
    _derivatives[0] = -0.5;
    _derivatives[1] =  0.5;
  }

}

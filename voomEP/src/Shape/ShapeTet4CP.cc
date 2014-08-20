#include "ShapeTet4CP.h"

namespace voom {  
  void ShapeTet4::compute(const CoordinateArray & s) {

    _functions[0] = 1.0-s(0)-s(1)-s(2);
    _functions[1] = s(0);
    _functions[2] = s(1);
    _functions[3] = s(2);

    _derivatives[0]=  -1.0,-1.0,-1.0;
    _derivatives[1]=  1.0, 0.0, 0.0;
    _derivatives[2]=  0.0, 1.0, 0.0;
    _derivatives[3]=  0.0, 0.0, 1.0;
   
  }
};


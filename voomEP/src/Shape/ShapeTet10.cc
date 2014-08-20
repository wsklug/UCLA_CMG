#include "ShapeTet10.h"

namespace voom {  
  void ShapeTet10::compute(const CoordinateArray & s) {

    // again, made consistent with the numbering of TetGen nodes

    _functions[0] = (1.0-s(0)-s(1)-s(2))*(1.0-2.0*s(0)-2.0*s(1)-2.0*s(2));
    _functions[1] = s(0)*(2.0*s(0)-1.0);
    _functions[2] = s(1)*(2.0*s(1)-1.0);
    _functions[3] = s(2)*(2.0*s(2)-1.0);

    _functions[4] = 4.0*s(0)*(1.0-s(0)-s(1)-s(2));
    _functions[5] = 4.0*s(0)*s(1);
    _functions[6] = 4.0*s(1)*(1.0-s(0)-s(1)-s(2));

    _functions[7] = 4.0*s(2)*(1.0-s(0)-s(1)-s(2));
    _functions[8] = 4.0*s(2)*s(0);
    _functions[9] = 4.0*s(1)*s(2);

    // old ordering from Cook book
    
    // _functions[5] = 4.0*s(1)*(1.0-s(0)-s(1)-s(2));
    // _functions[6] = 4.0*s(2)*(1.0-s(0)-s(1)-s(2));
    
    // _functions[7] = 4.0*s(0)*s(1);
    // _functions[8] = 4.0*s(1)*s(2);
    // _functions[9] = 4.0*s(2)*s(0);
    

    _derivatives[0]= -3.0+4.0*s(0)+4.0*s(1)+4.0*s(2), -3.0+4.0*s(0)+4.0*s(1)+4.0*s(2), -3.0+4.0*s(0)+4.0*s(1)+4.0*s(2);
    _derivatives[1] = 4.0*s(0)-1.0, 0.0, 0.0;
    _derivatives[2] = 0.0, 4.0*s(1)-1.0, 0.0;
    _derivatives[3] = 0.0, 0.0, 4.0*s(2)-1.0;

    _derivatives[4] = 4.0*(1.0-2.0*s(0)-s(1)-s(2)), -4.0*s(0), -4.0*s(0);
    _derivatives[5] = 4.0*s(1), 4.0*s(0), 0.0;
    _derivatives[6] = -4.0*s(1), 4.0*(1.0-s(0)-2.0*s(1)-s(2)), -4.0*s(1);

    _derivatives[7] = -4.0*s(2), -4.0*s(2), 4.0*(1.0-s(0)-s(1)-2.0*s(2));
    _derivatives[8] = 4.0*s(2), 0.0, 4.0*s(0);
    _derivatives[9] = 0.0, 4.0*s(2), 4.0*s(1);

    // old ordering from Cook book
    
    // _derivatives[5]= -4.0*s(1), 4.0*(1.0-s(0)-2.0*s(1)-s(2)), -4.0*s(1);
    // _derivatives[6]= -4.0*s(2), -4.0*s(2), 4.0*(1.0-s(0)-s(1)-2.0*s(2));
    
    // _derivatives[7]= 4.0*s(1), 4.0*s(0), 0.0;
    // _derivatives[8]= 0.0, 4.0*s(2), 4.0*s(1);
    // _derivatives[9]= 4.0*s(2), 0.0, 4.0*s(0);
    
  }
};


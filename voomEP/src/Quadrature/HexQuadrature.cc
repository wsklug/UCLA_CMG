
#include <iostream>
#include "HexQuadrature.h"

namespace voom {
  void HexQuadrature::_initialize(unsigned int order)
  {
      
    if( order == 1 ) {  // 1 point 
      _points.resize(1);

      _points[0].coords = 0.0, 0.0, 0.0;				
      _points[0].weight = 8.0;
	
    } else if ( order == 2 ) {  // 8 points
      _points.resize(8);
      
      _points[0].coords = -0.577350269, -0.577350269, -0.57735026;		
      _points[1].coords = 0.577350269, -0.577350269, -0.57735026;			
      _points[2].coords = -0.577350269, 0.577350269, -0.57735026;
      _points[3].coords = 0.577350269, 0.577350269, -0.57735026;
      _points[4].coords = -0.577350269, -0.577350269, 0.57735026;
      _points[5].coords = 0.577350269, -0.577350269, 0.57735026;
      _points[6].coords = -0.577350269, 0.577350269, 0.57735026;
      _points[7].coords = 0.577350269, 0.577350269, 0.57735026;

      _points[0].weight = 1.0;
      _points[1].weight = 1.0;
      _points[2].weight = 1.0;
      _points[3].weight = 1.0;
      _points[4].weight = 1.0;
      _points[5].weight = 1.0;
      _points[6].weight = 1.0;
      _points[7].weight = 1.0;
      
    } else {
      
      std::cout << "HexQuadrature::_initialize(): No quadrature rule is implmented for order " 
		<< order << ".  Quadrature uninitialized." << std::endl;
    }
    return;
  }

} // namespace voom

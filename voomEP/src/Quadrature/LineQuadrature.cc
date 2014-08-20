// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2008 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file QuadQuadrature.cc

  \brief Quadrature rule for a 1D element.

*/

#include <iostream>
#include "LineQuadrature.h"

namespace voom {
  void LineQuadrature::_initialize(unsigned int order)
  {
      
    if( order == 1 ) {  // 1 point 
      _points.resize(1);

      _points[0].coords = 0.0;				
      _points[0].weight = 2.0;
	
    } else if ( order == 3 ) {  // 2 points
      _points.resize(2);

      _points[0].coords = -0.577350269189626;				
      _points[1].coords =  0.577350269189626;
      _points[0].weight = _points[1].weight = 1.0;
	
    } else if ( order == 5 ) { // 3 points
      _points.resize(3);
      
      _points[0].coords = -0.774596669241483;
      _points[1].coords =  0.0;
      _points[2].coords =  0.774596669241483;
 
      _points[0].weight = 
	_points[2].weight = 5.0/9.0;

      _points[1].weight = 8.0/9.0;
    } else if (order == 7 ) { // 4 points
      _points.resize(4);
      _points[0].coords =  0.339981043584856;
      _points[1].coords = -0.339981043584856;
      _points[2].coords =  0.861136311594053;
      _points[3].coords = -0.861136311594053;
     
      _points[0].weight = _points[1].weight = 0.652145154862546;
      _points[2].weight = _points[3].weight = 0.347854845137454;
    } else if (order == 9) { // 5points
      _points.resize(5);
      _points[0].coords = 0.;
      _points[1].coords =  0.538469310105683;
      _points[2].coords = -0.538469310105683;
      _points[3].coords =  0.906179854938664;
      _points[4].coords = -0.906179854938664;

      _points[0].weight = 0.568888888888889;
      _points[1].weight = _points[2].weight = 0.478628670499366;
      _points[3].weight = _points[4].weight = 0.236926885056189;
    } else if (order == 11) { // 6 points
      _points.resize(6);
      _points[0].coords =  0.238619186083197;
      _points[1].coords = -0.238619186083197;
      _points[2].coords =  0.661209386466265;
      _points[3].coords = -0.661209386466265;
      _points[4].coords =  0.932469514203152;
      _points[5].coords = -0.932469514203152;

      _points[0].weight = _points[1].weight = 0.467913934572691;
      _points[2].weight = _points[3].weight = 0.360761573048139;
      _points[4].weight = _points[5].weight = 0.171324492379170;
    } else if (order == 13) { // 7 points
      _points.resize(7);
      _points[0].coords =  0.;
      _points[1].coords =  0.405845151377397;
      _points[2].coords = -0.405845151377397;
      _points[3].coords =  0.741531185599394;
      _points[4].coords = -0.741531185599394;
      _points[5].coords =  0.949107912342759;
      _points[6].coords = -0.949107912342759;

      _points[0].weight = 0.417959183673469;
      _points[1].weight = _points[2].weight = 0.381830050505119;
      _points[3].weight = _points[4].weight = 0.279705391489277;
      _points[5].weight = _points[6].weight = 0.129484966168870;
    } else if (order == 15) { // 8 points
      _points.resize(8);
      _points[0].coords =  0.183434642495690;
      _points[1].coords = -0.183434642495690;
      _points[2].coords =  0.525532409916329;
      _points[3].coords = -0.525532409916329;
      _points[4].coords =  0.796666477413627;
      _points[5].coords = -0.796666477413627;
      _points[6].coords =  0.960289856497536;
      _points[7].coords = -0.960289856497536;

      _points[0].weight = _points[1].weight = 0.362683783378362;
      _points[2].weight = _points[3].weight = 0.313706645877887;
      _points[4].weight = _points[5].weight = 0.222381034453374;
      _points[6].weight = _points[7].weight = 0.101228536290376;
    } else if (order == 17 ) { // 9 points
      _points.resize(9);
      _points[0].coords = 0.;
      _points[1].coords =  0.836031107326636;
      _points[2].coords = -0.836031107326636;
      _points[3].coords =  0.968160239507626;
      _points[4].coords = -0.968160239507626;
      _points[5].coords =  0.3242534234038089;
      _points[6].coords = -0.3242534234038089;
      _points[7].coords =  0.613371432700590;
      _points[8].coords = -0.613371432700590;

      _points[0].weight = 0.330239355001260;
      _points[1].weight = _points[2].weight = 0.180648160694857;
      _points[3].weight = _points[4].weight = 0.081274388361574;
      _points[5].weight = _points[6].weight = 0.312347077040003;
      _points[7].weight = _points[8].weight = 0.260610696402935;
    } else if (order == 19 ) { // 10 points
      _points.resize(10);
      _points[0].coords =  0.148874338981631;
      _points[1].coords = -0.148874338981631;
      _points[2].coords =  0.433395394129247;
      _points[3].coords = -0.433395394129247;
      _points[4].coords =  0.679409568299024;
      _points[5].coords = -0.679409568299024;
      _points[6].coords =  0.865063366688985;
      _points[7].coords = -0.865063366688985;
      _points[8].coords =  0.973906528517172;
      _points[9].coords = -0.973906528517172;
      
      _points[0].weight = _points[1].weight = 0.295524224714753;
      _points[2].weight = _points[3].weight = 0.269266719309996;
      _points[4].weight = _points[5].weight = 0.219086362515982;
      _points[6].weight = _points[7].weight = 0.149451349150581;
      _points[8].weight = _points[9].weight = 0.066671344308688 ;
    } else if (order == 21) { // 11 points
      _points.resize(11);
      _points[0].coords =  0.;
      _points[1].coords =  0.269543155952345;
      _points[2].coords = -0.269543155952345;
      _points[3].coords =  0.519096129110681;
      _points[4].coords = -0.519096129110681;
      _points[5].coords =  0.730152005574049;
      _points[6].coords = -0.730152005574049;
      _points[7].coords =  0.887062599768095;
      _points[8].coords = -0.887062599768095;
      _points[9].coords =  0.978228658146057;
      _points[10].coords = -0.978228658146057;

      _points[0].weight = 0.272925086777901;
      _points[1].weight = _points[2].weight = 0.262804544510247;
      _points[3].weight = _points[4].weight = 0.233193764591990;
      _points[5].weight = _points[6].weight = 0.186290210927734;
      _points[7].weight = _points[8].weight = 0.125580369464905;
      _points[9].weight = _points[10].weight = 0.055668567116174; 
    } else {
      
      std::cout << "LineQuadrature::_initialize(): No quadrature rule is implmented for order " 
		<< order << ".  Quadrature uninitialized." << std::endl;
    }
    return;
  }

  bool LineQuadrature::check(unsigned int d) const {

    // create a random polynomial of degree=d, i.e., 
    // \sum_{i+j<=d} a_{ij} s_1^i s_2^j
    srand(time(0));
    std::vector<double> a;
    for(int i=0; i<=d; i++) {
      a.push_back( 100.0*(static_cast<double>(rand())/RAND_MAX - 0.5) );
    }
    
    double I_exact=0.0;
    for(int i=0; i<=d; i++) {
      I_exact += (a[i]/(i+1))*( 1 + std::pow(-1.0,i) );
    }    

    // Integrate the polynomial P(\xi) numerically
    double I_numerical=0.0;
    for(LineQuadrature::ConstPointIterator p=this->begin(); p!=this->end(); p++) {
      double xi = p->coords(0);
      for(int i=0; i<=d; i++) {
	I_numerical += a[i]*pow(xi,i)*(p->weight);
      }
    }
    
    std::cout << "LineQuadrature::check("<<d<<"):"<<std::endl
	      << "I_exact     = " << I_exact << std::endl
	      << "I_numerical = " << I_exact << std::endl
	      << "Error       = " << std::abs(I_exact-I_numerical)/std::abs(I_exact) 
	      << std::endl;

    double tol = 1.0e-8;
    if ( std::abs(I_numerical-I_exact) <= tol*std::abs(I_exact) ) {  
      std::cout << "LineQuadrature::check("<<d<<") PASSED!"
		<<std::endl;
      return true;
    }
    std::cout << "LineQuadrature::check("<<d<<") FAILED!"
	      <<std::endl;
    
    return false;
  }
} // namespace voom

// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2005 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.2  2005/04/11 04:58:17  klug
// Now log info is added directly to source file.
//
//
//----------------------------------------------------------------------

/*! 
  \file QuadQuadrature.cc

  \brief Quadrature rule for a quadrilateral element.

*/

#include <iostream>
#include "QuadQuadrature.h"

namespace voom {
  void QuadQuadrature::_initialize(unsigned int order)
  {
      
    if( order == 1 ) {  // 1 point 
      _points.resize(1);

      _points[0].coords = 0.0, 0.0;				
      _points[0].weight = 4.0;
	
    } else if ( order == 2 ) {  // 4 points
      _points.resize(4);
      const Real val = 1./sqrt(3);
      _points[0].coords = -val, -val;				
      _points[1].coords = val, -val;				
      _points[2].coords = val, val;
      _points[3].coords = -val, val;
      _points[0].weight = 
	_points[1].weight =
	_points[2].weight =
	_points[3].weight = 1.0;
	
    } else if ( order == 3 ) { // 9 points
      _points.resize(9);   
      const Real val = sqrt(0.6);
      _points[0].coords = -val,-val;
      _points[1].coords =  0.0        ,-val;
      _points[2].coords =  val,-val;
      _points[3].coords = -val, 0.0;
      _points[4].coords =  0.0        , 0.0;
      _points[5].coords =  val, 0.0;
      _points[6].coords = -val, val;
      _points[7].coords =  0.0        , val;
      _points[8].coords =  val, val;
      _points[4].weight = 64.0/81.0;
 
      _points[0].weight =
	_points[2].weight = 
	_points[6].weight = 
	_points[8].weight = 25.0/81.0;

      _points[1].weight =
	_points[3].weight = 
	_points[5].weight = 
	_points[7].weight = 40.0/81.0;
    } else if (order == 5) {
      _points.resize(25);
      const Real r = sqrt(10./7.);
      const Real p1 = 0., p2 = (1./3.)*sqrt( 5. - 2.*r),
	p3 = (1./3.)*sqrt( 5. + 2.*r );
      const Real w1 = 128./225., w2 = (322. + 13.*sqrt(70))/900.,
	w3 = (322. - 13.*sqrt(70))/900.;
      Real p[] = {p1, p2, -p2, p3, -p3};
      Real w[] = {w1, w2, w2, w3, w3};
      for(unsigned int i = 0; i < order; i++) 
	for(unsigned int j = 0; j < order; j++) {
	_points[i*order + j].coords = p[i], p[j];
	_points[i*order + j].weight = w[i]*w[j];
      }
    } else if (order == 6 ) {
      _points.resize(36);
      Real p1 = 0.238619186083197, p2 = 0.661209386466265, 
	p3 = 0.932469514203152;
      Real w1 = 0.467913934572691, w2 = 0.360761573048139,
	w3 = 0.171324492379170;

      Real p[] = {p1, -p1, p2, -p2, p3, -p3};
      Real w[] = {w1, w1, w2, w2, w3, w3};
      for(unsigned int i = 0; i < order; i++) 
	for(unsigned int j = 0; j < order; j++) {
	_points[i*order + j].coords = p[i], p[j];
	_points[i*order + j].weight = w[i]*w[j];
      }
    } else if (order == 7 ) {
      _points.resize(49);
      Real p[] = {0., 0.4058451513, -0.4058451513, 0.7415311855, -0.7415311855,
		  0.9491079123, -0.9491079123};
      Real w[] = {0.4179591836, 0.3818300505, 0.3818300505,
		  0.2797053914, 0.2797053914, 0.1294849661, 0.1294849661};
      for(unsigned int i = 0; i < order; i++) 
	for(unsigned int j = 0; j < order; j++) {
	_points[i*order + j].coords = p[i], p[j];
	_points[i*order + j].weight = w[i]*w[j];
      }
    } else if (order == 8 ) {
      _points.resize(64);
      Real p[] = {0.1834346424,-0.1834346424,0.5255324099,-0.5255324099,
		  0.7966664774,-0.7966664774,0.9602898564,-0.9602898564};
      Real w[] = {0.3626837833,0.3626837833,0.3137066458,0.3137066458,
		  0.2223810344,0.2223810344,0.1012285362,0.1012285362};
      for(unsigned int i = 0; i < order; i++) 
	for(unsigned int j = 0; j < order; j++) {
	_points[i*order + j].coords = p[i], p[j];
	_points[i*order + j].weight = w[i]*w[j];
      }
    } else if (order == 9) {
      _points.resize(81);
      Real p[] = {0.,0.3242534234,-0.3242534234,0.6133714327,-0.6133714327,
		  0.8360311073,-0.8360311073,0.9681602395,-0.9681602395};
      Real w[] = {0.3302393550,0.3123470770,0.3123470770,0.2606106964,
		  0.2606106964,0.1806481606,0.1806481606,0.0812743884,
		  0.0812743883};
      for(unsigned int i = 0; i < order; i++) 
	for(unsigned int j = 0; j < order; j++) {
	_points[i*order + j].coords = p[i], p[j];
	_points[i*order + j].weight = w[i]*w[j];
      }
    } else if (order == 10) {
      _points.resize(100);
      Real p[] = {0.1488743389, -0.1488743389, 0.4333953941, -0.4333953941,
		  0.6794095682, -0.6794095682, 0.8650633666, -0.8650633666,
		  0.9739065285, -0.9739065285};
      Real w[] = {0.2955242247, 0.2955242247, 0.2692667193, 0.2692667193,
		  0.2190863625, 0.2190863625, 0.1494513491, 0.1494513491,
		  0.0666713443, 0.0666713443};
      for(unsigned int i = 0; i < order; i++) 
	for(unsigned int j = 0; j < order; j++) {
	_points[i*order + j].coords = p[i], p[j];
	_points[i*order + j].weight = w[i]*w[j];
      }
    } else if (order == 12) {
      _points.resize(144);
      Real p[] = {0.1252334085,-0.1252334085, 0.3678314989,-0.3678314989, 
		  0.5873179542,-0.5873179542, 0.7699026741, -0.7699026741, 
		  0.9041172563, -0.9041172563, 0.9815606342, -0.9815606342};
      Real w[] = {0.2491470458, 0.2491470458, 0.2334925365, 0.2334925365, 
		  0.2031674267, 0.2031674267, 0.1600783285, 0.1600783285, 
		  0.1069393259, 0.1069393259, 0.0471753363, 0.0471753363};
      for(unsigned int i = 0; i < order; i++) 
	for(unsigned int j = 0; j < order; j++) {
	_points[i*order + j].coords = p[i], p[j];
	_points[i*order + j].weight = w[i]*w[j];
      }
    } else if (order  == 16) {
      _points.resize(order*order);
      Real p[] = {-0.0950125098376374,0.0950125098376374,-0.2816035507792589,
		  0.2816035507792589,-0.4580167776572274,0.4580167776572274,
		  -0.6178762444026438,0.6178762444026438,-0.7554044083550030,
		  0.7554044083550030,-0.8656312023878318,0.8656312023878318,
		  -0.9445750230732326,0.9445750230732326,-0.9894009349916499,
		  0.9894009349916499};
      Real w[] = {0.1894506104550685,0.1894506104550685,0.1826034150449236,
		  0.1826034150449236,0.1691565193950025,0.1691565193950025,
		  0.1495959888165767,0.1495959888165767,0.1246289712555339,
		  0.1246289712555339,0.0951585116824928,0.0951585116824928,
		  0.0622535239386479,0.0622535239386479,0.0271524594117541,
		  0.0271524594117541};
      for(unsigned int i = 0; i < order; i++) 
	for(unsigned int j = 0; j < order; j++) {
	_points[i*order + j].coords = p[i], p[j];
	_points[i*order + j].weight = w[i]*w[j];
      }
    } else {
      
      std::cout << "QuadQuadrature::_initialize(): No quadrature rule is implmented for order " 
		<< order << ".  Quadrature uninitialized." << std::endl;
    }
    return;
  }

  bool QuadQuadrature::check(unsigned int d) const {

    // create a random polynomial of degree=d, i.e., 
    // \sum_{i+j<=d} a_{ij} s_1^i s_2^j
    srand(time(0));
    std::vector<double> a;
    for(int i=0; i<=d; i++) {
      for(int j=0; i+j<=d; j++) {
	a.push_back( 100.0*(static_cast<double>(rand())/RAND_MAX - 0.33) );
      }
    }
    
    // Integrate the polynomial exactly using the formula from Cook's text:
    //
    //   \int_A s_1^i s_2^j (1-s_1-s_2)^k dA = 2A\frac{i!j!m!}{(2+i+j+k)!}
    // 
    // or with k=0
    //
    //   \int_A s_1^i s_2^j dA = 2A\frac{i!j!}{(2+i+j)!}
    //
    // For the standard quadrilateral A=1.0
    
    double I_exact=0.0;
    for(int i=0, k=0; i<=d; i++) 
      for(int j=0; i+j<=d; j++, k++) {
	if ((i%2 != 0) || (j%2 !=0)) continue;
	I_exact += a[k]*(2./(i+1))*(2./(j+1));
      }
    

    // Integrate the polynomial P(s) numerically
    //  \sum_p P(s_p) w_p A
    double I_numerical=0.0;
    for(QuadQuadrature::ConstPointIterator p=this->begin(); p!=this->end(); p++) {
      double s1 = p->coords(0);
      double s2 = p->coords(1);
      for(int i=0, k=0; i<=d; i++) {
	for(int j=0; i+j<=d; j++, k++) {
	  I_numerical += a[k]*pow(s1,i)*pow(s2,j)*(p->weight);
	}
      }
    }
    
    std::cout << "QuadQuadrature::check("<<d<<"):"<<std::endl
	      << "I_exact     = " << I_exact << std::endl
	      << "I_numerical = " << I_numerical << std::endl
	      << "Error       = " << std::abs(I_exact-I_numerical)/std::abs(I_exact) 
	      << std::endl;

    double tol = 1.0e-8;
    if ( std::abs(I_numerical-I_exact) <= tol*std::abs(I_exact) ) {  
      std::cout << "QuadQuadrature::check("<<d<<") PASSED!"
		<<std::endl;
      return true;
    }
    std::cout << "QuadQuadrature::check("<<d<<") FAILED!"
	      <<std::endl;
    
    return false;
  }
} // namespace voom

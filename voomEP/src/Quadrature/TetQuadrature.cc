#include <iostream>
#include "TetQuadrature.h"

namespace voom {
  void TetQuadrature::_initialize(unsigned int order)
  {
    
    if( order == 1 ) {  // 1 point 
      _points.resize(1);

      _points[0].coords = 1.0/4.0, 1.0/4.0, 1.0/4.0;				
      _points[0].weight = 1.0/6.0;
	
    } else 

    if( order == 2) { // 4 points
      _points.resize(4);

      double alpha = 0.58541020;
      double beta = 0.13819660;

      _points[0].coords = alpha, beta,  beta;
      _points[1].coords = beta,  alpha, beta;
      _points[2].coords = beta,  beta,  alpha;
      _points[3].coords = beta,  beta,  beta;
      _points[0].weight = 1.0/4.0/6.0;
      _points[1].weight = 1.0/4.0/6.0;
      _points[2].weight = 1.0/4.0/6.0;
      _points[3].weight = 1.0/4.0/6.0;

    }  else 

    if( order == 3) { // 5 points
      _points.resize(5);

      _points[0].coords = 1.0/4.0, 1.0/4.0, 1.0/4.0;
      _points[1].coords = 1.0/2.0, 1.0/6.0, 1.0/6.0;
      _points[2].coords = 1.0/6.0, 1.0/2.0, 1.0/6.0;
      _points[3].coords = 1.0/6.0, 1.0/6.0, 1.0/2.0;
      _points[4].coords = 1.0/6.0, 1.0/6.0, 1.0/6.0;
      _points[0].weight = -4.0/5.0/6.0;
      _points[1].weight = 9.0/20.0/6.0;
      _points[2].weight = 9.0/20.0/6.0;
      _points[3].weight = 9.0/20.0/6.0;
      _points[4].weight = 9.0/20.0/6.0;

    } 
    return;
  }

  bool TetQuadrature::check(unsigned int d) const {

    TetQuadrature quad(*this);

    // create a random polynomial of degree=d, i.e., 
    // \sum_{i+j+k<=d} a_{ijk} s_1^i s_2^j s_3^k
    srand(time(0));
    std::vector<double> a;
    for(int i=0, q=0; i<=d; i++) {
      for(int j=0; i+j<=d; j++) {
	for(int k=0; i+j+k<=d; k++, q++) {
	  a.push_back( 100.0*(static_cast<double>(rand())/RAND_MAX - 0.33) );
	  // std::cout << "a(i=" << i << ",j=" << j << ",k=" << k << ")=" << a[q] << "with q=" << q << std::endl;
	}
      }
    }
    
    double I_exact=0.0;
    for(int i=0, q=0; i<=d; i++) {
      for(int j=0; i+j<=d; j++) {
        for(int k=0; i+j+k<=d; k++, q++) {
	  I_exact += a[q]*_factorial(i)*_factorial(j)*_factorial(k)/_factorial(3+i+j+k);
	  // std::cout << "contribution to I_exact with (ijk) at (" << i << "," << j << "," << k << ") is: " << a[q]*_factorial(i)*_factorial(j)*_factorial(k)/_factorial(3+i+j+k) << std::endl;
	}
      }
    }    
    // Integrate the polynomial P(s) numerically
    //  \sum_p P(s_p) w_p A
    double I_numerical=0.0;
    for(TetQuadrature::ConstPointIterator p=quad.begin(); p!=quad.end(); p++) {
      double s1 = p->coords(0);
      double s2 = p->coords(1);
      double s3 = p->coords(2);
      for(int i=0, q=0; i<=d; i++) {
	for(int j=0; i+j<=d; j++) {
	  for(int k=0; i+j+k<=d; k++, q++) {
	    I_numerical += a[q]*pow(s1,i)*pow(s2,j)*pow(s3,k)*(p->weight);
	    // std::cout << "contr to I_num w/ijk " << i << "," << j << "," << k << ") is:" << a[q]*pow(s1,i)*pow(s2,j)*pow(s3,k)*(p->weight)/6.0 << "q=" << q << " a=" << a[q]  << std::endl;
	  }
        }
      }
    }
      I_numerical *= 1.0/6.0;
      // NOTE: this was I_numerical *= 0.5 for the triangle.
    
      std::cout << "TetQuadrature::check("<<d<<"):"<<std::endl
		<< "I_exact     = " << I_exact << std::endl
		<< "I_numerical = " << I_numerical << std::endl
		<< "Error       = " << std::abs(I_exact-I_numerical)/std::abs(I_exact) 
		<< std::endl;
	       
      double tol = 1.0e-7;
      if ( std::abs(I_numerical-I_exact) <= tol*std::abs(I_exact) ) {  
	std::cout << "TetQuadrature::check("<<d<<") PASSED!"
		  <<std::endl;
	return true;
      }
      std::cout << "TetQuadrature::check("<<d<<") FAILED!"
		<<std::endl;
    
      return false;
    }
  //  }
  
} // namespace voom




#if !defined(__HexQuadrature_h__)
#define __HexQuadrature_h__

#include "Quadrature.h"

namespace voom
{
  class HexQuadrature
    :public Quadrature<3>
  {

  public:

    //! default constructor
    HexQuadrature() { _initialize(1); }
		
    HexQuadrature(unsigned int order) {_initialize(order);}
    
    //! assignment operator
    HexQuadrature & operator = (const HexQuadrature & q) {
      if( this != &q ) {
	Quadrature<3>::operator=(q);
      }
      return *this;
    }

    //! destructor
    virtual ~HexQuadrature() {}
    bool check(unsigned int order) const;

  private:
    //! initialize gauss points for different rules
    void _initialize(unsigned int order);  

    //! recursive funtion to compute factorial
    static int _factorial(int n) {
      int temp;
      if(n <= 1) return 1;
      temp = n * _factorial(n - 1);
      return temp;
    }
  };
	
} // namespace voom

#endif 

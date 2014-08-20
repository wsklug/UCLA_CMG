// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2008 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file LineQuadrature.h

  \brief Quadrature rule for a line element.

*/

#if !defined(__LineQuadrature_h__)
#define __LineQuadrature_h__

#include "Quadrature.h"

namespace voom
{
  //! Concrete Class representing 1-D Gaussian quadrature over a linear domain.
  /*! This class derives from <tt> Quadrature </tt> defining the
    method <tt>_initialize()</tt> for several quadrature orders.  In
    addition, this class provides a method <tt>check()</tt>
    which integrates exactly a random polynomial of the appropriate
    order and compares the result to that of quadrature.
  */
  class LineQuadrature
    :public Quadrature<1>
  {

  public:

    //! default constructor
    LineQuadrature() { _initialize(1); }
		
    LineQuadrature(unsigned int order) {_initialize(order);}
    
    //! assignment operator
    LineQuadrature & operator = (const LineQuadrature & q) {
      if( this != &q ) {
	Quadrature<1>::operator=(q);
      }
      return *this;
    }

    //! destructor
    virtual ~LineQuadrature() {}
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

#endif // __LineQuadrature_h__

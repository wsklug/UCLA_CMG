// -*- C++ -*-
//----------------------------------------------------------------------
//
//                 William S. Klug, Melissa M. Gibbons
//                University of California Los Angeles
//                 (C) 2004-2006 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
//
//----------------------------------------------------------------------

/*! 
  \file TetQuadrature.h

  \brief Quadrature rule for a tetrahedral element.

*/

#if !defined(__TetQuadrature_h__)
#define __TetQuadrature_h__

#include "Quadrature.h"

namespace voom
{
  //! Concrete Class representing 3-D Gaussian quadrature over a tetrahedral domain.
  /*! This class derives from <tt> Quadrature </tt> defining the
    method <tt>_initialize()</tt> for several quadrature orders.  In
    addition, this class provides a method <tt>check()</tt>
    which integrates exactly a random polynomial of the appropriate
    order and compares the result to that of quadrature.
  */
  class TetQuadrature
    :public Quadrature<3>
  {

  public:

    //! default constructor
    TetQuadrature() { _initialize(1); }
		
    //! Construct from a user-supplied quadrature order (1, 2, or 3)
    TetQuadrature(unsigned int order) {_initialize(order);}
    
    //! assignment operator
    TetQuadrature & operator = (const TetQuadrature & q) {
      if( this != &q ) {
	Quadrature<3>::operator=(q);
      }
      return *this;
    }

    //! destructor
    virtual ~TetQuadrature() {}
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

#endif // __TetQuadrature_h__

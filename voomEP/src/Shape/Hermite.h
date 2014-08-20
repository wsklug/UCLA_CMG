// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2008 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file Hermite.h

  \brief Class for cubic Hermitian shape functions

*/

#if !defined(__Hermite_h__)
#define __Hermite_h__

#include <vector>
#include "voom.h"


namespace voom
{
  //! 1-D Cubic Hermitian shape functions.
  /*! These cubic Hermitian shape functions are to be used for
    C0-continuous 1-D beam-like problems.  The single parametric
    coordinate has the range -1 <= xi <= 1.
  */
  class Hermite {
    
  public:

    typedef std::vector<double> FunctionArray;
    typedef std::vector<double> DerivativeArray;

    //! constructor
    Hermite(double xi) {
      const unsigned int n = 4;
      _functions.resize(n);
      _derivatives.resize(n);
      _2derivatives.resize(n);

      compute(xi);
    }

    //! return shape functions
    const FunctionArray & functions() const {return _functions;}
	
    //! return derivatives of shape functions
    const DerivativeArray & derivatives() const {return _derivatives;}

    //! return derivatives of shape functions
    const DerivativeArray & secondDerivatives() const {return _2derivatives;}

    //! compute the shape functions and derivatives
    void compute(double xi);

  protected:

    //! shape functions
    FunctionArray _functions;

    //! derivatives of shape functions 
    DerivativeArray _derivatives;

    //! second derivatives of shape functions 
    DerivativeArray _2derivatives;

  };
}


#endif //#define __Hermite_h__

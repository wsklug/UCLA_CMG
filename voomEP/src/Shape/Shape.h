// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file Shape.h

  \brief Templated base class for isoparametric Finite Element shape functions.

*/

#if !defined(__Shape_h__)
#define __Shape_h__

#include <vector>
#include "voom.h" 


namespace voom
{
  //! Templated base class for isoparametric Finite Element shape functions.
  /*! Objects derived from the Shape class represent isoparametric
    shape functions for \f$C^0\f$-conforming finite elements.  It is
    assumed that the mapping from parametric coordinates to physical
    coordinates (and its jacobian) will be computed by the client,
    thus the Shape object needs not deal with the embedding space.
    The parametric dimension is specified by the first template
    parameter.  

    The Shape template specifies the interface through which all
    isoparametric shape function objects are to be manipulated.
    Objects can be constructed by passing the parametric coordinates
    of the evaluation point.  The constructor immediately calls the
    purely virtual method <tt>compute()</tt>, which must be
    implemented by concrete derived classes to evaluate and store the
    shape functions and their parametric derivatives

    Note that specializations of the Shape template are purely virtual
    because of the purely virtual method <tt>compute()</tt>.  This
    public method can also be used directly by clients to reset the
    evaluation point and recompute the shape functions and
    derivatives.  The shape functions and their derivatives are stored
    in arrays which can be accessed by the <tt>functions()</tt> and
    <tt>derivatives()</tt> methods.  

    The Shape template class also provides two verification tests for
    derived classes.  The <tt>checkConsistency()</tt> method computed
    derivatives numerically and compares these to the directly
    computed values.  The <tt>checkC0Completeness()</tt> tests that
    the interpolation functions can exactly represent any linear
    polynomial in the parametric dimension.
  */    
  template<std::size_t PARAM_DIM>
  class Shape {
    
  public:
    typedef std::vector<double> FunctionContainer;
    typedef typename tvmet::Vector<double, PARAM_DIM> CoordinateArray;
    typedef std::vector<CoordinateArray> DerivativeContainer;
    typedef std::vector<CoordinateArray> PositionContainer;

    virtual ~Shape() {}

    //! return shape functions
    const FunctionContainer & functions() const {return _functions;}
	
    //! return derivatives of shape functions
    const DerivativeContainer & derivatives() const {return _derivatives;}

    //! compute the shape functions and derivatives
    virtual void compute(const CoordinateArray & s) = 0;

    //! check the consistency of the shape function derivatives

  protected:

    //! shape functions (size Node_n)
    FunctionContainer _functions;

    //! derivatives of shape functions 
    DerivativeContainer _derivatives;

    //! nodal positions corresponding to shape functions 
    PositionContainer _positions;

    //! jacobian of the isoparametric mapping
    double _jacobian;

  public:

    //! Verify linear completeness of shape functions.
    bool checkC0Completeness();

    //! Verify consistency of shape function derivatives
    bool checkConsistency(CoordinateArray s);

    virtual PositionContainer nodalCoordinates() = 0; 
    
  };
}

#include "Shape.icc"

#endif //#define __Shape_h__

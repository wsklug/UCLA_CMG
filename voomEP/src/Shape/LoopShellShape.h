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
// Revision 1.7  2005/04/21 01:49:12  klug
// Changed and simplified interface.
//
//----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////
//
// 
//         Loop thin shell shape functions
//         for triangle elements
//         derived from Shape class
//
//
////////////////////////////////////////////////////////////////////////////////////////

#ifndef _LOOPSHELLSHAPE_
#define _LOOPSHELLSHAPE_

#include <iostream>
#include "voom.h"
#include "VoomMath.h"

namespace voom
{
  //
  //! class defination
  class LoopShellShape {
  public:

    typedef tvmet::Vector<double, 2> CoordinateArray;

    typedef blitz::Array<double,1> FunctionArray;
    typedef blitz::Array<double,2> DerivativeArray;
    typedef blitz::Array<double,3> SecondDerivativeArray;

    typedef blitz::Array<double,2> SubdivisionMatrix;

    typedef tvmet::Vector<unsigned int, 3> CornerValences;
	
  protected:
    //! number of nodes
    int _nodes;

    //! shape functions
    FunctionArray  _functions;
	  
    //! derivatives of shape functions
    DerivativeArray  _derivatives;

    //! the second derivatives of shape functions
    SecondDerivativeArray _secondDerivatives;

   //! parametric coords
    CoordinateArray _coords;	
	
 
  public:

    /*! default constructor function for 12 nodes regular element. 2D
     *  parametric coord (s1 and s2) assumes barycentric gauss point
     */    
    LoopShellShape() {
      _nodes = 12;
      _functions.resize(_nodes);
      _derivatives.resize(_nodes);
      _secondDerivatives.resize(_nodes);
      CornerValences V(6, 6, 6);
      _initialize(12, V);
    };

    //! constructor assumes barycentric gauss point (1.0/3.0, 1.0/3.0)
    LoopShellShape(const int nodes,
		   const CornerValences& V) {
      _nodes = nodes;
      _functions.resize(_nodes);
      _derivatives.resize(_nodes,2);
      _secondDerivatives.resize(_nodes,2,2);
      
      _initialize(nodes, V);
    }

    //! constructor function with parrametric coords
    LoopShellShape(const int nodes,
		   const CornerValences & V,
		   const CoordinateArray & paraCoords) {
      _nodes = nodes;
      _functions.resize(_nodes);
      _derivatives.resize(_nodes,2);
      _secondDerivatives.resize(_nodes,2,2);
      _initialize(nodes, V, paraCoords);
    }
	
    //! default destructor function
    virtual ~LoopShellShape() {;}	

  private:
    //! pure virtual function:
    //! return shape functions
    void _computeFunctions( const SubdivisionMatrix & S);

    //! reutrn 1st derivatives of shape functions
    void _computeDerivatives( const SubdivisionMatrix & S);

    //! return 2nd derivatives of shape functions
    void _computeSecondDerivatives( const SubdivisionMatrix & S);
	
    //! compute sub-division matrix
    void _computeSubdivisionMatrix( const CornerValences& V, 
				   SubdivisionMatrix & S, 
				   bool needSubdivide );

    //! compute the parametric coords for sub-patch (irregular elements)
    void _convertParaCoords();

    //!
    void _initialize(const int nodes, const CornerValences & V, 
		    const CoordinateArray & paraCoords);
    void _initialize(const int nodes, const CornerValences & V) {
      CoordinateArray paraCoords;
      paraCoords = 1.0/3.0, 1.0/3.0;
      _initialize(nodes, V, paraCoords);
    }

  public:

    //! return shape functions
    const FunctionArray & functions() const {
      return _functions;
    }
		
    //! return derivatives of shape functions
    const DerivativeArray & derivatives() const {
      return _derivatives;
    }

    //! return 2nd derivatives of shape functions
    const SecondDerivativeArray & secondDerivatives() const {
      return _secondDerivatives;
    }

    //! verify the correction of the shape functions
    static void checkShapeFunctions();

    //! check shape functions numerically
    void checkShapeNumerically(int n,tvmet::Vector<unsigned,3> V);
	  
  };

}

#endif  //  _LOOPSHELLSHAPE_

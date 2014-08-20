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
// Revision 1.4  2005/06/27 03:49:06  klug
// Added const to some methods that needed it.
//
// Revision 1.3  2005/04/11 04:56:04  klug
// Now log info is added directly to source file.
//
// Revision 1.2  2005/04/11 04:53:45  klug
// *** empty log message ***
//
//

/*! 
  \file Quadrature.h

  \brief Virtual Base Class for the concept of a quadrature rule.

*/

#if !defined(__Quadrature_h__)
#define __Quadrature_h__

#include <vector>
#include "voom.h" 

namespace voom
{
	
  //! Virtual Base Class for the concept of a quadrature rule.

  /*! <tt>Quadrature</tt> objects provide storage of and access to a
    set of quadrature points, each of which is represented by a
    <tt>Point</tt> structure, containing coordinates and a weight.
    This set of points is stored in an STL vector, which clients can
    access through the method <tt>points()</tt>.  For convenience,
    clients can iterate through the vector of points by accessing
    iterators to the begining and end of the vector with the
    <tt>begin()</tt> and <tt>end()</tt> methods.

    This class is purely virtual due to the absence of a default
    implementation of the protected <tt>_initialize()</tt> method.
    This method must be implemented by concrete derived classes, to
    set the coordinates and weights of the quadrature rule, for a
    requested quadrature order.  This order determines the highest
    degree polynomial which is integrated exactly by the quadrature
    rule.

    The one template parameter of this class determines the dimension
    of the parametric space, i.e., the number of independent
    coordinates determining the position of a quadrature point.
    
  */
  template< std::size_t dim_n >
  class Quadrature
  {

  public:

    //! structure containing the coordinates and weight of a quadrature point
    /*! Point is a structure, so its members are public meaning
      clients can access them directly.
    */
    struct Point 
    {
      typedef typename tvmet::Vector<double,dim_n> CoordinateArray;      
      CoordinateArray coords;
      double weight;
    };
    
    typedef typename std::vector<Point> PointContainer;
    typedef typename std::vector<Point>::const_iterator ConstPointIterator;

    //! default constructor
    Quadrature() {}

    //! copy constructor 
    Quadrature(const Quadrature& q) {
      if( this != &q ) {
	_copy(q);
      }
    }
    
    //! assignment operator
    Quadrature & operator=(const Quadrature & q) {
      if( this != &q ) {
	_copy(q);
      }
      return *this;
    }
    
    //! destructor
    virtual ~Quadrature() {;}
    
    //! access quadrature points
    const PointContainer& points() const { return _points; }
    
    //! get iterator to first quadrature point
    ConstPointIterator begin() const { return _points.begin(); }
    
    //! get iterator to last quadrature point
    ConstPointIterator end() const { return _points.end(); }

    //! Size of quadrature points
    int size() const { return _points.size(); }

  protected:

    //! the quadrature points.
    PointContainer _points;

    //! initialize gauss points for different rules
    virtual void _initialize(unsigned int order) = 0;

    virtual void _copy(const Quadrature & q) {
      _points.clear();
      for(ConstPointIterator p=q._points.begin(); p!=q._points.end(); p++ )
	_points.push_back( *p );
    }

  };
	
} // namespace voom

#endif // __QuadratureRule_h__

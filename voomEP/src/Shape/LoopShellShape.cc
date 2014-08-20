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
// Revision 1.11  2005/04/21 01:49:12  klug
// Changed and simplified interface.
//
//----------------------------------------------------------------------

#include "LoopShellShape.h"

//#define DEBUG_SUB

namespace voom
{
  //! constructor function with parrametric coords
  void LoopShellShape::_initialize(const int nodes,
				  const CornerValences & V,
				  const CoordinateArray & paraCoords)
  {
    //
    // need subdivision
    bool needSubdivide = true;
    if ( V(0) == 6 && V(1)==6 && V(2) == 6 ) needSubdivide = false;

    //
    // debug
#ifdef DEUBG_SUB
    needSubdivide = true;
#endif
    // end of debug
    //
		
    _coords = paraCoords;

    // allocate memory for 2ndDerivs and set all components equal to 0
    _functions = 0.0;
    _derivatives = 0.0;
    _secondDerivatives = 0.0;

    // sub-division matrix
    SubdivisionMatrix  S(12, nodes);
    _computeSubdivisionMatrix(V, S, needSubdivide);
    //! convert the parametric coords for sub-patch (irregular patches)
    if( needSubdivide ) _convertParaCoords();
		
    // compute shape functions
    _computeFunctions( S );

    // compute the 1st derivatives of shape functions
    _computeDerivatives( S );
    // if subidivide, apply chain rule
    if ( needSubdivide ) _derivatives *= -2.0;
		
    // compute the 2nd derivatives of shape functions
    _computeSecondDerivatives( S );
    if ( needSubdivide ) _secondDerivatives *= 4.0;
	
  }


  void LoopShellShape::_computeFunctions( const Array2D& SDMatrix )
  {
    // set the parametric coords
    const double v = _coords(0);
    const double w = _coords(1);
    const double u = 1 - v - w;

    FunctionArray boxSplines(12);
    // computing ...
    // 12 shape functions for the regular patch element
    //
    boxSplines(0)=(u*u*u*u + 2.0*u*u*u*v)/12.0;

    boxSplines(1)=(u*u*u*u + 2.0*u*u*u*w)/12.0;

    boxSplines(2)=(u*u*u*u + 2.0*u*u*u*w + 6.0*u*u*u*v + 6.0*u*u*v*w +
		   12.0*u*u*v*v + 6.0*u*v*v*w + 6.0*u*v*v*v + 2.0*v*v*v*w +
		   v*v*v*v)/12.0;

    boxSplines(3)=(6.0*u*u*u*u + 24.0*u*u*u*w + 24.0*u*u*w*w + 8.0*u*w*w*w + 
		   w*w*w*w + 24.0*u*u*u*v + 60.0*u*u*v*w + 36.0*u*v*w*w + 
		   6.0*v*w*w*w + 24.0*u*u*v*v + 36.0*u*v*v*w + 12.0*v*v*w*w + 
		   8.0*u*v*v*v + 6.0*v*v*v*w + v*v*v*v)/12.0;

    boxSplines(4)=(u*u*u*u + 6.0*u*u*u*w + 12.0*u*u*w*w + 6.0*u*w*w*w + 
		   w*w*w*w + 2.0*u*u*u*v + 6.0*u*u*v*w + 6.0*u*v*w*w + 
		   2.0*v*w*w*w)/12.0;

    boxSplines(5)=(2.0*u*v*v*v + v*v*v*v)/12.0;

    boxSplines(6)=(u*u*u*u + 6.0*u*u*u*w + 12.0*u*u*w*w + 
		   6.0*u*w*w*w + w*w*w*w + 8.0*u*u*u*v + 36.0*u*u*v*w + 
		   36.0*u*v*w*w + 8.0*v*w*w*w + 24.0*u*u*v*v + 60.0*u*v*v*w + 
		   24.0*v*v*w*w + 24.0*u*v*v*v + 24.0*v*v*v*w + 
		   6.0*v*v*v*v)/12.0;

    boxSplines(7)=(u*u*u*u + 8.0*u*u*u*w + 24.0*u*u*w*w + 24.0*u*w*w*w + 
		   6.0*w*w*w*w + 6.0*u*u*u*v + 36.0*u*u*v*w + 60.0*u*v*w*w + 
		   24.0*v*w*w*w + 12.0*u*u*v*v + 36.0*u*v*v*w + 
		   24.0*v*v*w*w + 6.0*u*v*v*v + 8.0*v*v*v*w + v*v*v*v)/12.0;

    boxSplines(8) =(2.0*u*w*w*w + w*w*w*w)/12.0; 
    boxSplines(9)=(2.0*v*v*v*w + v*v*v*v)/12.0;


    boxSplines(10)=(2.0*u*w*w*w + w*w*w*w + 6.0*u*v*w*w + 6.0*v*w*w*w + 
		    6.0*u*v*v*w + 12.0*v*v*w*w + 2.0*u*v*v*v + 
		    6.0*v*v*v*w + v*v*v*v)/12.0;

    boxSplines(11)=(w*w*w*w + 2.0*v*w*w*w)/12.0;

    // computing the shape functions for irregular patches
    blitz::firstIndex i; 
    blitz::secondIndex j;  
    _functions = sum (boxSplines(j) * SDMatrix(j,i), j);
		
  }



  void LoopShellShape::_computeDerivatives( const Array2D& SDMatrix )
  {
    // set the parametric coords
    const double v = _coords(0);
    const double w = _coords(1);
    const double u = 1 - v - w;

    Array2D bsDerivatives(12,2);
    // computing ....
    //
    // 12 * 2 components of the 1st derivatives of the shape function for the
    // regular element
    bsDerivatives(0,0) = (-6.0*v*u*u - 2.0*u*u*u)/12.0;
    bsDerivatives(0,1) = (-6.0*v*u*u - 4.0*u*u*u)/12.0;

    bsDerivatives(1,0) = (-4.0*u*u*u-6.0*u*u*w)/12.0;
    bsDerivatives(1,1) = (-2.0*u*u*u-6.0*u*u*w)/12.0;

    bsDerivatives(2,0) = (-2.0*v*v*v-6.0*v*v*u
			  + 6.0*v*u*u+2.0*u*u*u)/12.0;
    bsDerivatives(2,1) = (-4.0*v*v*v-18.0*v*v*u
			  - 12.0*v*u*u-2.0*u*u*u
			  - 6.0*v*v*w-12.0*v*u*w
			  - 6.0*u*u*w)/12.0;

    bsDerivatives(3,0) = (-4.0*v*v*v-24.0*v*v*u
			  - 24.0*v*u*u-18.0*v*v*w 
			  - 48.0*v*u*w-12.0*u*u*w
			  - 12.0*v*w*w - 12.0*u*w*w
			  - 2.0*w*w*w)/12.0;

    bsDerivatives(3,1) = (-2.0*v*v*v-12.0*v*v*u
			  - 12.0*v*u*u-12.0*v*v*w
			  - 48.0*v*u*w-24.0*u*u*w
			  - 18.0*v*w*w-24.0*u*w*w
			  - 4.0*w*w*w)/12.0;

    bsDerivatives(4,0) = (-6.0*v*u*u-2.0*u*u*u
			  - 12.0*v*u*w-12.0*u*u*w
			  - 6.0*v*w*w-18.0*u*w*w
			  - 4.0*w*w*w)/12.0;

    bsDerivatives(4,1) = (2.0*u*u*u+6.0*u*u*w
			  - 6.0*u*w*w-2.0*w*w*w)/12.0;

    bsDerivatives(5,0) = (2.0*v*v*v+6.0*v*v*u)/12.0;
    bsDerivatives(5,1) = -v*v*v/6.0;

    bsDerivatives(6,0) = (24.0*v*v*u+24.0*v*u*u
			  + 4.0*u*u*u+12.0*v*v*w
			  + 48.0*v*u*w+18.0*u*u*w
			  + 12.0*v*w*w+12.0*u*w*w
			  + 2.0*w*w*w)/12.0;
  
    bsDerivatives(6,1) = (12.0*v*v*u+12.0*v*u*u
			  + 2.0*u*u*u-12.0*v*v*w
			  + 6.0*u*u*w-12.0*v*w*w
			  - 6.0*u*w*w-2.0*w*w*w)/12.0;

    bsDerivatives(7,0) = (-2.0*v*v*v-6.0*v*v*u
			  + 6.0*v*u*u+2.0*u*u*u
			  - 12.0*v*v*w+12.0*u*u*w
			  - 12.0*v*w*w+12.0*u*w*w)/12.0;

    bsDerivatives(7,1) = (2.0*v*v*v+12.0*v*v*u
			  + 18.0*v*u*u+4.0*u*u*u
			  + 12.0*v*v*w+48.0*v*u*w
			  + 24.0*u*u*w+12.0*v*w*w
			  + 24.0*u*w*w)/12.0;

    bsDerivatives(8,0) = -w*w*w/6.0;
    bsDerivatives(8,1) = (6.0*u*w*w+2.0*w*w*w)/12.0;

    bsDerivatives(9,0) = (4.0*v*v*v+6.0*v*v*w)/12.0;
    bsDerivatives(9,1) = v*v*v/6.0;

    bsDerivatives(10,0) = (2.0*v*v*v+6.0*v*v*u
			   + 12.0*v*v*w+12.0*v*u*w
			   + 18.0*v*w*w+6.0*u*w*w
			   + 4.0*w*w*w)/12.0;

    bsDerivatives(10,1)= (4.0*v*v*v+6.0*v*v*u
			  + 18.0*v*v*w+12.0*v*u*w
			  + 12.0*v*w*w+6.0*u*w*w
			  + 2.0*w*w*w)/12.0;

    bsDerivatives(11,0) = w*w*w/6.0;
    bsDerivatives(11,1) = (6.0*v*w*w+4.0*w*w*w)/12.0;
    //
    // computing the derivatives of shape function of the irregular pathc
    blitz::firstIndex i; 
    blitz::secondIndex j;
    blitz::thirdIndex k;
		
    _derivatives = sum (bsDerivatives(k, j) * SDMatrix(k,i), k);
    // Chain rule to compute 2nd derivatives w.r.t. the curvilinear
    // coords in original patch
    //_derivatives *= (-2.0);
    return;
		
  }



  void LoopShellShape::_computeSecondDerivatives( const Array2D& SDMatrix )
  {
    /* second order derivatives of the box spline shape functions */
    /* der( ,0) derivative with respect to vv                     */
    /* der( ,1) derivative with respect to ww                     */
    /* der( ,2) derivative with respect to vw                     */
    const double v = _coords(0);
    const double w = _coords(1);
    const double u = 1 - v - w;
  
    Array3D bs2ndDerivatives(12,2,2);

    bs2ndDerivatives(0,0,0) = v*u;
    bs2ndDerivatives(0,1,1) = v*u+u*u;
    bs2ndDerivatives(0,0,1) = (12.0*v*u+6.0*u*u)/12.0;

    bs2ndDerivatives(1,0,0) = u*u+u*w;
    bs2ndDerivatives(1,1,1) = u*w;
    bs2ndDerivatives(1,0,1) = (6.0*u*u+12.0*u*w)/12.0;
             
    bs2ndDerivatives(2,0,0) = -2.0*v*u;
    bs2ndDerivatives(2,1,1) = v*v+v*u+v*w+u*w;
    bs2ndDerivatives(2,0,1) = (6.0*v*v-12.0*v*u
			       -6.0*u*u)/12.0;
             
    bs2ndDerivatives(3,0,0) = v*v-2.0*u*u
      + v*w-2.0*u*w;
    bs2ndDerivatives(3,1,1) = -2.0*v*u-2.0*u*u
      + v*w+w*w;
    bs2ndDerivatives(3,0,1) = (6.0*v*v-12.0*u*u
			       + 24.0*v*w+6.0*w*w)/12.0;
             
    bs2ndDerivatives(4,0,0) = v*u+v*w+u*w+ w*w;
    bs2ndDerivatives(4,1,1) = -2.0*u*w;
    bs2ndDerivatives(4,0,1) = (-6.0*u*u-12.0*u*w 
			       + 6.0*w*w)/12.0;
             
    bs2ndDerivatives(5,0,0) = v*u;
    bs2ndDerivatives(5,1,1) = 0.0;
    bs2ndDerivatives(5,0,1) = -v*v/2.0;
             
    bs2ndDerivatives(6,0,0) = (-24.0*v*v+12.0*u*u-24.0*v*w
			       + 12.0*u*w)/12.0;
    bs2ndDerivatives(6,1,1) = (-24.0*v*v-24.0*v*u-24.0*v*w
			       - 24.0*u*w)/12.0;
    bs2ndDerivatives(6,0,1) = (-12.0*v*v+6.0*u*u-24.0*v*w
			       - 12.0*u*w-6.0*w*w)/12.0;
             
    bs2ndDerivatives(7,0,0) = -2.0*v*u-2.0*v*w-2.0*u*w- 2.0*w*w;
    bs2ndDerivatives(7,1,1) = v*u+u*u-2.0*v*w - 2.0*w*w;
    bs2ndDerivatives(7,0,1) = (-6.0*v*v-12.0*v*u+6.0*u*u 
			       - 24.0*v*w-12.0*w*w)/12.0;
             
    bs2ndDerivatives(8,0,0) = 0.0;
    bs2ndDerivatives(8,1,1) = u*w;
    bs2ndDerivatives(8,0,1) = -w*w/2.0; 
             
    bs2ndDerivatives(9,0,0) = (12.0*v*v+12.0*v*w)/12.0;
    bs2ndDerivatives(9,1,1) = 0.0;
    bs2ndDerivatives(9,0,1) = v*v/2.0;
             
    bs2ndDerivatives(10,0,0)= (12.0*v*u+12.0*v*w+12.0*u*w
			       + 12.0*w*w)/12.0;
    bs2ndDerivatives(10,1,1)= v*v+v*u+v*w+u*w;
    bs2ndDerivatives(10,0,1)= (6.0*v*v+12.0*v*u+24.0*v*w 
			       + 12.0*u*w+6.0*w*w)/12.0;
             
    bs2ndDerivatives(11,0,0)= 0.0;
    bs2ndDerivatives(11,1,1)= v*w+w*w;
    bs2ndDerivatives(11,0,1)= w*w/2.0;

    for(int i=0; i<12; i++) bs2ndDerivatives(i,1,0) = bs2ndDerivatives(i,0,1);
    //
    // computing the 2nd derivatives of the shape functions
    blitz::firstIndex i; 
    blitz::secondIndex j;
    blitz::thirdIndex k;
    blitz::fourthIndex l;
		
    _secondDerivatives = sum (bs2ndDerivatives(l, j, k) * SDMatrix(l,i), l);
    // Chain rule to compute 2nd derivatives w.r.t. the curvilinear
    // coords in original patch
    //_2ndDerivs *= 4.0;

	
    return;	
  }


  void LoopShellShape::_computeSubdivisionMatrix( const CornerValences & V, 
						 SubdivisionMatrix & S, 
						 bool needSubdivide)
  {
    // zeroize subdivision matrix		
    S = 0.0;
    //
    //  to understand this part, need the formula of Transformation Matrix 
    //  in the document.
    //
    if ( needSubdivide ){
      int N0 = V(0)-2,  N1 = V(1)-2,  N2 = V(2)-2;
  
      if(false){
	std::cout << "N0 = " << N0 << std::endl
		  << "N1 = " << N1 << std::endl
		  << "N2 = " << N2 << std::endl;
      }
      ///////////////////////////////////////////////////////////////////////
      //
      //  if the vertices of elements surrounding one element are overlapped,
      //  need new method to construct the Transformation Matrix.
      //
      assert( N0 != 1 || N1 != 1 || N2 != 1);
      //
      ///////////////////////////////////////////////////////////////////////

      //double w0 = ( 0.625 - sqr( ( 0.375 + 0.25*cos(2.0*M_PI/(N0+2)) ) ) )/(N0+2);
      //double w1 = ( 0.625 - sqr( ( 0.375 + 0.25*cos(2.0*M_PI/(N1+2)) ) ) )/(N1+2);
      //double w2 = ( 0.625 - sqr( ( 0.375 + 0.25*cos(2.0*M_PI/(N2+2)) ) ) )/(N2+2);
      //
      // warren way
      double w0 = 0.375/(N0+2);
      double w1 = 0.375/(N1+2);
      double w2 = 0.375/(N2+2);

      double oneMinusNw0 = 1.0 - ( N0 + 2) * w0;
      double oneMinusNw1 = 1.0 - ( N1 + 2) * w1;
      double oneMinusNw2 = 1.0 - ( N2 + 2) * w2;

      if(false){
	std::cout <<"w0 = " << w0 << std::endl
		  <<"w1 = " << w1 << std::endl
		  <<"w2 = " << w2 << std::endl;
	std::cout <<"oneMinusNw0 = " << oneMinusNw0 << std::endl
		  <<"oneMinusNw1 = " << oneMinusNw1 << std::endl
		  <<"oneMinusNw2 = " << oneMinusNw2 << std::endl;
      }
      const int n = static_cast<int>(_nodes);
      assert(n > 0);
 
      double oneEighth = 1.0/8.0;
      double threeEighths = 3.0/8.0;  

      //
      //  first row
      S(0, 0) = oneEighth;
      S(0, 1) = threeEighths;
      S(0, 3) = threeEighths;
      S(0, 4) = oneEighth;
      // 
      // second row
      S(1, 0) = threeEighths;
      S(1, 1) = oneEighth;
      S(1, 3) = threeEighths;
      S(1, n-1) = oneEighth;
      //
      // third row
      S(2, 0) = w1;
      S(2, 1) = oneMinusNw1;
      S(2, blitz::Range(2,3+N1-1) ) = w1;
      //
      // fourth row
      S(3, 0) = threeEighths;  
      S(3, 1) = threeEighths;  
      S(3, 2) = oneEighth;  
      S(3, 3) = oneEighth;  
      //
      // fifth row
      S(4, 0) = oneMinusNw0;
      S(4, 1) = w0;
      S(4, 2) = w0;
      S(4, 3) = w0;
      S(4, blitz::Range(n-N0+1, n-1) ) = w0;
      //
      // sixth row
      S(5, 1) = threeEighths;
      S(5, 2) = oneEighth;
      S(5, 3+N1-2) = oneEighth;
      S(5, 3+N1-1) = threeEighths;
      //
      // seventh row
      S(6, 0) = oneEighth;
      S(6, 1) = threeEighths;
      S(6, 2) = threeEighths;
      S(6, 3+N1-1) = oneEighth;
      //
      // eighth row
      S(7, 0) = threeEighths;
      S(7, 1) = oneEighth;
      S(7, 2) = threeEighths;
      S(7, 3+N1+N2-2) = oneEighth;
      //
      // ninth row
      S(8, 0) = threeEighths;
      S(8, 2) = oneEighth;
      S(8, 3+N1+N2-2) = threeEighths;
      if ( 3+N1+N2-1 == N1 + N2 + N0  )
	S(8,3) = 1.0/8.0;
      else if ( 3+N1+N2-1 < N1 + N2 + N0  )		
	S(8, 3+N1+N2-1) = oneEighth;
      else
	std::cout << "never thought about this case ..." << std::endl;
      //
      // tenth row
      S(9, 1) = oneEighth;
      S(9, 2) = threeEighths;
      S(9, 3+N1-1) = threeEighths;
      S(9, 3+N1) = oneEighth;
      //
      // eleventh row
      S(10, 0) = w2;
      S(10, 1) = w2;
      S(10, 2) = oneMinusNw2;
      S(10, blitz::Range( 3+N1-1, 3+N1+N2-2 ) ) = w2;
      //
      // twelfth row
      S(11, 0) = oneEighth;
      S(11, 2) = threeEighths;
      S(11, 3+N1+N2-3) = oneEighth;
      S(11, 3+N1+N2-2) = threeEighths;
      //
    }
    else{
      S(3 ,0 ) = 1.0;
      S(6 ,1 ) = 1.0;
      S(7 ,2 ) = 1.0;
      S(2 ,3 ) = 1.0;
      S(5 ,4 ) = 1.0;
      S(9 ,5 ) = 1.0;
      S(10 ,6) = 1.0;
      S(11 ,7) = 1.0;
      S(8 ,8 ) = 1.0;
      S(4 ,9 ) = 1.0;
      S(1, 10) = 1.0;
      S(0 ,11) = 1.0;						
    }

    const bool Output_Flag = false;
		
    if( Output_Flag ) {
//       std::cout.precision(4);
//       std::cout.setf(std::ios_base::scientific,std::ios_base::floatfield);
      std::cout << "***** Begin Subdivision Matrix *****"<< std::endl;
      //	      << S << std::endl
      for(int i=0; i<S.rows(); i++) {
	std::cout << "    ";
	for(int j=0; j<S.cols(); j++) {
	  std::cout << S(i,j) << "    "; 
	}
	std::cout << "}" << std::endl;
      }
      std::cout << "*****  End Subdivision Matrix  *****"<< std::endl;
    }	
  }



  void LoopShellShape::_convertParaCoords()
  {
    _coords(0) = 1.0 - 2.0 * _coords(0);
    _coords(1) = 1.0 - 2.0 * _coords(1);	
    //! now, the parametric coords are related to the sub-patch
  }


  void LoopShellShape::checkShapeFunctions()
  {
    unsigned Min = 4, Max = 12;
    for(unsigned i=Min; i<=Max; i++) {
      for(unsigned j=Min; j<=Max; j++) {
	for(unsigned k=Min; k<=Max; k++) {
	  tvmet::Vector<unsigned,3> V;
	  V = i, j, k;
	  unsigned N = V(0) + V(1) + V(2) - 6;
	  CoordinateArray s;
	  s = 1.0/3.0, 1.0/3.0;
	  srand(time(0));
	  s(0) += 0.001*(double)(rand())/RAND_MAX; 
	  s(1) += 0.001*(double)(rand())/RAND_MAX; 
	  LoopShellShape shape(N, V);
	  std::cout <<"N="<<N<<std::endl
		    <<"V="<<V[0]<<" "<<V[1]<<" "<<V[2]<<std::endl;
	  shape.checkShapeNumerically(N, V);
	}
      }
    }
    std::cout << std::endl;
    return;
  }

  void LoopShellShape::checkShapeNumerically(int n, tvmet::Vector<unsigned,3> V)
  {
    //
    //  check sum
    //
    const bool Output_Flag = true;
    //
    ////////////////////////////////////////////////////////////////////////////////////////
    //
    // debug info
    //
    // output shape function of regular patch
    double sum=0.0;
    //
    // output shape functions of irregular patch
    sum=0.0;
    for(int  i=0; i<_functions.size(); i++) {
      //std::cout << "n["<<i<<"] =\t" << _functions(i) << std::endl;
      sum += _functions(i);
    }
    if(Output_Flag){
      std::cout << std::endl;
      std::cout << "*************************** New shape Functions ********************************" << std::endl;
      std::cout << "sum = "<< sum << std::endl;
      std::cout << " *******************************  END   **************************************" << std::endl;
    }
    assert( fabs(sum-1.0) < 1.0e-10 );


    //
    ////////////////////////////////////////////////////////////////////////////////////////
    //
    // output for regular patch
    double sum0=0.0, sum1=0.0;
    //
    // output for irregular patch	
    sum0 = sum1 = 0.0;
    for(int i=0; i<_derivatives.rows(); i++) {
      //std::cout << "i="<<i<<"\t" << _derivatives(i,0)<< "\t" <<_derivatives(i, 1) 
      //	  <<std::endl;
      sum0 += _derivatives(i, 0); sum1 += _derivatives(i, 1);
    }
    if(Output_Flag){
      std::cout << std::endl;
      std::cout << "***********  The new first derivatives of  shape function  ******" 
		<< std::endl;
      std::cout << "sum0 = "<<sum0<<std::endl
		<< "sum1 = "<<sum1<<std::endl;
      std::cout << "***********  sum0 and sum1 should equal to 0.0  **************"  
		<< std::endl;
      std::cout << " *************************  END  *****************************" 
		<< std::endl;
    }
    assert( fabs(sum0) < 1.0e-10 );
    assert( fabs(sum1) < 1.0e-10 );
      
    //
    ////////////////////////////////////////////////////////////////////////////////////////
    //
    // output for regular patch		
    double sum01=0.0, sum10=0.0;
    //
    // output for irregular patch		
    sum0 = sum1 = sum01 = sum10 = 0.0;
    for(int i=0; i<_secondDerivatives.rows(); i++) {
      //       std::cout << "i="<<i
      // 		<< "\t" << _secondDerivatives(i,0,0)
      // 		<< "\t" << _secondDerivatives(i,1,1)
      // 		<< "\t" << _secondDerivatives(i,0,1)
      // 		<< "\t" << _secondDerivatives(i,1,0) <<std::endl;
      sum0  += _secondDerivatives(i,0,0); 
      sum1  += _secondDerivatives(i,1,1); 
      sum01 += _secondDerivatives(i,0,1);
      sum10 += _secondDerivatives(i,1,0);
    }
    if(Output_Flag){
      std::cout << std::endl;
      std::cout << "*******  The new second derivatives of  shape function   ******" 
		<< std::endl;
      std::cout << " sum0  = " <<sum0<<std::endl
		<< " sum1  = " <<sum1<<std::endl
		<< " sum01 = " <<sum01<<std::endl
		<< " sum10 = " <<sum10<<std::endl;
      std::cout << "********** sums should be equal to 0.0  ***********" << std::endl;

      std::cout << "************************  END  ****************************" 
		<< std::endl;
    }
    assert( fabs(sum0) < 1.0e-10 );
    assert( fabs(sum1) < 1.0e-10 );
    assert( fabs(sum01) < 1.0e-10 );
    assert( fabs(sum10) < 1.0e-10 );
    //
    // end of checking sum
    // 

    //
    //  check consistency between the numerical and analytical results
    //
    const double eps = 1.0e-8;
    CoordinateArray para;
    FunctionArray fns( _functions.shape() );
    DerivativeArray derivs( _derivatives.shape() );
	
    for (int i = 0; i < 2; i ++){
      // set to current parametric coords
      para = _coords;
      // disturb +
      para(i) = _coords(i) + eps;
      LoopShellShape shapePos(n, V, para);
      fns = shapePos.functions();
      derivs = shapePos.derivatives();
      // disturb -
      para(i) = _coords(i) - eps;
      LoopShellShape shapeNeg(n, V, para);
      fns -= shapeNeg.functions();
      derivs -= shapeNeg.derivatives();

      //
      // compute average
      //
      // Deriv w.r.t ith component
      fns /= 2.0*eps;
      //
      // 2nd Deriv w.r.t ith component
      derivs /= 2.0*eps;

      //
      // compare with analytical results
      for(int j = 0; j < fns.size(); j++){
	const double fstD  = fns(j) - derivatives()(j,i);
	const double sndD0 = derivs(j, 0) - secondDerivatives()(j,i,0);
	const double sndD1 = derivs(j, 1) - secondDerivatives()(j,i,1);
	if ( fstD >= 1.0e-6 || sndD0 >= 1.0e-6 || sndD1 >= 1.0e-6 ){
	  std::cout << "Comparision between the Numerical and analytical results did NOT pass!" << std::endl;
	  std::cout << j << " th Node "
		    << i << " th direction"
		    << std::endl;
	  std::cout << fns(j) << ", " << derivatives()(j,i) 
		    << std::endl;
	  std::cout << derivs(j, 0) << ", " << secondDerivatives()(j,i,0) 
		    << std::endl;
	  std::cout << derivs(j, 1) << ", " << secondDerivatives()(j,i,1) 
		    << std::endl;		
	  exit(0);
	}
						
      }
    }
    std::cout << "Comparision between the Numerical and analytical results passed!" << std::endl;
  }
}

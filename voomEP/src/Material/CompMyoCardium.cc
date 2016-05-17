// -*- C++ -*-
//----------------------------------------------------------------------
//
//              Luigi Perotti & Shankarjee Krishnamoorthi
//                University of California Los Angeles
//                 (C) 2004-2007 All Rights Reserved
//
//----------------------------------------------------------------------
//

#include <iostream>
#include "VoomMath.h"
#include "CompMyoCardium.h"

//using namespace std;
namespace voom {

  // Construction/Destruction
  CompMyoCardium::CompMyoCardium(const CompMyoCardium &Input)
  {
    if (this == &Input) return;

    _F = Input._F;

    _rho = Input._rho;
   
    _lambda =  Input._lambda;
    _a   = Input._a;
    _b   = Input._b;
    _af  = Input._af;
    _bf  = Input._bf;
    _as  = Input._as;
    _bs  = Input._bs;
    _afs = Input._afs;
    _bfs = Input._bfs;
  }

	
  void CompMyoCardium::_init(double rho, double lambda, double a, double b, 
			     double af, double bf,  double as, double bs,  
			     double afs, double bfs)
  {
    _rho = rho;

    _lambda  = lambda;
    _a   = a;
    _b   = b;
    _af  = af;
    _bf  = bf;
    _as  = as;
    _bs  = bs;
    _afs = afs;
    _bfs = bfs;

    // deformation gradient tensor
    _F =
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0;
    //
    // strain energy density
    _W = 0.0;
    //
    // first Piola stress tensor
    _P = 
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0;
  }

  // Accessors/mutators

  double CompMyoCardium::massDensity()
  {
    return _rho;
  }

  // Operators


  // General methods
        
  void CompMyoCardium::updateState(bool fl0, bool fl1, bool fl2,
				   const Vector3D f, const Vector3D s)
  {

    double I1 = 0.0, I4f = 0.0, I4s = 0.0, I8fs = 0.0;
    Tensor3D C(0.0), Uff(0.0), Uss(0.0), Vfs(0.0), Delta(0.0), FUff(0.0), 
      FUss(0.0), FVfs(0.0), invF(0.);
    C = (tvmet::trans(_F))*_F;
    double jac = determinant(_F);
    invert( _F, invF );
    _I1 = tvmet::trace( C );
    _I2 = 0.;
    _I3 = jac*jac;

    if (fl0 || fl1 || fl2)
    {
      // Compute invariants
      I1 = tvmet::trace(C);
      _I4f = I4f =  tvmet::dot(f,(C*f));
      _I4s = I4s =  tvmet::dot(s,(C*s));
      I8fs = tvmet::dot(f,(C*s));
    }
    
    // Energy density
    if (fl0) {
      _W = 0.5*( (_a/_b)*exp(_b*(I1-3.0)) +
		 (_af/_bf)*(exp(_bf*pow(I4f-1.0,2.0))-1.0) +
		 (_as/_bs)*(exp(_bs*pow(I4s-1.0,2.0))-1.0) +
		 (_afs/_bfs)*(exp(_bfs*pow(I8fs,2.0))-1.0) +
		 _lambda*pow(log(jac),2.0) - 2.0*_a*log(jac) );
    }
      
    if (fl1 || fl2)
    {
      for (unsigned int i = 0; i < 3; i++) {
	for (unsigned int j = 0; j < 3; j++) {
	  Uff(i,j) = f(i)*f(j); 
	  Uss(i,j) = s(i)*s(j); 
	  Vfs(i,j) = f(i)*s(j) + f(j)*s(i); 
	}
      }
    }

    // First piola Kirchhoff stress tensor
    if (fl1) {
      _P =  _a*exp(_b*(I1-3.0))*_F + 
	2.0*_F*(_af*(I4f-1.0)*exp(_bf*pow(I4f-1.0,2.0))*Uff +
		_as*(I4s-1.0)*exp(_bs*pow(I4s-1.0,2.0))*Uss ) +
	_afs*I8fs*exp(_bfs*pow(I8fs,2.0))*_F*Vfs + 
	(_lambda*log(jac) - _a)*(tvmet::trans(invF));
      _S = _P*_F/jac;
      // std::cout << _P << "\n";
    }

    // Tangent moduli tensor
    if (fl2) {
      Delta(0,0) = 1.0; Delta(1,1) = 1.0; Delta(2,2) = 1.0;
      FUff = _F*Uff;
      FUss = _F*Uss;
      FVfs = _F*Vfs;
      for(int i = 0; i < 3; i++)
	for(int J = 0; J < 3; J++)
	  for(int k = 0; k < 3; k++)
	    for(int L = 0; L < 3; L++)
	      _C(i,J,k,L) = _a*exp(_b*(I1-3.0))*(2.0*_b*_F(i,J)*_F(k,L) + 
						 Delta(i,k)*Delta(J,L)) +
		2.0*_af*exp(_bf*pow(I4f-1.0,2.0))*(2.*FUff(i,J)*FUff(k,L)*(1.0+2.0*_bf*pow(I4f-1.0,2.0)) + (I4f-1.0)*Delta(i,k)*Uff(L,J) ) + 
		2.0*_as*exp(_bs*pow(I4s-1.0,2.0))*(2.*FUss(i,J)*FUss(k,L)*(1.0+2.0*_bs*pow(I4s-1.0,2.0)) + (I4s-1.0)*Delta(i,k)*Uss(L,J) ) +
		_afs*exp(_bfs*pow(I8fs,2.0))*(FVfs(i,J)*FVfs(k,L)*(1.0+2.0*_bfs*pow(I8fs,2.0)) + I8fs*Delta(i,k)*Vfs(L,J) ) +
		_lambda*invF(J,i)*invF(L,k) -
		(_lambda*log(jac) - _a)*invF(J,k)*invF(L,i);
    }

    return;
  }


  void CompMyoCardium::ConsistencyTest()
  {
    std::cout << "checking consistency of 1st Piola stress tensor" << std::endl;
    updateState(false, true, false);
    Tensor3D PAna, PNum;  PAna = _P;          // current 1st Piola stress
    const double eps = 1.0e-8 * max(_F);      // perturbation
    for(int i = 0; i < PAna.rows(); i ++){
      for( int j = 0; j < PAna.cols(); j++){
	_F(i,j) += eps;
	updateState(true, false, false);
	double W = _W;
	_F(i,j) -= 2*eps;
	updateState(true, false, false);
	W -= _W;
	_F(i,j) += eps;       // restore value
	PNum(i,j) = W/2.0/eps;
      }
    }
    std::cout << "Analytical value of 1st Piola stress:" << std::endl;
    std::cout << PAna << std::endl;
    std::cout << "Numerical value of 1st Piola stress:" << std::endl;
    std::cout << PNum << std::endl;
    std::cout << "Difference between 1st Piola stresses:" << std::endl;
    PNum -= PAna;
    std::cout << PNum << std::endl;
		
    std::cout << "\n\n\n\n" << std::endl;
    return;
  }

  // Get material invariants
  Real CompMyoCardium::getInvariant(std::string request) const {
    if ( request == "I1") return _I1;
    else if (request == "I2") return _I2;
    else if (request == "I3") return _I3;
    else if (request == "I4s") return _I4s;
    else if (request == "I4f") return _I4f;	     
    else return 0.;    
  }

} // namespace voom

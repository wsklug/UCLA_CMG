// -*- C++ -*-
//----------------------------------------------------------------------
//
//                 William S. Klug, Melissa M. Gibbons
//                University of California Los Angeles
//                 (C) 2004-2007 All Rights Reserved
//
//----------------------------------------------------------------------
//

#include <iostream>
#include "VoomMath.h"
#include "NeoHookean.h"

//using namespace std;
namespace voom {

  // Construction/Destruction
  NeoHookean::NeoHookean(const NeoHookean &Input)
  {
    if (this == &Input) return;

    _F = Input._F;

    _rho = Input._rho;
    _C0 = Input._D0;
    _D0 = Input._D0;
  }

	
  void NeoHookean::_init(double rho, double C0, double D0)
  {
    _rho = rho;
    _C0 = C0;
    _D0 = D0;

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

  double NeoHookean::massDensity()
  {
    return _rho;
  }

  double NeoHookean::longitudinalWaveSpeed()
  {
    return 0.;
  }

  // Operators


  // General methods
        
  void NeoHookean::updateState(bool fl0, bool fl1, bool fl2,
				   const Vector3D f, const Vector3D s)
  {
    bool debug=false;
    // invF:  _F^-1
    Tensor3D invF;
    double jac = determinant(_F);
    invert( _F, invF );
    _I1 = tvmet::trace( tvmet::trans(_F)*_F );
    _I3 = jac*jac;

    // 3-D Elasticity
    // first PK stress = (lame*ln(J) - shear)*inverse(_F) + shear*_F
    if (fl1 || fl0) {
      Real c1 = _C0 * pow( jac, -2./3.);
      Real c2 = 2./_D0*(jac-1.)*jac - 2*_C0*_I1*pow(jac, -2./3.);
      _P = 2.*c1*_F + c2* tvmet::trans( invF );
      _S = _F * _P/jac;
    }

    if (fl0) {
      // C:  right Cauchy-green strain tensor
      _W = 0.0;
      _W = _C0*(_I1*pow(jac, -2./3.) - 3.0) + 1./_D0*(jac - 1)*(jac - 1.);
    }
  
    if (fl2) {
      // Computing the Lagrangian Moduli
      _C = 0.; // Initialize the Moduli
      Tensor3D delta(0.);
      Real J23 = pow(jac, -2./3.);
      for(int i = 0; i < 3; i++) delta(i,i) = 1.;

      for(int i = 0; i < 3; i++)
	for(int J = 0; J < 3; J++)
	  for(int k = 0; k < 3; k++)
	    for(int L = 0; L < 3; L++)
	      _C(i,J,k,L) = 2*_C0*J23*( -2./3.*invF(L,k)*_F(i,J) +     
					delta(i,k)*delta(J,L) - 
					2./3.*_F(k,L)*invF(J,i) +
					2./9.*_I1*invF(L,k)*invF(J,i) +	
					_I1/3.*invF(J,k)*invF(L,i) ) +	
		2*_D0*( (2*jac-1)*jac*invF(L,k)*invF(J,i) -	
			(jac*jac-jac)*invF(J,k)*invF(L,i) );
    }


    return;
  }


  void NeoHookean::ConsistencyTest()
  {
    std::cout << "checking consistency of 1st Piola stress tensor" << std::endl;
    updateState(false, true, false);
    Tensor3D PAna, PNum;  PAna = _P;                  // current 1st Piola stress
    const double eps = 1.0e-8 * max(_F);   // perturbation
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
    std::cout << "Numberical value of 1st Piola stress:" << std::endl;
    std::cout << PNum << std::endl;
    std::cout << "Difference between 1st Piola stresses:" << std::endl;
    PNum -= PAna; PNum /= max(PAna);
    std::cout << PNum << std::endl;
		
    std::cout << "\n\n\n\n" << std::endl;
    return;
  }
  
  // Get material invariants
  Real NeoHookean::getInvariant(std::string request) const{
    if ( request == "I1") return _I1;
    else if (request == "I2") return _I2;
    else if (request == "I3") return _I3;
    else return 0.;
    
  }

} // namespace voom

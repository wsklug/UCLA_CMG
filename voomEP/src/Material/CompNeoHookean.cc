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
#include "CompNeoHookean.h"

//using namespace std;
namespace voom {

  // Construction/Destruction
  CompNeoHookean::CompNeoHookean(const CompNeoHookean &Input)
  {
    if (this == &Input) return;

    _F = Input._F;

    _rho = Input._rho;
    _E = Input._E;
    _nu = Input._nu;
  }

	
  void CompNeoHookean::_init(double rho, double E, double nu)
  {
    _rho = rho;
    _E = E;
    _nu = nu;

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

  double CompNeoHookean::massDensity()
  {
    return _rho;
  }

  double CompNeoHookean::longitudinalWaveSpeed()
  {
    // lame: lame constant
    // shear: shear modulus
    const double lame = _nu*_E/((1.0+_nu)*(1.0-2.0*_nu)); 
    const double shear = 0.5*_E/(1.0+_nu);
    //                       ______________________
    //                      /
    //                     /  lame + 2.0 * shear
    //  P wave speed  ==  /  ____________________
    //                   /
    //                 \/            rho
    //
    return sqrt((lame + 2.0 * shear)/_rho);
  }

  // Operators


  // General methods
        
  void CompNeoHookean::updateState(bool fl0, bool fl1, bool fl2,
				   const Vector3D f, const Vector3D s)
  {
    bool debug=false;
    
    // lame: lame constant
    // shear: shear modulus
    const double lame = _nu*_E/((1.0+_nu)*(1.0-2.0*_nu)); 
    const double shear = 0.5*_E/(1.0+_nu);
        
    // invF:  _F^-1
    Tensor3D invF;
    double jac = determinant(_F);
    invert( _F, invF );
    _I1 = tvmet::trace( tvmet::trans(_F)*_F );
    _I3 = jac*jac;

    // 3-D Elasticity
    // first PK stress = (lame*ln(J) - shear)*inverse(_F) + shear*_F
    if (fl1 || fl0) {
      _P = (lame*log(jac) - shear)*(tvmet::trans(invF)) + shear*_F;
      _S = _F * _P/jac;
    }

    if (fl0) {
      // C:  right Cauchy-green strain tensor
      Tensor3D C;
      C = (tvmet::trans(_F))*_F;

      _W = 0.0;
      _W += 0.5*lame*(log(jac))*(log(jac)) - shear*log(jac) + 0.5*shear*(tvmet::trace(C)-3.0);
    }
  
    if (fl2) {
      // Computing the Lagrangian Moduli
      _C = 0.; // Initialize the Moduli
      Tensor3D delta(0.);
      for(int i = 0; i < 3; i++) delta(i,i) = 1.;

      for(int i = 0; i < 3; i++)
	for(int J = 0; J < 3; J++)
	  for(int k = 0; k < 3; k++)
	    for(int L = 0; L < 3; L++)
	      _C(i,J,k,L)=lame*invF(J,i)*invF(L,k) + 
		shear*delta(i,k)*delta(J,L) -
		(lame*log(jac) - shear)*invF(J,k)*invF(L,i);
    }


    return;
  }


  void CompNeoHookean::ConsistencyTest()
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
  Real CompNeoHookean::getInvariant(std::string request) const{
    if ( request == "I1") return _I1;
    else if (request == "I2") return _I2;
    else if (request == "I3") return _I3;
    else return 0.;
    
  }

} // namespace voom

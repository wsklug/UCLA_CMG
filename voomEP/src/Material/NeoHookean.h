// -*- C++ -*-
//----------------------------------------------------------------------
//
//                 William S. Klug, Melissa M. Gibbons
//                University of California Los Angeles
//                 (C) 2004-2007 All Rights Reserved
//
//----------------------------------------------------------------------
//
/*! 
  \file NeoHookean.h

  \brief Interface for a compressible Neo-Hookean finite deformation
  hyperelasticity model.

*/

#ifndef _NEOHOOKEAN_H_
#define _NEOHOOKEAN_H_

#include "Material.h"
#include "VoomMath.h"

namespace voom {
  class NeoHookean : public Material
  {
    
  public:
    
    // Constructors/destructors:
    //! Default constructor
    NeoHookean() { 
      _init(0.0,0.0,0.0); 
      _C.resize(3,3,3,3);    
    }
    
    //! Construct material from density, C0 and D0
    NeoHookean(Real data[]) {
      _init(data[0], data[1], data[2]);
      _C.resize(3,3,3,3);    
    }
    //! Destructor
    virtual ~NeoHookean() {}
    //! Copy constructor
    NeoHookean(const NeoHookean &);
    
    // Accessors/mutators:
    //! Returns density
    inline double massDensity();
    //! Calculates longitudinal wave speed
    inline double longitudinalWaveSpeed();
    
    //! Get an invariant from the model
    Real getInvariant(std::string request) const;
    
    // General methods:
    /*! 
      Based on new deformation gradient tensor, F, calculates state of 
      material (strain energy density, first Piola-Kirchhoff stress tensor, 
      Lagrangian Moduli)
    */
    void updateState(bool f0, bool f1, bool f2, 
		     const Vector3D f=Vector3D(0.,0.,1.), 
		     const Vector3D s=Vector3D(1.,0.,0.));
    
    // Tests:
    //! Consistency test
    void ConsistencyTest();
    static void MFITest();
    static void IsotropyTest();
    
    
  private:
    
    // Members:
    //! Density
    double _rho;
    //! Young's Modulus
    double _C0;
    //! Poisson's ratio
    double _D0;

    //! Initializes material properties and zeros out deformation/stress tensors
    void _init(double rho, double C0, double D0);
  };
  
}
#endif // _NEOHOOKEAN_H_

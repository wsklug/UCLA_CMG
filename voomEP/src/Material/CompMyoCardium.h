// -*- C++ -*-
//----------------------------------------------------------------------
//
//              Luigi Perotti & Shankarjee Krishnamoorthi
//                University of California Los Angeles
//                 (C) 2004-2007 All Rights Reserved
//
//----------------------------------------------------------------------
//
/*! 
  \file MyoCardium.h

  \brief Interface for passive myocardium material model according to:
   Holzapfel G.A. amd Ogden R.W., 2009:"Constitutive modelling of passive 
   myocardium: a structurally based framework for material characterization". 
   Phylosophical transactions of the royal society A, 367, 3445-3475.
   +
   compressible part
*/

#ifndef _COMPMYOCARDIUM_H_
#define _COMPMYOCARDIUM_H_

#include "Material.h"
#include "VoomMath.h"

namespace voom {
  
  class CompMyoCardium : public Material
  {
  public:
    // Constructors/destructors:
    //! Default constructor
    CompMyoCardium() { 
      _init(1.053e-3, 1.0e3, 59.0, 8.023, 18472.0, 16.026, 2481.0, 11.120, 
	    216.0, 11.436); 
      // Standard initialization using values reported in Holzapfel and 
      // Ogden 2009
      // (rho[g/mm^3], lambda, a[Pa], b, af[Pa], bf, as[Pa], bs, afs[Pa], bfs)
      _C.resize(3,3,3,3);    
    }
    
    //! Construct material from assigned material parameters
    CompMyoCardium(Real data[]) {
      _init(data[0], data[1], data[2], data[3], data[4], data[5], data[6],
	    data[7], data[8], data[9]);
      // std::cout << "f: " << f << "\n";
      // std::cout << "s: " << s << "\n";
      _C.resize(3,3,3,3);    
    }
    //! Destructor
    virtual ~CompMyoCardium() {}
    //! Copy constructor
    CompMyoCardium(const CompMyoCardium &);
    
    // Accessors/mutators:
    //! Returns density
    inline double massDensity();
    
    // Operators
    
    // General methods:
    /*! 
      Based on new deformation gradient tensor, F, calculates state of 
      material (strain energy density, first Piola-Kirchhoff stress tensor, 
      Lagrangian Moduli)
    */
    void updateState(bool f0, bool f1, bool f2,
		     const Vector3D f=Vector3D(0.,0.,1.),
		     const Vector3D s=Vector3D(1.,0.,0.) );
    
    // Tests:
    //! Consistency test
    void ConsistencyTest();
    static void MFITest();
    static void IsotropyTest();

    //! Get an invariant from the model
    Real getInvariant(std::string request) const;    
    
  private:
    
    // Members:
    //! Density
    double _rho;
   
    double _lambda;
    double _a;
    double _b;
    double _af;
    double _bf;
    double _as;
    double _bs;
    double _afs;
    double _bfs;
    double _I4f, _I4s;

    //! Initializes material properties and zeros out deformation/stress tensors
    void _init(double rho, double lambda, double a, double b, double af, 
	       double bf,  double as, double bs,  double afs, double bfs);
  };
  
}
#endif // _COMPMYOCARDIUM_H_

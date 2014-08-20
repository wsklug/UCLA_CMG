//-*-C++-*-
/*! \brief
  Tusscher et al ionic model for human ventricular tissue. Developed based on 
  the paper "Alternans and spiral breakup in a human ventricular tissue model"
  by K. H. W. J. ten Tusscher and A. V. Panfilov, Am J Physiol Heart Circ 
  Physiol 291:H1088-H1100, 2006. First published 24 March 2006
  - Implemented methods described in "Efficient Numerical Technique for the 
  solution of the monodomain and Bidomain Equations" by Whiteleley. We can
  take larger time steps for Ionic Solve
*/
#ifndef _Tusscher_h_
#define _Tusscher_h_

#include "IonicModel.h"

namespace voom{
  class Tusscher: public IonicModel {
  private:
    //! Constants for the model
    Real *_Constants;
    
    //! State Variables
    std::vector<Real>  _State;

    //! Rate Variables
    Real  _Rates[19];

    //! Algebraic Variables
    Real _Algebraic[70];

    //! Adaptive Update of State Variables
    void UpdateStateVariables(const Real dt);

  public:
    //! Constructor
    Tusscher(Real *constants){
      _Constants = constants;
      _State.resize( 19 );
      /* As in Tusscher Paper
      // State Variables
      _State[0] = -85.423;    _State[1] = 138.52;    _State[2] = 10.132;
      _State[3] = 0.000153;   _State[4] = 0.0165;    _State[5] = 0.473;
      _State[6] = 0.0174;     _State[7] = 0.00165;   _State[8] = 0.749;
      _State[9] = 0.6788;     _State[10] = 0.00042;  _State[11] = 3.288e-5;
      _State[12] = 0.7026;    _State[13] = 0.9526;   _State[14] = 0.9942;
      _State[15] = 0.999998;  _State[16] = 2.347e-8; _State[17] = 4.272;
      _State[18] = 0.8978;
      */

      // Modified State Variables as requested by KCL
      _State[0] = -85.423; // Resting Voltage
      _State[1] = 136.89;  // Intracellular_K
      _State[2] = 8.604;   // Intracellular_Ni
      _State[3] = 0.000126;// Intracellular Calcium
      _State[4] = 0.00621; // Xr1
      _State[5] = 0.4712;  // Xr2
      _State[6] = 0.0095;  // Xs
      _State[7] = 0.00172; // m
      _State[8] = 0.7444;  // h
      _State[9] = 0.7045;  // j
      _State[10] = 0.00036; // Subspace Calcium
      _State[11] = 3.373e-5; // L_Type d
      _State[12] = 0.7888;   // L_Type f
      _State[13] = 0.9755;   // L_type f2
      _State[14] = 0.9953;   // L_Type fCass
      _State[15] = 0.999998; // Transient outward current s gate
      _State[16] = 2.42e-8;  // Transient outward current r gate
      _State[17] = 3.64;    // Sarcoplasmic Reticulum Calcium
      _State[18] = 0.9073;   // R_prime
      
      //! Surface Area to Volume Ratio from the paper mentioned above
      _Xi        = 2000.;
    }

    //! Destructor
    ~Tusscher() {;}

    //! Compute Rates
    void ComputeRates();

    //! Compute Ionic Current
    Real Compute_Ion(Real Xi, bool userXi, Real C_m, Real dt, Real volt,
		     Real istim);

    //! Get gamma
    Real getGamma();

    // get internal variables
    const std::vector<Real>& getInternalParameters(int& nData) const { 
      nData = 19; return _State; }

    // Set internal variables
    void setInternalParameters(const std::vector<Real>& data) {
      for(int i = 0; i < 19; i++) _State[i] = data[i];
    }      
  };
}
#endif



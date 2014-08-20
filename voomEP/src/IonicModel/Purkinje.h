//-*-C++-*-
/*! \brief
  Ionic model for Rabbit Purkinje structure. Developed based on the paper by
  Corrias et al titled "Ionic mechanisms of electrophysiological properties 
  and repolarization abnormalities in rabbit Purkinje fibers", American
  Journal of HeartCirculation Physiology, Feb 2011
  - Added Faster computation as suggested by Whiteley. We do Explicit 
  update on Linear ODE and Newton-Raphson on Non-Linear ODE's.

*/
#ifndef _Purkinje_h_
#define _Purkinje_h_

#include "IonicModel.h"

namespace voom {
  class Purkinje: public IonicModel {
  private:
    //! Constants for the model
    Real  *_Constants;

    //! State variables
    std::vector<Real>  _State;

    //! Rate Variables
    Real  _Rates[22];

    //! Algebraic variable array
    Real  _Algebraic[69];

    //! Surface Area of the Cell
    Real  _SurfaceArea;

    //! Forward Euler Adaptive Update of State Variables
    void UpdateStateVariables(const Real dt);

    //!
  public:
    //! Constructor
    Purkinje(Real *constants){
      _Constants = constants;
      _State.resize( 22 );
      // Initialize States
      _State[0] = -88.34;      _State[1] = 0.00001;      _State[2] = 0.000032;
      _State[3] = 0.17;        _State[4] = 6.7;          _State[5] = 140;    
      _State[6] = 0.001337;    _State[7] = 0.01;         _State[8] = 0.000003;
      _State[9] = 0.1;         _State[10] = 0.7;         _State[11] = 0.0;
      _State[12] = 0.7;        _State[13] = 0.000007;    _State[14] = 0.978861;
      _State[15] = 0.000012;   _State[16] = 0.864489;    _State[17] = 0.25;
      _State[18] = 1.0;        _State[19] = 0.0;         _State[20] = 0.011099;
      _State[21] = 0.011099;
      // Surface Area to Volume ratio from the paper cm^-1
      _Xi        =  2500.;
      // Surface Area of the Cell from the paper 4612 um^2. Converted to cm^2
      _SurfaceArea = 4.612E-5 ;// cm^2
    }

    //! Destructor
    ~Purkinje() {;}
    
    //! Compute Rates
    void ComputeRates();    

    //! Compute Variables
    void ComputeVariables();

    //! Compute Ionic Current
    Real Compute_Ion(Real Xi, bool userXi, Real C_m, Real dt, Real volt,
		     Real istim);

    //! Get gamma
    Real getGamma();

    // get internal variables
    const std::vector<Real>& getInternalParameters(int &nData) const { 
      nData = 22; return _State; }

    // Return voltage
    Real getVoltage() {
      return _State[0];
    }

    // Set internal variables
    void setInternalParameters(const std::vector<Real>& data) {
      for(int i = 0; i < 22; i++) _State[i] = data[i];
    }
  };  
}

#endif


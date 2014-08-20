//-*-C++-*-
#ifndef __Mahajan_fail_h__
#define __Mahajan_fail_h__

/*! \file ucla.c
 *  \brief Source code for the UCLA rabbit ventricle model.
 *
 *  This code is an implementation of the UCLA rabbit ventricular action
 *  potential model as described in Mahajan et al. (2008).
 *
 *  Implementation:
 *  Daisuke Sato 
 *  Yohannes Shiferaw  
 *  Enno de Lange
 *
 * (c) 2011, UCLA Dept. of Medicine (Cardiology)
 * Implemented in VOOM
 */

#include "IonicModel.h"

namespace voom {
  class Mahajan_fail: public IonicModel {
  private:
    //! Number of Variables
    int _nVar;
    
    //! Internal State Variables
    std::vector<Real> _State;

    //! Internal Rate Variables
    Real _Rates[26];

    //! Constants
    Real *_Constants;

    //! Compute I_na
    Real comp_ina(const Real dt);

    //! Compute I_kr
    Real comp_ikr(const Real dt);

    //! Compute I_ks
    Real comp_iks(const Real dt);

    //! Compute I_k1
    Real comp_ik1();

    //! Compute I_total
    Real comp_ito(const Real dt);

    //! Compute I_NaK
    Real comp_inak();

    //! Compute I_NaCa
    Real comp_inaca();

    //! Compute I_CalPo
    Real comp_icalpo(const Real dt);

    //! Compute I_up
    Real comp_iup();

    //! Compute I_leak
    Real comp_ileak();

    //! Compute Inst_Buffer
    Real comp_inst_buffer(Real cs);

    //! Compute rxa
    Real comp_rxa();

    //! Compute Q
    Real comp_Q();

    //! Copute dir
    Real comp_dir(Real po, Real Qr, Real rxa, Real dcj);

    //! Compute dcp
    Real comp_dcp(Real po, Real Qr, Real rxa);

    //! Compute Insca
    Real comp_insca();

    //! Forward Euler Update of State Variables
    void UpdateStateVariables(const Real dt, const Real i_stim);
  public:
    //! Constructor. Pos defines Epi/Myo/Endo/Apex/Center/base 0 - 8 index
    Mahajan_fail(Real *Constants, int);

    //! Destructor
    ~Mahajan_fail() {;}
    
    /*!
      Euler Time stepping: Compute the rhs of the UCLA cell model for given
      dt. Stimulus current is in pA/pF. 
    */
    void ucla_rhsfun(const Real h, const Real istim);

    //! Get default Parameters
    void getDefaultParameters();

    /*! 
    Compute Ion Function. The stimulus current is converted to pA/pF inside
    the compute routine. The body class need not perform this
    */
    Real Compute_Ion(Real Xi, bool userXi, Real C_m, Real dt, Real Volt,
		     Real istim);

    //! get Gamma
    Real getGamma() {
	return 0.;
    }

    // Return voltage
    Real getVoltage() {
      return _State[0];
    }

    // get internal variables
    const std::vector<Real>& getInternalParameters(int& nData) const {
      nData = 26; return _State; }

    // Set internal variables
    void setInternalParameters(const std::vector<Real>& data) {
      for(int i = 0; i < 26; i++) _State[i] = data[i];
    }
  };
} 

#endif

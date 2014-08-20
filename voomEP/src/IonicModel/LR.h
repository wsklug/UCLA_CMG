//-*-C++-*-
// Luo Rudy ionic model class 
/*
Luo, C. H., & Rudy, Y. (1994). A dynamic model of the cardiac ventricular action potential. I. Simulations of ionic currents and concentration changes. Circulation research, 74(6), 1071-1096.
*/

#if !defined(__LR_h__)
#define __LR_h__
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "voom.h"
#include "LRGateTable.h"
#include "IonicModel.h"


using namespace std;
namespace voom {  
  class LuoRudy: public IonicModel{    
  public:    
    void reinitialize(Real y[]);  
    Real Compute_Ion(Real Xi, bool userXi, Real C_m, Real dt, Real volt,
		     Real istim);
    void make_parameter1(Real v,Real tab[]);
    void initialize();    
    void setGateTable(LuoRudyGateTable * t) { gate = t;}
    void setGna(Real val);
    void setGsi(Real val);     

    //! default constructor
    LuoRudy(bool useGate, LuoRudyGateTable* Gate){
      _useGate = useGate;
      gate = Gate;
      _in.resize(8);
      initialize();
    }

    // get internal parameters
    const std::vector<Real>& getInternalParameters(int& nData) const {
      nData = 8; return _in; 
    }

    // set internal parameters
    void setInternalParameters(const std::vector<Real>& data) {
      for(int i = 0; i < 8; i++) _in[i] = data[i];
    }

    std::vector<Real> _in;
    Real Gna,Gsi,c_K0,GK;
    Real c_Na0,c_Nai,c_Ki,RT_F;
    Real E_Na, E_K,E_K1,E_kp,G_K1;
    LuoRudyGateTable *gate; 

    //! Get gamma
    Real getGamma();

  private:
    //! Use LuoRudy Gate Table or Not
    bool _useGate; 
  };
  
}  // namespace voom

#endif

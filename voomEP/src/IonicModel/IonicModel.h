//-*-C++-*-
/*!\brief
  A base Ionic class. This will be used by Cardiac Body to typecast to the 
  derived class and this pointer will be provided to the Body class
*/

#ifndef _Ionic_Model_h_
#define _Ionic_Model_h_
#include "voom.h"
#include <vector>

// Calling LAPACK Ax=b Solver
extern "C" void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV,
		       double *B, int *LDB, int *INFO);
namespace voom {
  class IonicModel{
  protected:
    //! Surface Area to Volume Ratio Xi
    Real        _Xi;
  public:
    //! Constructor
    IonicModel(){;}

    //! Destructor
    ~IonicModel(){;}

    //! get SurfaceArea to Volume Ratio
    Real getSurfaceAreaToVolumeRatio() const {
	return _Xi;
    }

    //! Compute Ion virtual function. Return dV/dt to body class
    virtual Real Compute_Ion(Real Xi, bool userXi, Real C_m, Real dt, 
			     Real Volt, Real istim) = 0;

    /*! 
      Get gamma. Magnitude of active deformation. This depends on the 
      [Ca] concentration.
    */
    virtual Real getGamma() = 0;
    
    virtual Real getVoltage() { return 0.0; };

    //! Get all internal parameters
    virtual const std::vector<Real>& getInternalParameters(int& nData) const 
    = 0;

    //! Set all internal parameters
    virtual void setInternalParameters(const std::vector<Real>& data) = 0;
  };
}
#endif

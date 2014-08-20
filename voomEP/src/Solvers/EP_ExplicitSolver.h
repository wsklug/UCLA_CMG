//-*-C++-*-
#ifndef _EP_ExplicitSolver_h_
#define _EP_ExplicitSolver_h_
#include "EPSolver.h"

namespace voom{
  /*! \brief
    Forward Euler Solver for the diffusion part. We solve the PDE
    \f[
    \chi  C\frac{\partial V}{\partial t} = -\nabla \cdot (\boldsymbol \sigma \,
    \nabla V)
    \f]
    using Forward Euler Scheme as
    \f[
    \chi\mathbf{C} \frac{V^{n+1}-V^n}{\Delta t} = \boldsymbol \sigma  V^n
    \f]
    There is no Matrix solution involved here. We model the capacitance as
    lumped. Hence we get
    \f[
    \chi C V^{n+1} = \boldsymbol \sigma V^n \,\Delta t + \chi C V^n
    \f]
    Here \f$\sigma\f$ is the conductivity matrix of the entire system.
    In the Ionic Solve we have
    \f[
    \int_{\Omega}C_{ij}N_i N_j dV \dot{V}= \int_{\Omega} I_{ion} N_i N_j dV
    \f]
    We lump the capacitance matrix and hence we get
    \f[
    C_m \dot{V} \left( \int_{\Omega} N_i N_j dV \right) = I_{ion} \left( \int_{
    \Omega} N_i N_j dV \right)
    \f]
    This yields
    \f[
    C_m \dot{V} = I_{ion}
    \f]
  */
  class EP_ExplicitSolver : public EPSolver {
  private:
    //! Conductivity Matrix
    Epetra_FECrsMatrix *_D;

    //! Check if Diffusion matrix has been initialized
    bool                _check;

    /*! 
      Initialize some internal variables. In this case the terms in the 
      conductivity matrices are filled up.
    */
    void Initialize();
  public:
    //! Constructor
    EP_ExplicitSolver(Epetra_MpiComm *comm, const MeshContainer& mesh, 
		      const BodyContainer& cardiacBody,
		      const SmoothCoeff& coeff):
      EPSolver(comm, mesh, cardiacBody, coeff) {
      _D = new Epetra_FECrsMatrix(Copy, *_targetMap, 10);
      // Call the compute in body for diffusive stiffness calculations
      for(int i = 0 ;i < _myBody.size(); i++) 	  
	_myBody[i]->Compute();
      _check = true;
    }

    //! Destructor
    ~EP_ExplicitSolver() {delete _D;}


    //! Virtual function from base class overwritten.
    void Solve(const Real dtDiff = 0.05);
  };
}

#endif

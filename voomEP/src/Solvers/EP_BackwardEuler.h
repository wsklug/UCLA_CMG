//-*-C++-*-
#ifndef _EP_BackwardEuler_h_
#define _EP_BackwardEuler_h_
#include "EPSolver.h"

#include <Epetra_FECrsMatrix.h>
namespace voom{
  /*! \brief
    Backward Euler for the diffusion part. This is the basic implicit solver.
    \f[
    \chi  C\frac{\partial V}{\partial t} = -\nabla \cdot (\boldsymbol \sigma \,
    \nabla V)
    \f]
    Using backward Euler we have
    \f[
    \chi \mathbf{C} \frac{V^{n+1} - V^n}{\Delta t} = \boldsymbol \sigma V^{n+1}
    \f]
    Using lumped Capacitance and rewriting the equation we get
    \f[
    ( \chi C - \boldsymbol \sigma \Delta t ) V^{n+1} = \chi C V^n
    \f]
    We have a problem of the form \f$ Ax=b \f$. This is solved using the 
    Conjugate Gradient Solver in Trilinos. We use the Jacobi Preconditioner.
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

  class EP_BackwardEuler : public EPSolver {
  private:
    //! Stiffness Matrix
    Epetra_FECrsMatrix *_stiffness;

    //! Using AZTEC Solver
    Epetra_LinearProblem* _problem;
    AztecOO*              _solver;

    //! Check if Diffusion matrix has been built
    bool                  _check;

    //! Verbose Output
    bool                  _verbose;

    //! Initialize
    void Initialize();

  public:
    //! Constructor
    EP_BackwardEuler(Epetra_MpiComm *comm, const MeshContainer& mesh, 
		     const BodyContainer& cardiacBody,
		     const SmoothCoeff& coeff, 
		     const Real dT = 1.):
      EPSolver(comm, mesh, cardiacBody, coeff) {
      // Call the compute in body for diffusive stiffness calculations
      for(int i = 0 ;i < _myBody.size(); i++) 	  
	_myBody[i]->Compute();
      // Create the stiffness matrix
      _stiffness = new Epetra_FECrsMatrix(Copy, *_targetMap, 2);
      _dT = dT;
      _check = true;
      _verbose = false;
    }


    //! Destructor
    ~EP_BackwardEuler() {
      delete _stiffness; delete _problem; delete _solver;
    }
    
    //! Verbose version
    void setVerbose(bool cond = true) {
      _verbose = cond;
    }

    //! Defining pur virtual function in base class
    void Solve(const Real dtDiff = 0.1);
  };
}

#endif

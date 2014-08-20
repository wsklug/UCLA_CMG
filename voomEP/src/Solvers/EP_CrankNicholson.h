//-*-C++-*-
#ifndef _EP_CrankNicholson_h_
#define _EP_CrankNicholson_h_
#include "EPSolver.h"

namespace voom{
  /*! \brief
    CrankNicholson Implicit solver. Gives us temporal accuracy of 
    \f$O(\Delta t^2)\f$. 
    \f[
    \chi  C\frac{\partial V}{\partial t} = -\nabla \cdot (\boldsymbol \sigma \,
    \nabla V)
    \f]
    We discretize this as
    \f[
    \chi C \frac{V^{n+1}-V^n}{\Delta t} = \frac{1}{2} \left( \boldsymbol 
    \sigma V^{n+1} + \boldsymbol \sigma V^n \right)
    \f]
    Rearranging terms we get
    \f[
    \left ( \chi C - \boldsymbol \sigma \frac{\Delta t}{2} \right) V^{n+1} = 
    \left( \frac{ \boldsymbol \sigma \Delta t}{2} + \chi C \right) V^n
    \f]
    This is solved using Conjugate Gradient solver. We use operator 
    splitting to get use different time steps. In the Ionic Solve we have
    \f[
    \int_{\Omega}C_{ij}N_i N_j dV \dot{V} = \int_{\Omega} I_{ion} N_i N_j dV
    \f]
    To improve the solve time we use lumping preferentially. We lump 
    on the LHS and not on RHS. Let \f$\int_{\Omega}N_iN_j dV = M_{ij}\f$. We 
    get
    \f[
    C_L \dot{V} = M I_{ion}
    \f]
    where \f$C_L \f$ is the lumped capacitance matrix.
  */
  
  class EP_CrankNicholson : public EPSolver {
  private:
    //! Stiffness Matrix
    Epetra_FECrsMatrix *_stiffness;

    //! Diffusion Matrix
    Epetra_FECrsMatrix *_D;

    //! Capacitance Matrix
    Epetra_FECrsMatrix *_C, *_CL;

    //! ICI
    Epetra_Vector *_ICI;

    //! AMESOS problem variables
    Epetra_LinearProblem* _problem;
    AztecOO*              _solver;

    //! Ionic Problem Variables
    Epetra_LinearProblem* _ionicProblem;
    AztecOO*              _ionicSolver;
    
    //! Check if Diffusion matrix has been built
    bool                  _check;

    //! Initialize
    void Initialize();

  public:
    //! Constructor
    EP_CrankNicholson(Epetra_MpiComm *comm, const MeshContainer& mesh, 
		      const BodyContainer& cardiacBody,
		      const SmoothCoeff& coeff, 
		      const Real dT = 1.):
      EPSolver(comm, mesh, cardiacBody, coeff) {
      // Call the compute in body for diffusive stiffness calculations
      for(int i = 0 ;i < _myBody.size(); i++) 	  
	_myBody[i]->Compute();
      // Create the Diffusion  matrix. 10 here is just a guess
      _stiffness = new Epetra_FECrsMatrix(Copy, *_targetMap, 10);
      _D = new Epetra_FECrsMatrix(Copy, *_targetMap, 10);
      _C = new Epetra_FECrsMatrix(Copy, *_targetMap, 20);
      _CL = new Epetra_FECrsMatrix(Copy, *_targetMap, 1);
      _ICI =  new Epetra_Vector(*_targetMap);
      _dT = dT;
      _check = true;
    }

    //! Destructor
    ~EP_CrankNicholson() {
      delete _stiffness; delete _problem; delete _solver; delete _D;
      delete _C; delete _ionicProblem; delete _ionicSolver; delete _ICI;
      delete _CL;
    }

    //! Defining pure virtual function in base class
    void Solve(const Real dtDiff = 1.);
  };
}

#endif

#include "EP_CrankNicholson.h"

namespace voom{
  
  void EP_CrankNicholson::Initialize() {
    Real *capacitance;
    _capacitance->ExtractView(&capacitance);    
    if (_comm->MyPID() == 0)
      std::cout << "** Begin Constructing System Matrices " << std::endl;
    /* 
       Get Body to fill in entries first. Note we solve diffusion for dT/2. 
       Hence for entries in A we need $\Delta t/2$ which is _dT/4
    */
    // Fill in C_lumped in relevant matrices
    for(int i = 0; i < _nLocalNodes; i++) {
      Real Value = capacitance[i];
      int gRow = _localNodesID[i];
      _CL->InsertGlobalValues(1, &gRow, 1, &gRow, &Value);
    }
    // Assemble it over all Processors and set up the solver
    _CL->GlobalAssemble();
    _myBody[0]->getDiffusionMatrix(_stiffness, _dT/4.);
      
    // Get _D assembled globally. Getting DV is faster with Epetra Calls
    _myBody[0]->getDiffusionMatrix(_D, -_dT/4.);

    // For FLL add capacitance to sitffness
    _myBody[0]->getCapacitanceMatrix( _stiffness );
    _myBody[0]->getCapacitanceMatrix( _D );
    // For smoothing we add alpha*D to capacitance matrix
    if (_myBody[0]->useDiffusiveSmoothing() ) {
      const Real factor = _myBody[0]->getSmoothingCoeff()*4./_dT;
      _myBody[0]->getDiffusionMatrix( _stiffness, factor );
      _myBody[0]->getDiffusionMatrix(_D, -factor);
      _myBody[0]->getDiffusionMatrix(_C, factor*_dT/4.);
    }

    _stiffness->GlobalAssemble();
    _D->GlobalAssemble();

    // Create Capacitance Matrix
    _myBody[0]->getCapacitanceMatrix(_C);
    // Assemble C over processors
    _C->GlobalAssemble();
    if (_myBody[0]->useDiffusiveSmoothing() ) {
      Real CLNorm = _CL->NormFrobenius();
      Real CFNorm = _C->NormFrobenius();
      if (_comm->MyPID() == 0) {
	std::cout << "** Lumped capacitance Norm    : " << CLNorm << "\n";
	std::cout << "** Consistent capacitance Norm: " << CFNorm << "\n";
      }
    }

    // Set up Solver Objects
    _problem = new Epetra_LinearProblem(_stiffness, _Voltage, 
					_diffusiveCurrent);
    
    _solver = new AztecOO(*_problem);
    _solver->SetAztecOption(AZ_precond, AZ_Jacobi);
    _solver->SetAztecOption(AZ_solver, AZ_cg);
    _solver->SetAztecOption(AZ_diagnostics, AZ_none);
    _solver->SetAztecOption(AZ_output, 0);

    // Ionic Problem Solver
    _ionicProblem = new Epetra_LinearProblem(_CL, _Voltage, 
					_diffusiveCurrent);
    
    _ionicSolver = new AztecOO(*_ionicProblem);
    _ionicSolver->SetAztecOption(AZ_precond, AZ_Jacobi);
    _ionicSolver->SetAztecOption(AZ_solver, AZ_cg);
    _ionicSolver->SetAztecOption(AZ_diagnostics, AZ_none);
    _ionicSolver->SetAztecOption(AZ_output, 0);
    _check = false;
    if (_comm->MyPID() == 0)
      std::cout << "** Finished Constructing System Matrices" << std::endl;
  }

  void EP_CrankNicholson::Solve(const Real dtDiff) {
    /* \brief
      The steps performed are
      - Solve for diffusive current for a time dT/2
      - Solve for ionic current for time dT taking _ionicTimestep as 
         incremental time step
      - Solve for diffusive current for a time dT/2
    */    
    Real *voltage, *capacitance, *diffCurrent, *ionicCurrent, *ici;
    const int* indexMap = _myMesh[0]->getIndexMap();    

    // Check if Matrix has been built
    if (_check) Initialize();
    
    // Extract view to do stuff faster
    _Voltage->ExtractView(&voltage);
    _capacitance->ExtractView(&capacitance);
    _diffusiveCurrent->ExtractView(&diffCurrent);
    _ionicCurrent->ExtractView(&ionicCurrent);
    _ICI->ExtractView(&ici);

    // Step 1: Compute the diffusive Current from body
    _diffusiveCurrent->PutScalar( 0. );

    // Get Diffusion Current from Body. This returns -DV_n
    _D->Multiply(false, *_Voltage, *_diffusiveCurrent);
    _solver->Iterate(1000,1E-7);

    if( _verbose)
      if (_comm->MyPID() == 0)
	std::cout << "Iterations: " << _solver->NumIters() << " Residual: "
		  << _solver->TrueResidual() << std::endl;
    // Update Ghost Voltages
    _ghostVoltage->Import(*_Voltage, *_exporter, Insert);
    //    std::cout << "Diff Solve 1\n" << *_Voltage << "\n";
    // Step 2: Get Voltage from Ionic current only first body only
    Real total = 0.;
    Real dt = (_ionicTimeStep < _dT) ? _ionicTimeStep : _dT;
    while ( total < _dT ){
      // get Ionic Current
      _myBody[0]->getIonicCurrent(_ionicCurrent, _Voltage, _ghostVoltage, dt); 
      /* 
	 We are trying to solve M_L dV/dt = M I (M_L - > Lumped Mass)
	 Instead we are doing C_lumped dV/dt = C I
	 We have C and C_lumped calculated. We are using them instead.
      */
      _CL->Multiply(false, *_ionicCurrent, *_ICI);

      // Get CL*V^n. Ionic Model retruns dv/dt
      _CL->Multiply(false, *_Voltage, *_diffusiveCurrent);

      for(int i = 0; i < _nLocalNodes; i++)
	diffCurrent[i] += dt* ici[i];

      // Solve
      _ionicSolver->Iterate(1000, 1E-10);
      if (_verbose)
	if (_comm->MyPID() == 0)
	  std::cout << "Iterations: " << _ionicSolver->NumIters() << 
	    " Residual: " << _ionicSolver->TrueResidual() << std::endl;
      total += dt;
      dt = ( total + dt < _dT) ? dt : _dT - total;
      // Update Ghost Voltages
      _ghostVoltage->Import(*_Voltage, *_exporter, Insert);
      //      std::cout << "Ionic Solve\n" << *_Voltage << "\n";
    }

    // Step 3: Same as step 1
    _diffusiveCurrent->PutScalar( 0. );

    // Get Diffusion Current from Body. This returns -DV_n
    _D->Multiply(false, *_Voltage, *_diffusiveCurrent);
    _solver->Iterate(1000,1E-7);
    if( _verbose)
      if (_comm->MyPID() == 0)
	std::cout << "Iterations: " << _solver->NumIters() << " Residual: "
		  << _solver->TrueResidual() << std::endl;
    //    std::cout << "Diff Solve 2\n" << *_Voltage << "\n";
  } // End of Solve routine

} // namespace

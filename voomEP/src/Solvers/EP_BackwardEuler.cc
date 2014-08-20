#include "EP_BackwardEuler.h"

namespace voom{
  
  void EP_BackwardEuler::Initialize() {
    /*! \brief 
       Now add entries from Capacitance to _stiffness matrix. Capacitance
       is a diagonal matrix stored as a vector. We will add entries to the 
       diagonal of stiffness matrix.
    */
    if (_comm->MyPID() == 0)
      std::cout << "** Begin Constructing System Matrices " << std::endl;
    // Get Body to fill in entries first
    _myBody[0]->getDiffusionMatrix(_stiffness, _dT/2.);
    
    // Capacitance should be computed before this step
    Real *capacitance;
    const int *indexMap = _myMesh[0]->getIndexMap();
    _capacitance->ExtractView(&capacitance);
    

    // Do for Local Nodes First
    for(int i = 0; i < _nLocalNodes; i++) {
      int id = indexMap[ _localNodesID[i] ] - 1;
      _stiffness->InsertGlobalValues(1, &_localNodesID[i], 1, 
				     &_localNodesID[i], &capacitance[id]);
    }

    
    // Assemble it over all Processors and set up the solver
    _stiffness->GlobalAssemble();

    // Set up Solver Objects
    _problem = new Epetra_LinearProblem(_stiffness, _Voltage, 
					_diffusiveCurrent);
    
    _solver = new AztecOO(*_problem);
    _solver->SetAztecOption(AZ_precond, AZ_Jacobi);
    _solver->SetAztecOption(AZ_solver, AZ_cg);
    _solver->SetAztecOption(AZ_diagnostics, AZ_none);
    _solver->SetAztecOption(AZ_output, 0);
    _check = false;
    if (_comm->MyPID() == 0)
      std::cout << "** Finished Constructing System Matrices" << std::endl;
  }

  void EP_BackwardEuler::Solve(const Real dtDiff) {
    /*! \brief
      The steps performed are
      1. Solve for diffusive current for a time dT/2
      2. Solve for ionic current for time dT taking _ionicTimestep as 
         incremental time step
      3. Solve for diffusive current for a time dT/2
    */    
    Real *voltage, *capacitance, *diffCurrent;
    
    // Check if Matrix has been built
    if (_check) Initialize();
    
    // Extract view to do stuff faster
    _Voltage->ExtractView(&voltage);
    _capacitance->ExtractView(&capacitance);
    _diffusiveCurrent->ExtractView(&diffCurrent);
    
    // Step 1: Compute the diffusive Current from body
    for(int i = 0; i < _nLocalNodes; i++)
      diffCurrent[i] = capacitance[i] * voltage[i]; 
    _solver->Iterate(1000,1E-10);

    if( _verbose)
      std::cout << "Iterations: " << _solver->NumIters() << " Residual: "
		<< _solver->TrueResidual() << std::endl;

    // Update Ghost Voltages
    _ghostVoltage->Import(*_Voltage, *_exporter, Insert);    
    
    // Step 2: Get Voltage from Ionic current only first body only
    _myBody[0]->getVoltageFromIonicCurrent(_Voltage, _ghostVoltage, _dT, 
					   _ionicTimeStep, _nLocalNodes);

    // Step 3: Same as step 1
    for(int i = 0; i < _nLocalNodes; i++)
      diffCurrent[i] = capacitance[i] * voltage[i]; 

    _solver->Iterate(1000,1E-10);
    if( _verbose)
      std::cout << "Iterations: " << _solver->NumIters() << " Residual: "
		<< _solver->TrueResidual() << std::endl;


  } // End of Solve routine

} // namespace

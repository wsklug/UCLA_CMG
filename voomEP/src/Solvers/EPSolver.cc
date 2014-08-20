#include "EPSolver.h"
#include "EP_ExplicitSolver.h"
#include "EP_BackwardEuler.h"
#include "EP_CrankNicholson.h"

namespace voom{
  EPSolver::EPSolver(Epetra_MpiComm *comm, const MeshContainer& mesh, 
		     const BodyContainer& cardiacBody, 
		     const SmoothCoeff& coeff) {
    _comm             = comm;
    _myMesh           = mesh;
    _myBody           = cardiacBody;
    _coeff            = coeff;
    _verbose          = false;
    _displacementData = false;
    _maxIter          = 20;
    _maxCGIter        = 20000;

    // Use data from the first mesh. 
    _nLocalNodes    = _myMesh[0]->getNumLocalNodes();
    _nGhostNodes    = _myMesh[0]->getNumGhostNodes();
    _ghostNodesID   = _myMesh[0]->getGhostNodeID();
    _localNodesID   = _myMesh[0]->getLocalNodeID();
    _indexMap       = _myMesh[0]->getIndexMap();
    int nRest       = _myMesh[0]->getNumberOfRestraints();
    int nDof        = _myMesh[0]->getDOF();

    _dlocalNodeID    = _myMesh[0]->getDisplacementLocalNodeID();
    _dghostNodeID    = _myMesh[0]->getDisplacementGhostNodeID();
    
    // Define Epetra vectors
    _targetMap = new Epetra_Map(-1, _nLocalNodes, _localNodesID, 0, *comm);
    _sourceMap = new Epetra_Map(-1, _nGhostNodes, _ghostNodesID, 0, *comm);

    _exporter = new Epetra_Export(*_sourceMap, *_targetMap);

    _dispTargetMap = new Epetra_Map(-1, nDof*_nLocalNodes, _dlocalNodeID, 0, 
				    *comm);
    _dispSourceMap = new Epetra_Map(-1, nDof*_nGhostNodes, _dghostNodeID, 0, 
				    *comm);

    _dispExporter = new Epetra_Export(*_dispSourceMap, *_dispTargetMap);
    
    // Local Values
    _diffusiveCurrent = new Epetra_Vector(*_targetMap);
    _Voltage          = new Epetra_Vector(*_targetMap);
    _capacitance      = new Epetra_Vector(*_targetMap);
    _ionicCurrent     = new Epetra_Vector(*_targetMap);

    // _U holds increments in displacement solution.
    _U                = new Epetra_Vector(*_dispTargetMap);
    _force            = new Epetra_Vector(*_dispTargetMap);
    _fInternal        = new Epetra_Vector(*_dispTargetMap);
    _displacement     = new Epetra_Vector(*_dispTargetMap);

    // Ghost Values
    _ghostDiffusiveCurrent = new Epetra_Vector(*_sourceMap);
    _ghostVoltage          = new Epetra_Vector(*_sourceMap);
    _ghostCapacitance      = new Epetra_Vector(*_sourceMap);
    _ghostU                = new Epetra_Vector(*_dispSourceMap);
    _ghostDisplacement     = new Epetra_Vector(*_dispSourceMap);
    _ghostfInternal        = new Epetra_Vector(*_dispSourceMap);
    // tolerance for Mechanics solve
    _tolerance = 1E-5;

    // CG Solver Tolerance
    _cgSolverTolerance = 1E-7;

    // Print if model is Iso/Aniso
    if (_comm->MyPID() == 0){
	if ( _myMesh[0]->isIsotropic() )
	    std::cout << "** ISOTROPIC MODEL\n";
	else std::cout << "** ANISOTROPIC MODEL\n";
    }
  }

  //! Impose Essential Boundary conditions
  void EPSolver::imposeEssentialBoundaryConditions(Epetra_FECrsMatrix *K) {
    const std::vector<int>& nodalID     = _myMesh[0]->getRestraintNodes();
    const std::vector<int>& nodalDOF    = _myMesh[0]->getRestraintDOF();
    const std::vector<Real>& constraint = _myMesh[0]->getRestraintValues();
    const int DOF = _myMesh[0]->getDOF();
    const int *dindexMap = _myMesh[0]->getDisplacementIndexMap();
    Real *force;
    _force->ExtractView( &force );
    
    for(int i = 0; i < nodalID.size(); i++) {
      int gRow = (nodalID[i] - 1)*DOF + (nodalDOF[i] - 1);
      int id = dindexMap[ gRow ];
      if ( id == 0 ) continue;
      if ( id > 0 ) {
	int gCol = gRow;
	Real value = 1.; 
	K->InsertGlobalValues(1, &gRow, 1, &gCol, &value);
	force[id-1] += value * constraint[i];
      }
    }// i loop
  }

  //! Impose Natural Boundary Conditions
  void EPSolver::imposeNaturalBoundaryConditions() {
    const std::vector<int>& nodalID     = _myMesh[0]->getForceNodes();
    const std::vector<int>& nodalDOF    = _myMesh[0]->getForceDOF();
    const std::vector<Real>& values     = _myMesh[0]->getForceValues();
    const int DOF = _myMesh[0]->getDOF();
    const int *dindexMap = _myMesh[0]->getDisplacementIndexMap();
    Real *force;
    _force->ExtractView( &force );
    for(int i = 0; i < nodalID.size(); i++) {
      int id = dindexMap[ (nodalID[i] - 1)*DOF + (nodalDOF[i] - 1) ];
      if ( id > 0 )  force[id - 1] += values[i];
    }// i loop
  }

  //! Write Voltage Output
  void EPSolver::WriteOutput(const Real time, std::string path){
    char name[200], dname[200];
    Real *voltage, *displacement;
    const int DOF = _myMesh[0]->getDOF();
    // Format Path/Voltage_time.procNum
    sprintf(name,"%s/Voltage_time_%08.2f.%02.0f", path.c_str(), time, 
	    Real(_comm->MyPID()) );
    sprintf(dname,"%s/Displacement_time_%08.2f.%02.0f", path.c_str(), time, 
	    Real(_comm->MyPID()) );

    ofstream out(name, ios::out|ios::binary);

    if( !out.is_open() ){
      std::cerr << "Error Creating File : " << name << std::endl;
      exit(0);
    }
    _Voltage->ExtractView(&voltage);
    for(int i = 0 ; i < _nLocalNodes; i++) 
      out.write((char*)(&voltage[i]), sizeof(Real));

    out.close();
    // if Displacement data exists write them also
    if ( time == 0. ) _displacementData = true;
    if ( _displacementData ) {
      _displacementData = false;
      ofstream disp(dname, ios::out|ios::binary);
      if( !disp.is_open() ){
	std::cerr << "Error Creating File : " << dname << std::endl;
	exit(0);
      }
      _displacement->ExtractView(&displacement);
      for(int i = 0 ; i < _nLocalNodes; i++) 
	for(int j = 0; j < DOF; j++) 
	  disp.write((char*)(&displacement[i*DOF + j]), sizeof(Real));
      disp.close();
    }
  }

  // Write State Internal variables
  void EPSolver::writeRestart(const Real time, std::string path) {
    std::vector<int> numVar;
    std::vector< std::vector<Real> > data;
    _myBody[0]->getStateVariables( numVar, data);
    char name[200];
    sprintf(name,"%s/InterVariables_time_%08.2f.%02.0f", path.c_str(), time, 
	    Real(_comm->MyPID()) );
    ofstream out(name, ios::out|ios::binary);

    if( !out.is_open() ){
      std::cerr << "Error Oepning File : " << name << std::endl;
      exit(1);
    }
    
    // Write the data
    for(int i = 0 ; i < _nLocalNodes; i++) {
      out.write((char*)(&numVar[i]), sizeof(int));
      for(int j = 0; j < numVar[i]; j++)
	out.write((char*)(&data[i][j]), sizeof(Real));
    }
    out.close();
  }

  // Read restart output - Restart capability
  void EPSolver::readRestart(const Real time, std::string path) {
    char name[200], fname[200];
    Real *voltage, v;
    std::vector< std::vector<Real> > data;
    int nVar;
    // Format Path/Voltage_time.procNum
    sprintf(name,"%s/Voltage_time_%08.2f.%02.0f", path.c_str(), time, 
	    Real(_comm->MyPID()) );
    sprintf(fname,"%s/InterVariables_time_%08.2f.%02.0f", path.c_str(), time, 
	    Real(_comm->MyPID()) );
    ifstream out(name, ios::out|ios::binary);
    ifstream cell(fname, ios::out|ios::binary);
    if( !out.is_open() ){
      std::cerr << "Error Oepning File : " << name << std::endl;
      exit(1);
    }
    _Voltage->ExtractView(&voltage);
    data.resize( _nLocalNodes );
    for(int i = 0 ; i < _nLocalNodes; i++) {
      // Read Voltage
      out.read((char*)(&v), sizeof(Real));
      voltage[i] = Real(v);
      // Read Cell internal data
      cell.read((char*)(&nVar), sizeof(int));
      std::vector<Real> intData(nVar);
      for(int j = 0; j < nVar; j++) 
	cell.read((char*)(&intData[j]), sizeof(Real));
      data[i] = intData;
    }
    _myBody[0]->setStateVariables( data );
    out.close();
    cell.close();
  }

  //! Compute lumped capacitance
  void EPSolver::computeCapacitance() {
    _myBody[0]->computeLumpedCapacitance(_capacitance, _ghostCapacitance);
    // For ghost nodes we need to sum up values from all processors
    _capacitance->Export(*_ghostCapacitance, *_exporter, Add);
    _ghostCapacitance->Import(*_capacitance, *_exporter, Insert);
  }

  //! Perform Consistency checks
  void EPSolver::consistencyCheck() {
    bool result = _myBody[0]->consistencyCheck( _material );
    if (result)
      std::cout << "** Consistency Checks PASSED\n";
    else
      std::cout << "** Consistency Checks FAILED\n";
  }

  //! Solve Mechanics Portion
  int EPSolver::solveMechanics() {
    Real error = 100.;
    const int numNodes = _myMesh[0]->getNumLocalNodes() + 
      _myMesh[0]->getNumGhostNodes();
    std::vector<Real> errorValues;
    _force->PutScalar( 0. );
    this->imposeNaturalBoundaryConditions();
    // Update ghost voltages
    _ghostVoltage->Import(*_Voltage, *_exporter, Insert);

    // Set displacementdata to true
    _displacementData = true;
    // Get Gamma from Body
    std::vector<Real> nodalGamma;
    std::map<int, int> nodeID; 
    _myBody[0]->getActiveContraction( nodalGamma, nodeID, numNodes);
    // Interpolate Gamma to Integration Points
    std::vector<Real> gamma(_myMesh[0]->interpolateScalar(nodalGamma, nodeID) );
    // Get Fa Inverse
    std::vector<Tensor3D> 
      FaInv(_myBody[0]->computeActiveDeformationInverse( gamma ));
    int ctr = 0;
    int maxCGIter = _maxCGIter;
    while (error > _tolerance) {
      // Mechanics related initializations
      Epetra_FECrsMatrix *K = new Epetra_FECrsMatrix(Copy, *_dispTargetMap, 5);
      // Set internal force to 0
      _fInternal->PutScalar( 0. );
      _ghostfInternal->PutScalar( 0. );
      
      // Compute deformation Gradient
      _myMesh[0]->computeDeformationGradient();
      // Get Internal force vector
      _myBody[0]->getInternalForce(_material, _fInternal, 
				   _ghostfInternal, FaInv);
      _fInternal->Export( *_ghostfInternal, *_dispExporter, Add);
      _ghostfInternal->Import( *_fInternal, *_dispExporter, Insert);   
      // fint = fext - fint
      _fInternal->Update( 1., *_force, -1.);
      // Get Stiffness matrix. Body calls mesh to compute Fe
      _myBody[0]->getStiffnessMatrix(K, _material, FaInv);
      // Impose Essential Boundary Conditions
      this->imposeEssentialBoundaryConditions( K );   
      // Assemble over all processors
      K->GlobalAssemble();
      // Solve
      _U->PutScalar( 0. ); _ghostU->PutScalar( 0. );
      Real fNorm; _fInternal->Norm2( &fNorm );

      Real Tol = ( fNorm*_cgSolverTolerance > _cgSolverTolerance) ? 
	fNorm*_cgSolverTolerance : _cgSolverTolerance;

      // Aztec Data
      Epetra_LinearProblem *mechanicsProblem = 
	new Epetra_LinearProblem(K, _U, _fInternal);
      // Iterative solver using AztecOO
      AztecOO *mechanicsSolver = new AztecOO(*mechanicsProblem);

      // mechanicsSolver->SetAztecOption(AZ_solver, AZ_cg);
      // mechanicsSolver->SetAztecOption(AZ_solver, AZ_gmres);
      mechanicsSolver->SetAztecOption(AZ_solver, AZ_bicgstab);
      // mechanicsSolver->SetAztecOption(AZ_solver, AZ_cgs);
      // mechanicsSolver->SetAztecOption(AZ_solver, AZ_tfqmr);
      // mechanicsSolver->SetAztecOption(AZ_precond, AZ_dom_decomp);
      mechanicsSolver->SetAztecOption(AZ_precond, AZ_Neumann);
      // mechanicsSolver->SetAztecOption(AZ_precond, AZ_Jacobi);
      mechanicsSolver->SetAztecOption(AZ_overlap,0);
      mechanicsSolver->SetAztecOption(AZ_subdomain_solve, AZ_ilut);
      
      mechanicsSolver->SetAztecOption(AZ_diagnostics, AZ_none);
      mechanicsSolver->SetAztecOption(AZ_output, AZ_warnings);
      mechanicsSolver->SetAztecOption(AZ_conv, AZ_noscaled);

      mechanicsSolver->SetAztecOption(AZ_drop, Tol);
      // Iterative Solver
      mechanicsSolver->Iterate(maxCGIter, Tol);    
      const Real *status = mechanicsSolver->GetAztecStatus();

      if ( status[AZ_why] == AZ_ill_cond ) {
	std::cout << "@@ WARNING: Stiffness matrix is Ill Conditioned\n";
      }
      if (status[AZ_r] > 0.1*error && ctr > 2) 	maxCGIter *= 1.2;
      delete K; delete mechanicsProblem; delete mechanicsSolver;
      K = NULL; mechanicsProblem = NULL; mechanicsSolver = NULL;

      _fInternal->Norm2( &error );
      _U->Export( *_ghostU, *_dispExporter, Add);
      _ghostU->Import( *_U, *_dispExporter, Insert);
      // Update displacement EpetraVector  _disp = _disp + _U;
      _displacement->Update(1., *_U, 1.);
      _myMesh[0]->updateNodalPosition( _U , _ghostU );
      ctr++;
      if (_comm->MyPID() == 0)
	std::cout << "** NR Iteration: " << ctr << " Error: " << error << "\n";
      if (ctr == _maxIter) {
	error = 0.;
	if (_comm->MyPID() == 0) {
	  std::cout << "@@ Newton Raphson needs more than maxIter (" << _maxIter
		    << ") steps\n";
	  std::cout << "@@ Try reducing force increment per step\n";
	}
      }
      errorValues.push_back( error );
    } // while loop
    return 0;
  }

  //! Gateway to Solver Class
  EPSolver* EPSolver::New(Epetra_MpiComm *comm, const MeshContainer& mesh,
			  const BodyContainer& cardiacBody, 
			  const SmoothCoeff& coeff, Method SolverType) {
    EPSolver* mySolver;
    switch(SolverType)
      {
      case 0:
	mySolver = new EP_ExplicitSolver(comm, mesh, cardiacBody, coeff);
	break;
      case 1:
	mySolver = new EP_BackwardEuler(comm, mesh, cardiacBody, coeff);
	break;
      case 2:
	mySolver = new EP_CrankNicholson(comm, mesh, cardiacBody, coeff);
	break;
      default:
	{
	  std::cout << "Unknown Solver specified\n";
	  exit(0);
	}
      }
    return mySolver;
  }

}// End of namespace

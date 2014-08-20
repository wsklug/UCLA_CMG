/*!
  EPSimulator - VoomEP based application for solving CardiacEP and Mechanics.
  Modified version of Heart3D. Heart3D uses command line version whereas
  EPSimulator uses an input file and hence is more flexbile.

  Shankarjee Krishnamoorthi
 */
// VoomEP Header Files
#include "CardiacBody.h"
#include "Mesh.h"
#include "EPSolver.h"
#include "Material.h"
#include "QuadQuadrature.h"

// Standard C++ Header Files
#include <iostream>
#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>


struct stimData {
  int stimType;
  Real center[3], radius, normal[3];
  int stimAxis;
  Real limitValue;
  string nodeSetFile;
  Real t_start, t_stop, i_stim;
};
using namespace boost;
using namespace std;

void pLine() {
  for(int i = 0; i < 59;i++) cout << "=";
  cout << endl;
}

//! A Small helper function
MeshType getMeshType(string meshType) {
  to_upper( meshType );
  if ( meshType == "QUAD" ) return QUADMESH;
  if ( meshType == "TRIA" ) return TRI3MESH;
  if ( meshType == "HEX" ) return HEXMESH;
  if ( meshType == "TET" ) return TETMESH;
  if ( meshType == "MIX3D" ) return MIXEDMESH3D;
  if ( meshType == "MF2D" ) return MESHFREE2D;
  if ( meshType == "MF3D" ) return MESHFREE3D;
  
  cerr << "Unknown Mesh Type provided. Allowed Types: ";
  cerr << "QUAD, TRIA, HEX, TET, MF2D, MF3D,  MIX3D (Case Sensitive)\n";
  exit(1);
}

//! Another helper for Cell model type
CellModelType getCellModelType(string cellModel){
  to_upper( cellModel );
  if ( cellModel == "LUORUDY" ) return LUORUDY;
  if ( cellModel == "LUORUDYWITHGATE" ) return LUORUDYWITHGATE;
  if ( cellModel == "TUSSCHER" ) return TUSSCHER;
  if ( cellModel == "MAHAJAN") return MAHAJAN;
  if ( cellModel == "MAHAJAN_FAIL") return MAHAJAN_FAIL;

  cerr << "Unknown Cell Model Type provided. Allowed Types: ";
  cerr << "LUORUDY, LUORUDYWITHGATE, TUSSCHER, MAHAJAN (Case Sensitive)\n";
  exit(1);
}

// Helper Function for Solver
Method getSolverType(string solverType) {
  if ( solverType == "EXPLICIT" ) return EXPLICIT;
  if ( solverType == "BACKWARDEULER" ) return BACKWARDEULER;
  if ( solverType == "CRANKNICKOLSON" ) return CRANKNICHOLSON;
  
  cerr << "Unknown Solver type provided. Allowed Type: ";
  cerr << "EXPLICIT, BACKWARDEULER and CRANKNICKOLSON (Case Sensitive)" <<
    endl;
  exit(1);
}

// Central Processing
int main(int argc, char** argv) {
  // Input Parameters
  string modelName, outputDir = "Output/", cellModel = "MAHAJAN", 
    solverType = "CRANKNICKOLSON", meshType = "MIX3D", matType = "HOLZAPFEL", 
    initialStates;
  int quadOrder = 2;
  Real dt_Ionic = 0.1, dt_diffusion = 0.05, simulation_time = 300., dt = 0.1;
  Real frequency = 1., twoD_Radius = 0.0275, dtMech = 9999.;
  Real mechNRToler = 1E-5;
  bool consistencyCheck = false;
  std::vector< struct stimData > stimulusSet;
  std::vector< Real > purkinjeStates;
  std::vector< Real > myocardiumStates;
  Real matData[] = {1.053e-3, 6.0e4, 59.0e0, 8.023, 18.472e0, 16.026, 2.481e0, 
		    11.120,0.216e0, 11.436e0};
  // Real D[] = { 0.001, 0.0005, 0.00025 };     // Normal Heart
  // Real D[] = {0.0005, 0.00025, 0.00025};	// 2:1:1 (Modeling Connexin43)  
  // Real D[] = {0.00025, 0.00025, 0.00025};	// 1:1:1	  
  Real D[] = { 0.00075, 0.000375, 0.0001875 };  // VF with slower conduction velocity (~75% of D)
  // Real D[] = { 0.0005, 0.0005, 0.0005 };  // Isotropic with cross-fiber velocity
  Real Dp[] = {0.0032, 0.0032, 0.0032}; // 2 times DpOrig
  Real DpLeft[] = { 0.0032, 0.0032, 0.0032 }; // 4 times DpOrig
  // Real Dp[] = {0.0064, 0.0064, 0.0064}; // 2 times DpOrig
  // Real DpLeft[] = { 0.0064, 0.0064, 0.0064 }; // 4 times DpOrig
  ifstream inp;
  int linSolveMaxIteration = 20000;
  Real linSolveTolerance = 1E-6;
  int maxNRIteration = 100, convRateCheck = 8;
  bool readRestart = false, writeRestart = false;
  Real restartReadTime = 99999., restartWriteTime = 99999.;
  bool smoothing = false;
  Real alpha = 0.1;

  // Mpi Related Stuff
  MPI_Init(&argc, &argv);
  Epetra_MpiComm* mpicomm = new Epetra_MpiComm(MPI_COMM_WORLD);
  if (mpicomm->MyPID() == 0 ) {
    pLine();
    cout << "File     : " << __FILE__ << endl;
    cout << "Compiled : " << __DATE__ << " " << __TIME__ << endl;
    pLine(); 
  }

  if (argc == 1) {
    if ( mpicomm->MyPID() == 0 ) 
      std::cerr << "Usage : " << argv[0] << " inputFile \n";
    MPI_Finalize();
    return 1;
  }

  inp.open( argv[1] );
  if (!inp.is_open() ) {
    if ( mpicomm->MyPID() == 0 )
      std::cerr << "Error Opening input file: " << argv[1] << 
	". Exiting " << std::endl;
    MPI_Finalize();
    return 1;
  }

  string line;
  // Parse Input File
  while (getline(inp, line) ) {
    trim(line);
    if ( find_first(line,"#") || find_first(line, "$") || find_first(line,"*"))
      continue; // * or $ or # is comment
    vector<string> strs;
    split( strs, line, is_any_of("\t "), token_compress_on);
    to_upper( strs[0] );
    if (strs[0] == "MODEL") modelName = strs[1];
    if (strs[0] == "OUTPUTDIR") outputDir = strs[1];
    if (strs[0] == "CELLMODEL") cellModel = strs[1];
    if (strs[0] == "MESHTYPE") meshType = strs[1];
    if (strs[0] == "QUADORDER") quadOrder = atoi(strs[1].c_str());
    if (strs[0] == "IONICDT") dt_Ionic = atof(strs[1].c_str());
    if (strs[0] == "DT") dt = atof( strs[1].c_str() );
    // Read Restart information
    if (strs[0] == "READRESTART") {
      readRestart = true;
      restartReadTime = atof( strs[1].c_str() );
    }
    // Write Restart Information
    if (strs[0] == "WRITERESTART") {
      writeRestart = true;
      restartWriteTime = atof( strs[1].c_str() );
    }
    if (strs[0] == "CONSISTENCYCHECK") {
      to_upper(strs[1]);
      if (strs[1] == "YES") consistencyCheck = true;
    }
    if (strs[0] == "DIFFUSIONDT") dt_diffusion = atof(strs[1].c_str());
    if (strs[0] == "SIMULATIONTIME") simulation_time = atof(strs[1].c_str());
    if (strs[0] == "FREQUENCY") frequency = atof(strs[1].c_str());
    if (strs[0] == "TWODRADIUS") twoD_Radius = atof(strs[1].c_str());
    // Tissue Diffusion
    if (strs[0] == "TISSUEDIFFUSION")
      for (int j = 0; j <3; j++) D[j] = atof( strs[j+1].c_str() );
    // Purkinje Diffusion
    if (strs[0] == "PURKINJEDIFFUSION")
      for (int j = 0; j <3; j++) Dp[j] = atof( strs[j+1].c_str() );
    // Purkinje LV Diffusion
    if (strs[0] == "LVPURKINJEDIFFUSION")
      for (int j = 0; j <3; j++) DpLeft[j] = atof( strs[j+1].c_str() );
    // Mechanics DT
    if (strs[0] == "MECHDT") dtMech = atof(strs[1].c_str());
    // Mechanics NR Iteration
    if (strs[0] == "MAXNRITER") maxNRIteration = atoi( strs[1].c_str() );
    // Mechanics Lin Solve Tolerance
    if (strs[0] == "LINSOLVETOLER")
      linSolveTolerance = atof( strs[1].c_str() );
    // Mechanics Lin Solve Iteration
    if (strs[0] == "MAXLINSOLVEITER") 
      linSolveMaxIteration = atoi( strs[1].c_str() );
    // Convergence rate check iteration
    if (strs[0] == "CONVERGENCERATECHECK")
      convRateCheck = atoi( strs[1].c_str() );
    // Material Data
    if (strs[0] == "MATERIAL") {
      matType = strs[1]; to_upper(matType);
      for( int j = 2; j < strs.size(); j++) 
	matData[j-2] = atof(strs[j].c_str() );
    }
    // Diffusion smoothiing
    if (strs[0] == "SMOOTHING") {
      smoothing = true;
      if (strs.size() == 2 ) alpha = atof( strs[1].c_str());
    }

    // Initial States
    if (strs[0] == "INITIALSTATES") {
      initialStates = strs[1];
    }
    if (strs[0] == "MECHTOLERANCE") mechNRToler = atof( strs[1].c_str() );
    // Read a Stimulus input
    if (strs[0] == "STIMULUS") {
      struct stimData stimulus;
      map< string, int> Type, Axis;
      stimulus.i_stim = 50000.;
      Type["NODAL" ]=0; Type["CIRC"]=1; Type["ZONAL"]=2; Type["PLANE"]=3;
      Axis["X"] = 0; Axis["Y"] = 1; Axis["Z"] = 2; to_upper(strs[1]);
      stimulus.stimType = Type[ strs[1] ]; // Stim Type
      switch ( stimulus.stimType ) {
      case 0: // Nodal Stimulus
	stimulus.nodeSetFile = strs[2];
	stimulus.t_start = atof(strs[3].c_str() );
	stimulus.t_stop  = atof(strs[4].c_str() );
	if (strs.size() == 6) stimulus.i_stim  = atof(strs[5].c_str() );
	break;
      case 1: // Circular Stimulus
	stimulus.center[0] = atof( strs[2].c_str() );
	stimulus.center[1] = atof( strs[3].c_str() );
	stimulus.center[2] = atof( strs[4].c_str() );
	stimulus.radius    = atof( strs[5].c_str() );
	stimulus.t_start = atof(strs[6].c_str() );
	stimulus.t_stop  = atof(strs[7].c_str() );
	if (strs.size() == 9) stimulus.i_stim  = atof(strs[8].c_str() );	
	break;
      case 2: // Zonal Stim
	stimulus.stimAxis = Axis[ strs[2] ];
	stimulus.limitValue = atof( strs[3].c_str() );
	stimulus.t_start = atof( strs[4].c_str() );
	stimulus.t_stop  = atof( strs[5].c_str() );
	if (strs.size() == 7) stimulus.i_stim  = atof(strs[6].c_str() );	
	break;
      case 3: // Planar Stimulus
	stimulus.center[0] = atof( strs[2].c_str() );
	stimulus.center[1] = atof( strs[3].c_str() );
	stimulus.center[2] = atof( strs[4].c_str() );
	stimulus.normal[0] = atof( strs[5].c_str() );
	stimulus.normal[1] = atof( strs[6].c_str() );
	stimulus.normal[2] = atof( strs[7].c_str() );
	stimulus.t_start = atof(strs[8].c_str() );
	stimulus.t_stop  = atof(strs[9].c_str() );
	if (strs.size() == 11) stimulus.i_stim  = atof(strs[10].c_str() );	
	break;
      default:
	std::cerr << "Unknown Stimulus. Exiting...\n";
	MPI_Finalize();
	return 1;
      } // Switch statement
      stimulusSet.push_back( stimulus );
    }
      
  }
  inp.close();

  // Echoing User input parameters for the model
  if (mpicomm->MyPID() == 0 ) {
    cout << "MODEL PAREMETERS " << endl;
    cout << "***********************************************************\n";
    pLine();
    cout << "MESH SETTINGS" << endl; pLine();
    cout << "ModelName        : " << modelName << "\n";
    cout << "Output Directory : " << outputDir << "\n";
    cout << "Quadrature order : " << quadOrder << "\n";
    cout << "Mesh Type        : " << meshType << "\n\n"; pLine();
    cout << "ELECTROPHYSIOLOGY SETTINGS" << endl; pLine();
    cout << "CellModel(Tissue): " << cellModel << "\n";
    if (initialStates != "") {
      cout << "Initial State File Provided " << endl; }
    cout << "Time Step        : " << dt << "(ms)" << endl;
    cout << "Tissue D         : " << D[0] << " " << D[1] << " "
	 << D[2] << endl;
    cout << "Purkinje D       : " << Dp[0] << " " << Dp[1] << " " 
	 << Dp[2] << endl;
    cout << "LVPurkinje D     : " << DpLeft[0] << " " << DpLeft[1] << " " 
	 << DpLeft[2] << endl;
    cout << "Diffusion dt     : " << dt_diffusion << "(ms)" << endl;
    cout << "Ionic dt         : " << dt_Ionic << "(ms)"<< endl;
    cout << "Simulation time  : " << simulation_time << "(ms)" 
	 << endl << endl; pLine();
    cout << "STIMULUS STEPS" << endl; pLine();
    cout << "nSteps           : " << stimulusSet.size() << endl;
    for(int i = 0; i < stimulusSet.size(); i++ ) {
      cout << "Step             : " << i + 1 << endl;
      switch (stimulusSet[i].stimType) {
      case 0:
	cout << "Node Set Stimulus: " << stimulusSet[i].nodeSetFile << "\n";
	break;
      case 1:
	{
	  cout << "Circular Stimulus " << endl;
	  cout << "Center           : " << scientific << 
	    stimulusSet[i].center[0] << " " << stimulusSet[i].center[1] << " "  
	       << stimulusSet[i].center[2] << endl;
	  cout << "Radius           : " << scientific << stimulusSet[i].radius 
	       << endl;
	  break;
	}
      case 2: 
	{
	  map<int, string> Axis;
	  Axis[0] = "X"; Axis[1] = "Y"; Axis[2] = "Z";
	  cout << "Zonal Stimulus" << endl;
	  cout << "Axis             : " << Axis[stimulusSet[i].stimAxis] <<endl;
	  cout << "Limit Value      : " << stimulusSet[i].limitValue << 
	    "(cm)" << endl;
	  break;
	}
      case 3:
	{
	  cout << "Planar Stimulus " << endl;
	  cout << "Center           : " << scientific << 
	    stimulusSet[i].center[0] << " " << stimulusSet[i].center[1] << " "  
	       << stimulusSet[i].center[2] << endl;
	  cout << "Normal           : " << scientific << 
	    stimulusSet[i].normal[0] << " " << stimulusSet[i].normal[1] << " "  
	       << stimulusSet[i].normal[2] << endl;
	  break;
	}
      } // Switch statement
      cout << resetiosflags( ios::floatfield ) ;
      cout << "Stimulus Start   : " << stimulusSet[i].t_start 
	   << "(ms)" << endl;
      cout << "Stimulus Stop    : " << stimulusSet[i].t_stop
	   << "(ms)" << endl;
      cout << "Stimulus Current : " << stimulusSet[i].i_stim 
	   << " (uA/cc)" << endl << endl;
    } // if loop

    if (smoothing) {
      pLine(); cout << "DIFFUSIVE SMOOTHING\n"; pLine();
      cout << "Smoothing factor : " << alpha << endl;
    }

    // Restart Information
    pLine(); cout << "RESTART INFORMATION" << endl; pLine();
    cout << "Restart Read     : " << restartReadTime << "(ms)" << endl;
    cout << "Restart Write    : " << restartWriteTime << "(ms)" << endl 
	 << endl;

    pLine(); cout << "EPSOLVER PARAMETERS" << endl; pLine();
    cout << "Output Frequency : " << frequency << "(ms)" << endl;
    cout << "2D ElementRadius : " << twoD_Radius << "(cm)" << endl;
    cout << "Solver Type      : " << solverType << endl << endl;

    pLine(); cout << "MECHANICS PARAMETERS" << endl; pLine();
    cout << "Material Type    : " << matType << endl;
    cout << "Mechanics dT     : " << dtMech << "(ms)" << endl;
    cout << "Convergence Check: " << convRateCheck << endl;
    cout << "NR Tolerance     : " << mechNRToler << endl;
    cout << "Max NR Iteration : " << maxNRIteration << endl;
    cout << "NR Solver Toler  : " << linSolveTolerance << endl;
    cout << "Max LinSolve Iter: " << linSolveMaxIteration << endl << endl;

    cout << "***********************************************************\n";
  }
  cout << resetiosflags( ios::floatfield ) ;

  // Create Mesh Object
  MeshType mType = getMeshType(meshType);
  MeshContainer meshes;
  Mesh *myMesh = Mesh::New( mpicomm, modelName, quadOrder, mType );
  myMesh->setRadius( twoD_Radius );
  meshes.push_back(myMesh);
  if (mpicomm->MyPID() == 0 ) cout << "Created Mesh Object\n";

  // Real Initial States
  if (initialStates != "") {
    string cellTypeFlag;
    inp.open(initialStates.c_str());
    while(getline(inp,line)) {
      to_upper( line );
      if (line.empty()){
	cellTypeFlag = "";
	continue;
      }
      else if (line.find("PURKINJE") != string::npos) {
	cellTypeFlag  = "PURKINJE";
	continue;
      }
      else if (line.find("CELL") != string::npos) {
	cellTypeFlag = "CELL";
	continue;
      }

      if (cellTypeFlag == "PURKINJE")
	purkinjeStates.push_back(atof(line.c_str()));
      else if (cellTypeFlag == "CELL")
	myocardiumStates.push_back(atof(line.c_str()));
    }
    inp.close();
    
  }
  
  // Create a Cardiac Body
  BodyContainer bodies;
  // If initial states are provided use this constructor.
  CardiacBody *myBody;
  if (initialStates != "") {
    myBody = CardiacBody::New(myMesh, getCellModelType(cellModel), purkinjeStates, myocardiumStates); }
  else {
    myBody = CardiacBody::New(myMesh, getCellModelType(cellModel)); }
  bodies.push_back(myBody);

  // Set Xi for Body
  if (mpicomm->MyPID() == 0 ) cout << "Created Body Object\n";
  
  myBody->setDiffusion( D );
  myBody->setPurkinjeDiffusion( Dp );
  myBody->setPurkinjeDiffusionInLV( Dp );
  if (smoothing) 
    myBody->setDiffusiveSmoothing( true, alpha);
  
  // Smoothing Coeffs
  SmoothCoeff coeff;
  coeff.push_back(1.);

  // Create Material 
  Material *Tissue;
  // Real compData[] = {100., 1e4, 0.20 };
  Real compData[] = {100., 1e3, 0.30 };
  if (matType == "COMPNEOHOOKEAN") 
  Tissue = Material::New( compData, COMPNEOHOOKEAN);
  if (matType == "HOLZAPFEL") 
    Tissue = Material::New( matData, COMPMYOCARDIUM );

  Real fiberData[] = {100., 1E4, 0.};
  Material* Fiber  = Material::New( fiberData, COMPNEOHOOKEAN);

  const long int steps_per_frame = int(frequency/dt);
  long int disp_steps_per_frame = int(dtMech/dt);
  const long int max_disp_spf = disp_steps_per_frame;

  // Compute Capacitance for body
  const Real C_m = 1.;// Units are uF/cm^2.
  
  myBody->setNodalCapacitance( C_m );
  // Create Solver Object
  EPSolver *mySolver = EPSolver::New(mpicomm, meshes, bodies, coeff,
				     getSolverType(solverType));
  // const Real Initial_Voltage = -85.23; // Normal
  const Real Initial_Voltage = -81.1751;  // Rapid Pacing
  Real simtime = 0.;
  
  if (mpicomm->MyPID() == 0) 
    cout << "** Number of EssentialBC : " << myMesh->getNumberOfRestraints()
	 << endl;
  // Set initial Voltage
  //mySolver->InitializeVoltage( Initial_Voltage );
  
  // Use this InitializeVoltage to set Voltage based on Paced State[0]:
  mySolver->InitializeVoltage(myBody->getInitialVoltage(), myBody->getInitialGhostVoltage());
  mySolver->setIonicTimeStep(dt_Ionic);
  mySolver->computeCapacitance();
  mySolver->setLinearSolverTolerance( linSolveTolerance );
  mySolver->setMaterial( Tissue );

  // Mechanics Solver Settings
  mySolver->setMaxIteration( maxNRIteration );
  mySolver->setMaxLinSolveIteration( linSolveMaxIteration );
  mySolver->setTolerance( mechNRToler ); // 1 milliDyne force
  
  mySolver->setSolveTime( dt );
  if (mpicomm->MyPID() == 0 ) cout << "Created Solver Object\n";

  if (consistencyCheck)
    mySolver->consistencyCheck();

  long int count = 0, startStimCtr = 0, stopStimCtr = 0;
  long int pulse_start, pulse_stop, restart_write;
  // If some stimulus is given
  if ( stimulusSet.size() >= 1) {
    pulse_start = long( stimulusSet[0].t_start/dt );
    pulse_stop  = long( stimulusSet[0].t_stop/dt );
  } 
  // In case of restart no stimulus might be given
  else {
    pulse_start = -1;
    pulse_stop  = -1;
  }
  restart_write = long(restartWriteTime/dt);

  // Read Restart if needed
  if (readRestart) {
    if (mpicomm->MyPID() == 0)
      cout << "@@ Reading Restart information " << endl;
    simtime = restartReadTime;
    mySolver->readRestart( restartReadTime, outputDir );
    count = long(simtime/dt);
  } else 
    // Write output at 0 ms
    mySolver->WriteOutput( 0.0, outputDir);


  // Let the game begin
  while (simtime < simulation_time) {
    // Start Stimulus
    if ( count == pulse_start ) {
      vector<int> nodeSet;
      int data;
      switch (stimulusSet[startStimCtr].stimType) {
      case 0:
	inp.open( stimulusSet[startStimCtr].nodeSetFile.c_str() );
	while (inp >> data) nodeSet.push_back( data );
	if (mpicomm->MyPID() == 0)
	  cout << "** Will Stimulate " << nodeSet.size() << " Node(s)\n";
	inp.close();
	myBody->setStimulusCurrent( nodeSet, stimulusSet[startStimCtr].i_stim);
	break;
      case 1:
	myBody->setStimulusCurrent( stimulusSet[startStimCtr].center,
				    stimulusSet[startStimCtr].radius,
				    stimulusSet[startStimCtr].i_stim);
	break;
      case 2:
	myBody->setStimulusCurrent( stimulusSet[startStimCtr].stimAxis,
				    stimulusSet[startStimCtr].limitValue,
				    stimulusSet[startStimCtr].i_stim);
	break;
      case 3:
	myBody->setStimulusCurrent( stimulusSet[startStimCtr].center,
				    stimulusSet[startStimCtr].normal,
				    stimulusSet[startStimCtr].i_stim);     
      } // Switch Statement
      if (mpicomm->MyPID() == 0 ) 
	cout << "@@ Starting Stimulus " << startStimCtr + 1 << endl;
      startStimCtr++;
      if ( startStimCtr == stimulusSet.size() ) startStimCtr--;
      pulse_start = int( stimulusSet[startStimCtr].t_start/dt );
    }
    // Stop Stimulus
    if ( count == pulse_stop ) {
      vector<int> nodeSet;
      int data;
      switch (stimulusSet[stopStimCtr].stimType) {
      case 0:
	inp.open( stimulusSet[stopStimCtr].nodeSetFile.c_str() );
	while (inp >> data) nodeSet.push_back( data );
	inp.close();
	myBody->setStimulusCurrent( nodeSet );
	break;
      case 1:
	myBody->setStimulusCurrent( stimulusSet[stopStimCtr].center,
				    stimulusSet[stopStimCtr].radius );
	break;
      case 2:
	myBody->setStimulusCurrent( stimulusSet[stopStimCtr].stimAxis,
				    stimulusSet[stopStimCtr].limitValue );
      case 3:
	myBody->setStimulusCurrent( stimulusSet[startStimCtr].center,
				    stimulusSet[startStimCtr].normal);
	break;
      } // Switch Statement
      if (mpicomm->MyPID() == 0 ) 
	cout << "@@ Stopping Stimulus " << stopStimCtr + 1 << endl;
      stopStimCtr++;
      if ( stopStimCtr == stimulusSet.size() ) stopStimCtr--;
      pulse_stop  = int( stimulusSet[stopStimCtr].t_stop/dt );
    }
    // Call Solver
    mySolver->Solve();
    simtime += dt;

    // Solve Mechanics
    if ( (count+1)%disp_steps_per_frame == 0) {
      if (mpicomm->MyPID() == 0) 
	cout << "** Mechanics Solve at " << simtime << "(ms)" << endl;
      mySolver->solveMechanics();
    } // check it time to do mechanics

    count ++;      
    // Writing Output
    if ( count%steps_per_frame == 0) {
      mySolver->WriteOutput(simtime, outputDir);
      if (mpicomm->MyPID() == 0 ) 
	cout << "** Writing Output at Time: " << simtime << " (ms)\n";
    }
    // Write output if needed
    if (count%restart_write == 0) {
      mySolver->writeRestart(simtime, outputDir);
      if (mpicomm->MyPID() == 0 ) 
	cout << "@@ Writing Restart at Time: " << simtime << " (ms)\n";
    } 

  } // While time stepping loop
  delete Tissue; delete Fiber;
  delete myMesh; delete myBody; delete mySolver; delete mpicomm;
  MPI_Finalize();
  return 0;

}
